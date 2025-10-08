#!/bin/bash
# 
# exp-add-sims.sh: Generate a SLURM batch script for an Experiment/Family
#
# Given an Experiment or Family consisting of multiple EXOSIMS scripts, and a
# seed or list of seeds, makes a batch script ready for SLURM submission.
# Features:
#   - Submits an array job over all Scripts in the Family/Experiment
#   - With -j N, allows N-way parallelism over seeds in each Script
#   - With -m TARGET, allows postprocessing with "make reduce", etc.
#
# Simple Usage:
#   exp-add-sims.sh -0 -j 16 -m reduce,html-only Scripts/example.fam Experiment/seed100.txt
#   exp-add-sims.sh -0 -j 1  -m reduce,html-only Scripts/example.fam =777
# 
# Complete Usage:
#   exp-add-sims.sh [-0] [-j N] [-/ S] [[-m target,target,...] ...] [-q|-v] Experiment Seeds
#
# where arguments are:
#   Experiment: a directory name within Scripts/ of an Experiment or Family
#   Seeds: a file of integer seeds, or =S to use a single integer seed S.
#
# and options are:
#   -0: Warm caches before ensemble seed-by-seed runs
#   -j N: Ensemble is created with N-way parallelism
#   -/ S: Divide the "M" scripts in the Experiment into "S" batches, 
#         producing "S" scripts. Default 1. Needed for large M.
#   -m TARGET: Run "make S==... TARGET postprocessing after ensemble
#
# and less-used options:
#   -q: less chatty
#   -v: chatty about job invocation
#   -h: print this help
#
# Note: comma-separated postprocessing steps are done in one invovation
# of "make", and repeated "-m" options cause repeated invocations of
# "make". So -m reduce,obs-timeline-2 -m html-only will:
#   make S=... reduce obs-timeline-2
#   make S=... html-only
# within each script directory. The comma-separated form is quicker 
# to type, but the trailing "-m html-only" ensures that a final HTML 
# page is generated with the timelines in it.
#
# Context:
#   - Your working directory is the /scratch or /scratch-edge Sandbox
#   - Your scripts are in a Family/Experiment in same Sandbox
#   - You should be in the Python VENV you wish to use
#
# The script submits to the "slurm" or "edge" cluster, depending on 
# your working directory (/scratch vs. /scratch-edge). 
#
# Example:
# 
# # we are in /scratch, in a VENV
# $ pwd; echo $VIRTUAL_ENV
# /scratch/exo-yield/Sandbox/hwo
# /scratch-jpl/exo-yield/Sandbox/Python-venvs/exosims-202509-mjt
# 
# # create the script
# $ exp-add-sims.sh -0 -j16 -m reduce,obs-timeline-2 -m html-only Scripts/Test.exp Experiments/seed100.txt 
# (...output...)
# exp-add-sims.sh: Batch script placed in: Scripts/Test.exp/Batch/Run.sh
# 
# # submit the generated file
# $ sbatch Scripts/Test.exp/Batch/Run.sh
# Submitted batch job 4987542
#   
##

# fail-on-error, fail-on-undefined
set -eu

PROGNAME=$(basename $0)

# group-write of new files and dirs
umask 002

############################################################
#
# Argument parsing and validation
#
############################################################

# argument defaults
verbosity=1
jobs=1
split=1
bonus_args=
declare -a make_targets
make_targets=()

while getopts "0vqhj:/:m:" opt; do
    case $opt in
        j)
            # parallel jobs
            jobs="$OPTARG"
            [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: job count set to: $jobs"
            ;;
        /)
            # split the job array
            split="$OPTARG"
            [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: Splitting into batches: $split"
            ;;
        m)
            # make post-processing
            make_targets+=("$(echo "$OPTARG" | sed 's/,/ /g')")
            [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: make postprocessing: $OPTARG"
            ;;
        0)
            # place -0 in add-sims args
            bonus_args="$bonus_args -0"
            ;;
        q)
            verbosity=0
            ;;
        v)
            verbosity=2
            ;;
        h)
            # help text
            sed 's/^# \?//' "$(which "$0")" | awk '/^#$/{exit};NR>1{print}'
            exit 2
            ;;
        \?)
            echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
            exit 2
            ;;
    esac
done
shift $((OPTIND-1))


# enforce 2 arguments
if [ $# -ne 2 ]; then
   echo "${PROGNAME}: Error: Need exactly 2 arguments, use -h for help." >&2
   exit 2
fi

# seeds is a count for now
script_dir="$1"
seeds="$2"

# need an existing Scripts/.../*.{fam,exp}
if [[ "$script_dir" != Scripts/* ]]; then
    echo "${PROGNAME}: Error: first argument must be in Scripts/ " >&2
    exit 1
fi

if [[ "$script_dir" != *.fam && "$script_dir" != *.exp ]]; then
    echo "${PROGNAME}: Error: first argument must be .fam or .exp" >&2
    exit 1
fi

if [ ! -d "$script_dir" ]; then
    echo "${PROGNAME}: Error: first argument must be an existing family/experiment directory" >&2
    exit 1
fi

# validate seeds
if [[ "$seeds" == =* ]]; then
    # single literal seed like =777 or =0
    n_seed=1
    if [[ ${#make_targets[@]} -eq 0 ]]; then
        jobs=1
    fi
    seed_generator="echo ${seeds/=/}"
else
    if [ ! -r "$seeds" ]; then
        echo "${PROGNAME}: Error: Could not read seed list file." >&2
        exit 1
    fi
    n_seed=$(wc -l < "$seeds")
    seed_generator="cat $seeds"
fi
[ "$verbosity" -gt 0 ] && echo "${PROGNAME}: Ensembles are of size $n_seed."

# require a VENV
#   (first construct must work if $VIRTUAL_ENV is undefined/unset)
if [ "${VIRTUAL_ENV+defined}" = defined ]; then
    [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: python venv \`$(basename $VIRTUAL_ENV)' active, will use it"
else 
    echo "${PROGNAME}: Error: Not running in a Python VENV." >&2
    exit 1
fi

############################################################
#
# Set up directories/files
#
############################################################

## # needs help -- commented for now
## tag=final
## if [ -z "$tag" ]; then
##  tag=$(date +%F_%H-%M-%S)
## fi

batch_src_dir="$script_dir/Batch"
batch_run_dir=$(echo "$script_dir" | sed 's|Scripts/|sims/|')
batch_log_dir="$batch_run_dir/Batch/log"

# ensure directories
mkdir -p $batch_src_dir $batch_run_dir $batch_log_dir

# list of all scenarios in the .fam/.exp
scenario_file=$batch_src_dir/Scenario.lst
rm -f $scenario_file

# we may or may not be able to use this index file
scenario_index="$script_dir/s_index.json"
if [ -r "$scenario_index" ]; then
    [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: Using existing scenario index: s_index.json"
    jq -r '.[].script_name' < $scenario_index | sed "s|^|$script_dir/|" | shuf --random-source Experiments/entropy.out > $scenario_file
else
    [ "$verbosity" -gt 0 ] && echo "${PROGNAME}: Generating scenario index"
    # get json's that look like EXOSIMS scripts -- ok for now
    grep -l '"modules"' $script_dir/*.json | shuf --random-source Experiments/entropy.out > $scenario_file
fi

nscript=$(wc -l < $scenario_file)

# we can submit to slurm (jpl) or edge or all, depending on 
# command line options or defaults


# Typical run parameters
# maximum wallclock minutes per single EXOSIMS run/cache
# (pessimistic even for 3.2)
#minutes_per_job=20
#minutes_per_cache=30
# (closer for 3.6.5)
minutes_per_job=8
minutes_per_cache=20

for s in $(seq $split); do
  # submittable batch script file
  if [[ $split -gt 1 ]]; then
    batch_file=$batch_src_dir/Run-$s.sh
  else
    batch_file=$batch_src_dir/Run.sh
  fi
  # ensure batchfile is always executable
  rm -f $batch_file
  touch $batch_file
  chmod +x $batch_file

############################################################
#
# Batch template text [start]
#
############################################################
cat << EOF >> $batch_file
#!/bin/bash
#SBATCH -A exo-yield
#SBATCH -J Exo-AddSims
#SBATCH -o $batch_log_dir/%A/stdout-%a.log
#SBATCH -e $batch_log_dir/%A/stderr-%a.log
#SBATCH -p compute     # --partition
#SBATCH -N 1
#SBATCH -n $jobs
#SBATCH --mem=$((4 + $jobs * 2 / 3))G      # no decimals (2025/09: 16 jobs ~= 5.4GB)
#SBATCH --time=$(($minutes_per_cache + $minutes_per_job * $n_seed / $jobs))   # [minutes]
#SBATCH --mail-type=all
#SBATCH --array=$s-$nscript:$split
: No-op signals end of SBATCH directives

############################
# THIS IS A GENERATED FILE #
############################

# 
# Note: Provided you have "$script_dir" (including "$batch_src_dir")
# in your cluster scratch directory, this batch file is cluster-agnostic:
#  - You may submit this file to a specific cluster with sbatch -M
#  - Without -M, sbatch will submit to the default cluster for that node
#

# fail-on-error, fail-on-undefined
set -eu

# name of this script at runtime (not creation time)
PROGNAME=\$(basename \$0)

# Parameters defined at creation time
#   Batch script generated by $PROGNAME for $USER
#   Date: $(date)
#   VENV: $VIRTUAL_ENV
#   Scenario: $script_dir
#   Scenario size: $nscript
#   Ensemble: $seeds
#   Parallelism: $jobs

echo -n "\${PROGNAME}: *** Begin "
date

# Output useful parameters to the logfile
#   parallelism is also in the "jobs" variable
echo "\${PROGNAME}: TASK_ID within array: \$SLURM_ARRAY_TASK_ID"
echo "\${PROGNAME}: Running on cluster: \$SLURM_CLUSTER_NAME"
echo "\${PROGNAME}: Running on node(s): \$SLURM_JOB_NODELIST"
echo "\${PROGNAME}: Parallelism here: \$SLURM_NTASKS"


############################################################
# Configure the runtime environment
############################################################

# add gnu parallel to path
module load parallel

# which cluster are we are on
if [[ "\$SLURM_CLUSTER_NAME" == "slurm" ]]; then
    cluster=scratch
elif [[ "\$SLURM_CLUSTER_NAME" == "edge" ]]; then
    cluster=scratch-edge
else
    echo "\${PROGNAME}: Error: Unknown cluster" >&2
    exit 1
fi
# set up the correct EXOSIMS_PARAMS so param files are on local /scratch
export EXOSIMS_PARAMS=/\$cluster/exo-yield/Sandbox/Parameters/EXOSIMS_external_files 

# enter the correct VENV
# this allows the VENV to be dynamic at runtime based on cluster
cluster_venv=/\$cluster/exo-yield/Sandbox/Python-venvs/$(basename $VIRTUAL_ENV)
if [ ! -f \$cluster_venv/bin/activate ]; then
    echo "\${PROGNAME}: Error: venv \$cluster_venv not present on \$cluster" >&2
    echo "\${PROGNAME}: Error: Was it installed?" >&2
    exit 1
fi
source \$cluster_venv/bin/activate

# turn off progress bars for batch submission
export TQDM_DISABLE=1

############################################################
# Add simulations with add-sims
############################################################

if [[ $verbosity -gt 1 ]]; then
    set -x
fi

# this_script is ID'th line in the scenario list
this_script=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $scenario_file)

# -b: batch logging mode, even if just one seed
add-sims.sh -b $bonus_args -j $jobs \$this_script $seeds

# link the correct entry from the "pile" of log files generated by slurm to 
# conventional/findable locations under this_sim/log/slurm
this_sim=\$(echo "\$this_script" | sed -e 's|Scripts/|sims/|' -e 's|[.]json||')
this_sim_log_dir="\${this_sim}/log/slurm"
mkdir -p "\$this_sim_log_dir"
ln -f "$batch_log_dir/\${SLURM_ARRAY_JOB_ID}/stdout-\${SLURM_ARRAY_TASK_ID}.log" "\$this_sim_log_dir"
ln -f "$batch_log_dir/\${SLURM_ARRAY_JOB_ID}/stderr-\${SLURM_ARRAY_TASK_ID}.log" "\$this_sim_log_dir"

############################################################
# "make" postprocessing
############################################################

# re-package the arguments back into a literal array
# (below can be modified after the fact, just adjust make_n below)
slurm_make_targets=($(printf "'%s' " "${make_targets[@]}"))
if [[ ${#make_targets[@]} -gt 0 ]]; then
    # 0-based indexing
    make_n=$((${#make_targets[@]} - 1))
    for make_i in \$(seq 0 \$make_n); do
        make_target="\${slurm_make_targets[\$make_i]}"
        echo "\${PROGNAME}: Running make of target(s): \$make_target"
        # make_target below is unquoted
        make -j $jobs S=\$this_script \$make_target
        echo "\${PROGNAME}: Finished make of target(s): \$make_target"
    done
fi

echo -n "\${PROGNAME}: *** End "
date
# (end of generated file)
EOF
############################################################
#
# Batch template text [end]
#
############################################################
echo "${PROGNAME}: Batch script placed in: $batch_file"
done


############################################################

[ "$verbosity" -gt 0 ] && echo "${PROGNAME}: Exit OK."

