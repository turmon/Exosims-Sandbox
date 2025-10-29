#!/usr/bin/env bash
#
# add-sims: Perform EXOSIMS simulations, accumulating runs into subdirectories
#
# This is a wrapper around the python driver, sandbox_driver.py.
#
# Typical Usage:
#   add-sims.sh [-S | -Z] SCRIPT SEEDS
# 
# Most General Usage:
#   add-sims.sh [-q] [-A | -S | -Z] [-j JOBS] [-v VERB] [-x SCRIPT] [-O OPTS] SCRIPT SEEDS
#
# Uses the JSON script SCRIPT and performs a series of parallel runs given by SEEDS.
# * The runs are done on the local computer, from a pool that depends on machine
# (four less than the number of cores, so for example: mustang2/3 = 68).
# * If just one SEED is given, Exosims output is sent directly to standard output - 
# useful for checking validity.  If more than one SEED, the Exosims output is sent
# to a set of log files (sims/SCRIPT/log/*), and a summary of current job status
# is sent to standard output.
# * To repeat an error-causing run interactively, paste the command-line printed by 
# this script (the one that calls sandbox_driver.py) into your terminal,
# and add "--interactive" to the options.
#
# Typical Usage:
#   (a) Run a machine-dependent number of jobs spanning 64 seeds on the local machine
#     add-sims.sh Scripts/ExampleScript.json Experiments/seed64.txt 
#   (b) Run 100 jobs on mustang* and aftac*, across 100 deterministic seeds
#     add-sims.sh -Z Scripts/ExampleScript.json Experiments/seed100.txt 
#   (c) Do a test run for a single arbitrary seed (we use 777):
#     add-sims.sh Scripts/ExampleScript.json =777
#   (d) Write cache files for a new script, but do not run a sim:
#     add-sims.sh Scripts/ExampleScript.json =0
#
# Arguments:
#   SCRIPT -- a JSON script suitable for EXOSIMS
#   SEEDS  -- either:
#             (1) if a plain integer, that number of randomly-chosen initial seeds.
#             (2) if of the form =SEED, where SEED is an integer, that single integer 
#                 is the seed. If SEED is 0, only cache warming is done (no run_sim).
#             (3) else, it is a filename giving a list of integer seeds, one per line.
#
# Options:
#   -h        => show this help message and exit.
#   -0        => generate caches first, then do the requested SEEDS
#   -S        => run using Speedy's (mustang2/3/4 -- 24 + 24 + 16 jobs = 64 jobs)
#   -Z        => run using all (mustang2/3/4 + aftac1/2/3 -- total of 100 jobs)
#   -c        => chatty console output (for debugging) even if #SEEDS > 1
#
# Less-used Options:
#   -j JOBS   => runs only JOBS parallel jobs (not used with -A)
#   -v VERB   => set verbosity to VERB (0=quiet or 1=verbose)
#   -q        => quiet object creation
#   -P        => do not show the gnu parallel progress bar (good for scripted batch jobs)
#   -x SCRIPT => an extra scenario-specific script loaded on top of the argument SCRIPT
#   -p PATH   => EXOSIMS path is PATH instead of the default in your environment.
#                If you are in a Python venv, this will be detected and you should not 
#                need to specify -p.
#                Using PATH=@ abbreviates the command-line default (what you get from
#                "python -c import EXOSIMS").
#   -O OPTS   => output options, a string of comma-separated tags telling
#       which EXOSIMS variables should be written, and to which files.
#       See the notes in the python driver file for how it works.
#       Default: OPTS = 'drm:pkl,spc:spc'  -- this indicates the DRM
#       is stored as a pickle with extension .pkl, and the planet parameters
#       are stored as a pickle with extension .spc.
# 
# turmon oct 2017, 2018, 2020, 2023
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# files/dirs created by any child commands will be group-writable
umask 0002

# EXOSIMS driver program
# see also <EXOSIMS>/run/run_ipcluster_ensemble.py
DRIVER=Local/sandbox_driver.py

# option processing
# 0: EXOSIMS path
# EXO_PATH=$(pwd)/EXOSIMS # (formerly)
EXO_PATH=@
EXO_PATH_SET=no
# 1: driver options
OUT_OPT=
SEED_OPT=
# 2: generic options to the driver
DRIVER_OPT=
# 3: "dispatcher" that provides shell-level parallelism
DISPATCHER=parallel
# by default, use 4 less than the number of cores
DISPATCHER_JOBS="--jobs -4"
DISPATCHER_BATCH_PROGRESS="--progress"
BATCH=0
CHATTY=0
ECHO_CMD=0
ALL=0
CACHE=0
while getopts "0AHSZDhbcj:p:x:qv:PO:" opt; do
    case $opt in
	0)
	    # generate cache files first
	    CACHE=1
	    ;;
	A)
	    # use remote machines
	    ALL=1
	    ;;
	S)
	    # use only "speedy" remote machines
	    ALL=2
	    ;;
	Z)
	    # use even more remote machines
	    ALL=3
	    ;;
	D)
	    # fewer jobs on all remote machines (for Dakota)
	    ALL=4
	    ;;
	H)
	    # fewer jobs on mustang2/3 (for Dakota, temp. 2024-05)
	    ALL=5
	    ;;
	b)
	    # batch mode output to log-files
	    BATCH=1
	    ;;
	c)
	    # chatty output to console
	    CHATTY=1
	    ;;
	j)
	    # explicit job number
	    DISPATCHER_JOBS="--jobs $OPTARG"
	    ;;
	p)
	    # alternate EXOSIMS path
	    EXO_PATH="$OPTARG"
            EXO_PATH_SET=yes
	    ;;
	x)
	    # x-tra specs - ensure it is quoted
	    DRIVER_OPT="${DRIVER_OPT} --xspecs $OPTARG"
	    ;;
	q)
	    # quiet object creation
	    DRIVER_OPT="${DRIVER_OPT} --quiet"
	    ;;
	v)
	    # verbosity 0/1
	    DRIVER_OPT="${DRIVER_OPT} --verbose $OPTARG"
	    ;;
	P)
	    # suppress progress bar
            DISPATCHER_BATCH_PROGRESS=""
	    ;;
	O)
	    # output options
	    OUT_OPT="--outopts $OPT_ARG"
	    ;;
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
            echo "${PROGNAME}: Note: By default, using python at" $(which python)
	    exit 2
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
# save the original args for possible -0
orig_args=( "$@" )
shift $((OPTIND-1))

# enforce 2 arguments
if [ $# -ne 2 ]; then
   echo "${PROGNAME}: Error: Need exactly two arguments" >&2
   exit 1
fi

SCRIPT="$1"
SEEDS="$2"

##
## Recursive call if -0
##
if [[ $CACHE -eq 1 && "$SEEDS" != "=0" ]]; then
    # the && above prevents infinite recursion
    echo "${PROGNAME}: Generating caches (=0) before main run ($SEEDS)."
    echo "${PROGNAME}: { ---------------"
    # array construct: "all but the last argument" (the seed)
    ${PROGNAME} "${orig_args[@]::${#orig_args[@]}-1}" =0
    echo "${PROGNAME}: } ---------------"
    echo "${PROGNAME}: Returned OK from generating caches."
fi

##
## Handle SEEDS
##

# handle the 3 cases for the SEEDS argument
if expr "$SEEDS" : '^=[0-9][0-9]*$' > /dev/null; then
    # literal seed -- echoed into xargs at invocation time
    SEED_USED=$(echo $SEEDS | sed 's/=//')
    SEED_INPUT="echo $SEED_USED"
    NRUN=1
elif expr "$SEEDS" : '^[0-9][0-9]*$' > /dev/null; then
    # it is a number-of-runs
    # draw $SEEDS 9-digit random integers from awk
    SEED_INPUT="seq $SEEDS | awk 'BEGIN {srand()}; {print int(rand()*999999999+1)}'"
    NRUN=$SEEDS
else
    if [ ! -r "$SEEDS" ]; then
	echo "${PROGNAME}: Error: SEEDS argument looks like a file, but is not readable." >&2
	exit 1
    fi
    # seeds are in file
    SEED_INPUT="cat $SEEDS"
    NRUN=$(wc -l < $SEEDS)
fi

echo "${PROGNAME}: Performing $NRUN EXOSIMS run(s)."

# collapse NRUN, BATCH, CHATTY into just BATCH
if [ $NRUN -gt 1 ]; then
    BATCH=1
fi
# allow over-ride by -c
if [ $CHATTY == 1 ]; then
    BATCH=0
fi

##
## EXOSIMS path
##

# detect a VENV, and if so, set up to use the correct EXO_PATH
#   (first construct must work if $VIRTUAL_ENV is undefined/unset)
if [ "${VIRTUAL_ENV+defined}" = defined -a "$EXO_PATH_SET" = no ]; then
    EXO_PATH=@
    echo "${PROGNAME}: Note: python venv \`$(basename $VIRTUAL_ENV)' active, will use its EXOSIMS"
else 
    echo "${PROGNAME}: WARNING: Not running in a Python VENV." >&2
fi
# ensure EXOSIMS path
# tricky:
#  - "." is always first in sys.path, and (alas!) ./EXOSIMS exists, so "import EXOSIMS" can pick it up,
#    we import EXOSIMS.Prototypes to keep this from happening
#  - mustang* default python also has its own outdated EXOSIMS (stuck with this)
#  - ~/.local/lib/.../site-packages can have or point to an EXOSIMS
# So, we echo it back here, and *force* it via PYTHONPATH later on
if [ $EXO_PATH = @ ]; then
    # useful when using venv's: this specifies the venv's EXOSIMS as the one to use
    EXO_PATH=$(python -c 'import EXOSIMS.Prototypes; import os.path as p; print(p.dirname(p.dirname(p.dirname(EXOSIMS.Prototypes.__file__))))')
    echo "${PROGNAME}: Note: EXOSIMS path set to \`${EXO_PATH}'."
fi
if [ ! -r $EXO_PATH/EXOSIMS/__init__.py ]; then
   echo "${PROGNAME}: Error: EXOSIMS path seems invalid" >&2
   exit 1
else
   echo "${PROGNAME}: Note: EXOSIMS is found at \`${EXO_PATH}'."
fi
# Detect and warn if not in a git repo
make_diff=0
if git -C "$EXO_PATH" status >& /dev/null; then
    echo "${PROGNAME}: Note: EXOSIMS is in a repository, branch: $(git -C "$EXO_PATH" rev-parse --abbrev-ref HEAD)"
    if git -C "$EXO_PATH" status | grep -q 'working tree clean'; then
	echo "${PROGNAME}: Note: EXOSIMS working tree is clean."
    else
	echo "${PROGNAME}: WARNING: EXOSIMS has uncommitted changes in working tree." >&2
	echo "${PROGNAME}: Note: Will place \`git diff' in sim log directory." >&2
	make_diff=1
    fi
else
   echo "${PROGNAME}: WARNING: EXOSIMS is not in a git repository" >&2
fi
    

##
## Script/sim name
##

if [ ! -r $SCRIPT ]; then
    echo "${PROGNAME}: Error: Could not read the script file \`${SCRIPT}'"
    echo "${PROGNAME}: Fatal: Could not read the script file \`${SCRIPT}'" >&2
    exit 1
fi
# basename for sim results is derived from script name: if it's in
# the "special" Scripts/... dir, just strip Scripts/ and .json out, 
# else use the basename less .json
if [[ "$SCRIPT" =~ ^Scripts/.* ]]; then
    sim_base=$(echo "$SCRIPT" | sed -e 's:Scripts/::' -e 's:\.json$::')
else
    # (not the conventional Sandbox way)
    sim_base=$(basename $SCRIPT .json)
fi
sim_base=sims/$sim_base
logdir=$sim_base/log/console
pardir=$sim_base/log/gnu_parallel
mkdir -p "$logdir" "$pardir"

# optionally put "git diff" in the logs
if [ $make_diff == 1 ]; then
    diffdir=$sim_base/log/git-diff/$(date -Iminutes)
    mkdir -p "$diffdir"
    git -C "$EXO_PATH" status    > "$diffdir/status.out"
    git -C "$EXO_PATH" diff      > "$diffdir/diff.out"
    echo git -C "$EXO_PATH" diff > "$diffdir/diff.cmd"
fi

##
## Prepare for run
##

# sanity check: DIR/Module/__init__.py should exist, for each
# colon-separated DIR in PYTHONPATH, if you want to use "import Module".
export PYTHONPATH=${EXO_PATH}:$(pwd)/Local

# 6/2019: sometimes we run out of threads, if many parallel instances?
export OPENBLAS_NUM_THREADS=8

# 10/2025: Avoid error thrown by (PlanetPhysicalModel/Forecaster.py) when opening
# an HDF5 data file, during initialization, when many jobs are started on gattaca2:
# Unable to synchronously open file (unable to lock file, errno = 37, error message = 'No locks available'
export HDF5_USE_FILE_LOCKING="FALSE"

# options for gnu parallel - rather complex
# 1: multiple machines
#   cores reported: aftac1, 32; aftac2, 48; aftac3, 48.
#   Other considerations: mustang5 exists, but python not installed there
# 2: stdout and logging
#   --files ==> put chatter on stdout/stderr into log files instead of showing on stdout
#   --progress ==> show a progress indication
# 3: environment variables
#   PATH: export it, else ssh gets a vanilla PATH without matlab, our python, etc.
#   EXOSIMS_PARAMS: variable is used in JSON scripts for .fits files, etc.
#   TQDM_DISABLE: allow to kill the python "tqdm" package progress bar
#   note: parallel does not object if you give --env FOO when FOO is unset
par_env_var_names=(PATH PYTHONPATH VIRTUAL_ENV EXOSIMS_PARAMS OPENBLAS_NUM_THREADS HDF5_USE_FILE_LOCKING TQDM_DISABLE)
declare -a PAR_ENV_OPTS
for name in "${par_env_var_names[@]}"; do
    PAR_ENV_OPTS+=("--env" "$name")
done
# file logging
#   use specific % abbreviations to prevent locale issues
date_time=$(date +%Y-%m-%d_%H:%M:%S)
joblog_file=${pardir}/JobLog_${date_time}.txt
if [ $BATCH == 1 ]; then
    PAR_LOG_OPTS="--joblog $joblog_file --results ${logdir}/{}.out $DISPATCHER_BATCH_PROGRESS --files"
    # no progress bars in logfiles
    export TQDM_DISABLE=1
else
    # --line-buffer => allow lines of jobs to intermix
    # no --progress or --results: allow regular chatter to show progress
    PAR_LOG_OPTS="--joblog $joblog_file --line-buffer"
    #PAR_LOG_OPTS="--line-buffer"
fi
#   --sshdelay ==> required because sshd can't accept connections super-fast
#   -S ==> remote hosts, N/aftac1 means at most 2 jobs to aftac1, etc.
#   The special --workdir value . uses the current working dir
#   2024-jan: mustangs seem to not like --controlmaster
par_ssh_remote_base="--workdir ."
if [ $ALL == 1 ]; then
    PAR_SSH_OPTS="${par_ssh_remote_base} --sshdelay 0.06 -S 18/aftac1,23/aftac2,23/aftac3"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 2 ]; then
    PAR_SSH_OPTS="${par_ssh_remote_base} --sshdelay 0.12 -S 26/mustang2,26/mustang3,18/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 3 ]; then
    PAR_SSH_OPTS="${par_ssh_remote_base} --sshdelay 0.12 -S 12/aftac1,12/aftac2,12/aftac3,26/mustang2,24/mustang3,18/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 4 ]; then
    PAR_SSH_OPTS="${par_ssh_remote_base} --sshdelay 0.12 -S 6/aftac1,6/aftac2,6/aftac3,12/mustang2,12/mustang3,8/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 5 ]; then
    PAR_SSH_OPTS="${par_ssh_remote_base} --sshdelay 0.12 -S 12/mustang2,12/mustang3"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
else
    # local host only
    PAR_SSH_OPTS="-S :"
fi

# this can enable debugging a failing invocation
if [ $CHATTY == 1 ]; then
    set -x
fi
# Taking this apart, the most critical elements are:
#   SEED_INPUT | DISPATCHER [opts] python DRIVER [opts] SCRIPT
# where:
#   SEED_INPUT = input stream of newline-separated seeds, for parallel {}
#   DISPATCHER [opts] = gnu parallel (which is given appropriate options)
#   python DRIVER [opts] SCRIPT = python sim-runner, with opts and 1 arg ($SCRIPT)
eval $SEED_INPUT | $DISPATCHER $DISPATCHER_JOBS "${PAR_ENV_OPTS[@]}" $PAR_LOG_OPTS $PAR_SSH_OPTS \
python $DRIVER $DRIVER_OPT --seed {} $OUT_OPT "$SCRIPT"
