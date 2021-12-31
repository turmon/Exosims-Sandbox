#!/usr/bin/env bash
#
# add-to-exp.sh: perform for multiple "Experiment" runs
#
# This performs batch runs for experiments.
#
# Usage:
#   add-to-exp.sh EXP SEEDS
#
# where EXP is an experiment name, and SEEDS gives the seeds to use (file of seeds,
# or single number giving one seed explicitly).
# In this case, Scripts/EXP must contain several json scripts.
# One of them is an index, and the others will be used for ensemble runs.
# 
# turmon may 2020
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

while getopts "h" opt; do
    case $opt in
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
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
   echo "${PROGNAME}: Error: Need exactly two arguments." >&2
   exit 1
fi

# this dir is of the form Scripts/FOO for an Experiment named FOO
expdir="$1"
# seed specifier, as in the -S option to add-sims.sh
seeds="$2"

# check for existence of the given EXP directory
if [ ! -d "$expdir" ]; then
    echo "${PROGNAME}: Error: Given Experiment dir \`$expdir' not readable, exiting." >&2
    exit 1
fi    

# Get the number of scripts as a sanity mechanism
expnum=$(ls $expdir/*.json | grep -v _index.json | wc -l)

if [ $expnum -eq 0 ]; then
   echo "${PROGNAME}: Error: Could not find scripts within $expdir" >&2
   exit 1
fi

echo "${PROGNAME}: Found $expnum scripts in \`$expdir'."

# transform Scripts/foo.exp -> sims/foo.exp
# logs will go into this directory
simdir=$(echo $expdir | sed 's:Scripts/:sims/:')

# check for existence of simulation results directory
if [ ! -d "$simdir" ]; then
    echo "${PROGNAME}: Could not find Experiment results dir \`$simdir', did you do the preflight?" >&2
    echo "${PROGNAME}: Error: Could not find Experiment results dir \`$simdir', exiting." >&2
    exit 1
fi    

# make a list of the script base names, of the form:
#  tag
#  s_c1=0.01_c3=0.99
#  s_c1=0.05_c3=0.95
#  s_c1=0.10_c3=0.90
#  [...]
# these are used to name the log files as well as a script argument
index_list=$expdir/index.csv
(cd $expdir && ls *.json) | awk 'BEGIN {print "tag"}; /.*_index.json/ {next}; {sub(".json", ""); print $1}' > $index_list

# place for the logs
logdir=$simdir/log
mkdir -m 775 -p $logdir

# options for gnu parallel - rather complex
#   cores reported: aftac1, 32; aftac2, 48; aftac3, 48.
#   need PATH exported because otherwise ssh will just get a vanilla PATH
#   no current driver for PYTHONPATH being exported, but do it anyway
#   --header : --colsep , ==> format of the index file (it currently has no columns,
#                             and in fact multiple columns make the logging really ugly)
#   --files ==> put chatter on stdout/stderr into log files instead of showing on stdout
#   --progress ==> show a progress indication
P_OPTS="--env PATH --env PYTHONPATH --joblog $logdir/Joblog --header : --colsep , --results $logdir --files --progress"
#   --sshdelay ==> required because sshd can't accept connections super-fast
#   -S ==> remote hosts, N/aftac1 means at most 2 jobs to aftac1, etc.
R_OPTS="--workdir . --sshdelay 0.04 -S 2/aftac1,2/aftac2,2/aftac3"

echo "${PROGNAME}: Running $expnum scripts beginning now ($(date))."
parallel $P_OPTS $R_OPTS add-sims.sh -S "$seeds" -P 12 "$expdir/{tag}.json" 1 < $index_list

echo "${PROGNAME}: Finished running ($(date))."

# count the number of 0 exit status messages
exit_ok=$(awk 'NR>1 {ok += ($7 == 0)}; END {print ok}' $logdir/Joblog)

echo "${PROGNAME}: $exit_ok of $expnum scripts finished OK."
echo "${PROGNAME}: Job summary in \`$logdir/Joblog'."

