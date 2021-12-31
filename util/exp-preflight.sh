#!/usr/bin/env bash
#
# exp-preflight: Helper script to prepare an "Experiment" for runs
#
# This wraps a couple of commands that warm the caches for large
# multi-ensemble runs.  This serves multiple purposes:
#  -- Ensures the script really works, lowering the chance of mass failures
#     in a batch run.
#  -- Sets up cache files so that all subsequent runs will find cached
#     results, e.g. for keepout.  We do not want multiple batch runs
#     to simultaneously attempt to write cache files.
#  -- Give a time estimate for a single run so that a time projection
#     for a full run can be made.
# Note: If cache files were generated, the time estimate may be off.  If
# you want to use the time estimate in this case, just re-run the script.
#
# Usage:
#   exp-preflight.sh EXP
#
# where EXP is an experiment name.  By convention, EXP is in the Scripts/ directory
# and ends in '.exp'.  The directory named by EXP must contain several json scripts.
# One of them is an index, and the others will be used for ensemble runs.
#
# Simplest usage:
#   $ util/exp-preflight.sh Scripts/ExampleExp.exp
# 
# turmon may 2020
#
## [end comment block]

# exit-on-error
# (pipefail left out, causes issues when the pipe is closed by "head -1")
set -eu

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

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly one argument" >&2
   exit 1
fi

# this dir is of the form sims/FOO for an Experiment named FOO
expdir="$1"

# check for existence of EXP directory
if [ ! -d "$expdir" ]; then
    echo "${PROGNAME}: Given Experiment dir \`$expdir' not readable, exiting." >&2
    exit 1
fi    

# any integer seed is OK
SEED=777

# find a single script file within EXP
expnum=$(   ls $expdir/*.json | sed -e '/index\.json/d' -e '/\/base\.json/d' | wc -l)
expsingle=$(ls $expdir/*.json | sed -e '/index\.json/d' -e '/\/base\.json/d' | head -1)
expindex=$( ls $expdir/*index.json)

if [ $expnum -eq 0 ]; then
   echo "${PROGNAME}: Error: Could not find scripts within \`$expdir'." >&2
   exit 1
fi
if [ ! -r $expindex ]; then
   echo "${PROGNAME}: Error: No index script (..._index.json) within \`$expdir'." >&2
   exit 1
fi
if [ ! -r $expsingle ]; then
   echo "${PROGNAME}: Error: Found script \`$expsingle' that is not readable." >&2
   exit 1
fi

echo "${PROGNAME}: Found $expnum scripts in $expdir."
echo "${PROGNAME}: Warming cache with $expsingle"

# the one simulation we will add
cmdline="./add-sims.sh $expsingle =$SEED"

# for repeatability, say exactly what we are doing
echo "${PROGNAME}: Running:"
echo "$cmdline"

# run the command
# note: eval propagates error status, so script will exit on error
cmd_start=`date +%s`
echo "----------- preflight run begins -----------"
eval $cmdline
echo "----------- preflight run ends   -----------"
cmd_end=`date +%s`

# single-sim runtime
runtime=$(echo $cmd_end - $cmd_start | bc -l)
# naively forecast total runtime
ens_n=4 # example number of ensemble members
all_sec=$(echo $runtime \* $expnum \* $ens_n | bc -l)
all_min=$(echo $all_sec / 60.0 | bc -l)
all_hour=$(echo $all_sec / 3600.0 | bc -l)
all_day=$(echo $all_sec / 3600.0 / 24.0 | bc -l)
summary=$(printf "%.0f sec = %.2f min = %.2f hours = %.2f days" $all_sec $all_min $all_hour $all_day)

echo "${PROGNAME}: One sim took $runtime seconds."
echo "${PROGNAME}: $ens_n serially-running sims X $expnum scripts projects to: $summary."

# transform script filename -> result dirname
# e.g., Scripts/ExampleExp.exp/s_c1=0.01_c3=0.99.json -> sims/ExampleExp.exp/s_c1=0.01_c3=0.99
temp=$(dirname $expsingle)/$(basename $expsingle .json)
result_dir=$(echo $temp | sed -e 's:Scripts/:sims/:')

# Deliberately not doing rm -r here, so we won't delete useful stuff
echo "${PROGNAME}: Removing generated DRM/SPC files in $result_dir"
rm -f "$result_dir/drm/$SEED.pkl"
rm -f "$result_dir/spc/$SEED.spc"

