#!/usr/bin/env bash
#
# summarize-int-times.sh: summarize (characterization) integration times
#
# This extracts characterization integration times from a given set of
# DRMs, and finds means and other scalar indexes for them
#
# Usage:
#   summarize-int-times.sh SIM
#
# where SIM is a simulation directory, for example:
#    SIM = sims/SPIE_2023_coro.fam/H6C_CO_DulzE_omniNUV_ideal_20231229
# In this case, SIM/drm must contain several DRM's as .pkl files
#
# This assumes the "mlr" CSV utility is in your PATH. It is at:
#    /proj/exep/rhonda/Sandbox/Tools/bin/mlr
# so you can:
#    export PATH=${PATH}:/proj/exep/rhonda/Sandbox/Tools/bin
# or place the above in your .bash_profile
#
# 
# turmon jan 2024
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
	    sed 's/^# \?//' $(which $0) | awk '/^#/{exit};NR>1{print}'
	    exit 2
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
shift $((OPTIND-1))

# enforce # of arguments
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly one argument." >&2
   exit 1
fi

# this dir is of the form sims/FOO for an Experiment named FOO
simdir="$1"

# enforce mlr presence
command -v mlr > /dev/null || { echo "${PROGNAME}: Fatal: No mlr executable found, run with -h" >&2; exit 1; }

# detect a VENV, and if so, set up to use the correct EXO_PATH
#   (first construct must work if $VIRTUAL_ENV is undefined/unset)
if [ "${VIRTUAL_ENV-undefined}" = undefined ]; then
    echo "${PROGNAME}: Note: python venv not detected, hope that is ok" >&2
fi


# check for existence of the given sim directory
if [ ! -d "$simdir" ]; then
    echo "${PROGNAME}: Error: Given Simulation dir \`$simdir' not readable, exiting." >&2
    exit 1
fi    

drmdir="$simdir/drm"
if [ ! -d "$drmdir" ]; then
    echo "${PROGNAME}: Error: Given DRM dir \`$drmdir' not readable, exiting." >&2
    exit 1
fi    

# Get the number of drms as a sanity mechanism
drmnum=$(ls $drmdir/*.pkl | wc -l)

if [ $drmnum -eq 0 ]; then
   echo "${PROGNAME}: Error: Could not find DRMs within $drmdir" >&2
   exit 1
fi

# output this to stderr so we can still redirect a CSV on stdout
echo "${PROGNAME}: Found $drmnum scripts in \`$drmdir'." >& 2

# DRM tabulation incantation
#  -ns1: output obs-number, seed number, and header
#  -m char_info: insist that it's a char
#  -a Tint:char_info.0.char_time: output char_info[0]['char_time'] as 
#     a csv field, Tint
#  Not present: -a char_info.0.char_status (don't need it)
# mlr incantations:
#   (1) find count, sum, mean, etc., *grouped by seed* 
#   (2) find the mean *across seed*
#   (3) rename fields for readability


if [ 1 -eq 1 ]; then
  util/drm-tabulate.py -sn1 -m char_info -a arrival_time -a Tint:char_info.0.char_time $drmdir/*.pkl | mlr --csv cut -f seed,Tint | mlr --csv stats1 -a count,sum,mean,median -f Tint -g seed | mlr --csv stats1 -a mean --fx 'seed' | mlr --csv --ofmt '%.3f' rename -r '_mean$,_ens,Tint_count,char_count' 
else
  : # dead code, leave for reference
  # util/drm-tabulate.py -ns1 -m char_info -a arrival_time -a Tint:char_info.0.char_time $drmdir/*.pkl | mlr --csv cut -f seed,Tint | mlr --csv stats1 -a mean,stddev,p25,p50,p75 -f Tint -g seed | mlr --csv stats1 -a mean --fx 'seed' | mlr --csv put '$Tint_iqr_mean = $Tint_p75_mean - $Tint_p25_mean' | mlr --csv rename -g -r '_mean$,_ens'
fi




