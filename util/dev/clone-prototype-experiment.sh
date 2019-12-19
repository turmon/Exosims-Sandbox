#!/bin/sh
#
# clone-prototype-experiment: Copy a skeleton directory structure to create a new experiment
#
# Usage:
#  clone-prototype-experiment NEW_EXPERIMENT
#
# to make a skeleton directory structure named NEW_EXPERIMENT, suitable for
# Exosims ensemble runs.
#
# turmon jan 2018
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# option processing - not much used now
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

EXP="$1"

# enforce run directory
if [ ! -d Prototype ]; then
   echo "${PROGNAME}: Error: Need to run from correct directory" >&2
   exit 1
fi

# don't overwrite
if [ -e $EXP ]; then
    echo "${PROGNAME}: Error: Experiment \`${EXP}' already exists."
    exit 1
fi

############## actual beginning

# files/dirs created by any child commands will be group-writable
umask 0002

# clone Prototype into the given $EXP
rsync -av Prototype/ $EXP

# make $EXP and dirs beneath it setgid; its group will be 'exosims'
find $EXP -type d -exec chmod g+sw {} \;

# give group-write to files below $EXP
# chmod g+w $EXP/README
find $EXP -type f -exec chmod g+rw {} \;


