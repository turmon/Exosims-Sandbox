#!/usr/bin/env bash
# 
# json-script-takeout.sh: make a tarfile of the external files referenced by a script
#
# usage:
#   json-script-takeout.sh [-o outfile] script
# 
# where:
#   script: JSON script
#       e.g.: Scripts/HabEx_4m_TSDD_pop100DD_revwt.json
# optionally:
#   outfile: the destination filename for the tar archive
#            (if not supplied, a ".tgz" is added to the script basename)
#
# We make the assumption that the files are referenced by the $EXOSIMS_PARAMS
# variable.
#
# For example,
#
#  json-script-takeout.sh Scripts/HabEx_4m_TSDD_pop100DD_revwt.json
#
# See also: json-attr-edit
#
# Michael Turmon, JPL, Dec 2023

# exit-on-error
set -euo pipefail

progname=`basename $0`
USAGE="Usage: $progname [ -o outfile ] SCRIPT"

# get options
outfile=
while getopts "o:" opt; do
    case "$o" in
	o)    outfile="$opt";;
	[?])  echo "$USAGE" 1>&2
	      exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 1 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi
# ok to get the args
script="$1"

##############################################

# tar output file
if [ -z "$outfile" ]; then
    outfile=$(basename $script .json).tgz
    echo "${progname}: Output to: ${outfile}"
fi

# make created files group-writable
umask 0002

# tempfile
fn=/tmp/takeout-$$.txt

# if command error from here down, hint to the user to check the tempfile
trap "echo ${progname}: Check the list-of-files in: $fn" 0

# assumes files are in the $EXOSIMS_PARAMS lines 
# this line will error-out if some file is not present (due to the ls)
# the "gron" standardizes the file contents
gron $script | grep '$EXOSIMS_PARAMS/' | awk -F= '{print $2}' | tr -d '"; ' | sed 's/$EXOSIMS_PARAMS/./' > $fn

# within a subshell, check the files
(cd $EXOSIMS_PARAMS && ls -l $(cat $fn) > /dev/null)

# (if still here, the above "ls" worked)

# create the tarfile
tar czf "$outfile"  -C $EXOSIMS_PARAMS --files-from="$fn"

# echo
file_count=$(wc -l < $fn)
echo "${progname}: Files archived: $file_count"

# remove the tempfile
rm -f "$fn"

# remove the trap before exit
trap "" 0

echo "${progname}: Done."
