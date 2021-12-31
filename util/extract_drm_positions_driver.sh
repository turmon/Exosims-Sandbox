#!/usr/bin/env bash
#
# extract_drm_positions_driver.sh: extract observer (etc) positions from DRM
#
# This is a wrapper around the python script, extract_drm_positions.py
# It mainly handles annoying filename transformations, so you only need to
# supply the ensemble root directory.
#
# Usage:
#   extract_drm_positions_driver.sh [-D] DRM
#
# to summarize the paths in the given directory DRMS.  Note, the star-planet config
# file (SPC) corresponding to a DRM is also needed, as is the JSON script.
# These filenames are generated from the DRM filename, see "Conventions" below.
#
# where:
#  -D   signals to run python with the debugger turned on
#
# Conventions:
#   The script file (.json) is deduced from the DRM filename (directory component).
#   The output CSV name is generated from the DRM.
#   E.g.:
#       drm  = sims/HabEx_4m_TS_dmag26p0_20180206/drm/NNNN.pkl
#       json = Scripts/HabEx_4m_TS_dmag26p0_20180206.json
#       csv  = sims/HabEx_4m_TS_dmag26p0_20180206/pos/NNNN-position.csv (output)
# 
# turmon oct 2018
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

# ep_opts: extract positions program options (none so far)
# py_opts: python driver options
ep_opts=
py_opts=
while getopts "hD" opt; do
    case $opt in
	D)
	    # python debugging
	    py_opts="-m pdb"
	    ;;
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

drm="$1"

# check for existence of drm
if [ ! -r "$drm" ]; then
    echo "${PROGNAME}: Given DRM \`$drm' not readable, exiting." >&2
    exit 1
fi    

# find parts within input - fragile manipulations
# e.g., consider drm = ./sims/HabEx_4m_TS_dmag26p0_20180206/drm/NNNN.pkl
script=$(basename $(dirname $(dirname $drm))) # the "HabEx_4m_TS_dmag26p0_20180206" part

# synthesize script name
script_fn=Scripts/$script.json
if [ ! -r "$script_fn" ]; then
    echo "${PROGNAME}: Inferred script-file \`$script_fn' not readable, exiting." >&2
    exit 1
fi    

# output file template
output=sims/$script/pos/%s-position.csv
outdir=$(dirname $output)

# ensure the directory
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir" 
    chmod g+w "$outdir"
fi

# need to import EXOSIMS and local modules
export PYTHONPATH=$(pwd)/EXOSIMS:$(pwd)/Local

# $ep_opts is unquoted so as to expand the various words within it
# $py_opts, same (typically is empty)
cmd_line="util/extract_drm_positions.py $ep_opts -p $output $script_fn $drm"

# for repeatability, say exactly what we are doing
echo "${PROGNAME}: Running:"
echo "$cmd_line"

# do it
python $py_opts $cmd_line

