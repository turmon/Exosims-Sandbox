#!/usr/bin/env bash
#
# ens-path-summary: Convert DRM ensemble to plot set summarizing paths
#
# This is a wrapper around the python driver, `ens-path-graphics`.
# It mainly handles annoying filename transformations, so you only need to
# supply the ensemble root directory.
#
# ## Usage:
#   `ens-path-summary.sh [-a] [-D] DRMS`
#
# to summarize the paths in the given directory DRMS.  Note, the star-planet config
# file (SPC) corresponding to a DRM is also needed, as is the JSON script.
# These filenames are generated from the DRM filename, see "Note" below.
#
# ### Options:
# 
# * `-a` signals to omit the animation of the 3d path summary.
# * `-D` signals to run python with the debugger on
#
# ### Conventions:
#   The SPC file will be deduced from the DRM filename(s) by `ens-path-graphics`.
#   The script file (`.json`) is deduced from the DRM filename (directory component).
#   The image output name is generated from the DRM.
#   For example,
# ``` shell
#       drms  = sims/HabEx_4m_TS_dmag26p0_20180206/drm
#       spc   = sims/HabEx_4m_TS_dmag26p0_20180206/spc/???.spc (one will be selected)
#       json  = Scripts/HabEx_4m_TS_dmag26p0_20180206.json
#       movie = sims/HabEx_4m_TS_dmag26p0_20180206/path-ens/*.* (output)
# ```
# 
# turmon feb 2018
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

# ep_opts: ensemble path program options
# py_opts: python driver options
ep_opts=
py_opts=
while getopts "hDa" opt; do
    case $opt in
	D)
	    # python debugging
	    py_opts="-m pdb"
	    ;;
	a)
	    # skip animation
	    ep_opts="-a"
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

drmdir="$1"

# check for existence of drm directory
if [ ! -d "$drmdir" ]; then
    echo "${PROGNAME}: Given DRM-dir \`$drmdir' not readable, exiting." >&2
    exit 1
fi    

# find parts within input - fragile manipulations
# e.g., consider drmdir = ./sims/HabEx_4m_TS_dmag26p0_20180206/drm
#                drmdir = ./sims/Habex_X.exp/s_c1_c2_c3/drm
script=$(echo $(dirname $drmdir) | sed 's:.*sims/::') # the "HabEx_4m_TS_dmag26p0_20180206" part

# synthesize script name
script_fn=Scripts/$script.json
if [ ! -r "$script_fn" ]; then
    echo "${PROGNAME}: Inferred script-file \`$script_fn' not readable, exiting." >&2
    exit 1
fi    

# output file template
output=sims/$script/path-ens/path-%s.%s
outdir=$(dirname $output)

# ensure the directory
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir" 
    chmod g+w "$outdir"
fi

# 
export PYTHONPATH=$(pwd)/EXOSIMS:$(pwd)/Local

# enforce that no cluster needs to be started
#xspecs='!{"ensemble_mode": "init-only"}'

# $ep_opts is unquoted so as to expand the various words within it
# $py_opts, same (typically is empty)
cmd_line="util/ens-path-graphics.py $ep_opts --outfile $output $script_fn $drmdir/*.pkl"

# for repeatability, say exactly what we are doing
echo "${PROGNAME}: Running:"
echo "$cmd_line"

# do it
python $py_opts $cmd_line

