#!/usr/bin/env bash
#
# plot-obs-timelines.sh: Wrapper script, makes plots of observations from a DRM
#
# Usage:
#   plot-obs-timelines.sh [-o TEMPLATE] [-j JSON] DRM
#
# where:
#  -o TEMPLATE  gives the explicit output file template
#               must contain two %s's, such as timelines-%s.%s
#               By default, this is deduced from Sandbox conventions
#  -j JSON      gives the JSON script name
#               By default, this is deduced from Sandbox conventions
#  -D           signals to run python with the debugger on
#
# This is a wrapper around (currently) two python plotter-scripts:
#   plot-timeline.py
#   plot-star-obs-trace.py
# The above scripts make *per-run* plots (one-plot-per-drm). Such plots are
# somewhat analogous to the movies, and to the keepout plots. But both
# of those take more time to run, so they are broken out separately.
#
# The reason for wrapping these is to make a one-call interface for the Makefile
# to use (make obs-timeline is the target). Because we create more than one output
# file, a sentinel text file is created by this script upon successful exit,
# signaling that all files were made OK.
#
# The star-planet config file (SPC) corresponding to a DRM is also needed,
# as is the JSON script. By default, these filenames are generated from the
# DRM filename (see "Conventions" below), but they can also be supplied.
#
# Conventions:
#   The SPC file will be deduced from the DRM filename by the called routines (not
#   by this script).
#   The script file (.json) is deduced from the DRM filename (directory component),
#   but can be supplied explicitly with -j
#   The image output name is generated from the DRM, but can be explicitly given.
#   E.g.:
#       drm  = sims/HabEx_4m_TS_dmag26p0_20180206/drm/777.pkl
#       spc  = sims/HabEx_4m_TS_dmag26p0_20180206/spc/777.spc
#       json = Scripts/HabEx_4m_TS_dmag26p0_20180206.json
# 
# turmon jul 2023
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

echo "${PROGNAME}: Called via:"
echo "$PROGNAME" "$@"

# attempt to give group-write to created files
umask 002

# ep_opts: ensemble path program options
# py_opts: python driver options
py_opts=
out_opt=
script_opt=
while getopts "hDj:o:" opt; do
    case $opt in
	j)
	    # JSON script
	    script_opt="$OPTARG"
	    ;;
	o)
	    # output filename template
	    out_opt="$OPTARG"
	    ;;
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

##############################################################
##
## DRM argument
##

# enforce 1 argument
if [ $# -ne 1 ]; then
    echo "${PROGNAME}: Error: Need exactly one argument (DRM)" >&2
    exit 1
fi

# get the argument
drm="$1"

# check for existence of drm
if [ ! -r "$drm" ]; then
    echo "${PROGNAME}: Given DRM \`$drm' not readable, exiting." >&2
    exit 1
fi    

echo "${PROGNAME}: Timeline graphics: $(basename $drm)"

##############################################################
##
## script argument
##

# if not given: guess script name from DRM
if [ -z "$script_opt" ]; then
    # e.g., consider:
    #    drm = ./sims/HabEx_4m_TS_dmag26p0_20180206/drm/777.pkl
    #    drm = ./sims/Habex_X.exp/s_c1_c2_c3/drm/777.pkl
    script_opt=$(echo "$drm" | sed -e 's|/drm/.*|.json|' -e 's|.*sims/|Scripts/|')
    echo "${PROGNAME}: From DRM, guess script is \`$script_opt'."
fi
if [ ! -r "$script_opt" ]; then
    echo "${PROGNAME}: Inferred script-file \`$script_opt' not readable, exiting." >&2
    exit 1
fi    

##############################################################
##
## output argument
##

# if not given: guess output-template name from DRM
# enforce "out_opt" argument
if [ -z "$out_opt" ]; then
    out_opt=$(echo "$drm" | sed -e 's|.pkl$|-%s.%s|' -e 's|/drm/|/path/|')
    echo "${PROGNAME}: From DRM, guess outfile template is \`$out_opt'."
fi
# ensure exactly 2 "%s" in output file pattern (portable)
if [ $(echo "ZZZ${out_opt}ZZZ" | grep -Fo "%s" | wc -l) -ne 2 ]; then
   echo "${PROGNAME}: Error: -o template must contain two %s markers" >&2
   exit 1
fi

# make a sentinel filename
sentinel=$(printf "$out_opt" obs-timelines txt)
outdir=$(dirname $sentinel)

# ensure the directory for the sentinel file
# (it will also contain the outputs)
if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir" 
    chmod g+w "$outdir"
fi

# ensure we can write the sentinel
if [ ! -w "$outdir" ]; then
    echo "${PROGNAME}: Error: Cannot create files in \`$(outdir)'." >&2
    exit 1
fi


##############################################################
##
## run timeline command
##

cmd_line="util/plot-timeline.py -d -o $out_opt $script_opt $drm"

# say exactly what we are doing
echo "${PROGNAME}: Running:"
echo "$cmd_line"

PYTHONPATH=Local python $py_opts $cmd_line

##############################################################
##
## run obs-trace command
##

cmd_line="util/plot-star-obs-trace.py -o $out_opt $drm"

# say exactly what we are doing
echo "${PROGNAME}: Running:"
echo "$cmd_line"

python $py_opts $cmd_line

##############################################################
##
## signal success
##

# note, umask above makes $sentinel group-writable
cat > "$sentinel" <<EOF
Sentinel file
${PROGNAME}: completed OK by ${USER} on `date`
Outputs in: $out_opt
EOF

