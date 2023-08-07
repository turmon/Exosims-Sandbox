#!/usr/bin/env bash
#
# drm-to-movie: Convert DRM to movie
#
# This is a wrapper around the python driver, keepout_path_graphics.py.
# It mainly handles annoying filename transformations, so you only need to
# supply the DRM file name.
#
# Usage:
#   drm-to-movie.sh [-F] [-l N] [-d D] DRM
#
# to make a movie from the given DRM file.  Note, the star-planet config
# file (SPC) corresponding to the DRM is also needed, as is the JSON script.
# These filenames are generated from the DRM filename, see "Conventions" below.
#
# where:
#  -F   signals to make only the final frame.  Don't use -l with -F.
#  -l N is the movie length in years.  Use -l 0 to use the "missionLife" in the script.
#       Default is the keepout_path_graphics default, 0.5 years.
#  -d D is the delta-t between frames in days.
#       Default is the keepout_path_graphics default, 5 days.
#  -c   signals to output certain cumulative graphics.
#  -f   signals to output frame-by-frame graphics.  Don't use with -F.
#  -D   means to load the python debugger while running -- developer use only
#
# Conventions:
#   The SPC file is deduced from the DRM filename.
#   The script file (.json) is deduced from the DRM filename (directory component).
#   The movie output name is generated from the DRM.
#   E.g.:
#       drm   = sims/HabEx_4m_TS_dmag26p0_20180206/drm/163191934.pkl
#       spc   = sims/HabEx_4m_TS_dmag26p0_20180206/spc/163191934.spc
#       json  = Scripts/HabEx_4m_TS_dmag26p0_20180206.json
#       movie = sims/HabEx_4m_TS_dmag26p0_20180206/path/163191934.mp4 (output)
# 
# turmon feb 2018
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# attempt to give group-write to created files
umask 002

# ko_opts: keepout program options
# py_opts: python driver options
ko_opts=
py_opts=
cume_out=0
frame_out=0
while getopts "hDFfcd:l:" opt; do
    case $opt in
	D)
	    # python debugging
	    py_opts="-m pdb"
	    ;;
	F)
	    # final-frame-only
	    # --start 1 means:
	    #  start a fraction of 1.0 through the movie, i.e., at the end
	    ko_opts="$ko_opts --start 1 --length 0"
	    ;;
	f)
	    # frame-by-frame output, 0/1
	    frame_out=1
	    ;;
	c)
	    # cumulative output, 0/1
	    cume_out=1
	    ;;
	d)
	    # delta between frames in days
	    ko_opts="$ko_opts --delta $OPTARG"
	    ;;
	l)
	    # movie length in years
	    ko_opts="$ko_opts --length $OPTARG"
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

# check for existence of input
if [ ! -r "$drm" ]; then
    echo "${PROGNAME}: Given DRM-file \`$drm' not readable, exiting." >&2
    exit 1
fi    

# find parts within input - fragile manipulations
# e.g., consider drm of:
#   ./sims/HabEx_4m_TS_dmag26p0_20180206/drm/163191934.pkl
#   sims/HabEx_4m_TSDD_ETHZnoDD_20190628.exp/s_c1=0.85_c3=0.15/drm/128389321.pkl
drm_num=$(basename $drm .pkl) # the 163191934 part
# script=$(basename $(dirname $(dirname $drm))) # formerly
script=$(echo $drm | sed 's|.*sims/\(.*\)/drm.*|\1|')

# synthesize script name
script_fn=Scripts/$script.json
if [ ! -r "$script_fn" ]; then
    echo "${PROGNAME}: Inferred script-file \`$script_fn' not readable, exiting." >&2
    exit 1
fi    

# synthesize SPC name
spc=sims/$script/spc/$drm_num.spc
if [ ! -r "$spc" ]; then
    echo "${PROGNAME}: Inferred SPC-file \`$spc' not readable, exiting." >&2
    exit 1
fi    

### Output files and directories

# output files
movie=sims/$script/path/$drm_num.mp4
movie_dir=$(dirname $movie)

# ensure the directory
if [ ! -d "$movie_dir" ]; then
    mkdir -p "$movie_dir" 
    chmod g+w "$movie_dir"
fi

# support cumulative-graphics output
if [ $cume_out = 1 ]; then
    cume_dir=sims/$script/path/${drm_num}-cume
    if [ ! -d "$cume_dir" ]; then
	mkdir -p "$cume_dir" 
	chmod g+w "$cume_dir"
    fi
    ko_opts="$ko_opts -c $cume_dir"
fi

# support frame-by-frame output
if [ $frame_out = 1 ]; then
    frame_dir=sims/$script/path/${drm_num}-frames
    if [ ! -d "$frame_dir" ]; then
	mkdir -p "$frame_dir" 
	chmod g+w "$frame_dir"
    fi
    ko_opts="$ko_opts -f $frame_dir"
fi

### Invocation of the command


# formerly:
# pick up Exosims from its default location
#export PYTHONPATH=$(pwd)/EXOSIMS:$(pwd)/Local
# we typically do not need Local to be on PYTHONPATH,
# because keepout_path_graphics sets the SurveyEnsemble
# to the Prototype. Commenting out.
#export PYTHONPATH=Local

# ensure we can at least import EXOSIMS
if python -c 'import EXOSIMS' ; then
    source=$(python -c 'import os, EXOSIMS; print(os.path.dirname(EXOSIMS.__file__))')
    echo "${PROGNAME}: Using EXOSIMS from $source"
    preamble=
else
    echo "${PROGNAME}: Cannot import EXOSIMS from your python. Setting PYTHONPATH to use Sandbox EXOSIMS"
    preamble="PYTHONPATH=EXOSIMS"
fi    

# FIXME 2023-08: delete the below: keepout_path_graphics uses the Prototype SurveyEnsemble
# enforce that no ipyparallel cluster needs to be started
xspecs='!{"ensemble_mode": "init-only"}'

# $ko_opts is unquoted so as to separate the various words within it
# $py_opts, same (typically is empty)
$preamble python $py_opts util/keepout_path_graphics.py -e $ko_opts -x "$xspecs" \
	--drm $drm --spc $spc --movie "$movie" "$script_fn"
