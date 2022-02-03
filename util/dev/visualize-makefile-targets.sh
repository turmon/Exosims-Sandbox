#!/bin/bash
#
# visualize-makefile-targets: generate graphical summaries of some Makefile targets
#
# Usage:
#  visualize-makefile-targets [ -S SCRIPT ] [ -t TARGET ] ROOTDIR
#
# to make some plots and place them in the given directory.
#
# The file Scripts/SCRIPT.json (and more importantly, the directory sims/SCRIPT/drm/,
# and a small ensemble of runs) must exist.
#
# By default, if -S is not given, a certain SCRIPT will be used. Other convenience
# script abbreviations:
#   SCRIPT == BASE -> same default script
#   SCRIPT == FAM  -> one ensemble from a 'family' of ensembles
#   SCRIPT == EXP  -> an experiment
#
# A default list of Makefile targets is summarized. To add more, use one or more
# `-t TARGET' options -- a separate summary is generated for each -t you give.
# If you want two targets to appear on the same graph output, use quotes:
# `-t "target1 target2 target3" '
#
# If you give -T, the default target list is cleared and more can be added
# with -t. This is for Experiments which don't have the usual targets.
#
# Typical Usage:
#   util/dev/visualize-makefile-targets.sh ~/basic
#   util/dev/visualize-makefile-targets.sh -S FAM ~/family
#   util/dev/visualize-makefile-targets.sh -S EXP -T -t exp-reduce -t exp-html ~/experiment
#
# turmon jan 2022
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# make2graph program
#  a short C program from github with only standard dependencies
MAKE2GRAPH=~turmon/unix/bin/x86_64/make2graph

# default script names
# the .exp is special, its ensembles have no = signs
DEF_BASE_SCRIPT=ATSA_CODulzE_base_O2_EZimp_20220118
DEF_FAML_SCRIPT=test.fam/Subdirectory_Test
DEF_EXP_SCRIPT=HabEx_4m_coro_ljsDDPC_X_20181120.exp
SCRIPT=

# initial target list
#   the quoted group will all appear in the same output graph
TARGETS=( path-ensemble "reduce graphics html" keepout-2 path-movie-2 obs-timeline-2 )

# option processing
while getopts "hS:t:T" opt; do
    case $opt in
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
	    exit 2
	    ;;
	S)
	    # supply the script
	    SCRIPT="$OPTARG"
	    ;;
	T)
	    # zero the targets
	    TARGETS=()
	    ;;
	t)
	    # another target for the target-list
	    TARGETS+=("$OPTARG")
	    ;;
	\?)
	    echo "${PROGNAME}: Invalid option, exiting.  Try -h." >&2
	    exit 2
	    ;;
    esac
done
shift $((OPTIND-1))

# enforce the run directory
if [ ! -r Makefile ]; then
    echo "${PROGNAME}: Error: Did not find \`Makefile'" >&2
    exit 1
fi

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly one argument" >&2
   exit 1
fi
DIR="$1"

############## process/massage arguments

# script choice
if [ -z "$SCRIPT" ]; then
    SCRIPT=$DEF_BASE_SCRIPT
elif [ BASE = "$SCRIPT" ]; then
    SCRIPT=$DEF_BASE_SCRIPT
elif [ FAM = "$SCRIPT" ]; then
    SCRIPT=$DEF_FAML_SCRIPT
elif [ EXP = "$SCRIPT" ]; then
    SCRIPT=$DEF_EXP_SCRIPT
fi
echo "${PROGNAME}: Making ${#TARGETS[@]} graphs for the script \`$SCRIPT'."

# enforce that data is there, otherwise the make usually is nonsense
if [ ! -d sims/$SCRIPT/drm ]; then
    echo "${PROGNAME}: Warning: Did not find drms in \`sims/$SCRIPT/drm'" >&2
    #exit 1
fi

# files/dirs created by any child commands will be group-writable
umask 0002

# output dir
if [ ! -d $DIR ]; then
    echo "${PROGNAME}: Making destination directory \`${DIR}'"
    mkdir -p $DIR
else
    echo "${PROGNAME}: Destination \`${DIR}' exists, will overwrite some contents."
fi

############## actual beginning

# used to substitute the (arbitrary) script name with the string Script
# in the summary output -- but only at the bottom (ensemble) level
SCRIPT_base=$(basename $SCRIPT)
# different sub depending on whether the experiment or the script should
# be tidied up
if [[ $SCRIPT =~ .*\.exp ]]; then
    SCRIPT_SUB=example.exp
else
    SCRIPT_SUB=Script
fi

# dpi=300          -- needed to get good resolution on png for text boxes
# size=50 (inches) -- sets the maximum size, the graph will be scaled down
#                     if needed to fit in this square. LaTeX/graphicx does not
#                     like graphics bigger than about 400 inches, but dot
#                     does not want the size constraint much bigger than this.
DOT_ARGS="-Nshape=rect -Gfontsize=20 -Gfontname="courier:bold" -Gdpi=300 -Gsize=50 -Edir=back"

for t in "${TARGETS[@]}"; do
    # short filename
    name=$(echo $t | awk '{print $1}')
    # output file root
    outfile=$DIR/dependencies-$name
    # progress indication
    echo "${PROGNAME}: Images for \`make S=$SCRIPT $t' -> \`${outfile}*'"

    # produce the DAG in .dot format (text)
    make -Bnd S=$SCRIPT $t | \
	$MAKE2GRAPH | \
	grep -v Makefile | \
	sed "s|$SCRIPT_base|$SCRIPT_SUB|g" > $outfile-orig.dot

    # mangle the colors (still a text file)
    #   if sims/ is present, is a file target (implying certain colors);
    #   otherwise, is a .PHONY target (other colors)
    < $outfile-orig.dot \
      awk '/sims\// {sub("\"red\"", "\"blue\""); sub("\"green\"", "\"darkgreen\""); print $0; next}; {sub("\"red\"", "\"gold\""); print}' \
      > $outfile-mod.dot

    # produce the graphical DAG
    dot -Tpng -Glabel="make $t" $DOT_ARGS < $outfile-mod.dot > $outfile.png
    dot -Tpdf -Glabel="make $t" $DOT_ARGS < $outfile-mod.dot > $outfile.pdf
done
