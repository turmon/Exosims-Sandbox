#!/bin/sh
#
# add-sims: Perform EXOSIMS simulations, accumulating runs into subdirectories
#
# This is a wrapper around the python driver, ipcluster_ensemble_jpl.py.
#
# Usage:
#   add-sims [-3] [-q] [-v VERB] [-S SEEDS] [-P PAR] [-x SCRIPT] [-e] [-E ADDR] [-O OPTS] SCRIPT N
#
# to use the JSON script SCRIPT (a file) and perform N (an integer) runs.
#
# Options:
#   -h        => show this help message and exit.
#   -q        => quiet object creation
#   -3        => run the EXOSIMS sim-runner with python3
#   -p PATH   => EXOSIMS path is PATH instead of the default EXOSIMS/EXOSIMS
#   -v VERB   => set verbosity to VERB (0=quiet or 1=verbose)
#   -P PAR    => run without ipython parallel, using PAR independent jobs
#   -S SEEDS  => perform one sim for each integer seed, one per line,
#                in the file SEEDS.  This over-rides "N".
#                Exception: As a shortcut, if SEEDS is an integer, that integer 
#                is taken as the literal seed, not a filename.
#   -x SCRIPT => an extra scenario-specific script loaded on top of the argument SCRIPT
#   -e        => send email to current user at localhost on completion
#   -E ADDR   => send email to an arbitrary `addr' on completion
#   -O OPTS   => output options, a string of comma-separated tags telling
#       which EXOSIMS variables should be written, and to which files.
#       See the notes in the python driver file for how it works.
#       Default: OPTS = 'drm:pkl,spc:spc'  -- this indicates the DRM
#       is stored as a pickle with extension .pkl, and the planet parameters
#       are stored as a pickle with extension .spc.
# 
# turmon oct 2017, apr 2018
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)
TMP_ROOT=/tmp/$USER.addsims.$$ # tempfile root: matching rm at end of file

# assumed ipython parallel setup directory - see also the Makefile
# of the form "ipyparallel/aftac1-rhonda", etc.
IPYDIR=ipyparallel/${USER}/$(hostname -s)
IPYCLI=$IPYDIR/security/ipcontroller-client.json

# ipyparallel ensemble driver program
# DRIVER=EXOSIMS/run/run_ipcluster_ensemble.py
# DRIVER=Local/exo_s_ipcluster_ensemble.py
DRIVER_OLD=Local/ipcluster_ensemble_jpl_old.py
DRIVER_NEW=Local/ipcluster_ensemble_jpl_driver.py
# this is the driver we will use
DRIVER=$DRIVER_NEW

# option processing
# 0: EXOSIMS path
EXO_PATH=$(pwd)/EXOSIMS
# 1: driver options that are or were special
EMAIL_OPT=
OUT_OPT=drm:pkl,spc:spc
STAND_OPT=
SEED_OPT=
# 2: generic options to the driver
DRIVER_OPT=
# 3: "dispatcher" that provides shell-level parallelism
DISPATCH_INPUT=":" # no-op, by default
DISPATCHER=
DISPAT_OPT=
# 4: "executive" that choses python2/python3
PYTHON_EXECUTIVE_2=python
PYTHON_EXECUTIVE_3=/usr/local/anaconda3/envs/cornell/bin/python3.7
PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_2
while getopts "h23p:x:P:S:neqv:E:O:" opt; do
    case $opt in
	3)
	    # set python3 executive
	    PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_3
	    ;;
	2)
	    # set python2 executive (currently, the default)
	    PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_2
	    ;;
	n)
	    # alternate driver program - developer-only
	    DRIVER=$DRIVER_OLD
	    ;;
	p)
	    # alternate EXOSIMS path
	    EXO_PATH=$OPTARG
	    ;;
	P)
	    # no ipyparallel => pass stand-alone mode down to python level
	    DRIVER_OPT="${DRIVER_OPT} --standalone"
	    DISPAT_OPT="-P $OPTARG"
	    # set up xargs as dispatcher, unless it was already set elsewhere
	    if [ -z "$DISPATCHER" ]; then
		DISPATCHER="xargs -I {}"
	    fi
	    ;;
	S)
	    # multi-seeds -- assume standalone, use xargs as dispatcher
	    DRIVER_OPT="${DRIVER_OPT} --standalone"
	    if expr "$OPTARG" : '[0-9 ]' > /dev/null; then
		# literal seed -- echoed into xargs at invocation time
		DISPATCH_INPUT="echo $OPTARG"
		DISPATCHER="xargs -I{}"
	    else
		# file
		DISPATCHER="xargs -I{} -a $OPTARG"
	    fi
            # (removed this, -P defaults to 1 anyway, b/c did not want to make -S override -P)
	    # DISPAT_OPT="-P 8" # can over-ride with subsequent -P
	    SEED_OPT="--seed {}"
	    ;;
	x)
	    # x-tra specs - ensure it is quoted
	    DRIVER_OPT="${DRIVER_OPT} --xspecs $OPTARG"
	    ;;
	q)
	    # quiet object creation
	    DRIVER_OPT="${DRIVER_OPT} --quiet"
	    ;;
	v)
	    # verbosity 0/1
	    DRIVER_OPT="${DRIVER_OPT} --verbose $OPTARG"
	    ;;
	e)
	    # send e-mail to $USER
	    EMAIL_OPT="--email $USER"
	    ;;
	E)
	    # send e-mail to named arg
	    EMAIL_OPT="--email $OPTARG"
	    ;;
	O)
	    # output options
	    OUT_OPT="$OPTARG"
	    ;;
	h)
	    # help text
	    sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
            echo "Note: Using python at" $(which python)
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
   echo "${PROGNAME}: Error: Need exactly two arguments" >&2
   exit 1
fi

SCRIPT="$1"
NRUN=$2

# ensure EXOSIMS path
if [[ ! -r $EXO_PATH/EXOSIMS/__init__.py ]]; then
   echo "${PROGNAME}: Error: EXOSIMS path seems invalid" >&2
   exit 1
else
   echo "${PROGNAME}: Using EXOSIMS at \`${EXO_PATH}'."
fi

# handle ordinary multi-runs by sending $NRUN lines to the
# standalone dispatcher (xargs) - condition is:
# DISPATCHER is set, but SEED_OPT is not set
if [ -n "$DISPATCHER" -a -z "$SEED_OPT" ]; then
    seq ${NRUN} > ${TMP_ROOT}.seq
    DISPATCHER="${DISPATCHER} -a ${TMP_ROOT}.seq"
    NRUN=1
fi
# force NRUN=1 if SEED_OPT is set: multi-runs are handled by
# the xargs mechanism in this case
if [ -n "$SEED_OPT" ]; then
    NRUN=1
fi

# files/dirs created by any child commands will be group-writable
umask 0002

if [ ! -r $SCRIPT ]; then
    echo "${PROGNAME}: Error: Could not read the file \`${SCRIPT}'"
    exit 1
fi

# basename for sim results is derived from script name: if it's in
# the "special" Scripts/... dir, just script Scripts/ out, else
# use the basename less .json
if [[ "$SCRIPT" =~ ^Scripts/.* ]]; then
    sim_base=$(echo "$SCRIPT" | sed -e 's:Scripts/::' -e 's:\.json$::')
else
    sim_base=$(basename $SCRIPT .json)
fi

# (1) sanity check: DIR/Module/__init__.py should exist, for each
#     colon-separated DIR in PYTHONPATH, if you want to use "import Module".
# (2) this PYTHONPATH should also be used for the remote engines,
#     and the engine PYTHONPATH is set at engine-creation time using 
#     the Local/00-path.py startup file, which is copied to the
#     $IPYDIR/startup directory.
# formerly:
# export PYTHONPATH=$(pwd)/EXOSIMS:$(pwd)/Local
export PYTHONPATH=${EXO_PATH}:$(pwd)/Local

# 6/2019: sometimes we run out of threads, if many parallel instances?
export OPENBLAS_NUM_THREADS=8

# Taking this apart, the most critical elements are:
#   DISPATCH_INPUT | DISPATCHER [opts] PYTHON_EXECUTIVE DRIVER [opts] SCRIPT NRUN
# where:
#   DISPATCH_INPUT = optional input stream of newline-separated seeds, for xargs
#   DISPATCHER [opts] = xargs (which is given appropriate options)
#   PYTHON_EXECUTIVE = python2 or python3
#   DRIVER [opts] SCRIPT NRUN = python sim-runner, with opts and 2 args
$DISPATCH_INPUT | $DISPATCHER $DISPAT_OPT \
$PYTHON_EXECUTIVE $DRIVER $EMAIL_OPT $DRIVER_OPT $SEED_OPT \
  --outpath "sims/$sim_base" --outopts "$OUT_OPT" --controller $IPYCLI $SCRIPT $NRUN

# remove tempfiles, if any - double-checks that the variable is set first
[ -n "$TMP_ROOT" ] && rm -f ${TMP_ROOT}.*

