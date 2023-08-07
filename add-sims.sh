#!/usr/bin/env bash
#
# add-sims: Perform EXOSIMS simulations, accumulating runs into subdirectories
#
# This is a wrapper around the python driver, ipcluster_ensemble_jpl.py.
#
# Usage:
#   add-sims.sh [-3] [-q] [-A | -S | -Z] [-j JOBS] [-v VERB] [-x SCRIPT] [-e] [-E ADDR] [-O OPTS] SCRIPT SEEDS
#
# Uses the JSON script SCRIPT and performs a series of parallel runs given by SEEDS.
# * The runs are done on the local computer, from a pool that depends on machine
# (aftac1 = 28, aftac2/3 = 44).
# * If just one SEED is given, Exosims output is send directly to the screen, which
# is useful for checking validity.  If more than one SEED, the Exosims output is sent
# to a set of log files (sims/SCRIPT/log_sim/1/*), and a summary of current job status
# is sent to the screen.
# * To repeat an error-causing run interactively, paste the command-line printed by 
# this script (the one that calls ipcluster_ensemble_jpl.py) into your terminal,
# and add "--interactive" to the arguments.
#
# Typical Usage:
#   (a) Run 100 jobs on mustang* and aftac*, across 100 deterministic seeds
#     add-sims.sh -Z Scripts/ExampleScript.json Experiments/seed100.txt 
#   (b) Run a machine-dependent number of jobs spanning 64 seeds on the local machine
#     add-sims.sh Scripts/ExampleScript.json Experiments/seed64.txt 
#   (c) Do a test run for a single arbitrary seed (we use 777):
#     add-sims.sh Scripts/ExampleScript.json =777
#   (d) Write cache files for a new script, but do not run a sim:
#     add-sims.sh Scripts/ExampleScript.json =0
#
# Arguments:
#   SCRIPT -- a JSON script suitable for EXOSIMS
#   SEEDS  -- either:
#             (1) if a plain integer, that number of randomly-chosen initial seeds.
#             (2) if of the form =SEED, where SEED is an integer, that single integer 
#                 is the seed. If SEED is 0, only cache warming is done (no run_sim).
#             (3) otherwise, it is a filename giving a list of integer seeds, one per line.
#
# Options:
#   -h        => show this help message and exit.
#   -A        => run using Afta's (aftac1/2/3 -- 18 + 23 + 23 jobs = 64 jobs)
#   -S        => run using Speedy's (mustang2/3/4 -- 24 + 24 + 16 jobs = 64 jobs)
#   -Z        => run using all (mustang2/3/4 + aftac1/2/3 -- total of 100 jobs)
#   -3        => run the EXOSIMS sim-runner with python3; by default, "python" is used.
#   -p PATH   => EXOSIMS path is PATH instead of the default EXOSIMS/EXOSIMS
#                PATH=@ abbreviates the command-line default (what you get from
#                "python -c import EXOSIMS"), which may not be the same as
#                EXOSIMS/EXOSIMS (i.e., Sandbox EXOSIMS). Useful for venv's.
#
# Less-used Options:
#   -j JOBS   => runs only JOBS parallel jobs (not used with -A)
#   -v VERB   => set verbosity to VERB (0=quiet or 1=verbose)
#   -q        => quiet object creation
#   -x SCRIPT => an extra scenario-specific script loaded on top of the argument SCRIPT
#   -c        => chatty console output to aid debugging
#   -e        => send email to current user at localhost on completion
#   -E ADDR   => send email to an arbitrary `addr' on completion
#   -O OPTS   => output options, a string of comma-separated tags telling
#       which EXOSIMS variables should be written, and to which files.
#       See the notes in the python driver file for how it works.
#       Default: OPTS = 'drm:pkl,spc:spc'  -- this indicates the DRM
#       is stored as a pickle with extension .pkl, and the planet parameters
#       are stored as a pickle with extension .spc.
# 
# turmon oct 2017, apr 2018, may 2020
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# Assumed ipython parallel setup directory - see also the Makefile
# of the form "ipyparallel/aftac1-rhonda", etc.
# This option is no longer used.
IPYDIR=ipyparallel/${USER}/$(hostname -s)
IPYCLI=$IPYDIR/security/ipcontroller-client.json

# ensemble driver program
# see also EXOSIMS/run/run_ipcluster_ensemble.py
DRIVER=Local/ipcluster_ensemble_jpl_driver.py

# option processing
# 0: EXOSIMS path
EXO_PATH=$(pwd)/EXOSIMS
EXO_PATH_SET=no
# 1: driver options that are or were special
EMAIL_OPT=
OUT_OPT=drm:pkl,spc:spc
STAND_OPT=
SEED_OPT=
# 2: generic options to the driver
DRIVER_OPT=
# 3: "dispatcher" that provides shell-level parallelism
DISPATCHER=parallel
# by default, use 4 less than the number of cores
DISPATCHER_JOBS="--jobs -4"
BATCH=0
CHATTY=0
ECHO_CMD=0
ALL=0
# 4: "executive" that choses python2/python3
PYTHON_EXECUTIVE_2=python
PYTHON_EXECUTIVE_3=/usr/local/anaconda3/envs/cornell/bin/python3.7
PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_2
while getopts "ASZDh23bcj:p:x:qv:eE:O:" opt; do
    case $opt in
	A)
	    # use remote machines
	    ALL=1
	    ;;
	S)
	    # use only "speedy" remote machines
	    ALL=2
	    ;;
	Z)
	    # use even more remote machines
	    ALL=3
	    ;;
	D)
	    # fewer jobs on all remote machines (for Dakota)
	    ALL=4
	    ;;
	3)
	    # set python3 executive
	    PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_3
	    ;;
	2)
	    # set python2 executive (currently, the default)
	    PYTHON_EXECUTIVE=$PYTHON_EXECUTIVE_2
	    ;;
	b)
	    # batch mode output to log-files
	    BATCH=1
	    ;;
	c)
	    # chatty output to console
	    CHATTY=1
	    ;;
	j)
	    # explicit job number
	    DISPATCHER_JOBS="--jobs $OPTARG"
	    ;;
	p)
	    # alternate EXOSIMS path
	    EXO_PATH=$OPTARG
            EXO_PATH_SET=yes
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
            echo "Note: By default, using python at" $(which python)
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
SEEDS=$2

# files/dirs created by any child commands will be group-writable
umask 0002

# handle the 3 cases for the SEEDS argument
if expr "$SEEDS" : '^=[0-9][0-9]*$' > /dev/null; then
    # literal seed -- echoed into xargs at invocation time
    SEED_USED=$(echo $SEEDS | sed 's/=//')
    SEED_INPUT="echo $SEED_USED"
    NRUN=1
elif expr "$SEEDS" : '^[0-9][0-9]*$' > /dev/null; then
    # it is a number-of-runs
    # draw $SEEDS 9-digit random integers from awk
    # SEED_INPUT="awk 'BEGIN {srand(); for (i=0; i<$SEEDS; i++) print int(rand()*1000000000+1)}'"
    SEED_INPUT="seq $SEEDS | awk 'BEGIN {srand()}; {print int(rand()*999999999+1)}'"
    NRUN=$SEEDS
else
    if [ ! -r "$SEEDS" ]; then
	echo "${PROGNAME}: Error: SEEDS argument looks like a file, but is not readable." >&2
	exit 1
    fi
    # seeds are in file
    SEED_INPUT="cat $SEEDS"
    NRUN=$(wc -l < $SEEDS)
fi

echo "${PROGNAME}: Performing $NRUN EXOSIMS run(s)."

# collapse NRUN, BATCH, CHATTY into just BATCH
if [ $NRUN -gt 1 ]; then
    BATCH=1
fi
# allow over-ride
if [ $CHATTY == 1 ]; then
    BATCH=0
fi

# ensure EXOSIMS path
if [ $EXO_PATH = @ ]; then
    # useful when using venv's: this specifies the venv's EXOSIMS as the one to use
    EXO_PATH=$(python -c 'import EXOSIMS; import os; print(os.path.dirname(os.path.dirname(EXOSIMS.__file__)))')
    echo "${PROGNAME}: EXOSIMS path set to \`${EXO_PATH}'."
fi
# first construct must work if $VIRTUAL_ENV is undefined/unset
if [ "${VIRTUAL_ENV+defined}" = defined -a "$EXO_PATH_SET" = no ]; then
    echo "${PROGNAME}: Warning: In a VENV, but still using default path \`${EXO_PATH}'."
fi
if [ ! -r $EXO_PATH/EXOSIMS/__init__.py ]; then
   echo "${PROGNAME}: Error: EXOSIMS path seems invalid" >&2
   exit 1
else
   echo "${PROGNAME}: Using EXOSIMS found at \`${EXO_PATH}'."
fi

if [ ! -r $SCRIPT ]; then
    echo "${PROGNAME}: Error: Could not read the file \`${SCRIPT}'"
    exit 1
fi

# basename for sim results is derived from script name: if it's in
# the "special" Scripts/... dir, just strip Scripts/ and .json out, 
# else use the basename less .json
if [[ "$SCRIPT" =~ ^Scripts/.* ]]; then
    sim_base=$(echo "$SCRIPT" | sed -e 's:Scripts/::' -e 's:\.json$::')
else
    sim_base=$(basename $SCRIPT .json)
fi
sim_base=sims/$sim_base
logdir=$sim_base/log_sim
mkdir -p $logdir

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

# options for gnu parallel - rather complex
#   cores reported: aftac1, 32; aftac2, 48; aftac3, 48.
#   need PATH exported because otherwise ssh will just get a vanilla PATH
#   --files ==> put chatter on stdout/stderr into log files instead of showing on stdout
#   --progress ==> show a progress indication
# (parallel does not object if the --env FOO is not set)
PAR_ENV_OPTS="--env PATH --env PYTHONPATH --env OPENBLAS_NUM_THREADS --env EXOSIMS_PARAMS"
if [ $BATCH == 1 ]; then
    PAR_LOG_OPTS="--joblog $logdir/Joblog --results $logdir --files --progress"
else
    # --line-buffer => allow lines of jobs to intermix
    PAR_LOG_OPTS="--joblog $logdir/Joblog --line-buffer"
fi
#   --sshdelay ==> required because sshd can't accept connections super-fast
#   -S ==> remote hosts, N/aftac1 means at most 2 jobs to aftac1, etc.
if [ $ALL == 1 ]; then
    PAR_SSH_OPTS="--workdir . --sshdelay 0.06 -S 18/aftac1,23/aftac2,23/aftac3"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 2 ]; then
    PAR_SSH_OPTS="--workdir . --sshdelay 0.10 -S 24/mustang2,24/mustang3,16/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 3 ]; then
    PAR_SSH_OPTS="--workdir . --sshdelay 0.12 -S 14/aftac1,11/aftac2,11/aftac3,26/mustang2,24/mustang3,18/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
elif [ $ALL == 4 ]; then
    PAR_SSH_OPTS="--workdir . --sshdelay 0.10 -S 6/aftac1,6/aftac2,6/aftac3,12/mustang2,12/mustang3,8/mustang4"
    # remove the --jobs option, which interacts with the above
    DISPATCHER_JOBS=""
else
    # local host only
    PAR_SSH_OPTS="-S :"
fi

# this can enable debugging a failing invocation
if [ $CHATTY == 1 ]; then
    set -x
fi
# Taking this apart, the most critical elements are:
#   DISPATCH_INPUT | DISPATCHER [opts] PYTHON_EXECUTIVE DRIVER [opts] SCRIPT NRUN
# where:
#   DISPATCH_INPUT = optional input stream of newline-separated seeds, for xargs
#   DISPATCHER [opts] = xargs (which is given appropriate options)
#   PYTHON_EXECUTIVE = python2 or python3
#   DRIVER [opts] SCRIPT NRUN = python sim-runner, with opts and 2 args
eval $SEED_INPUT | $DISPATCHER $DISPATCHER_JOBS $PAR_ENV_OPTS $PAR_LOG_OPTS $PAR_SSH_OPTS \
$PYTHON_EXECUTIVE $DRIVER $EMAIL_OPT $DRIVER_OPT --standalone --seed {} \
 --outpath "$sim_base" --outopts "$OUT_OPT" --controller $IPYCLI "$SCRIPT" 1

