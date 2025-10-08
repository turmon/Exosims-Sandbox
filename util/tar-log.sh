#!/usr/bin/env bash
#
# tar-log.sh: archive log files for scenario(s)
#
# Make an archive of log files to save inodes on gattaca.
#
# Usage:
#   `tar-log.sh [-v] [-r] SCENARIO`
#
# (1) If SCENARIO is an ensemble directory, then typically SCENARIO/drm, 
# etc., exist within it.  The log directory
#   SCENARIO/log
# is compiled into a tarfile indexed with the current date, and the
# *log directory is removed*.
# (2) If SCENARIO is a .fam/exp directory, then the log directory:
#    SCENARIO/Batch/log
# is compiled into a tarfile as above, and the *log directory is removed*.
# Also, in this case, if -r is given, then we descend into the child
# directories (.fam/.exp/ensemble), and compress their log directories also.
# 
# Give -v to increase verbosity.
#
# Simplest usage:
# ```
#   $ util/tar-log.sh sims/Example.fam/coroOnlyScenario
# ```
#
# turmon 2025-oct
#
## [end comment block]

# exit-on-error
# (pipefail left out, causes issues when the pipe is closed by "head -1")
set -eu

PROGNAME=$(basename $0)

START_TIME=$SECONDS

# attempt to give group-write to created files
umask 002

# arguments to recursive invocation, if any
args_prog=
verbose=0
recurse=0
indent=
while getopts "hrvI:" opt; do
    case $opt in
        I)
            # indent string
            indent="$OPTARG"
            ;;
        r)
            # recursion OK?
            recurse=1
            args_prog="$args_prog -r"
            ;;
        v)
            # verbosity
            verbose=1
            args_prog="$args_prog -v"
            ;;
        h)
            # help text
            sed 's/^#//' $(which $0) | awk '/^#/{exit};NR>1{print}'
            exit 2
            ;;
        \?)
            echo "${PROGNAME}:${indent} Invalid option, exiting.  Try -h." >&2
            exit 2
            ;;
    esac
done
shift $((OPTIND-1))

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}:${indent} Error: Need exactly 1 argument, use -h for help." >&2
   exit 2
fi

# scenario (simulation) directory -- strip trailing slash if any
sim_dir="${1%/}"

####################
#
# validation for input
#
####################

# For now, insist that we're doing the tar from sims/... 
# so we can unpack from sims/ eventually
if [[ $sim_dir != sims/* ]]; then
    echo "${PROGNAME}:${indent} Error: given Scenario must be in sims/" >&2
    exit 1
fi

if [ ! -d "$sim_dir" ]; then
    echo "${PROGNAME}:${indent} Error: given Scenario must be an existing directory in sims/ " >&2
    exit 1
fi

# Divide into mutually-exclusive cases
#   logic: [-r] may be turned on even though we're just
#   in a scenario or an unknown directory, so we cannot let
#   $recurse determine what we do
if [[ $sim_dir == sims/*.exp || $sim_dir == sims/*.fam ]]; then
    # ExpFam: always compress, and then recurse if flag set
    arg_type=ExpFam
elif [ -d $sim_dir/drm ]; then
    # Scenario: always compress, and never recurse (even if recurse=1)
    arg_type=Scenario
else
    # Unknown: just exit, do not ever recurse
    arg_type=Unknown
fi


####################
#
# Utilities
#
####################

# exit with logging
early_out () {
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "${PROGNAME}:${indent} Done (${ELAPSED_TIME} seconds)"
    exit $1
}

# ensure we can create a new file here -- will error out if we cannot
release_canary() {
    # we're doing rm -f, so be careful
    if [ -z "$1" ]; then
        echo "${PROGNAME}:${indent} Fatal: Zero-length canary: should not happen." >&2
        early_out 1
    fi
    # fly, birdie
    rm -f "$1"
    if ! touch "$1"; then
        echo "${PROGNAME}:${indent} Error: Cannot create a new test file here (inode quota issue?)" >&2
        early_out 1
    fi
    rm -f "$1"
}


####################
#
# Dispatch by cases 
#   (known to be mutually exclusive)
#
####################

# a universal tag so we can repeat the operation
datenow=$(date -I)

# 1: Unknown
if [[ $arg_type == Unknown ]]; then
    if [ $verbose -eq 1 ]; then
        echo "${PROGNAME}:${indent} Skipping: Unrecognized dir ($sim_dir)" >&2
    fi
    early_out 0
fi

# 2: ExpFam
if [[ $arg_type == ExpFam ]]; then

    # 2A: tar the experiment-wide logfiles, if any
    if [ ! -d "$sim_dir/Batch/log" ]; then
        # (skip, but no early-out)
        echo "${PROGNAME}:${indent} Experiment has no overall Batch/log area. Skipping that part."
    else
        # ensure we can create a new file here -- error out if we cannot
        release_canary "$sim_dir/.tar_file_canary"

        echo "${PROGNAME}:${indent} Consolidating the Batch/log files"
        # make the logfiles match a pattern
        reps=0
        tarfile="$sim_dir/Batch/log-${datenow}_${reps}.tgz"
        while [ -r $tarfile ]; do
            # echo "${PROGNAME}:${indent} Note: Batch tarfile would over-write ($tarfile)" >&2
            reps=$(( reps + 1 ))
            tarfile="$sim_dir/Batch/log-${datenow}_${reps}.tgz"
        done

        tar czf "$tarfile" "$sim_dir/Batch/log"
        echo "${PROGNAME}:${indent} Batch: $tarfile"
        echo "${PROGNAME}:${indent} Removing the Batch/log files"
        rm -r "$sim_dir/Batch/log"
    fi

    # if no [-r], do not go downward
    if [[ $recurse == 0 ]]; then
        if [ $verbose -eq 1 ]; then
            echo "${PROGNAME}:${indent} No recursive descent into $sim_dir" >&2
        fi
        early_out 0
    fi

    # 2B: tar Exp/Fam/Scenarios below
    #     (below here, we do go downward)
    if [ $verbose -eq 1 ]; then
        echo "${PROGNAME}:${indent} Descending into: $sim_dir" >&2
    fi
    find "$sim_dir" -mindepth 1 -maxdepth 1 -type d -print0 | \
    while read -d $'\0' d; do
        # (will always be a directory)
        # (don't usually do this, it's too much)
        if [ $verbose -gt 1 ]; then
            echo "${PROGNAME}:${indent} Downward: $d"
        fi
        # the recursive call, with arguments including -v and -r
        #   add spaces to indentation
        $0 -I "$indent  " $args_prog "$d"
    done

    # 2C: Exit now
    early_out 0
fi

####################
#
# 3: Scenario

if [ ! -d "$sim_dir/log" ]; then
    if [ $verbose -eq 1 ]; then
        echo "${PROGNAME}:${indent} Scenario has no log directory. Skipping."
    fi
    # not an error
    early_out 0
fi

# just the script name part
s_name=$(basename $sim_dir)
echo "${PROGNAME}:${indent} Ensemble: $s_name"

# tar the pooled logfiles

# ensure we can create a new file here -- error out if we cannot
release_canary "$sim_dir/.tar_file_canary"

# echo "${PROGNAME}:${indent} Consolidating log files"

# make the logfiles match a pattern
reps=0
tarfile="$sim_dir/log-${datenow}_${reps}.tgz"
while [ -r $tarfile ]; do
    # echo "${PROGNAME}:${indent} Note: Batch tarfile would over-write ($tarfile)" >&2
    reps=$(( reps + 1 ))
    tarfile="$sim_dir/log-${datenow}_${reps}.tgz"
done

tar czf "$tarfile" "$sim_dir/log"
echo "${PROGNAME}:${indent} + $tarfile"
[ $verbose -gt 0 ] && echo "${PROGNAME}:${indent} Removing the log files"
rm -r $sim_dir/log

# exit OK
early_out 0

