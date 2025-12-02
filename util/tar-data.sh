#!/usr/bin/env bash
#
# tar-data.sh: archive data files for scenario(s)
#
# Make an archive of data files to save inodes on gattaca.
#
# Usage:
#   `tar-data.sh [-v] [-d] [-l] SCENARIO`
#
# (1) If SCENARIO is an ensemble directory, then typically SCENARIO/drm, 
# etc., exist within it.  The data directories
#   SCENARIO/{drm,spc}
# are compiled into a tarfile indexed with the current date, and the
# *drm/spc directories are removed*.
# (2) If SCENARIO is a .fam/exp directory, we descend into the child
# directories (but just down *one level* into ensembles there), and
# tar-and-delete their data directories also.
# 
# Options:
#   -d: dry run. No tar-files created, no data deleted. Informational scenario counts are given.
#   -t THRESH: only archive data for scenarios with yield <= THRESH. Implies -l.
#   -l: only archive data for low-yield scenarios; shortcut for -t 4.
#   -v to increase verbosity.
#
# Simplest usage:
# ```
#   $ util/tar-data.sh sims/Example.fam/coroOnlyScenario
# ```
#
# turmon 2025-dec
# note:
#   - There may be interactions with other components that use existence of
#     SCENARIO/drm to validate that it's a real Exosims scenario.
#   - Best with bash 4.2+
#
## [end comment block]

# exit-on-error
# (pipefail left out, causes issues when the pipe is closed by "head -1")
set -eu

PROGNAME=$(basename $0)

# run the while loop below in the current shell (bash 4.2+)
if [[ $BASH_VERSION == 3.* || $BASH_VERSION == 4.0* || $BASH_VERSION == 4.1* ]]; then
    echo "${PROGNAME}: Old bash, works, but counts will be incorrect."
else
    shopt -s lastpipe
fi

START_TIME=$SECONDS

# attempt to give group-write to created files
umask 002

verbose=0
do_low_only=0
low_thresh=4
dryrun=0
while getopts "dlt:hv" opt; do
    case $opt in
        d)
            # dryrun -- no tar, no rm
            dryrun=1
            ;;
        l)
            # tar only some low-yield scenarios
            do_low_only=1
            ;;
        t)
            # threshold for earth-chars of what "low" means
            low_thresh=$OPTARG
	    do_low_only=1
            ;;
        v)
            # verbosity
            verbose=1
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
   echo "${PROGNAME}: Error: Need exactly 1 argument, use -h for help." >&2
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
    echo "${PROGNAME}: Error: given Scenario must be in sims/" >&2
    exit 1
fi

if [ ! -d "$sim_dir" ]; then
    echo "${PROGNAME}: Error: given Scenario must be an existing directory in sims/ " >&2
    exit 1
fi


####################
#
# Utilities
#
####################

# exit with logging
early_out () {
    ELAPSED_TIME=$(($SECONDS - $START_TIME))
    echo "${PROGNAME}: Done (${ELAPSED_TIME} seconds)"
    exit $1
}

# ensure we can create a new file here -- will error out if we cannot
release_canary() {
    # we're doing rm -f, so be careful
    if [ -z "$1" ]; then
        echo "${PROGNAME}: Fatal: Zero-length canary: should not happen." >&2
        early_out 1
    fi
    # fly, birdie
    rm -f "$1"
    if ! touch "$1"; then
        echo "${PROGNAME}: Error: Cannot create a new test file here (inode quota issue?)" >&2
        early_out 1
    fi
    rm -f "$1"
}


####################
#
# Utility: make data tarfile for one scenario
#
####################

tar_scenario() {

    local sim_dir="$1"
    
    # find data directories, if any
    local -a datas
    datas=()

    if [ -d "$sim_dir/drm" ]; then
	datas+=("$sim_dir/drm")
    fi
    if [ -d "$sim_dir/spc" ]; then
	datas+=("$sim_dir/spc")
    fi

    if [[ ${#datas[@]} -eq 0 ]]; then
	if [ $verbose -eq 1 ]; then
            echo "${PROGNAME}: Scenario has no data directories. Skipping."
	fi
	# not an error
	return 0
    fi

    # just the script name part
    local s_name=$(basename $sim_dir)
    if [ $verbose -eq 1 ]; then
	echo "${PROGNAME}: Ensemble: $s_name"
    fi

    # tar the data directories

    # ensure we can create a new file here -- will error out if we cannot
    release_canary "$sim_dir/.data_file_canary"

    # echo "${PROGNAME}: Consolidating data files"

    # make the logfiles match a pattern
    reps=0
    tarfile="$sim_dir/data-${datenow}_${reps}.tgz"
    while [ -r $tarfile ]; do
	# echo "${PROGNAME}: Note: Batch tarfile would over-write ($tarfile)" >&2
	reps=$(( reps + 1 ))
	tarfile="$sim_dir/data-${datenow}_${reps}.tgz"
    done

    if [ $dryrun -eq 1 ]; then
	# no-op
	precmd=":"
    else
	precmd=
    fi

    $precmd tar czf "$tarfile" "${datas[@]}" 
    echo "${PROGNAME}: + $tarfile"
    [ $verbose -gt 0 ] && [ $dryrun -eq 0 ]] && echo "${PROGNAME}: Removing data directories"
    $precmd rm -r "${datas[@]}"

    # global tarfile count
    n_tarfile=$((n_tarfile + 1))

    return 0
}

####################
#
# Main program
#
####################

# a universal tag so we can repeat the operation
datenow=$(date -I)

# examine the given Scenario
if [[ $sim_dir == sims/*.exp || $sim_dir == sims/*.fam ]]; then
    # .exp/.fam
    arg_type=ExpFam
elif [ -d $sim_dir/drm -o -d $sim_dir/spc ]; then
    # plain Scenario
    arg_type=Scenario
else
    arg_type=Unknown
fi

if [ $verbose -eq 1 ]; then
    echo "${PROGNAME}: Checking into: $sim_dir" >&2
fi

if [ $do_low_only -eq 1 ]; then
    # TODO: allow setting the threshold manually
    args="-t $low_thresh"
else
    args=""
fi

# number of scenarios found
n_scenario=0
# number of new tarfiles created
n_tarfile=0

# always-used arguments:
#   -S => sandbox conventions for CSV-file naming
#   -0 => output line termination with "\0"
#   -n T => output all matching lines
#   bottom => sort from low-to-high (although order does not matter)
#   -o experiment => output the experiment (scenario) name
#   -k chars_earth_unique => sort/threshold on earth chars
util/select_ensembles.py $args -S -0 -n T -o experiment -k chars_earth_unique bottom "$sim_dir" | \
    while read -d $'\0' d; do
	# skip nested dirs
	# (could make this into a recursive call, but no need right now)
	if [[ $d == *.exp || $d == *.fam ]]; then
            if [ $verbose -gt 0 ]; then
		echo "${PROGNAME}: Skipping nested dir: $d"
            fi
	    continue;
	fi
	# number of matching scenarios considered
	n_scenario=$((n_scenario + 1))
        # (will always be a directory)
        # (don't usually do this, it's too much)
        if [ $verbose -gt 1 ]; then
            echo "${PROGNAME}: Downward: $d"
        fi
	if [[ $arg_type == ExpFam ]]; then
	    # .exp/.fam: need to add the scenario name
	    target_dir="$sim_dir/$d"
	else
	    # otherwise, the name is complete as given
	    target_dir="$sim_dir"
	fi
        # pack the given directory
        tar_scenario "$target_dir"
    done

echo "${PROGNAME}: Scenarios examined: $n_scenario"
echo "${PROGNAME}: New tarfiles created: $n_tarfile"
if [ $dryrun -eq 1 ]; then
    echo "${PROGNAME}: NOTE: dry run: No actions taken"
fi
    

early_out 0
# (above line exits)


