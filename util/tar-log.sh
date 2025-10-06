#!/usr/bin/env bash
#
# tar-log.sh: archive log files for a scenario
#
# Make an archive of log files to save inodes on gattaca.
#
# Usage:
#   `tar-log.sh SCENARIO`
#
# where SCENARIO is a simulation directory (typically SCENARIO/drm, etc.,
# exist within it).  The directory:
#   SCENARIO/log
# is compiled into a tarfile indexed with the current date, and the
# *directory is removed*.
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

# attempt to give group-write to created files
umask 002

verbose=0
while getopts "hv" opt; do
    case $opt in
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
# validation on input
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

if [ ! -d "$sim_dir/log" ]; then
    echo "${PROGNAME}: Note: Scenario has no log directory. Skipping."
    # not an error
    exit 0
fi

# just the script name part
s_name=$(basename $sim_dir)
echo "${PROGNAME}: Processing $s_name"

# make a tag, so we can repeat the operation later
datenow=$(date -I)

# tar the pooled logfiles

# ensure we can create a new file here -- will error out if we cannot
canary="$sim_dir/.tar_file_canary"
rm -f "$canary"
if ! touch "$canary"; then
    echo "${PROGNAME}: Error: Cannot create a new test file here (inode quota issue?)" >&2
    exit 1
fi

rm -f "$canary"

# echo "${PROGNAME}: Consolidating log files"

# make the logfiles match a pattern
reps=0
tarfile="$sim_dir/log-${datenow}_${reps}.tgz"
while [ -r $tarfile ]; do
    # echo "${PROGNAME}: Note: Batch tarfile would over-write ($tarfile)" >&2
    tarfile="$sim_dir/log-${datenow}_${reps}.tgz"
    reps=$(( reps + 1 ))
done

tar czf $tarfile $sim_dir/log
echo "${PROGNAME}: $tarfile"
[ $verbose -gt 0 ] && echo "${PROGNAME}: Removing the log files"
rm -r $sim_dir/log

echo "${PROGNAME}: Done."
