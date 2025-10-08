#!/bin/bash
#
# ad hoc function to remove generated html
#

# safe mode
set -eu

PROGNAME=$(basename $0)

# enforce 1 argument
if [ $# -ne 1 ]; then
   echo "${PROGNAME}: Error: Need exactly 1 argument, use -h for help." >&2
   exit 2
fi

# experiment/fam directory
exp_dir="$1"

# need an existing sims/.../*.{fam,exp}
if [[ "$exp_dir" != sims/* ]]; then
    echo "${PROGNAME}: Error: first argument must be in sims/ " >&2
    exit 1
fi

if [ ! -d "$exp_dir" ]; then
    echo "${PROGNAME}: Error: first argument must be an existing directory in sims/ " >&2
    exit 1
fi

# for a tag so we can repeat the operation
datenow=$(date -I)

# remove the "script-by-script" HTML
# not the best way to do this loop
for d in $(find "$exp_dir" -maxdepth 2 -type d -name html); do
    if [ ! -d $d/../drm ]; then
        echo "${PROGNAME}: Something wrong ($d)"
        exit 1
    fi
    bn=$(dirname $d)
    # just the script name part
    s_name=$(basename $bn)
    echo "${PROGNAME}:     Processing $s_name"
    # just the html
    rm -r $d
done

echo "${PROGNAME}: Done."
