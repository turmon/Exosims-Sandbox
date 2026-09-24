#!/usr/bin/env bash
#
# select_ensembles_exp.sh: record ensemble selection orderings as a Makefile fragment
#
# For an Experiment or Family (`.exp`/`.fam`) directory, write a Makefile
# fragment giving the full ordering of the ensembles within it, for each
# of the two selection modes used by the `exp-*` targets in the Makefile:
#
# * `EXP_ORDER_top`: by yield, highest first (key `chars_earth_unique`)
# * `EXP_ORDER_mix`: arbitrary but stable (MD5 hash of key `experiment`)
#
# Operational context/motivation: Selecting the first N elements of either list is then the same
# as `select_ensembles.py -n N`.
# The Makefile includes this fragment to generate dependencies for the 
# `exp-*` targets. The Makefile knows to remake it (with this script)
# after the experiment is reduced (altering the reduction summary file).
#
# ## Usage:
# ```
#   select_ensembles_exp.sh [-o OUTFILE] EXPDIR
# ```
#
# where EXPDIR is the experiment directory, e.g., `sims/Example.exp`.
# The selection is made from `EXPDIR/reduce-yield-plus.csv`.  If that is
# not present (e.g., nothing has been reduced), the lists are empty.
#
# ### Options:
#
# * `-o OUTFILE` gives the output file. The default is `EXPDIR/exp-select.mk`.
#
# Typical usage (the Makefile does this):
# ```
#   $ util/select_ensembles_exp.sh sims/Example.exp
# ```
#
## [end comment block]

# exit-on-error
set -euo pipefail

PROGNAME=$(basename $0)

# the command line, for the generated-file comment
cmd_line="util/$PROGNAME $*"

# attempt to give group-write to created files
umask 002

# selector program, and its common options
#   -q => quiet if not reduced already
#   -o experiment => output this key (ens. name)
SELECT_PROG="util/select_ensembles.py -q -o experiment -n T"

outfile=
while getopts "ho:" opt; do
    case $opt in
	o)
	    # output file
	    outfile="$OPTARG"
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

# experiment directory, without trailing /
expdir="${1%/}"

if [ ! -d "$expdir" ]; then
    echo "${PROGNAME}: Given experiment directory \`$expdir' not readable, exiting." >&2
    exit 1
fi

# default output file
if [ -z "$outfile" ]; then
    outfile="$expdir/exp-select.mk"
fi

csv="$expdir/reduce-yield-plus.csv"
tmpfile="$outfile.tmp"

rm -f "$tmpfile"
{
    echo "# Generated file -- do not edit."
    echo "# Made by: $cmd_line"
    echo "# On: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "# Full ordering of the ensembles within $expdir, from:"
    echo "#   $csv"
    $SELECT_PROG -k chars_earth_unique -M EXP_ORDER_top top "$csv"
    $SELECT_PROG -k experiment         -M EXP_ORDER_mix mix "$csv"
} > "$tmpfile"

chmod g+w "$tmpfile"
mv -f "$tmpfile" "$outfile"
