#!/usr/bin/env bash
# 
# rad-sma-rectangle-plot-driver.sh: sh driver for radius/luminosity box-format plots
#
# usage:
#   rad-sma-rectangle-plot-driver.sh [-d] [-o OP] in_tmpl
# where:
#
# +  in_tmpl: filename template for data input
#       e.g.: sims/ExoS_A_SAG13/reduce-%s.%s
# +  dest_tmpl: the destination template for graphical output
#       e.g.: sims/ExoS_A_SAG13/gfx/plot-%s.%s
# +  -q indicates to plot quantiles, otherwise, means are plotted.
#
# For example,
#
# +  `rad-sma-rectangle-plot-driver.sh self`
#
## [end comment block]

# Michael Turmon, JPL, Aug 2019

# exit-on-error
set -euo pipefail

PROGNAME=`basename $0`
USAGE="Usage: $PROGNAME [-q] [-o output_tmpl] in_tmpl"

# program that we're the driver of
PLOT_PROG=util/rad-sma-rectangle-bin-plot.py

# default field list: csv data
csv_fields=population,char_full,char_strict,char_tput_full,char_tput_strict
# default field list: 'self' data
self_fields=kopparapu,sag13,sdet

# get options
prog_arg=
output_opt=
fields_opt=
while getopts "qhf:o:" o; do
    case "$o" in
	q)
	    prog_arg="${prog_arg} --quantile";;
	f)
	    fields_opt="$OPTARG";;
	o)
	    output_opt="$OPTARG";;
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
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 1 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi

echo "${PROGNAME}: Beginning."

# ok to get the arg
in_tmpl="$1"

# infer the output location if it was not given
if [ -z "$output_opt" ]; then
    if [ $in_tmpl != self ]; then
	# usual case for scripted operation
	output_opt=$(dirname $in_tmpl)/gfx/det-rad-sma-%s.%s
    else
	# make default for builtin data from "self"
	output_opt=./rad-sma-rect-%s.%s
    fi
fi

# fill in default fields if not given
if [ -z "$fields_opt" ]; then
    if [ $in_tmpl != self ]; then
	fields_opt="$csv_fields"
    else
	fields_opt="$self_fields"
    fi
fi

# make created files group-writable
umask 0002

# try to create the destination dir
mkdir -p "$(dirname $output_opt)"
# attempt to make it group-writable (and ignore errors)
chmod g+sw "$(dirname $output_opt)" 2> /dev/null || true

# Run the command -- something like:
# util/radlum-rectangle-bin-plot.py --quantile -o /tmp/mjt/det-%s.%s ...
#    -f population,char_full,char_strict,char_tput_full,char_tput_strict ...
#    sims/HabEx_4m_TSDD_DulznoDD_TF14_20190707/reduce-%s.csv

$PLOT_PROG --eta --percent $prog_arg -o "$output_opt" -f "$fields_opt" "$in_tmpl"
