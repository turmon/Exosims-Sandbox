#!/bin/sh
# get det/char QOI stuff out of the JS source
# hacky, but works
grep 'fieldname:.*h_.*_mean' ../www-resources/star-target-plots.js | sed -e 's/[{ ][a-z]*://g' -e 's/  */ /g' -e "s/'//g" -e 's/[{}]//g' -e 's/h_star_//' | awk -F, 'BEGIN {OFS=" | "};{print "",$2, $4, "`" $1 "`", "TBD", ""}' | sed -e 's/` /`/' -e 's/^ //' -e 's/ $//'

