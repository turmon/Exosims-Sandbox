#!/bin/sh
#
# Put lines for index at the end of a designated header file.
#
# Typical usage:
#   util/compile-plot-index.sh
# or 
#   util/compile-plot-index.sh >> index.md
#
# The convention is that files beginning with "index" are not themselves 
# put in the index!
#

echo ""
grep '^Title: ' $(ls *.md | grep -v '^index') \
  | sed -e 's/Title: //' -e 's/[.]md/.html/' \
  | awk -F : '{print "* [" $2 "](" $1 ")"; print ""}'
