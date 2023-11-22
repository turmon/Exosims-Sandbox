#!/bin/sh
#
# stdin -> awk -> stdout filter to extract metadata to m4 format
#
# Typical usage:
#   util/meta-extract.sh < obs-timeline.md
#
# Extracts an initial block of metadata from the input, 
# assumed to be of the form:
#
# Title: A. U. Thor
# Date: 2023 B.C.
# Subject: Rocks
#
# that is, <tag>:<whitespace><value>
# terminated by a blank line.
#
# We extract it into a series of macro-definitions like:
# define({{Title}},A. U. Thor)dnl
# These definitions can be directly read by m4.

awk -F: 'BEGIN {print "dnl Contents: Metadata in M4 format"}; /^[A-Za-z]*:/ {gsub(/^[ \t]+/, "", $2); print "define({{META_" $1 "}}," $2 ")dnl"}; /^[ \t]*$/ {exit 0}'

