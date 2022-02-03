#!/bin/sh -f
# 
# generate some graphics capturing makefile dependencies
#
# these are used in the brief paper I wrote summarizing the 
# use of Make in the Sandbox
#
# Usage: Just run from the Sandbox, outputs placed in in ~/Media
#
# 2022-02-01 turmon

# plain script, standard targets
util/dev/visualize-makefile-targets.sh ~/Media/make-dep/basic
# plain script, other targets
util/dev/visualize-makefile-targets.sh -T -t keepout-5 -t path-movie-5 ~/Media/make-dep/basic

# family, standard targets
util/dev/visualize-makefile-targets.sh -S FAM ~/Media/make-dep/family

# Experiment, looping over several targets
util/dev/visualize-makefile-targets.sh -S EXP -T -t exp-reduce -t exp-html -t exp-path-ensemble -t exp-path-movie-5 -t exp-graphics ~/Media/make-dep/experiment
