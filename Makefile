# Makefile for EXOSIMS data reduction and plot generation.
#
# Usage:
#    $ make [ -B ] S=SCENARIO TARGET
#    $ make TARGET
# where TARGET is one of the below.
#
# * make -B forces the action, e.g., when the plot code changes.
# * The S=SCENARIO argument is required for data reduction, plots, and html
#   SCENARIO is the base name of the relevant script file, e.g., for
#   Scripts/HabEx_4m_TS_20180201.json, use S=HabEx_4m_TS_20180201
# * Alternately: set S to the json script, or the sim directory, so that
#   shell filename completion can fill in the name.
# * Simulations are added by the `add-sim.sh' command, separate from "make".
#
# Targets:
# (1) Data reduction and plotting
#   All these targets require a scenario name.
#   reduce:        reduce DRMs to tabulated CSV files for later plotting
#   graphics:      make detections-vs-time plots, and radius-luminosity bar plots.
#   html:          re-generate the index.html that summarizes the given scenario
#   html-only:     same as html, but do not re-reduce the data or remake graphics.
#   path-ensemble: make lon/lat plots of slews taken by an ensemble.
#   star-visits:   per-star tabulation of successful detection visits
#   path-movie-N:  make "N" path-movies and final frames 
#   path-final-N:  make "N" final frames, only
#   obs-timeline-N:make "N" observing-target timelines
#   keepout-N:     make "N" keepout-vs-time plots
#     (above targets make graphics for N arbitrary sims from the scenario, 
#      N = 1, 2, 5, 10, 20, 50, 100, or T, where T=all)
#   tar-log:       replace scenario/log with its "tar" archive to save space
#   status:        list the current contents of DRMs for this scenario (like "ls")
# (2) Multi-script reduction and plotting
#   Note, all these targets require an *experiment* name.
#   exp-reduce:      makes "reduce" for all ensembles within the experiment
#   exp-html-top-N:  makes html (inc. graphics) for the N top (by yield) ensembles
#   exp-html-mix-N:  makes html (inc. graphics) for N selected-arbitrarily ensembles
#   exp-html:        makes html for 10 top + 20 selected ensembles - can use make -j2
#   exp-path-ensemble* \   Same pattern as html above with -mix or -top, and a
#   exp-graphics*       \  number saying how many.  Also, can leave off -top-N
#   exp-html-only*      /  and just make 10 top + 20 selected.
#   exp-movie-M-*      /   Make M movies in each of (top/mix)-N ensembles.
#   exp-obs-timeline-M-*   Make M obs-timelines in each of (top/mix)-N ensembles.
# (3) Web-server
#   html-ensure: start Apache httpd web-server, if not running already
#   html-start: start Apache httpd web-server
#   html-stop: stop Apache httpd web-server
#   html-status: show running web-servers, if any
#
## turmon oct 2017, mar 2018, feb 2022

# set the default shell
SHELL:=/bin/sh

# clear builtin suffix rules (for .c, etc.)
.SUFFIXES:

# clear builtin pattern rules to get files out of source control
%: %,v
%: RCS/%,v
%: RCS/%
%: s.%
%: SCCS/s.%

## # debugging: echo rule chains
# $(warning Debug mechanism on)
# OLD_SHELL := $(SHELL)
# SHELL = $(warning Building $@$(if $<, (from $<))$(if $?, ($? newer)))$(OLD_SHELL)

## Normalize the S input from a (possible) file/dir name -> script name
# 1: ensure S is simply-expanded variable, so it can be redefined
ifdef S
 # S_COPY serves as an indicator that S was supplied
 S_COPY:=$(S)
 override undefine S
 S:=$(S_COPY)
else
 S:=ScriptNotDefined
endif
# 2: remove some pathname components, if present
#    recall: the Scenario contains neither sims/... nor Scripts/...
# 2a: pull off Scripts/ if given as a file/dir (S=Scripts/foo.exp -> S=foo.exp)
S := $(patsubst Scripts/%,%,$(S)) 
# 2b: get script basename if given as a JSON (S=foo.json -> S=foo)
S := $(patsubst %.json,%,$(S)) 
# 2c: pull off leading sims/, if present (S=sims/foo -> S=foo)
S := $(patsubst sims/%,%,$(S)) 
# 2d: strip trailing / which could be present (S=sims/foo/ -> S=foo)
S := $(patsubst %/,%,$(S))
# 2e: strip added space at the end of S
S := $(strip $(S))
# 3: repeat script value back, if supplied
ifdef S_COPY
 $(info Make: Scenario name: "$(S)")
endif

# needed to escape equal signs in some Experiments, alas
EQUAL:= =

# hostname
HOST_NAME=$(shell hostname -s)

# program used for data reduction
REDUCE_PROG=util/reduce_drms.py
# program used for reduction of multi-ensemble experiments
# (-E: expand the given "parent" dir to subdirs with Ensembles)
REDUCE_ENS_PROG=util/reduce_drm_sets.py -E
# programs used for matlab/python graphics
GRAPHICS_PROG=util/plot_drms.sh
GRAPHALT_PROG=util/plot_drms_driver.py
GRAPHYCS_PROG=util/rad-sma-rectangle-plot-driver.sh -q
# generates tables
TABLES_PROG=util/tabulate_csv.py -q
# make tarfile of data directories
TAR_DATA_PROG=util/tar-data.sh
# make tarfile of log directory (recursive)
TAR_LOG_PROG=util/tar-log.sh -r
# generates detection visits report
STAR_VISIT_PROG=util/star_visit_pmf_tabulate.py
# select a given number of runs from a single Scenario (e.g., for timelines)
SELECT_RUN_PROG=util/select_runs.py -q
# Programs analogous to the path-movie maker...
#   driver script for making (per-drm) observation timeline plots
#   (these run rather quickly)
TIMELINE_PROG=util/plot-obs-timelines.sh
#   program for making (per-drm) keepout + observation timelines
KEEPOUT_PROG=util/plot-keepout-and-obs.py
# drm path-movie maker, generating command lines like:
#   drm-to-movie.sh -c -l 0 sims/HabEx_4m_TS_dmag26p0_20180206/drm/297992296.pkl
#   -c specifies to make the cumulative plots as well (which does not
#   make sense for the FINAL mode).
PATH_PROG=util/drm-to-movie.sh -c -l 0
PATH_PROG_FINAL=util/drm-to-movie.sh -F
# path-ensemble-plot driver script
PATH_ENS_PROG=util/ens-path-summary.sh -a
# program used for html summary
#   -i: to regenerate the global index.html as well as that for $(S)
HTML_PROG=util/html-summary.py -i
# program to select a given number of ensembles within an experiment
#  also needs a sort key (-k) argument before use
SELECT_PROG=util/select_ensembles.py -q -o experiment
# Ensemble counts to select ensembles within an experiment
EXP_COUNTS:=1 2 5 10 20 50 100 T

# default target prints help text
default:
	@ awk '/^##/{exit};NR>0{print}' $(MAKEFILE_LIST) | sed 's/^#//' 

# External variable, passed in via command line
# S = "scenario" name = script file basename without extension
# check that the script file corresponding to $S exists
script-exists:
	@ [ -r Scripts/$(S).json -o -d Scripts/$(S) ] || \
		(echo "Require a script file \`Scripts/$(S).json' or experiment directory \`Scripts/$(S)'" && exit 1)
experiment-exists:
	@ [ -d Scripts/$(S) ] || \
		(echo "Require an experiment/family directory \`Scripts/$(S)'" && exit 1)

.PHONY: script-exists experiment-exists default

# do not try to remake Makefile
Makefile:;

########################################
## Simulation status and misc

# list drms that have been made
status: script-exists
	util/drm-ls.py -l sims/$(S)/drm

.PHONY: exp-preflight
exp-preflight:
	util/exp-preflight.sh Scripts/$(S)

# compress logfiles - don't require script-exists
.PHONY: tar-log
tar-log:
	$(TAR_LOG_PROG) sims/$(S)

# compress data-files - don't require script-exists
.PHONY: tar-data
tar-data:
	$(TAR_DATA_PROG) sims/$(S)

# compress low-yield (-l) data-files - don't require script-exists
.PHONY: tar-some-data
tar-some-data:
	$(TAR_DATA_PROG) -l sims/$(S)

########################################
## Data reductions
##
.PHONY: reduce reduce-only exp-reduce exp-reduce-only
# (TODO: do without for-loop. the loop forces remake of top-level
# index files for every scenario)
# the presence of sims/X/drm is the cue that X is an ensemble
exp-reduce: experiment-exists
	@ for d in sims/$(S)/*; do \
	        [ -d $$d/drm -a -d $$d/spc ] || continue; \
		d_prime=$$(echo $$d | sed 's:sims/::'); \
		echo $(MAKE) S=$$d_prime reduce; \
		$(MAKE) --no-print-directory S=$$d_prime reduce || exit $$?; \
	done
	@ echo "Make: Reducing overall experiment..."
	$(REDUCE_ENS_PROG) sims/$(S)

exp-reduce-only: experiment-exists
	@ echo "Make: Reducing ONLY overall experiment..."
	$(REDUCE_ENS_PROG) sims/$(S)

# 'make reduce' flows from sims/ down to $S/ through DIR/reduce-info.csv targets
# in all intermediate dirs: see the PROPAGATE_REDUCTION_UPWARD mechanism below
reduce: script-exists sims/reduce-info.csv

# 'make reduce-only' does not flow from sims -> $S: it just does the bottom level
reduce-only: script-exists sims/$(S)/reduce-info.csv

# dependence for a bottom-level reduction in sims/$S
sims/$(S)/reduce-info.csv: sims/$(S)/drm
	@ echo "Make: Reducing $< ..."
	$(REDUCE_PROG) -O sims/$(S)/reduce-%s.%s sims/$(S)/drm/*.pkl

# Below: a variable, a macro, and a foreach link the top-level reduce
# target (sims/reduce-info.csv) to the base-level one (sims/$S/reduce-info.csv)
# This allows subdirectories of scripts.

# awk one-liner -- transforms $S into a chain of intermediate directories:
#   S=HabExSample [drms in sims/HabExSample/drm/...] -->
#     HabExSample
#   S=a.fam/b.fam/sub [drms in sims/a.fam/b.fam/sub/drm/...] -->
#     a.fam
#     a.fam/b.fam
#     a.fam/b.fam/sub
REDUCE_CHAIN:=$(shell echo $S | awk -F/ '{for(i=1; i<=NF; i++) {for (j=1; j<i; j++) printf "%s/", $$j; printf "%s\n", $$i;}}')

# The chain of REDUCE_ENS_PROG invocations uses rules of the form:
#   ParentDir(LINK)/reduce-info.csv: LINK/reduce-info.csv
# where LINK (a.k.a. $1 here) is, for example: sims/HabExSample
# below, note $(dir $1) preserves the trailing slash on the enclosing dir
# [turmon 2025-02: previously the key line below was:
# [	$(REDUCE_ENS_PROG) -O $$(@D)/reduce-%s.%s $$(@D)/*/
define PROPAGATE_REDUCTION_UPWARD
$(dir $1)reduce-info.csv: $1/reduce-info.csv
	@ echo "Make: Reducing parent: $$(@D)"
	$(REDUCE_ENS_PROG) $$(@D)
endef

# Expand the above rule into one transformation for each intermediate dir.
#   for S=HabExSample, just one rule --
#     sims/reduce-info.csv: sims/HabExSample/reduce-info.csv
#   for S=a.fam/b.fam/sub, three rules, such as:
#     sims/a.fam/reduce-info.csv: sims/a.fam/b.fam/reduce-info.csv
# Some instances of LINK below can have embedded = signs (Experiments), and
# these will cause a make parse error in the define above. (The = will
# make the "targ=et: depende=ncy" look like a variable being set.)
# Thus, the reduction rule line will seem to not be connected to a rule.
# Solution: sub in $(EQUAL) for every literal = in the target and dependency.
$(foreach LINK,$(REDUCE_CHAIN),$(eval $(call PROPAGATE_REDUCTION_UPWARD,$(subst =,$$(EQUAL),sims/$(LINK)))))


########################################
## Graphics - detections, fuel use
## (original Matlab version)
##
.PHONY: graphics
# delegate to the graphics status file
graphics: script-exists sims/$(S)/gfx/det-info.txt

# just one ensemble's graphics
sims/$(S)/gfx/det-info.txt: sims/$(S)/reduce-info.csv
	@ echo "Make: Graphics plotting into $(@D) ..."
	@ rm -f sims/$(S)/gfx/det-*.*
	$(GRAPHICS_PROG) sims/$(S)/reduce-%s.%s sims/$(S)/gfx/det-%s.%s
	$(GRAPHYCS_PROG) sims/$(S)/reduce-%s.csv
########################################
## Graphics - detections, fuel use
## (matplotlib version)
##
.PHONY: graphics-alt
# delegate to an ALT graphics status file
graphics-alt: script-exists sims/$(S)/gfx/det-info-ALT.txt

# just one ensemble's graphics
# the cp moves the "ALT" sentinel file into place
sims/$(S)/gfx/det-info-ALT.txt: sims/$(S)/reduce-info.csv
	@ echo "Make: mpl Graphics plotting into $(@D) ..."
	@ rm -f sims/$(S)/gfx/det-*.*
	$(GRAPHALT_PROG) sims/$(S)/reduce-%s.%s sims/$(S)/gfx/det-%s.%s
	@ cp -p sims/$(S)/gfx/det-info.txt sims/$(S)/gfx/det-info-ALT.txt
	$(GRAPHYCS_PROG) sims/$(S)/reduce-%s.csv

########################################
## Tables - promotion funnel (more to come?)
##
.PHONY: tables
# delegate to the table status file
tables: script-exists sims/$(S)/tbl/table-status.txt

# just one ensemble's graphics
sims/$(S)/tbl/table-status.txt: sims/$(S)/reduce-info.csv
	@ echo "Make: Tables into $(@D) ..."
	@ rm -f sims/$(S)/tbl/table-*.*
	$(TABLES_PROG) -o sims/$(S)/tbl/table-%s.%s all sims/$(S)/reduce-%s.%s

########################################
## Detection visits tables - for scheduler analysis
##
.PHONY: star-visits
# delegate to the html document
star-visits: script-exists sims/$(S)/sched/detection-visits.html

# one ensemble's detection visit document
sims/$(S)/sched/detection-visits.html: sims/$(S)/drm
	@ echo "Make: Detection visits document into $(@D) ..."
	$(STAR_VISIT_PROG) sims/$(S)


########################################
## Path ensemble graphics - starshade slew map
##
# delegate to the 'path-ens' for the named script
path-ensemble: script-exists sims/$(S)/path-ens/path-map.png

# one ensemble's path plots - they depend on the DRM-set, not the reduction
sims/$(S)/path-ens/path-map.png: sims/$(S)/drm
	@ echo "Make: Making ensemble tour summary graphic in \`$(basename $@)'"
	$(PATH_ENS_PROG) sims/$(S)/drm

########################################
## Path movies
##   
#  target is: path-movie-N and path-final-N,
#  for N = 1, 2, 5, 10, 20, etc.

# Rule to make a single-drm path movie, given a SCENARIO
sims/$(S)/path/%.mp4: sims/$(S)/drm/%.pkl
	@ echo "Make: Path movie \`$@'"
	$(PATH_PROG) $<

# Rule to make a single-drm path final-frame, given a SCENARIO
sims/$(S)/path/%-final.png: sims/$(S)/drm/%.pkl
	@ echo "Make: Path final-frame \`$@'"
	$(PATH_PROG_FINAL) $<

# Rule to make a single-drm timeline plot-set, given a SCENARIO
sims/$(S)/path/%-obs-timelines.txt: sims/$(S)/drm/%.pkl
	@ echo "Make: Timeline \`$@'"
	$(TIMELINE_PROG) -o sims/$(S)/path/$(*)-%s.%s -j Scripts/$(S).json $<

# Rule to make a single-drm keepout map, given a SCENARIO
sims/$(S)/path/%-keepout-and-obs.png: sims/$(S)/drm/%.pkl
	@ echo "Make: Keepout \`$@'"
	$(KEEPOUT_PROG) -o sims/$(S)/path/$(*)-%s.%s Scripts/$(S).json $<

# Rule to make a group of movies, given a count.  The rule ends up looking like:
#  path-movie-N: sims/$(S)/path/SEED1.mp4 sims/SCRIPT/path/SEED2.mp4 ...
# where N is the number of movies (mp4's).  The number N is given to head
# as a ceiling on the number of path movie targets requested.  The .mp4 targets,
# in turn, are made by $(PATH_PROG) as shown above.
# The obs-timelines are not movies, but piggy-back on this same setup.
## Note: if sims/$(S)/drm/ does not exist (e.g., if S is an "experiment"),
## the $(... find ...) will be empty and the Makefile would be syntactically
## invalid (causing a hard error) *unless* some target is present. So,
## script-exists has a dual purpose - raise error if there's no script, and
## ensure the rule syntax is OK.

define MAKE_N_MOVIES
.PHONY: path-movie-$1 path-final-$1 obs-timeline-$1
path-movie-$1: script-exists $(shell $(SELECT_RUN_PROG) -n $1 sims/$(S) | sed -e 's:/drm/:/path/:' -e 's:\.pkl:.mp4:')
	@ echo "Make: Placed movies in \`sims/$(S)/path'."
path-final-$1: script-exists $(shell $(SELECT_RUN_PROG) -n $1 sims/$(S) | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-final.png:')
	@ echo "Make: Placed final-frames in \`sims/$(S)/path'."
obs-timeline-$1: script-exists $(shell $(SELECT_RUN_PROG) -n $1 sims/$(S) | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-obs-timelines.txt:')
	@ echo "Make: Placed obs-timelines in \`sims/$(S)/path'."
keepout-$1: script-exists $(shell $(SELECT_RUN_PROG) -n $1 sims/$(S) | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-keepout-and-obs.png:')
	@ echo "Make: Placed keepout in \`sims/$(S)/path'."
endef

# Allowable path-movie counts to construct targets for below
# gnu head accepts T as an abbreviation for Tera, so T is in effect an alias for everything
MOVIE_COUNTS:=1 2 5 10 20 50 100 T

# For each N in MOVIE_COUNTS, the foreach constructs targets of the form:
#   path-movie-N
#   path-final-N
$(foreach N,$(MOVIE_COUNTS),$(eval $(call MAKE_N_MOVIES,$N)))

########################################
## HTML indexes
##   
.PHONY: html html-all html-only

# delegate to the 'html/index.html' for the named script
html: script-exists reduce sims/$(S)/html/index.html;

# same as html, but omit re-making the graphics
#   for experiments, generates the top-level index only
html-only: script-exists
	@ echo "Make: HTML index (only) $@ ..."
	$(HTML_PROG) $(S)

# one ensemble's html summary
sims/$(S)/html/index.html: sims/$(S)/gfx/det-info.txt sims/$(S)/tbl/table-status.txt
	@ echo "Make: HTML index $@ ..."
	$(HTML_PROG) $(S)

# recursively regenerate all index.html's for all sims,
# and then regenerate the global index.html.  Does *not*
# imply re-reduction or graphics remake.
html-all:
	$(HTML_PROG) -r

########################################
## Experiments = dirs *containing* ensembles
##
## targets: exp-graphics*, exp-path-ensemble*, exp-html*
## Exceptions:
##   exp-reduce is handled separately (reduce must precede selection)
##   exp-html-only has a direct rule that regenerates everything, but
##       the sub-targets (exp-html-only-top-10, etc.) are also defined here.

# The make targets for Experiments are prefixed exp-*.
# Many are defined by macro expansion using the following rule.  It uses a
# selection helper routine to pick out some ensembles *within* the experiment,
# and runs a sub-Make on each such ensemble, using a shell for-loop.
# The "exp-reduce" target is *not* done that way, because the reduce
# needs to be made before selection makes sense.

# Rule to make a generic target for N ensembles within an experiment, given a count.
#   $1 = make target
#   $2 = top or mix; ensembles with top yield OR mixed selection of ensembles
#   $3 = number of ensembles N
#   $$$$d = subdirectory of the experiment directory
# Notes to this somewhat complex macro:
# * Change the eval(...) below to info(...) to debug this macro.
# * Text here is expanded twice by make.  The $$$$d construction is expanded the
# first time (within the eval) to $$d.  During the later invocation of the rule
# itself, make turns $$d into $d, that is, the shell variable "d".
# * "+" prepended to rule below transmits -j setting to sub-make.
# * There are two types of calls to the selector program:
#  one for "top" -- keying off # earth chars
#  one for "mix" -- keying off of the MD5 hash of the (string) experiment name

SELECT_PROG_top=$(SELECT_PROG) -k chars_earth_unique
SELECT_PROG_mix=$(SELECT_PROG) -k experiment

define MAKE_EXP_OPERATION
.PHONY: exp-$1-$2-$3
exp-$1-$2-$3: experiment-exists
	@+ for d in `$(SELECT_PROG_$2) -n $3 $2 sims/$(S)/reduce-yield-plus.csv`; do \
	        [ -d sims/$(S)/$$$$d/drm ] || continue; \
		echo $(MAKE) S=$(S)/$$$$d $1; \
		$(MAKE) --no-print-directory S=$(S)/$$$$d $1 || exit $$$$?; \
	done
endef


## path-ensemble
# targets: exp-path-ensemble-{top,mix}-N 
# such as: exp-path-ensemble-mix-10
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,path-ensemble,top,$N)))
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,path-ensemble,mix,$N)))
# default target for above -- supports 2-way parallelism
.PHONY: exp-path-ensemble
exp-path-ensemble: exp-path-ensemble-top-10 exp-path-ensemble-mix-20

## path-movie
# targets: exp-path-movie-M-{top,mix}-N
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,path-movie-$M,top,$N))))
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,path-movie-$M,mix,$N))))
# default target for above -- supports 2-way parallelism
# 5 movies each within the top-10 and mix-10 ensembles = 10*5 + 10*5 = 100 movies
.PHONY: exp-path-movie exp-path-movie-5
exp-path-movie: exp-path-movie-5;
exp-path-movie-5: exp-path-movie-5-top-10 exp-path-movie-5-mix-10

## keepout
# targets: exp-keepout-M-{top,mix}-N
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,keepout-$M,top,$N))))
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,keepout-$M,mix,$N))))
# default target for above -- supports 2-way parallelism
# 5 keepout-maps each within the top-10 and mix-10 ensembles = 10*5 + 10*5 = 100 maps
.PHONY: exp-keepout exp-keepout-5
exp-keepout: exp-keepout-5;
exp-keepout-5: exp-keepout-5-top-10 exp-keepout-5-mix-10

## obs-timeline
# targets: exp-obs-timeline-M-{top,mix}-N
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,obs-timeline-$M,top,$N))))
$(foreach M,$(MOVIE_COUNTS),$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,obs-timeline-$M,mix,$N))))
# default target for above -- supports 2-way parallelism
# 5 timelines each within the top-10 and mix-10 ensembles = 10*5 + 10*5 = 100 timelines
.PHONY: exp-obs-timeline exp-obs-timeline-5
exp-obs-timeline: exp-obs-timeline-5;
exp-obs-timeline-5: exp-obs-timeline-5-top-10 exp-obs-timeline-5-mix-10

## graphics
# targets: exp-graphics-{top,mix}-N
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,graphics,top,$N)))
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,graphics,mix,$N)))
# default target for above -- supports 2-way parallelism
.PHONY: exp-graphics
exp-graphics: exp-graphics-top-10 exp-graphics-mix-20

## html -- within the sub-make, it will require graphics
# targets: exp-html-{top,mix}-N
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,html,top,$N)))
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,html,mix,$N)))
# default target for above -- supports 2-way parallelism
.PHONY: exp-html
exp-html: exp-html-top-10 exp-html-mix-20

## html-only -- does not re-make graphics
# targets: exp-html-only-{top,mix}-N
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,html-only,top,$N)))
$(foreach N,$(EXP_COUNTS),$(eval $(call MAKE_EXP_OPERATION,html-only,mix,$N)))
# default target for above -- handled specially
#   re-generates all html for all the ensembles (-r option)
.PHONY: exp-html-only
exp-html-only: experiment-exists
	@ echo "Make: HTML index (full experiment) $@ ..."
	$(HTML_PROG) -r $(S)

########################################
## http server start, stop, status
##

.PHONY: html-ensure html-start html-stop html-status
html-ensure:
	util/html-serve.sh ensure
	@echo "Stop with \`make html-stop'."

html-start:
	util/html-serve.sh start
	@echo "Stop with \`make html-stop'."

html-stop:
	util/html-serve.sh stop

html-status:
	util/html-serve.sh status

