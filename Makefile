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
#   reduce:          reduce DRMs to tabulated CSV files for later plotting
#   graphics:        make detection/char plots, and radius-luminosity bar plots.
#   graphics-extra:  make even more detection/char plots
#   graphics-clean:  remove existing detection/char plots, to allow re-make
#   html:            re-generate the index.html that summarizes the given scenario
#   html-only:       same as html, but do not re-reduce the data or remake graphics.
#   path-ensemble:   make lon/lat plots of slews taken by an ensemble.
#   star-visits:     per-star tabulation of successful detection visits
#   path-movie-N:    make "N" path-movies and final frames 
#   path-final-N:    make "N" final frames, only
#   obs-timeline-N:  make "N" observing-target timelines
#   keepout-N:       make "N" keepout-vs-time plots
#                    (the -N targets choose N arbitrary sims, where 
#                    N = 1, 2, 5, 10, 20, 50, 100, or T, where T=all)
#   tar-log:         replace scenario/log with its "tar" archive to save space
#   status:          list the current contents of DRMs for this scenario (like "ls")
# (2) Multi-script reduction and plotting
#   All these targets require an *experiment* name.
#   exp-reduce:      makes "reduce" for all ensembles within the experiment
#   exp-html-top-N:  makes html (inc. graphics) for the N top (by yield) ensembles
#   exp-html-mix-N:  makes html (inc. graphics) for N selected-arbitrarily ensembles
#   exp-html:        makes html for 10 top + 20 selected ensembles - can use make -jN
#   exp-path-ensemble* \   Same pattern as html above with -mix or -top, and a
#   exp-graphics*       \  number saying how many.  Also, can leave off -top-N
#   exp-html-only*      /  and just make 10 top + 20 selected.
#   exp-path-movie-M-*  /  Make M movies in each of (top/mix)-N ensembles.
#   exp-keepout-M-*    /   Make M keepout maps in each of (top/mix)-N ensembles.
#   exp-obs-timeline-M-*   Make M obs-timelines in each of (top/mix)-N ensembles.
#   (The exp-* targets first reduce any ensembles that need it.)
# (3) Web-server
#   html-ensure: start Apache httpd web-server, if not running already
#   html-start: start Apache httpd web-server
#   html-stop: stop Apache httpd web-server
#   html-status: show running web-servers, if any
#
## turmon oct 2017, mar 2018, feb 2022

# set the default shell
SHELL:=/bin/bash

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
# 2c(i): pull off leading sims/, if present (S=sims/foo -> S=foo)
S := $(patsubst sims/%,%,$(S))
# 2c(ii): allow for just sims *only*, if present (S=sims -> S=)
S := $(patsubst sims,,$(S))
# 2d: strip trailing / which could be present (S=sims/foo/ -> S=foo)
S := $(patsubst %/,%,$(S))
# 2e: strip added space at the end of S
S := $(strip $(S))
# 3: repeat script value back, if supplied
#    (but not again if make restarts, after remaking an included makefile)
ifdef S_COPY
 ifndef MAKE_RESTARTS
  $(info Make: Scenario name: "$(S)")
 endif
endif

# needed to escape equal signs in some Experiments, alas
EQUAL:= =

# OS name (Darwin/Linux, typically)
# (alternatively: is $(cwd) on same partition as $(uv cache dir)?)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    # on HPC, mustang, etc. -- local cache, not ~/, for efficiency
    UV_PREFIX=UV_CACHEDIR=./.uv-cache/
else
    # on single-user laptop -- leave cache as-is
    UV_PREFIX=
endif

# program used for data reduction
REDUCE_PROG=util/reduce_drms.py
# program used for reduction of multi-ensemble experiments
# (-E: expand the given "parent" dir to subdirs with Ensembles)
REDUCE_ENS_PROG=util/reduce_drm_sets.py -E
# programs used for graphics
GRAPHICS_PROG=util/plot_drm_driver.py
GRAPHYCS_PROG=util/rad-sma-rectangle-plot-driver.sh -s
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
# html summary of ensembles/experiments
#   -i: to regenerate the global index.html as well as that for $(S)
HTML_PROG=util/html-summary.py -i
#   same, but without -i: index just the named sim(s)
HTML_PROG_NOINDEX=util/html-summary.py
# html summary of emulator Analysis/ results
EMU_HTML_PROG=$(UV_PREFIX) util/emulator_html_summary.py -R Local/www-resources -S ensemble-reports.css -J sorttable.js
# analysis/plots of experiment results
EMU_PLOT_PROG=$(UV_PREFIX) util/run_emulator_workflows.py
# program to record the orderings of ensembles within an experiment,
#  as a Makefile fragment (see EXP_SELECT_MK)
EXP_SELECT_PROG=util/select_ensembles_exp.sh
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

analysis-exists:
	@ [ -d sims/$(S)/Analysis ] || \
		(echo "Require an experiment/family Analysis directory \`sims/$(S)/Analysis'" && exit 1)


.PHONY: default script-exists experiment-exists analysis-exists

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
# 'make exp-reduce' flows through the same dependency graph as 'make reduce':
#   sims/reduce-info.csv <- ... <- sims/$(S)/reduce-info.csv <- each ensemble's
# Because this is one make process, every node is reduced exactly once, however
# many ensembles the experiment holds.  (This was formerly a shell for-loop
# running a sub-make per ensemble; each of those independently walked the
# PROPAGATE_REDUCTION_UPWARD chain up to sims/, so every ancestor directory was
# re-reduced once per ensemble.)
exp-reduce: experiment-exists sims/reduce-info.csv

exp-reduce-only: experiment-exists
	@ echo "Make: Reducing ONLY overall experiment..."
	$(REDUCE_ENS_PROG) sims/$(S)

# 'make reduce' flows from sims/ down to $S/ through DIR/reduce-info.csv targets
# in all intermediate dirs: see the PROPAGATE_REDUCTION_UPWARD mechanism below
reduce: script-exists sims/reduce-info.csv

# 'make reduce-only' does not flow from sims -> $S: it just does the bottom level
reduce-only: script-exists sims/$(S)/reduce-info.csv

# dependence for a bottom-level reduction: any directory holding a drm/ is an
# ensemble.  Writes DIR/reduce-info.csv, and many others.
# This is a pattern rule, not an explicit rule for $(S) alone, so that the
# ensembles *within* an experiment can also be built as prerequisites (see
# exp-reduce above).  It declines to match container directories, which have
# no drm/, so those fall through to PROPAGATE_REDUCTION_UPWARD below.
sims/%/reduce-info.csv: sims/%/drm
	@ echo "Make: Reducing $< ..."
	$(REDUCE_PROG) $<

# When $(S) is an experiment/family -- it has no drm/ of its own -- its
# reduction depends on the reduction of every ensemble it contains.  This rule
# is what replaces the old exp-reduce for-loop.  As before, the presence of
# both drm/ and spc/ is the cue that a subdirectory is an ensemble.
# An ensemble with no DRMs yet (e.g., runs in progress) is skipped: reducing
# it would succeed without writing reduce-info.csv, so it would be perpetually
# out of date -- and, via EXP_SELECT_MK below, make would restart forever.
# The ifeq guard matters: when $(S) is an ensemble this rule must not exist,
# or its recipe would shadow the pattern rule above.
EXP_ENSEMBLES     = $(patsubst %/drm,%,$(wildcard sims/$(S)/*/drm))
EXP_ENSEMBLE_CSVS = $(foreach d,$(EXP_ENSEMBLES),\
                      $(if $(and $(wildcard $d/spc),$(wildcard $d/drm/*.pkl)),$d/reduce-info.csv))
ifeq ($(wildcard sims/$(S)/drm),)
sims/$(S)/reduce-info.csv: $(EXP_ENSEMBLE_CSVS)
	@ echo "Make: Reducing overall experiment: $(@D) ..."
	$(REDUCE_ENS_PROG) $(@D)
endif

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
## Graphics
##

.PHONY: graphics graphics-clean graphics-extra
# distinguished sentinel file for make dependency chain
GRAPHICS_SENTINEL:=sims/$(S)/gfx/det-info.txt

# ** This is the main graphics target **
# delegate to the graphics sentinel file
graphics: script-exists $(GRAPHICS_SENTINEL)

# newer graphics - one ensemble
# This, and the other per-ensemble product rules below, are pattern rules
# (keyed on the ensemble directory, %) rather than explicit rules for $(S)
# alone.  That way the ensembles *within* an experiment can be built by the
# exp-* targets in this same make process, without a sub-make per ensemble.
sims/%/gfx/det-info.txt: sims/%/reduce-info.csv
	@ echo "Make: Graphics (new-format) into $(@D) ..."
	@ rm -f sims/$*/gfx/det-*.*
	$(GRAPHICS_PROG) sims/$*/reduce-%s.%s sims/$*/gfx/det-%s.%s
	$(GRAPHYCS_PROG) sims/$*/reduce-%s.csv

# imperatively remove existing graphics, allowing clean re-make
graphics-clean: script-exists
	@ echo "Make: Removing existing graphics in sims/$(S)/gfx ..."
	rm -f sims/$(S)/gfx/det-*.*

# extra (and normal) graphics - one ensemble
# this is imperative, not delegated to $(GRAPHICS_SENTINEL)
graphics-extra: script-exists sims/$(S)/reduce-info.csv
	@ echo "Make: Graphics (normal + extras) into sims/$(S)/gfx ..."
	@ rm -f sims/$(S)/gfx/det-*.*
	$(GRAPHICS_PROG) --mode_op + sims/$(S)/reduce-%s.%s sims/$(S)/gfx/det-%s.%s
	$(GRAPHYCS_PROG) sims/$(S)/reduce-%s.csv

########################################
## Tables - promotion funnel (more to come?)
##
.PHONY: tables
# delegate to the table status file
tables: script-exists sims/$(S)/tbl/table-status.txt

# just one ensemble's tables
sims/%/tbl/table-status.txt: sims/%/reduce-info.csv
	@ echo "Make: Tables into $(@D) ..."
	@ rm -f sims/$*/tbl/table-*.*
	$(TABLES_PROG) -o sims/$*/tbl/table-%s.%s all sims/$*/reduce-%s.%s

########################################
## Detection visits tables - for scheduler analysis
##
.PHONY: star-visits
# delegate to the html document
star-visits: script-exists sims/$(S)/sched/detection-visits.html

# one ensemble's detection visit document
sims/%/sched/detection-visits.html: sims/%/drm
	@ echo "Make: Detection visits document into $(@D) ..."
	$(STAR_VISIT_PROG) sims/$*


########################################
## Path ensemble graphics - starshade slew map
##
# delegate to the 'path-ens' for the named script
path-ensemble: script-exists sims/$(S)/path-ens/path-map.png

# one ensemble's path plots - they depend on the DRM-set, not the reduction
sims/%/path-ens/path-map.png: sims/%/drm
	@ echo "Make: Making ensemble tour summary graphic in \`$(basename $@)'"
	$(PATH_ENS_PROG) sims/$*/drm

########################################
## Path movies
##   
#  target is: path-movie-N and path-final-N,
#  for N = 1, 2, 5, 10, 20, etc.

# Enable deferred ("secondary") expansion of prerequisites for the rules
# defined below.  This does two jobs:
# (1) The per-DRM rules below map sims/ENS/path/SEED.* to sims/ENS/drm/SEED.pkl.
#     That needs two stems (ENS and SEED), but a pattern rule has only one.
#     So the stem is ENS/path/SEED, and $$(subst ...) computes the .pkl from it.
# (2) It makes the $(SELECT_RUN_PROG) calls, and the ensemble selection for
#     exp-* targets, lazy: they run only for a target actually asked for,
#     instead of once per (count x target-kind) combination on every single
#     invocation of make.
.SECONDEXPANSION:

# Script file for a per-DRM product: sims/ENS/path/FILE -> Scripts/ENS.json
PATH_TO_SCRIPT = $(patsubst sims/%/path/,Scripts/%.json,$(dir $1))

# Rule to make a single-drm path movie
sims/%.mp4: $$(subst /path/,/drm/,sims/$$*).pkl
	@ echo "Make: Path movie \`$@'"
	$(PATH_PROG) $<

# Rule to make a single-drm path final-frame
sims/%-final.png: $$(subst /path/,/drm/,sims/$$*).pkl
	@ echo "Make: Path final-frame \`$@'"
	$(PATH_PROG_FINAL) $<

# Rule to make a single-drm timeline plot-set
sims/%-obs-timelines.txt: $$(subst /path/,/drm/,sims/$$*).pkl
	@ echo "Make: Timeline \`$@'"
	$(TIMELINE_PROG) -o sims/$(*)-%s.%s -j $(call PATH_TO_SCRIPT,$@) $<

# Rule to make a single-drm keepout map
sims/%-keepout-and-obs.png: $$(subst /path/,/drm/,sims/$$*).pkl
	@ echo "Make: Keepout \`$@'"
	$(KEEPOUT_PROG) -o sims/$(*)-%s.%s $(call PATH_TO_SCRIPT,$@) $<

## Note: script-exists is the first prerequisite of each rule below.  It
## raises a clear error when S names neither a script nor an experiment
## (e.g., if S is an "experiment", sims/$(S)/drm/ does not exist and the
## selector below quietly returns nothing).

# Map a run-count onto the list of per-DRM products to build.
#   $1 = number of runs to select (or T for all)
#   $2 = suffix that replaces ".pkl" on each selected DRM
#   $3 = the ensemble directory, sims/...
# NOTE: the shell command lives here in a variable, rather than inline as
# $$(shell ...) in the prerequisite lists below.  Make scans a rule line for
# the target/prereq ":" before expanding anything, so a literal ":" inside
# $$(...) would be misread as a second rule separator ("*** multiple target
# patterns.  Stop.").  Keeping the command here makes the rule lines immune
# to the choice of sed delimiter.
SELECT_RUN_TARGETS = $(shell $(SELECT_RUN_PROG) -n $1 $3 | \
                       sed -e 's:/drm/:/path/:' -e 's:\.pkl:$2:')

# Targets to make a group of per-DRM products, given a count N, e.g.:
#   path-movie-5: sims/$(S)/path/SEED1.mp4 sims/$(S)/path/SEED2.mp4 ...
# N may be any non-negative integer, or T for all runs.  The individual
# products are made by the per-DRM rules above ($(PATH_PROG) and friends).
# The obs-timelines and keepout maps are not movies, but piggy-back on the
# same setup.
# NOTE: these targets must NOT be declared .PHONY.  Make skips pattern-rule
# search for phony targets, which would silently disable all four rules.
path-movie-%: script-exists $$(call SELECT_RUN_TARGETS,$$*,.mp4,sims/$(S))
	@ echo "Make: Placed movies in \`sims/$(S)/path'."

path-final-%: script-exists $$(call SELECT_RUN_TARGETS,$$*,-final.png,sims/$(S))
	@ echo "Make: Placed final-frames in \`sims/$(S)/path'."

obs-timeline-%: script-exists $$(call SELECT_RUN_TARGETS,$$*,-obs-timelines.txt,sims/$(S))
	@ echo "Make: Placed obs-timelines in \`sims/$(S)/path'."

keepout-%: script-exists $$(call SELECT_RUN_TARGETS,$$*,-keepout-and-obs.png,sims/$(S))
	@ echo "Make: Placed keepout in \`sims/$(S)/path'."

# Per-DRM counts used to construct the exp-* targets further below
# (the single-scenario -N targets above accept any N, so they do not use this).
# T is the abbreviation for "all runs".
MOVIE_COUNTS:=1 2 5 10 20 50 100 T

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
# When the ensemble is $(S) itself, -i (in HTML_PROG) also re-indexes the
# enclosing sims.  When it is an ensemble within experiment $(S) (see exp-html
# below), many of these can run at once under -j, and all their -i's would
# rewrite the same parent indexes.  So instead, flag the experiment's index
# as stale, and the exp-html target re-indexes once, at the end.
sims/%/html/index.html: sims/%/gfx/det-info.txt sims/%/tbl/table-status.txt
	@ echo "Make: HTML index $@ ..."
	$(if $(filter $*,$(S)),$(HTML_PROG),$(HTML_PROG_NOINDEX)) $*
	$(if $(filter $*,$(S)),,@ touch sims/$(S)/$(EXP_HTML_STALE))

# recursively regenerate all index.html's for all sims,
# and then regenerate the global index.html.  Does *not*
# imply re-reduction or graphics remake.
html-all:
	$(HTML_PROG) -r

########################################
## Experiments = dirs *containing* ensembles
##
## targets: exp-graphics*, exp-path-ensemble*, exp-html*, and per-DRM
##   products (exp-path-movie*, exp-keepout*, exp-obs-timeline*)
## Exceptions:
##   exp-reduce is handled separately (see "Data reductions")
##   exp-html-only has a direct rule that regenerates everything, but
##       the sub-targets (exp-html-only-top-10, etc.) are also defined here.

# The make targets for Experiments are prefixed exp-*.  Each uses a
# selection helper to pick out some ensembles *within* the experiment (by
# "top" yield, or by an arbitrary-but-stable "mix"), and then depends
# directly on the products (html, graphics, movies...) for those ensembles.
# These are built by the same pattern rules used for a single ensemble, all
# within this one make process: there is no sub-make per ensemble.
#
# The catch is that selection reads reduce-yield-plus.csv, so the reduction
# must be up to date *before* the selection is made.  (Otherwise, a newly-run
# ensemble not yet reduced would be invisible to the selector, and hence
# would never be reduced or processed.)  A plain $(shell ...) runs too early,
# when the Makefile is read.  So instead, the selection is recorded in a
# generated Makefile fragment, EXP_SELECT_MK, which is included below and
# depends on the full reduction.  GNU make first brings included makefiles up
# to date -- here, by reducing, and then regenerating the fragment -- and
# then restarts itself to read the new version.  This uses only features
# present in older GNU make (3.81+), so it does not need .WAIT (4.4+).
#
# EXP_SELECT_MK records the *full* ordering of ensembles, for "top" and for
# "mix", so one file serves every N.

# generated Makefile fragment defining EXP_ORDER_top and EXP_ORDER_mix
#   "top": by yield (# earth chars); "mix": by MD5 hash of the ensemble name
EXP_SELECT_MK:=sims/$(S)/exp-select.mk
# flag file (within the experiment): an ensemble html index was remade
EXP_HTML_STALE:=.exp-html-stale
# the goals (if any) that need EXP_SELECT_MK
EXP_SELECT_GOALS:=$(filter-out exp-reduce exp-reduce-only exp-preflight exp-html-only exp-analysis-%,\
                    $(filter exp-%,$(MAKECMDGOALS)))
# nonempty if this is a dry run (make -n)
DRY_RUN:=$(findstring n,$(firstword -$(MAKEFLAGS)))

# Only include (and thus possibly reduce) when an exp-* goal needs selection:
# included makefiles are always brought up to date, so an unguarded include
# would reduce the whole experiment for any goal at all.
ifneq ($(and $(EXP_SELECT_GOALS),$(wildcard Scripts/$(S)/.),$(wildcard sims/$(S)/.)),)
include $(EXP_SELECT_MK)

# Under make -n, GNU make still *really runs* the rule for an included
# makefile, and its prerequisites -- here, the full reduction.  For dry runs,
# drop that prerequisite, so the selection is made from existing CSVs.  (The
# reduction still appears in the dry run, as a prerequisite of each target.)
# Also drop it after a restart (MAKE_RESTARTS is set): if some reduction
# recipe succeeded without writing its target, that target would be remade,
# and EXP_SELECT_MK with it, on every restart -- an endless loop.  This
# ensures at most one restart.
$(EXP_SELECT_MK): $(if $(or $(DRY_RUN),$(MAKE_RESTARTS)),,sims/reduce-info.csv)
	@ echo "Make: Selecting ensembles within $(@D) ..."
	$(EXP_SELECT_PROG) -o $@ sims/$(S)
endif

# Ensemble directories selected within the experiment.
#   $1 = selection(s), each as MODE/N, e.g., top/10 or "top/10 mix/20"
#        where MODE is top or mix, and N is a count, or T for all
# Only directories with a drm/ qualify.  The sort removes duplicates, so an
# ensemble chosen by both top/10 and mix/20 is only listed once.
EXP_SELECT   = $(sort $(foreach s,$1,$(call EXP_SELECT_1,$(firstword $(subst /, ,$s)),$(lastword $(subst /, ,$s)))))
EXP_SELECT_1 = $(foreach d,$(if $(filter T,$2),$(EXP_ORDER_$1),$(wordlist 1,$2,$(EXP_ORDER_$1))),\
                 $(if $(wildcard sims/$(S)/$d/drm),sims/$(S)/$d))

# Products for one ensemble, for each kind of exp-* operation.
#   $1 = ensemble directory (sims/...)
#   $2 = number of per-DRM products (movies, etc.) within the ensemble
# Intermediate files in the chain (reduce-info.csv, det-info.txt, ...) are
# listed explicitly.  Otherwise make would regard them as "intermediate", and
# delete them after the build.
EXP_PRODUCTS_html          = $1/reduce-info.csv $1/gfx/det-info.txt $1/tbl/table-status.txt $1/html/index.html
EXP_PRODUCTS_html-only     =
EXP_PRODUCTS_graphics      = $1/reduce-info.csv $1/gfx/det-info.txt
EXP_PRODUCTS_path-ensemble = $1/path-ens/path-map.png
EXP_PRODUCTS_path-movie    = $(call SELECT_RUN_TARGETS,$2,.mp4,$1)
EXP_PRODUCTS_keepout       = $(call SELECT_RUN_TARGETS,$2,-keepout-and-obs.png,$1)
EXP_PRODUCTS_obs-timeline  = $(call SELECT_RUN_TARGETS,$2,-obs-timelines.txt,$1)

# Products for all selected ensembles.  $1 = kind, $2 = count, $3 = selection(s)
EXP_PRODUCTS = $(foreach d,$(call EXP_SELECT,$3),$(call EXP_PRODUCTS_$1,$d,$2))

# Recipes run after the products are made, for each kind (most have none).
#   $1 = selection(s)
# html: re-index the experiment and its enclosing sims, once, if an
# ensemble's index was remade (see the per-ensemble html rule)
EXP_RECIPE_html = @ if [ -e sims/$(S)/$(EXP_HTML_STALE) ]; then \
	  echo "Make: HTML index (experiment) $(S) ..."; \
	  $(HTML_PROG) $(S) && rm -f sims/$(S)/$(EXP_HTML_STALE); fi
# html-only: index the selected ensembles, and their enclosing sims
EXP_RECIPE_html-only = $(if $(call EXP_SELECT,$1),$(HTML_PROG) $(patsubst sims/%,%,$(call EXP_SELECT,$1)))

# Rule to make a generic target for ensembles within an experiment.
#   $1 = make target, e.g., exp-html-top-10
#   $2 = kind: html, html-only, graphics, path-ensemble,
#        path-movie, keepout, or obs-timeline
#   $3 = number of per-DRM products M (for path-movie etc.; else empty)
#   $4 = selection(s), as for EXP_SELECT
# Notes to this somewhat complex macro:
# * Change the eval(...) below to info(...) to debug this macro.
# * Text here is expanded by the call, by the eval, and then (for the
#   prerequisites) again by .SECONDEXPANSION, hence the $$$$.  The deferred
#   expansion means the selectors run only for the targets actually made.
# * sims/reduce-info.csv is already up to date (EXP_SELECT_MK depends on
#   it), except in a dry run, where it shows the reductions to be done.
define MAKE_EXP_OPERATION
.PHONY: $1
$1: experiment-exists sims/reduce-info.csv $$$$(call EXP_PRODUCTS,$2,$3,$4)
	$$(call EXP_RECIPE_$2,$4)
endef

# Define the families of targets, for each kind:
#   exp-KIND-{top,mix}-N, e.g., exp-html-mix-10
#   exp-KIND-M-{top,mix}-N, e.g., exp-path-movie-5-top-10 (per-DRM products)
$(foreach K,path-ensemble graphics html html-only,\
  $(foreach X,top mix,\
    $(foreach N,$(EXP_COUNTS),\
      $(eval $(call MAKE_EXP_OPERATION,exp-$K-$X-$N,$K,,$X/$N)))))
$(foreach K,path-movie keepout obs-timeline,\
  $(foreach M,$(MOVIE_COUNTS),\
    $(foreach X,top mix,\
      $(foreach N,$(EXP_COUNTS),\
        $(eval $(call MAKE_EXP_OPERATION,exp-$K-$M-$X-$N,$K,$M,$X/$N))))))

# Default targets for the above: 10 top + 20 mix ensembles (or 10 + 10 for
# per-DRM products).  These can use make -jN for N-way parallelism.
$(eval $(call MAKE_EXP_OPERATION,exp-path-ensemble,path-ensemble,,top/10 mix/20))
$(eval $(call MAKE_EXP_OPERATION,exp-graphics,graphics,,top/10 mix/20))
$(eval $(call MAKE_EXP_OPERATION,exp-html,html,,top/10 mix/20))
# 5 movies each within the top-10 and mix-10 ensembles = 10*5 + 10*5 = 100 movies
# (and likewise for keepout maps and obs-timelines)
$(eval $(call MAKE_EXP_OPERATION,exp-path-movie-5,path-movie,5,top/10 mix/10))
$(eval $(call MAKE_EXP_OPERATION,exp-keepout-5,keepout,5,top/10 mix/10))
$(eval $(call MAKE_EXP_OPERATION,exp-obs-timeline-5,obs-timeline,5,top/10 mix/10))
.PHONY: exp-path-movie exp-keepout exp-obs-timeline
exp-path-movie: exp-path-movie-5;
exp-keepout: exp-keepout-5;
exp-obs-timeline: exp-obs-timeline-5;

## html-only -- does not re-make graphics
# default target for above -- handled specially
#   re-generates all html for all the ensembles (-r option)
.PHONY: exp-html-only
exp-html-only: experiment-exists
	@ echo "Make: HTML index (full experiment) $@ ..."
	$(HTML_PROG) -r $(S)

## emulation -- in development
# make the HTML
.PHONY: exp-analysis-html
exp-analysis-html: analysis-exists exp-analysis-graphics
	@ echo "Make: Emulator analysis html (experiment) $@ ..."
	$(EMU_HTML_PROG) sims/$(S)
	$(HTML_PROG) $(S)

.PHONY: exp-analysis-graphics exp-analysis-graphics-all

# distinguished sentinel file for post-experiment analysis
# base workflow = analysis, graphics, tables
EXP_ANALYSIS_SENTINEL:=sims/$(S)/Analysis/base.wfl/graphics-info.txt

# delegate to the exp-analysis sentinel file
exp-analysis-graphics: analysis-exists $(EXP_ANALYSIS_SENTINEL);

# (Also will depend on one or more reduce-yield-plus.csv files,
# but the particular files are given in the workflow script)
$(EXP_ANALYSIS_SENTINEL): sims/$(S)/Analysis/workflow-base.json
	@ echo "Make: Emulator analysis (experiment) $^ ..."
	$(EMU_PLOT_PROG) -C "" sims/$(S)/Analysis/workflow-base.json

# don't remake this
sims/$(S)/Analysis/workflow-base.json: ;

# analysis, graphics, tables -- all available workflows
# this runs without checking dependencies
exp-analysis-graphics-all: analysis-exists
	@ echo "Make: Emulator analysis graphics (all) (experiment) $@ ..."
	shopt -s nullglob && $(EMU_PLOT_PROG) -C "" sims/$(S)/Analysis/workflow-*.json

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

########################################
## documentation (per-installation)
##
.PHONY: doc doc-lint doc-publish
doc:
	@echo "=============================="
	@echo ">>  Making plot documentation"
	@echo "=============================="
	cd Local/www-doc && make install
	@echo "=============================="
	@echo ">>  Making code documentation"
	@echo "=============================="
	cd util/doc_sandbox && make doc
	@echo "=============================="
	@echo ">>  Moving code docs into place"
	@echo "=============================="
	cd util/doc_sandbox && make export

# check doc formatting: docstring blocks in scripts
doc-lint:
	util/dev/doc-lint.py

# push the generated code documentation to github pages.
# NOTE: this pushes to the gh-pages branch. GH_REMOTE picks the site:
#   origin (default) => the JPL-internal Pages site
#   export           => the github.com  (turmon.github.io) site
# e.g. "make doc-publish GH_REMOTE=export" to update the public site.
doc-publish:
	cd util/doc_sandbox && make publish
