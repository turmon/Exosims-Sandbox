# Makefile for EXOSIMS data reduction and plot generation.
#
# Usage:
#    $ make <target>
#  or
#    $ make S=SCENARIO <target>
# where <target> is one of the below.  Use "make -B" to force a target to be
# regenerated, e.g., when the underlying plot code changes.
# In general, the S=SCENARIO argument is unused for ipython-parallel targets,
# but required for data reduction or plotting.  The SCENARIO is the base name
# of the relevant script file, i.e., for Scripts/HabEx_4m_TS_20180201.json,
# use S=HabEx_4m_TS_20180201
#
# Simulations are added using the `add-sim.sh' command, separate from "make".
#
# Targets:
# (1) Data reduction and plotting targets --
#   Note, all these targets require a scenario name.
#   status: list the current contents of DRMs for this scenario (like "ls")
#   reduce: reduce DRMs to tabulated CSV files for later plotting
#   graphics: make detections-vs-time plots, and radius-luminosity bar plots.
#   path-ensemble: make lon/lat plots of slews taken by an ensemble.
#   path-movie-N: make "N" path-movies and final frames 
#   path-final-N: make "N" final frames, only
#   obs-timeline-N: make "N" observing-target timelines
#   keepout-N: make "N" keepout-vs-time plots
#     (above targets make selected graphics for N sims chosen from those present,
#     for N = 1, 2, 5, 10, 20, 50, 100, or T, where T=all)
#   html: re-generate the index.html that summarizes the given scenario
#   html-only: same as html, but do not re-reduce the data or remake graphics.
# (2) Multi-script reduction and plotting targets --
#   Note, all these targets require an *experiment* name.
#   exp-reduce: makes "reduce" for all ensembles within the experiment
#   exp-html-top-N: makes html (inc. graphics) for the N top (by yield) ensembles
#   exp-html-mix-N: makes html (inc. graphics) for N selected-arbitrarily ensembles
#   exp-html: makes html for 10 top + 20 selected ensembles - can use make -j2
#   exp-path-ensemble* \   Same pattern as html above with -mix or -top, and a
#   exp-graphics*       \  number saying how many.  Also, can leave off -top-N
#   exp-html-only*      /  and just make 10 top + 20 selected.
#   exp-movie-M-*      /   For movie, can make M movies in each of N ensembles.
# (3) Ipython-parallel ("ipp") targets --
#   Note: targets apply only to the machine where the "make" is run (e.g., aftac1).
#   ipp-create: create an ipython-parallel profile for this ensemble (use just once).
#       Copies several files into the ipyparallel directory.
#       To undo, see ipp-nuke, below.
#   ipp-start: start the ipython-parallel controller + engines
#       Note: If EXOSIMS_ENGINE_N (an integer) is exported from the environment,
#       this many engines will be started up, otherwise, the system-default
#       will be used.  To use this, run, from the shell:
#         $ EXOSIMS_ENGINE_N=8 make ipp-start
#   ipp-stop: stop the above.  See also ipp-kill, below.
#   ipp-status: report status of the controller + engines
#       Note: This attempts to run trivial jobs on the remote engines, so it
#       will not work if the engines are busy with another job.  See ipp-ps, below.
#   ipp-ps: use the unix "ps" command to identify the number of running engines
#   ipp-kill: sometimes ipp-stop does not work and engines or controllers
#       are orphaned.  ipp-kill identifies these by process id, and kills them.
#   ipp-nuke: deletes your ipython-parallel profile.  The inverse of ipp-create.
#       (Note: attempts to "ipp-kill" first, so as to not leave engines running.)
# (4) Web-server targets:
#   html-start: start Apache httpd web-server
#   html-stop: stop Apache httpd web-server
#   html-status: show running web-servers, if any
#
## turmon oct 2017, mar 2018

# clear builtin suffix rules (e.g., for .c, etc.)
.SUFFIXES:

## # debugging: echo rule chains
# $(warning Debug mechanism on)
# OLD_SHELL := $(SHELL)
# SHELL = $(warning Building $@$(if $<, (from $<))$(if $?, ($? newer)))$(OLD_SHELL)

# options for ipcluster startup
IPCLUSTER_OPT:=--daemonize --clean-logs=True
# allow number of engines to be set by environment
ifdef EXOSIMS_ENGINE_N
  IPCLUSTER_OPT:=-n $(EXOSIMS_ENGINE_N) $(IPCLUSTER_OPT)
endif

# hostname
HOST_NAME=$(shell hostname -s)

# python interpreter -- if not in environment, it is just "python"
# Note, it's better to just ask the *user* to set up their shell PATH 
# so that the right python is run.  The construction below, if used
# and propagated elsewhere, will eventually cause confusion.
PYTHON?=python
# PYTHON=/usr/local/anaconda3/envs/cornell/bin/python

# program used for data reduction
REDUCE_PROG=util/reduce_drms.py
# program used for reduction of multi-ensemble experiments
REDUCE_ENS_PROG=util/reduce_drm_sets.py
# programs used for matlab/python graphics
GRAPHICS_PROG=util/plot_drms.sh
GRAPHYCS_PROG=util/rad-sma-rectangle-plot-driver.sh -q
# generates tables
TABLES_PROG=util/tabulate_csv.py -q
# Programs analogous to the path-movie maker...
#   program for making (per-drm) timelines
TIMELINE_PROG=util/plot-timeline.py -d
#   program for making (per-drm) keepout + observation timelines
KEEPOUT_PROG=PYTHONPATH=EXOSIMS:Local util/plot-keepout-and-obs.py
# drm path-movie maker, generating command lines like:
#  drm-to-movie.sh -c -l 0 sims/HabEx_4m_TS_dmag26p0_20180206/drm/297992296.pkl
# -c specifies to make the cumulative plots as well (which does not
# make sense for the FINAL mode).
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
EXP_COUNTS:=1 2 5 10 20 50 100

# ipython parallel setup directory
IPYDIR=ipyparallel/$(USER)/$(HOST_NAME)
IPYCLI=$(IPYDIR)/security/ipcontroller-client.json

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
	@ [ -d Scripts/$(S) ]

.PHONY: script-exists experiment-exists default

# do not try to remake Makefile
Makefile:;

# ##
# ## Playground: invoking python simulations
# ##  --- not currently used ---

# # NB, the -1 option does not exist now
# RUNSIMS_PROG=add-sims.sh -1

# NSIMS := 10
# NUMBERS := $(shell seq 1 ${NSIMS})
# SIMS := $(addprefix sim,${NUMBERS})
# .PHONY: all-sims ${SIMS}
# all-sims: ${SIMS}
# 	@ echo "$@ success: made all sims"

# ${SIMS}: sim%:
# 	@ echo Running: ${RUNSIMS_PROG} $*

# ##### (end playground)

########################################
## Simulation status

# list drms that have been made
status: script-exists
	util/drm-ls.py -l sims/$(S)/drm

########################################
## Data reductions
##
.PHONY: reduce exp-reduce
# (there may be a non-for-loop way to do this)
exp-reduce: experiment-exists
	@ for d in sims/$(S)/s_*; do \
	        [ -d $$d/drm ] || continue; \
		d_prime=$$(echo $$d | sed 's:sims/::'); \
		echo $(MAKE) S=$$d_prime reduce; \
		$(MAKE) --no-print-directory S=$$d_prime reduce || exit $$?; \
	done
	@ echo "Make: Reducing overall experiment..."
	$(REDUCE_ENS_PROG) -i Scripts/$(S)/s_index.json -O sims/$(S)/reduce-%s.%s sims/$(S)/*

# delegate to the 'reduce' for the named script
reduce: script-exists sims/$(S)/reduce-info.csv 

# perform one reduction
sims/$(S)/reduce-info.csv: sims/$(S)/drm
	@ echo "Make: Reducing $< ..."
	$(REDUCE_PROG) -O sims/$(S)/reduce-%s.%s sims/$(S)/drm/*.pkl

########################################
## Graphics - detections, fuel use
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

# Rule to make a single path movie, given a SCENARIO
sims/$(S)/path/%.mp4: sims/$(S)/drm/%.pkl
	@ echo "Make: Path movie \`$@'"
	$(PATH_PROG) $<

# Rule to make a single path final-frame, given a SCENARIO
sims/$(S)/path/%-final.png: sims/$(S)/drm/%.pkl
	@ echo "Make: Path final-frame \`$@'"
	$(PATH_PROG_FINAL) $<

# Rule to make a single timeline, given a SCENARIO
sims/$(S)/path/%-obs-timeline.png: sims/$(S)/drm/%.pkl
	@ echo "Make: Timeline \`$@'"
	$(TIMELINE_PROG) -o sims/$(S)/path/$(*)-%s.%s Scripts/$(S).json $<

# Rule to make a single keepout map, given a SCENARIO
sims/$(S)/path/%-keepout-and-obs.png: sims/$(S)/drm/%.pkl
	@ echo "Make: Keepout \`$@'"
	$(KEEPOUT_PROG) -o sims/$(S)/path/$(*)-%s.%s Scripts/$(S).json $<

# Rule to make a group of movies, given a count.  The rule ends up looking like:
#  path-movie-N: sims/$(S)/path/SEED1.mp4 sims/SCRIPT/path/SEED2.mp4 ...
# where N is the number of movies (mp4's).  The number N is given to head
# as a ceiling on the number of path movie targets requested.  The .mp4 targets,
# in turn, are made by $(PATH_PROG) as shown above.
# The obs-timelines are not movies, but piggy-back on this same setup.
## Note: if sims/$(S)/drm/ does not exist, e.g. for an "experiment", we get
## an error that seems to be caused by syntactic mis-construction of the make 
## rule in this case.  The extra "path-dummy" target fixes this.
path-dummy:;

define MAKE_N_MOVIES
.PHONY: path-movie-$1 path-final-$1 obs-timeline-$1
path-movie-$1: path-dummy $(shell find sims/$(S)/drm -maxdepth 1 -name '*.pkl' 2>/dev/null | head --lines=$1 | sed -e 's:/drm/:/path/:' -e 's:.pkl:.mp4:')
	@ echo "Make: Placed movies in \`sims/$(S)/path'."
path-final-$1: path-dummy $(shell find sims/$(S)/drm -maxdepth 1 -name '*.pkl' 2>/dev/null | head --lines=$1 | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-final.png:')
	@ echo "Make: Placed final-frames in \`sims/$(S)/path'."
obs-timeline-$1: path-dummy $(shell find sims/$(S)/drm -maxdepth 1 -name '*.pkl' 2>/dev/null | head --lines=$1 | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-obs-timeline.png:')
	@ echo "Make: Placed obs-timeline in \`sims/$(S)/path'."
keepout-$1: path-dummy $(shell find sims/$(S)/drm -maxdepth 1 -name '*.pkl' 2>/dev/null | head --lines=$1 | sed -e 's:/drm/:/path/:' -e 's:\.pkl:-keepout-and-obs.png:')
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
html: script-exists sims/$(S)/html/index.html;

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

.PHONY: html-start html-stop html-status
html-start:
	util/html-serve.sh start
	@echo "Stop with \`make html-stop'."

html-stop:
	util/html-serve.sh stop

html-status:
	util/html-serve.sh status

########################################
## IPython parallel targets
##
## Note: these should be delegated to a shell script as for html-* targets

# create profile - only do this once
ipp-create: $(IPYDIR)/ipcluster_config.py;
$(IPYDIR)/ipcluster_config.py:
	@mkdir -p $(IPYDIR)
	@chmod g+r $(dir $(IPYDIR))
	ipython profile create --profile-dir=$(IPYDIR) --parallel
	@echo 'Installing engine path setup file.'
	cp Local/00-path.py $(IPYDIR)/startup
	cp Local/10-threads.py $(IPYDIR)/startup

# start a controller + engines
# The short sleep gives a chance for the cluster to start fully before returning.
# can optionally put a "nice" in front...but this currently (2/2018) fails because
# other long jobs are running un-niced.
ipp-start: $(IPYCLI);
$(IPYCLI): FORCE
	ipcluster start --profile-dir=$(IPYDIR) $(IPCLUSTER_OPT)
	@ sleep 3

# stop a controller + engines using the normal interface
ipp-stop:
	ipcluster stop --profile-dir=$(IPYDIR)

# stop a controller + engines
#  redirect to /dev/null needed in case there is no controller/engines
ipp-kill:
	@echo "Killing controller and engines spawned from $(IPYDIR)"
	-kill $$(ps uxww | grep $(IPYDIR) | grep -v grep | awk '{print $$2}') 2> /dev/null || true
	@sleep 1
	@echo "Done, confirm with \`make ipp-ps'".

# remove a controller
ipp-nuke: ipp-kill ipp-ps
	@echo "Removing controller associated with $(IPYDIR)"
	-rm -r $(IPYDIR)
	@echo "Removed.  Can re-create with \`make ipp-create'".

# ensure controller is running by connecting to the engines
ipp-status:
	@if [ ! -d $(IPYDIR) ]; then\
	  echo "No controller has been created at" $(IPYDIR); exit 1;\
	fi
	@if [ ! -f $(IPYCLI) ]; then\
	  echo "No controller is running in" $(IPYDIR); exit 1;\
	fi
	util/ipcluster-check.py --file $(IPYCLI)

# see how many engines are running
#   might be useful to see their load usage
ipp-ps:
	@echo "Found" `ps uxw | grep $(IPYDIR) | grep -v grep | grep engine | wc -l` "engines running."
	@ps uxw | grep $(IPYDIR) | grep -v grep | grep engine || true

# targets not corresponding to files
.PHONY: ipp-create ipp-start ipp-stop ipp-kill ipp-nuke ipp-status ipp-ps FORCE
