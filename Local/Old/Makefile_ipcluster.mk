# Makefile fragment for ipcluster stuff (formerly part of root Makefile)
#
# (3) Ipython-parallel ("ipp")
#   Note: targets apply only to the machine where the "make" is run (e.g., aftac1).
#   ipp-create: create an ipython-parallel profile for this user (use just once).
#       To undo, see ipp-nuke, below.
#   ipp-start: start the ipython-parallel controller + engines
#       Note: If EXOSIMS_ENGINE_N is exported from the environment, this many
#       engines will be started, otherwise, the default is used. From shell:
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

# options for ipcluster startup
IPCLUSTER_OPT:=--daemonize --clean-logs=True
# allow number of engines to be set by environment
ifdef EXOSIMS_ENGINE_N
  IPCLUSTER_OPT:=-n $(EXOSIMS_ENGINE_N) $(IPCLUSTER_OPT)
endif

# ipython parallel setup directory
IPYDIR=ipyparallel/$(USER)/$(HOST_NAME)
IPYCLI=$(IPYDIR)/security/ipcontroller-client.json

########################################
## IPython parallel targets
##

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
