
## Obsolete
 
We stopped using iPython parallel due to fragility of the distributed
worker instances.

### iPython parallel support
 
As noted, simulations can also be run using iPython parallel, rather than
shell-level parallelism.
This mode starts up a given number of python processes ("engines"), which are
eventually given an Exosims function to run by `add-sims.sh`.
This creates extra state (a Python interpreter is held by each "engine"),
but avoids re-initialization of the Exosims object for each run.
See also the `SurveyEnsemble` documentation within Exosims.

To support creation, startup, and shutdown of these engines,
we added some iPython-parallel ("ipp") targets to the `Makefile`.
These targets operate independently of the simulation/results
infrastructure targets.
The targets (invoked like `make ipp-create`) are:

* `ipp-create`: create an ipython-parallel profile for this mission (use once
     per sandbox instance only).
     Copies several files into a per-user, per-machine ipyparallel directory.
     To undo, see `ipp-nuke`, below.
	 
* `ipp-start`: start the ipython-parallel controller and engines
     Note: If `EXOSIMS_ENGINE_N` (an integer) is exported from the environment,
     this many engines will be started up, otherwise, the system-default
     will be used.  To use this, run, from the shell:

       `$ EXOSIMS_ENGINE_N=8 make ipp-start`

* `ipp-stop`: stop the above processes.  See also `ipp-kill`, below.

* `ipp-status`: report status of the controller and engines
     Note: This attempts to run trivial jobs on the remote engines, so it
     will not work if the engines are busy with another job.  See `ipp-ps`, below.

* `ipp-ps`: use the unix `ps` command to identify the number of running
    engines. This works when the engines are busy, but is not informative
    about whether engines are responding to Python commands.

* `ipp-kill`: sometimes `ipp-stop` does not work and engines or controllers
     are orphaned.  `ipp-kill` identifies these by process id, and kills them.

* `ipp-nuke`: deletes your ipython-parallel profile.  The inverse of `ipp-create`.
     (Note: it attempts to `ipp-kill` first, so as to not leave engines running.)

As evidenced by the various status and termination commands, sometimes
using ipython parallel in this context can be annoying, because you have
to remember the state of the worker engines.
In particular, the engines will have to be restarted (`ipp-stop` followed by `ipp-start`)
when the underlying Exosims code changes, because the already-running
engines will hold a stale copy of the code.

