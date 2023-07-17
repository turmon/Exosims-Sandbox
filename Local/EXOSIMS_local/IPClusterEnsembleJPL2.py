from __future__ import print_function

from ipyparallel import Client
from EXOSIMS.Prototypes.SurveyEnsemble import SurveyEnsemble 
from EXOSIMS.util.get_module import get_module
import time
from IPython.core.display import clear_output
import sys
import os


class RedirectStdStreams(object):
    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

        
class IPClusterEnsembleJPL2(SurveyEnsemble):
    """Parallelized survey ensemble based on IPython parallel (ipcluster)
    
    """

    def __init__(self, ensemble_controller=None, ensemble_mode='', **specs):
        
        SurveyEnsemble.__init__(self, **specs)

        # set up attributes
        self.verb = specs.get('verbose', True)
        self.standalone = False
        self.interactive = False

        # allow bail-out
        if 'init-only' in ensemble_mode:
            self.vprint("SurveyEnsemble: initialize-only mode")
            bailout = True
        if 'standalone' in ensemble_mode:
            self.vprint("SurveyEnsemble: standalone mode: no ipyparallel")
            self.standalone = True
            bailout = True
        if 'interactive' in ensemble_mode:
            self.vprint("SurveyEnsemble: interactive mode: stdout not redirected to logfile")
            self.interactive = True
            bailout = True

        # No ipyparallel mechanics are needed
        if bailout:
            return

        # specify the cluster
        if ensemble_controller:
            if '.json' in ensemble_controller:
                arglist = dict(url_file=ensemble_controller)
            else:
                arglist = dict(profile=ensemble_controller)
        else:
            arglist = dict()
        # access the cluster
        self.rc = Client(**arglist)
        self.dview = self.rc[:]
        self.dview.block = True
        # these are the imports needed by the run_one()
        # (cannot use import ... as ... for this purpose)
        with self.dview.sync_imports():
            import EXOSIMS, EXOSIMS.util.get_module, EXOSIMS_local, \
                time, os, os.path, random, numpy, cPickle, gzip, traceback
        if 'logger' in specs:
            specs.pop('logger')
        # pop the seed from the specs to force re-seeding
        if 'seed' in specs:
            specs.pop('seed')
        # push the specs to the engines
        self.dview.push(dict(specs=specs))
        # instantiate a SurveySimulation in the global workspace on each engine
        res = self.dview.execute(
            "SS = EXOSIMS.util.get_module.get_module_from_specs" +
            "(specs, 'SurveySimulation')(**specs)"
            )
        self.vprint("Created SurveySimulation objects on %d engines."%len(self.rc.ids))
        # pull the seeds from each engine
        seeds = self.dview.pull('SS.seed', block=True)
        # print stdout/stderr of each engine's activity - this is likely to be captured
        # in the invoking function.  Note, we don't have access to the parent SS.seed here.
        if True:
            for row,erow,id,seed in zip(res.stdout, res.stderr, res.engine_id, seeds):
                print('==== Engine = %d, Seed = %d ====' % (id, seed))
                if erow:
                    msg = ''.join(['[#%d] Error: %s\n' % (id, line) for line in erow.split('\n') if line])
                    print(msg)
                    sys.stderr.write(msg)
                print(''.join(['[#%d] %s\n' % (id, line) for line in row.split('\n')]))
        # we will use the load-balanced view for cluster Exosims runs
        self.lview = self.rc.load_balanced_view()

    def run_ensemble(self, sim, nb_run_sim, run_one=None, genNewPlanets=True,
            rewindPlanets=True, kwargs={}):
        
        if self.standalone:
            return self.run_ensemble_stand(sim, nb_run_sim, run_one, genNewPlanets, rewindPlanets, kwargs)
        else:
            return self.run_ensemble_ipp(  sim, nb_run_sim, run_one, genNewPlanets, rewindPlanets, kwargs)

    def run_ensemble_stand(self, sim, nb_run_sim, run_one, genNewPlanets=True, rewindPlanets=True, kwargs={}):
        r'''Stand-alone simulation runner.'''
        t1 = time.time()
        res = []
        for j in range(nb_run_sim):
            if nb_run_sim > 1:
                print('Survey simulation: %s/%s' % (j + 1, int(nb_run_sim)))
            seed = sim.seed
            fn = os.path.join(kwargs['outpath'], 'log', 'log-%d.out' % (seed,))
            if self.interactive:
                logfile = None # e.g., allow breakpoint() within the run_one
            else:
                logfile = open(fn, 'w')
            with RedirectStdStreams(stdout=logfile):
                ar = run_one(genNewPlanets=genNewPlanets, rewindPlanets=rewindPlanets, **kwargs)
            res.append(ar)
        t2 = time.time()
        self.vprint("Completed %s simulation(s) in %d sec"%(int(nb_run_sim), t2 - t1))
        return res


    def run_ensemble_ipp(self, sim, nb_run_sim, run_one=None, genNewPlanets=True,
            rewindPlanets=True, kwargs={}):
        
        if not run_one:
            raise ValueError('Require a run_one function to be provided')
        t1 = time.time()
        async_res = []
        for j in range(nb_run_sim):
            ar = self.lview.apply_async(run_one, genNewPlanets=genNewPlanets,
                    rewindPlanets=rewindPlanets, **kwargs)
            async_res.append(ar)
        print("Submitted %d tasks." % len(async_res))
        ar = self.rc._asyncresult_from_jobs(async_res)
        # ad hoc status-reporting
        progress = 0
        while not ar.ready():
            ar.wait(10.)
            clear_output(wait=True)
            if ar.progress == 0:
                forecast = 'not yet able to forecast time remaining.'
            elif ar.progress > progress:
                # update forecast right after we learn more about job-completion rate,
                # otherwise, the accuracy of the rate is diminished
                progress = ar.progress
                timeleft = ar.elapsed/ar.progress * (nb_run_sim - ar.progress)
                if timeleft > 3600.:
                    timeleftstr = "%2.2f hours" % (timeleft/3600.)
                elif timeleft > 60.:
                    timeleftstr = "%2.2f minutes" % (timeleft/60.)
                else:
                    timeleftstr = "%2.2f seconds" % timeleft
                forecast = 'about ' + timeleftstr + ' to go.'

            print("%4i/%i tasks finished after %4i s -- %s" % (ar.progress, nb_run_sim, ar.elapsed, forecast), end="")
            sys.stdout.flush()
        #self.rc.wait(async_res)
        #self.rc.wait_interactive(async_res)
        t2 = time.time()
        print("\nCompleted in %d sec" % (t2 - t1))
        # output the ipp engine stdout's to log-files
        for j, ar1 in enumerate(async_res):
            # retrieve result - just the seed, actually
            seed1 = ar1.get()
            fn = os.path.join(kwargs['outpath'], 'log', 'log-%s.out' % seed1)
            with open(fn, 'w') as fp:
                for line in ar1.stdout:
                    fp.write(line)
            if ar1.stderr:
                fn = os.path.join(kwargs['outpath'], 'log', 'log-%s.err' % seed1)
                with open(fn, 'w') as fp:
                    for line in ar1.stderr:
                        fp.write(line)
        # return the list of seeds
        return [ar.get() for ar in async_res]
