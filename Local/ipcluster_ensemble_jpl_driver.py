#!/usr/bin/env python

r"""
Top-level script for EXOSIMS runs in the Sandbox environment.

Simple Usage:
  ipcluster_ensemble_jpl_driver.py [--outpath PATH] SCRIPT N_RUNS
  ipcluster_ensemble_jpl_driver.py --help 

Detailed Usage:
  ipcluster_ensemble_jpl_driver.py  [--outpath PATH] [--outopts OPTS]
                               [--controller CONTROLLER] [--email EMAIL]
                               [--toemail TOEMAIL]
                               SCRIPT N_RUNS
where:
  positional arguments:
    SCRIPT               Path to a json script with EXOSIMS parameters.
    N_RUNS               Number of ensemble runs (integer).

  optional arguments:
    -h, --help            show this help message and exit
    --outpath PATH        Path to output directory.  Created if not present.
                          Default: basename of scriptfile.
    --outopts OPTS        Output result-file mode. Default: 'drm'.  See below.
    --controller CONTROLLER
                          Controller name (as used in ipcluster) for ipython
                          engines.
    --xspecs SCRIPT       an extra scenario-specific script loaded on top of the argument SCRIPT
    --email EMAIL         Email address to notify when run is complete.a
    --toemail TOEMAIL     Additional email to notify when run is complete.

Notes:  
  * An ipcluster instance must be running and accessible in order
    to use this script.  If everything is already set up properly,
    this is usually a matter of executing, from the shell, a command like:
       $ ipcluster start
    In the Sandbox, this is:
       $ make ipp-start
  * The output OPTS control which result files are written.  It defaults to
    'drm', which writes the DRM to drm/SEED.pkl, where SEED is the large 
    integer random number seed. 
    In general, this is a comma-separated list of file-types-to-write.  If a 
    colon follows the file-type, the string following the colon is used instead 
    of the default '.pkl'.
    If the supplied extension ends in '.gz', it is gzipped.  Thus, 
    --outopts 'drm:pkl.gz,spc:spc.gz' writes the gzipped pickled DRM to 
    drm/SEED.pkl.gz, and the gzipped spc dictionary to spc/SEED.spc.gz.
  * The --email option sends to TO_EMAIL using the localhost mail server.
    Or, it will authenticate with Gmail's server if TO_EMAIL contains "@gmail".
  * If an output directory is reused, new simulation output files will be added.  
    Any generated errors will be appended to any existing log.err file.

"""

# turmon 2017, starting from a version by cdx, 2016-2017

import numpy as np
import numpy
import argparse
import sys
import time
import socket
import json
import smtplib
import getpass
import tempfile
import shutil
# imports needed by run_one, but not elsewhere in the file
import EXOSIMS
import EXOSIMS.MissionSim
import os
import os.path
import cPickle # 'import as' is not allowed by ipyparallel
import gzip
import traceback


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


def run_one(genNewPlanets=True, rewindPlanets=True, outpath='.', outopts='', reseedRunSim=False):
    r'''Worker function that is compatible with iPython parallel dispatch.

    The function signature here is different than the Prototype.
    The SS object is in the global scope.
    '''
    global SS

    def generate_output_path(outname, outpath, outopts, seed):
        r'''Helper function: return the file-name and file-writer to use for the given outname.

        Returns None, None if that particular outname is not wanted.'''
        # turn a string like 'drm:pkl.gz,spc' into ['drm:pkl.gz', 'spc']
        optlist = outopts.split(',')
        for o in optlist:
            if not o.startswith(outname):
                continue
            # it is a match: determine the filename and file-writer
            # file extension: .pkl, or what is after the : (e.g., drm:pkl.gz -> pkl.gz)
            ext = 'pkl' if ':' not in o else o.split(':')[1]
            # determine file-writer
            writer = gzip.open if ext.endswith('.gz') else open
            # determine path, e.g., 'outpath/drm/01234321.pkl.gz'
            path = os.path.join(outpath, outname, '%s.%s' % (str(seed), ext))
            return path, writer
        return None, None
                
    def ensure_permissions(fn):
        r'''Ensure correct permissions on the named data file.  
        We use rw-rw-r-- = 664 (octal), to allow group-write.'''
        try:
            os.chmod(fn, 0o664)
        except OSError:
            pass # e.g., don't own the file

    # wrap the run_sim in a try/except loop
    seed = SS.seed
    try:
        # run one survey simulation
        print '%%% Running simulation...'
        SS.run_sim()
        DRM = SS.DRM[:]
        systems = SS.SimulatedUniverse.dump_systems()
        print '%%% Running simulation: done'
    except Exception as e:
        # if anything goes wrong, log the error to save it persistently
        fn = os.path.join(outpath, 'log', 'log-%s.err' % str(seed))
        with open(fn, 'w') as f:
            f.write('===== Exception raised: %s\n' % time.ctime())
            f.write(repr(e))
            f.write('\n')
            f.write(traceback.format_exc())
            f.write('\n\n')
        ensure_permissions(fn)
        raise # re-raise the exception
    
    print '%%% Preparing to pickle simulation...'

    # SPC object: save out a set of star/planet characteristics
    # turmon 2018/02: added coords so that star lat/lon can be recovered
    # turmon 2018/02: added Spec, star spectral type
    # turmon 2019/01: added MsTrue, stellar mass (!= MsEst)
    # turmon 2019/01: added SS.promoted_stars
    # turmon 2019/07: added TL.comp0; SU.e, SU.I; SS.t_char_earths
    param_retain = [
        (SS.SimulatedUniverse,
             ('a', 'e', 'I', 'Rp', 'Mp', 'nPlans', 'p', 'plan2star', 's', 'sInds')),
        (SS.PlanetPopulation,
             ('SAG13coeffs', )),
        (SS.TargetList,
            ('L', 'dist', 'Name', 'coords', 'Spec', 'MsTrue', 'comp0')),
        (SS,
             ('promoted_stars', 't_char_earths'))]
    spc_params = {}
    for (module, fieldlist) in param_retain:
        for field in fieldlist:
            try:
                spc_params[field] = getattr(module, field)
            except AttributeError:
                pass # allow the field to be absent

    # omnibus dict with multiple keys
    combo_output = dict(DRM=DRM, sys=systems)
    # all possible outputs - we may write any of these
    possible_outputs = dict(drm=DRM, spc=spc_params, sys=systems, combo=combo_output)

    # this could be generalized further:
    #    allow use of a json-dumper as well as a pickler
    #    allow use of filenames besides <seed>.extension
    for outname in possible_outputs:
        path, writer = generate_output_path(outname, outpath, outopts, seed)
        if not path: continue # no output of that type
        #raise ValueError(' '.join((outname,outpath,outopts,path)))
        with writer(path, 'wb') as f:
            cPickle.dump(possible_outputs[outname], f)

    # reset simulation object AFTER the above files are written
    # note, reset_sim() will pop the seed from SS.specs, which causes 
    # the seed to be regenerated next time
    SS.reset_sim(genNewPlanets=genNewPlanets, rewindPlanets=rewindPlanets)
    
    # caller will have the seed, but not the sim
    return seed


def ensure_email(args):
    r'''Ensure email can be accessed, fetching password if needed.

    If args.passwd is None, no password will be asked-for, otherwise, it
    is retrieved and returned via args.passwd.
    '''
    if not args.email:
        return
    if args.passwd is not None:
        args.passwd = getpass.getpass("Password for %s:\n" % args.email)
    # raise an error now in case of a bad e-mail setup
    server = smtplib.SMTP(args.smtp)
    server.ehlo()
    if args.passwd is not None:
        server.starttls()
        server.login(args.email, args.passwd)
    server.quit()

def send_email(args, message):
    r'''Send completion email, if desired.'''
    if not args.email:
        return
    server = smtplib.SMTP(args.smtp)
    server.ehlo()
    if args.passwd is not None:
        server.starttls()
        server.login(args.email, args.passwd)

    if args.toemail is not None:
        to_email = [args.email, args.toemail]
    else:
        to_email = [args.email]

    msg = "\r\n".join([
      "From: %s" % args.email,
      "To: %s" % ";".join(to_email),
      "Subject: Run Completed",
      "",
      "Results for: %s\nStored in: %s\nTiming: %s\n\nCome see what I've done." % (
          args.scriptfile, args.outpath, message)
      ])
    server.sendmail(args.email, to_email, msg)
    server.quit()


def main(args, xpsecs):
    r'''Prepare and execute simulations.'''

    # allow setting up an SS object for the run_one to hook into
    # this is only used for standalone mode
    global SS
    
    # ensure various output paths
    outpath = args.outpath
    for d in ('run', 'log', 'drm', 'sys', 'spc'):
        outpath1 = os.path.join(outpath, d)
        if not os.path.exists(outpath1):
            print "Creating output path `%s'." % outpath1
            try:
                os.makedirs(outpath1)
            except OSError:
                pass # usually a race condition with parallel jobs
        if not os.access(outpath1, os.W_OK | os.X_OK):
            raise ValueError("Cannot write to outpath `%s'." % outpath1)
    # the run summary (outspec, seeds); the logs
    outpath_run = os.path.join(outpath, 'run')
    outpath_log = os.path.join(outpath, 'log')

    # recall xspecs = additional specs as distinct from script file
    ensemble_mode = []
    # ensemble_controller is picked up by the IPClusterEnsembleJPL __init__
    if args.controller:
        xspecs['ensemble_controller'] = args.controller
    # ensemble_mode is picked up by the IPClusterEnsembleJPL __init__
    if args.standalone:
        ensemble_mode.append('standalone')
        # [11/2018] delay to avoid parallel race conditions (does help)
        time.sleep((os.getpid() % 40)/40.0)
    if args.verbose is not None:
        if args.verbose == 0:
            xspecs['verbose'] = False
        else:
            xspecs['verbose'] = True
    # place an explicit seed, if given, in xspecs
    if args.seed is not None:
        xspecs['seed'] = int(args.seed)
    # don't put the key in if it will be empty
    if len(ensemble_mode) > 0:
        xspecs['ensemble_mode'] = ','.join(ensemble_mode)

    # if quiet-mode, put output from object init in a tempfile until we know the seed
    fp_log = tempfile.NamedTemporaryFile(delete=False) if args.quiet else None
        
    # set up a simulation object
    # - for standalone mode, this is the sim we will actually use
    # - for ipyparallel mode, this forces precomputed cache-files to be written, 
    #   but this sim is not the one active on the engines.
    with RedirectStdStreams(stdout=fp_log):
        sim = EXOSIMS.MissionSim.MissionSim(args.scriptfile, **xspecs)
    seed = sim.SurveySimulation.seed
    res = sim.genOutSpec(tofile = os.path.join(outpath_run, 'outspec_%d.json' % seed))
    # place the output in a properly-named file
    if fp_log:
        fn_log = os.path.join(outpath_log, 'log-%d.init' % seed)
        shutil.move(fp_log.name, fn_log)
        os.chmod(fn_log, 0o664) # ensure group-write

    subtime = time.ctime()
    print "Beginning run on: %s" % subtime
    # if standalone: manually set up the global SS variable for run_one to access
    # (if run under ipyparallel, this is done on the engines, via ipyparallel,
    # in the SurveyEnsemble __init__)
    if args.standalone:
        SS = sim.SurveySimulation
    kwargs = {'outpath': outpath, 'outopts': args.outopts} # kwargs are flowed to run_one
    res = sim.run_ensemble(int(args.numruns), run_one=run_one, kwargs=kwargs)

    # result is a list of seeds
    with open(os.path.join(outpath_run, 'outseed_%d.txt' % seed), 'w') as f:
        f.write('# EXOSIMS run list\n# user: %s \n# host: %s\n# time: %s\n# command: %s\n' %
                    (os.getenv('USER'), socket.gethostname(), subtime, args.command))
        for r in res:
            f.write('%s\n' % str(r))

    # basic message
    return 'Run started on %s, complete on %s.' % (subtime, time.ctime())


if __name__ == "__main__":
    # allow cut-and-paste reproducibility
    print 'Invoked as:'
    print ' '.join(sys.argv)

    parser = argparse.ArgumentParser(description='Run an EXOSIMS ensemble under ipyparallel.')
    parser.add_argument('scriptfile', type=str, metavar='SCRIPT', help='Path to scriptfile.')
    parser.add_argument('numruns', type=int, metavar='N_RUNS', help='Number of ensemble runs.')
    # options
    parser.add_argument('--verbose', type=int, default=None, help='Verbosity (0/1).')
    parser.add_argument('--seed',    type=int, default=None, help='Random number seed (int).')
    parser.add_argument('-q', '--quiet', default=False, action='store_true', help='Send object creation messages to a log file.')
    parser.add_argument('--standalone', default=False, action='store_true', help='Stand-alone mode (no ipyparallel).')
    parser.add_argument('--outpath', type=str, metavar='PATH',
            help='Path to output directory, created if not present. Default: basename of SCRIPT.')
    parser.add_argument('--outopts', type=str, metavar='OPTS', default='drm,spc',
            help='Output result-file mode.')
    parser.add_argument('-x', '--xspecs', help='extra json spec-file to add on top of SCRIPT',
                      dest='xspecs', metavar='FILE', default='')
    parser.add_argument('--controller', type=str, help='Controller name (as used in ipcluster) for ipython engines.')
    parser.add_argument('--email', type=str, help='Email address to notify when run is complete.')
    parser.add_argument('--toemail', type=str, help='Additional email to notify when run is complete.')

    args = parser.parse_args()

    ### 1: Process arguments
    # validate scriptfile
    if not os.path.isfile(args.scriptfile):
        raise ValueError("Supplied script `%s' not found" % args.scriptfile)
    # default output path
    if args.outpath is None:
        args.outpath = os.path.join(os.path.abspath('.'),
                                    os.path.splitext(os.path.basename(args.scriptfile))[0])
    # load extra specs, if any
    if args.xspecs and not args.xspecs.startswith('!'):
        assert os.path.isfile(args.xspecs), "%s is not a file." % args.xspecs
        try:
            script = open(args.xspecs).read()
            xspecs = json.loads(script)
        except ValueError as err:
            print "Error: Input xspec file `%s' improperly formatted." % (args.xspecs,)
            print "Error: JSON error was: ", err
            # re-raise here to suppress the rest of the backtrace.
            # it is only confusing details about the bowels of json.loads()
            raise ValueError(err)
        except:
            print "Error: %s", (sys.exc_info()[0],)
            raise
    elif args.xspecs and args.xspecs.startswith('!'):
        xspecs = json.loads(args.xspecs[1:])
    else:
        xspecs = {}
        
    # save the command somewhere
    args.command = ' '.join(sys.argv)
    if args.email and '@gmail' in args.email:
        args.smtp = 'smtp.gmail.com:587'
        args.passwd = 'TBD' # we will prompt for password and fill this in
    else:
        args.smtp = 'localhost'
        args.passwd = None # signals no password needed

    ### 2: Run EXOSIMS
    ensure_email(args)
    message = main(args, xspecs)
    send_email(args, message)

    sys.exit(0)
    
