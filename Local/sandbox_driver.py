#!/usr/bin/env python

r"""
Top-level script for EXOSIMS runs in the Sandbox environment.

Simple Usage:
  sandbox_driver.py SCRIPT

Detailed Usage:
  sandbox_driver.py  [--outpath PATH] [--outopts OPTS] [--interactive] SCRIPT
where:
  positional arguments:
    SCRIPT               Path to a json script with EXOSIMS parameters.

  optional arguments:
    -h, --help            show this help message and exit
    --interactive         Interactive mode: post-init stdout not redirected to logfile
    --outpath PATH        Path to output directory.  Created if not present.
                          The default value (sims/basename(SCRIPT)) follows the 
                          Sandbox convention, so you don't usually need to specify it.
    --outopts OPTS        Results-file writing scheme. 
                          To write no files: --outopts ''
                          The default ('drm:pkl,spc:spc') follows the Sandbox convention,
                          so you don't usually need to specify it.
    --xspecs SCRIPT       an extra scenario-specific script loaded on top of the argument SCRIPT

Notes:  
  * The output OPTS control which result files are written.  It defaults to
    'drm:pkl,spc:spc', which writes the DRM to drm/SEED.pkl, and the star-
    planet configuration to spc/SEED.spc, where SEED is the large 
    integer random number seed. 
    In general, this is a comma-separated list of file-types-to-write.  If a 
    colon follows the file-type, the string following the colon is used instead 
    of the default '.pkl' extension.
    If the supplied extension ends in '.gz', the file is gzipped.  Thus, 
    --outopts 'drm:pkl.gz,spc:spc.gz' writes the gzipped pickled DRM to 
    drm/SEED.pkl.gz, and the gzipped spc dictionary to spc/SEED.spc.gz.
  * If an output directory is reused, new simulation output files will be added
    to that directory.

"""

# turmon 2017-2023

import subprocess
import argparse
import sys
import time
import socket
import json
import tempfile
import datetime
import shutil
import numpy as np
import numpy
import astropy
# imports needed by run_one, but not elsewhere in the file
import EXOSIMS
import EXOSIMS.MissionSim
import os
import os.path
import six.moves.cPickle # 'import as' is not allowed by ipyparallel
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
        print('%%% Running simulation...')
        SS.run_sim()
        DRM = SS.DRM[:]
        systems = SS.SimulatedUniverse.dump_systems()
        print('%%% Running simulation: done')
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
    
    print('%%% Preparing to pickle simulation...')

    # SPC object: save out a set of star/planet characteristics
    # turmon 2018/02: added coords so that star lat/lon can be recovered
    # turmon 2018/02: added Spec, star spectral type
    # turmon 2019/01: added MsTrue, stellar mass (!= MsEst)
    # turmon 2019/01: added SS.promoted_stars
    # turmon 2019/07: added TL.comp0; SU.e, SU.I; SS.t_char_earths
    # turmon 2020/09: added SS.known_{stars,rocky,earths} -- all in SS proto
    # turmon 2023/07: added TL.int_comp for EXOSIMS v3+, replacing comp0
    param_retain = [
        (SS.SimulatedUniverse,
             ('a', 'e', 'I', 'Rp', 'Mp', 'nPlans', 'p', 'plan2star', 's', 'sInds')),
        (SS.PlanetPopulation,
             ('SAG13coeffs', )),
        (SS.TargetList,
            ('L', 'dist', 'Name', 'coords', 'Spec', 'MsTrue', 'comp0', 'int_comp')),
        (SS,
             ('promoted_stars', 't_char_earths', 'known_stars', 'known_rocky', 'known_earths'))]
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
            six.moves.cPickle.dump(possible_outputs[outname], f)
        # change modtime of enclosing folder so that "make reduce" is not fooled
        # when the file contents (but not the filename) change
        # that is, writing the file may not alter the dir modtime
        now = datetime.datetime.now()
        epoch = now.timestamp()
        out_dir = os.path.dirname(path)
        try:
            os.utime(out_dir, (epoch, epoch))
        except PermissionError:
            # need try/except: utime fails if not dir owner (it happens!)
            # the below fallback creates an empty file in the directory
            tf = tempfile.NamedTemporaryFile(dir=out_dir, prefix='.', delete=True)
            tf.write(bytes('Dummy file\n', 'utf-8'))
            # will automatically delete the above file
            tf.close()

    # reset simulation object AFTER the above files are written
    # note, reset_sim() will pop the seed from SS.specs, which causes 
    # the seed to be regenerated next time
    # turmon 2023/11: no sense in doing this, we're always in standalone mode
    #SS.reset_sim(genNewPlanets=genNewPlanets, rewindPlanets=rewindPlanets)
    
    # caller will have the seed, but not the sim
    return seed

def sandbox_git_id():
    try:
        p = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True)
        commit_id = p.stdout.decode().strip()
    except:
        print('Could not get sandbox git information, continuing.', file=sys.stderr)
        commit_id = 'None'
    return commit_id


def main(args, xpsecs):
    r'''Prepare and execute simulations.'''

    # allow setting up an SS object for the run_one to hook into
    # this is only used for standalone mode
    global SS
    
    # ensure various output paths
    # turmon 01/2024: no longer use sys, adding
    outpath = args.outpath
    for d in ('log', 'log/environ', 'log/outspec', 'log/console', 'drm', 'spc'):
        outpath1 = os.path.join(outpath, d)
        if not os.path.exists(outpath1):
            try:
                os.makedirs(outpath1)
            except OSError:
                pass # usually a race condition with parallel jobs
        if not os.access(outpath1, os.W_OK | os.X_OK):
            raise ValueError("Cannot write to outpath `%s'." % outpath1)
    # various logs
    outpath_envi = os.path.join(outpath, 'log/environ') # S/W env.
    outpath_spec = os.path.join(outpath, 'log/outspec') # outspec
    outpath_log  = os.path.join(outpath, 'log/console')

    # recall xspecs = additional specs as distinct from script file
    ensemble_mode = []
    # ensemble_mode is picked up by the IPClusterEnsembleJPL __init__
    if args.standalone:
        ensemble_mode.append('standalone')
        # [11/2018] delay to avoid parallel race conditions (does help)
        time.sleep((os.getpid() % 40)/60.0)
    # ensemble_mode is picked up by the IPClusterEnsembleJPL __init__
    if args.interactive:
        ensemble_mode.append('interactive')
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
    if args.seed is not None and args.seed == 0:
        print('Seed given as 0. Cache warming only. Instantiation complete and caches in place.')
        return 'Cache-warming complete. No run performed.'
    res = sim.genOutSpec(tofile = os.path.join(outpath_spec, '%d.json' % seed))
    # place the output in a properly-named file
    if fp_log:
        fn_log = os.path.join(outpath_log, '%d_init.out' % seed)
        shutil.move(fp_log.name, fn_log)
        os.chmod(fn_log, 0o664) # ensure group-write

    printable_time = lambda tm: time.strftime("%Y-%m-%d %H:%M:%S", tm)
    subtime = time.localtime()
    print(f'Beginning run: {printable_time(subtime)}')
    # if standalone: manually set up the global SS variable for run_one to access
    # (if run under ipyparallel, this is done on the engines, via ipyparallel,
    # in the SurveyEnsemble __init__)
    if args.standalone:
        SS = sim.SurveySimulation
    kwargs = {'outpath': outpath, 'outopts': args.outopts} # kwargs are flowed to run_one
    # res is a list of seeds (a singleton in current usage)
    numruns = 1
    res = sim.run_ensemble(numruns, run_one=run_one, kwargs=kwargs)

    with open(os.path.join(outpath_envi, '%d.txt' % seed), 'w') as f:
        file_params = [
            ('file_type', 'EXOSIMS runtime environment'),
            ('user', os.getenv('USER')), 
            ('host', socket.gethostname()),
            ('time', printable_time(subtime)),
            ('shell_command', args.command),
            ('python_interpreter', ' '.join(sys.version.split())),
            ('python_venv', os.getenv('VIRTUAL_ENV', 'None')), 
            ('numpy_version', np.__version__),
            ('astropy_version', astropy.__version__),
            ('EXOSIMS_version', EXOSIMS.__version__),
            ('EXOSIMS_path', EXOSIMS.__path__[0]),
            ('sandbox_git_commit', sandbox_git_id()),
            ('seed', ','.join((str(r) for r in res)))]
        for name, value in file_params:
            f.write('%s: %s\n' % (name, value))
        #for r in res:
        #    f.write('%s\n' % str(r))
    # basic message
    return f'Run started {printable_time(subtime)}, complete {printable_time(time.localtime())}'


if __name__ == "__main__":
    # allow cut-and-paste reproducibility
    print('Invoked as:')
    print(' '.join(sys.argv))

    parser = argparse.ArgumentParser(description='Run an EXOSIMS in the Sandbox.')
    parser.add_argument('scriptfile', type=str, metavar='SCRIPT', help='Path to scriptfile.')
    # options
    parser.add_argument('--verbose', type=int, default=None, help='Verbosity (0/1).')
    parser.add_argument('--seed',    type=int, default=None, help='Random number seed (int); 0 for cache-warming only.')
    parser.add_argument('-q', '--quiet', default=False, action='store_true', help='Send object creation messages to a log file.')
    parser.add_argument('--interactive', default=False, action='store_true', help='Interactive mode (stdout not redirected to logfile).')
    parser.add_argument('--outpath', type=str, metavar='PATH', default=None,
            help='Path to output directory, created if not present. Default: basename of SCRIPT.')
    parser.add_argument('--outopts', type=str, metavar='OPTS', default='drm:pkl,spc:spc', help='Output result-file mode.')
    parser.add_argument('-x', '--xspecs', help='extra json spec-file to add on top of SCRIPT',
                      dest='xspecs', metavar='FILE', default='')

    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    # force standalone - it's no longer an option
    args.standalone = True

    ### 1: Process arguments
    # validate scriptfile
    if not os.path.isfile(args.scriptfile):
        raise ValueError("Supplied script `%s' not found" % args.scriptfile)
    # default output path
    #   derive sim results dir from script name: if it's in
    #   the "special" Scripts/... dir, just strip Scripts/ and .json out, 
    #   else use the basename less .json
    if args.outpath is None:
        if args.scriptfile.startswith('Scripts/'):
            # in Scripts/: use the tail
            pre_outpath = args.scriptfile.replace('Scripts/', 'sims/').replace('.json', '')
        else:
            # not in Scripts: use the basename only
            pre_outpath = os.path.splitext(os.path.basename(args.scriptfile))[0]
        args.outpath = os.path.join(os.path.abspath('.'), pre_outpath)
    # load extra specs, if any
    if args.xspecs and not args.xspecs.startswith('!'):
        assert os.path.isfile(args.xspecs), "%s is not a file." % args.xspecs
        try:
            script = open(args.xspecs).read()
            xspecs = json.loads(script)
        except ValueError as err:
            print("Error: Input xspec file `%s' improperly formatted." % (args.xspecs,))
            print("Error: JSON error was: ", err)
            # re-raise here to suppress the rest of the backtrace.
            # it is only confusing details about the bowels of json.loads()
            raise ValueError(err)
        except:
            print("Error: %s", (sys.exc_info()[0],))
            raise
    elif args.xspecs and args.xspecs.startswith('!'):
        xspecs = json.loads(args.xspecs[1:])
    else:
        xspecs = {}
        
    # save the command somewhere
    args.command = ' '.join(sys.argv)

    ### 2: Run EXOSIMS
    message = main(args, xspecs)
    print(f'{args.progname}: {message}')
    print(f'{args.progname}: Done.')
    sys.exit(0)
    
