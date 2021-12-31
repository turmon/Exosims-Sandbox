#!/usr/bin/env python
#
# Check python and EXOSIMS readiness to run.
#
# Finds python and python-package version information,
# and optionally attempt a test instantiation and run_sim()
# of an EXOSIMS object.
#
# Usage:
#   ...
# Simplest usage:
#   python-version-check.py
#
# For more on usage, use the -h option.
# Some options may be described there but not documented here.
#

# author:
#  Michael Turmon, JPL, 2021
#

from __future__ import print_function
import argparse
import os
import time
import sys
import subprocess
import importlib
import warnings

PACKAGES = (
    'numpy',
    'matplotlib',
    'scipy',
    'scipy.optimize',
    'astropy',
    'h5py',
    'ipyparallel',
    'jplephem',
    )

# want MissionSim also because we use it in the test run
X_PACKAGES = (
    'EXOSIMS',
    'EXOSIMS.MissionSim',
    'EXOSIMS_local.IPClusterEnsembleJPL2',
    )

def verify_platform(args):
    print('Interpreter')
    # split+join removes embedded newlines
    print('{}python: {}'.format(args.prefix, ' '.join(sys.version.split())))
    print('{}interpreter path: {}'.format(args.prefix, sys.executable))
    return True

def verify_packages(args, pkgs):
    print('Packages')
    fails = 0
    for pkg in pkgs:
        with warnings.catch_warnings():
            # h5py has an annoying FutureWarning, not our focus here
            warnings.simplefilter("ignore")
            try:
                p = importlib.import_module(pkg)
            except:
                p = None
                # print('Failed to load "{}".'.format(pkg))
        if p and '__version__' in p.__dict__:
            version = p.__version__
        else:
            version = 'N/A'
        msg = '{}: {} (version {})'.format(pkg, ('ok' if p else 'FAILED'), version)
        print(args.prefix + msg)
        if not p: fails += 1
    return (fails == 0)

# get git commit info
def verify_commit(args):
    print('EXOSIMS Commit Information')
    import EXOSIMS
    src_dir = EXOSIMS.__path__[0]
    if not os.path.isdir(src_dir):
        print('{}EXOSIMS commit: could not find source directory {}'.format(args.prefix, src_dir))
        return False
    git_dir = src_dir + '/../.git'
    if not os.path.isdir(git_dir):
        print('{}EXOSIMS commit: could not access git info {}'.format(args.prefix, git_dir))
        return False
    cmd = 'git --git-dir={} log -1 | head -3'.format(git_dir)
    response = subprocess.check_output(cmd, shell=True)
    for s in str.split(response.decode(), '\n'):
        if not s: continue
        print('{}{}'.format(args.prefix, s))
    return True

def verify_script(args):
    r'''Script file is in args.script, or we use a default script file instead.'''
    # separator for object-related chatter
    sep_str = '#' * 12

    # Phase 1: Object initialization
    print('Object Initialization from Script')
    # we know these packages will load OK because the earlier step completed
    import EXOSIMS, EXOSIMS.MissionSim
    if not args.script:
        # scriptfile = 'EXOSIMS/EXOSIMS/Scripts/sampleScript_coron.json'
        scriptfile = os.path.join(EXOSIMS.__path__[0], 'Scripts', 'sampleScript_coron.json')
    else:
        scriptfile = args.script
    if not os.path.isfile(scriptfile):
        print('Script file "{}" was not found. Aborting.'.format(scriptfile))
        return False

    print('{} Object initialization follows...'.format(sep_str))
    print('{} [started at {}]'.format(sep_str, time.ctime()))
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile, cachedir=args.cachedir)
    print('{} Object initialization complete.'.format(sep_str))
    print('Object initialization: ok.')
    print('')

    # Phase 2: Test run (optional)
    if args.run:
        print('Script Test Run')
        print('{} Test Run follows...'.format(sep_str))
        print('{} [started at {}]'.format(sep_str, time.ctime()))
        sim.run_sim()
        print('{} Test Run complete.'.format(sep_str))
        print('Test Run: ok.')
        print('')
    return True
    


############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine.  Returns nonzero for trouble.'''
    if not verify_platform(args):
        return False
    if not verify_packages(args, PACKAGES):
        return False
    if args.exosims:
        if not verify_packages(args, X_PACKAGES):
            return False
        verify_commit(args) # ok for this to fail
    if args.init:
        if not verify_script(args):
            return False
    return True

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Check Python, Python package, and EXOSIMS version information. Optionally, instantiate an EXOSIMS object and run a default or specified script.",
                                     epilog='''Option strength: [-x] < [-i] < [-r], so giving -r will do everything. Typical fast usage: PYTHONPATH=EXOSIMS:Local util/python-version-check.py -x''')
    parser.add_argument('-x', '--exosims', action="store_true",
                            help='Check for EXOSIMS as well as python packages.')
    parser.add_argument('-i', '--init', action="store_true",
                            help='Try to initialize an EXOSIMS object with a basic script, "sampleScript_coron.json". Implies -x.')
    parser.add_argument('-r', '--run', action="store_true",
                            help='Run the simulation with run_sim() after initializing it. Implies -i and -x.')
    parser.add_argument('-s', '--script', type=str, default=None,
                            help='Script filename for object initialization and (optionally) run, to replace the default basic script.')
    parser.add_argument('-c', '--cachedir', type=str, default=None, 
                            help='Cache file directory to use instead of EXOSIMS default. Name any non-existing directory to force an EXOSIMS cache rebuild.')
    args = parser.parse_args()

    # run implies object init
    if args.run:
        args.init = True

    # init implies check-for-exosims
    if args.init:
        args.exosims = True

    # validate script filename
    if args.script:
        if not os.path.isfile(args.script):
            print('Given script filename "{}" does not exist. Exiting.'.format(args.script))
            sys.exit(2)

    # prefix for some stdout displays
    args.prefix = '    '

    ok = main(args)
    print('Overall: {}.'.format('OK' if ok else 'FAILED'))

    sys.exit(0 if ok else 1)


############################################################
#
# holding zone
#



'''
        # add in the SVN/Git revision^M
        path = os.path.split(inspect.getfile(self.__class__))[0]
        path = os.path.split(os.path.split(path)[0])[0]
        #handle case where EXOSIMS was imported from the working directory^M
        if path is '':
            path = os.getcwd()
        #comm = "git -C " + path + " log -1"^M
        comm = "git --git-dir=%s --work-tree=%s log -1"%(os.path.join(path,".git"),path)
        rev = subprocess.Popen(comm, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,shell=True)
        (gitRev, err) = rev.communicate()
        if sys.version_info[0] > 2:
            gitRev = gitRev.decode("utf-8")
        if isinstance(gitRev, basestring) & (len(gitRev) > 0):
            tmp = re.compile('\S*(commit [0-9a-fA-F]+)\n[\s\S]*Date: ([\S ]*)\n') \
                    .match(gitRev)
            if tmp:
                out['Revision'] = "Github " + tmp.groups()[0] + " " + tmp.groups()[1]

'''
