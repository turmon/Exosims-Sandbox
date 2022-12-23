#!/usr/bin/env python
#
# Check python and EXOSIMS readiness to run.
#
# Finds python and python-package version information,
# and optionally attempt a test instantiation and run_sim()
# of an EXOSIMS object.
#
# Usage:
#   python-version-check.py [-x] [-i] [-r] [-s SCRIPT] [-c CACHEDIR] [-R [sdet|conda|FILE|URL]]
#   
# Simplest usage:
#   python-version-check.py
#
# More complex usage:
#   python-version-check.py -R https://raw.githubusercontent.com/dsavransky/EXOSIMS/photometryUpdate/requirements.txt -x
#
# For more on usage, use the -h option.
# Some options are described there but not documented here.
#

# author:
#  Michael Turmon, JPL, 2021, 2022


from __future__ import print_function
import argparse
import os
import time
import sys
import subprocess
import shutil
import tempfile
import importlib
import pkg_resources
import warnings
import urllib.request


# dictionary mapping input "requirements flavor" keys to
# a list of requirements
# Note, this is almost in requirements.txt format as-is
PACKAGES = {
    'sdet': (
        'numpy',
        'matplotlib',
        'scipy',
        'scipy.optimize',
        'astropy',
        'h5py',
        'ipyparallel',
        'jplephem',
        ),
    'conda-v2': (
        'numpy>=1.2.0',
        'scipy>=1.7.2',
        'astropy>=4.3.1',
        'jplephem>=2.18',
        'ortools>=9.0',
        'h5py',
        'astroquery',
        'exodetbox', # NB: killed the dashes
        'tqdm',
        'pandas>=1.3',
        'MeanStars>=3.3.1',
        'synphot>=1.1.1',
        'ipyparallel>=8.0.0',
        )
    }

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


def get_package_list(args):
    r'''Use args.require to return either a package list, or a file-pointer to
    a requirements.txt file.'''
    req = args.require
    # try for old-style basic list first
    if req in PACKAGES:
        pack_list = [r.split('>=')[0] for r in PACKAGES[req]]
        return 'old', pack_list
    # otherwise, try for URL or file
    if req.startswith('http'):
        # there are other ways to do this that don't use a tempfile
        with urllib.request.urlopen(req) as response:
            with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                shutil.copyfileobj(response, tmp_file)
        pack_fp = open(tmp_file.name)
        args.tempfile = tmp_file.name
    else:
        pack_fp = open(req)
    return 'new', pack_fp
        

def verify_packages_new(args, pack_fp):
    r'''Packages are verified from a file pointer to a requirements.txt file.

    The file consists of lines like: numpy >= 1.2.0, and the pkg_resources
    checker is used to validate the version.'''
    # TODO: this was put together quickly, and might not be using pkg_resources
    # in the best-practices way
    print('Packages from <{}>'.format(args.require))
    fails = 0
    reqs = list(pkg_resources.parse_requirements(pack_fp))
    print('  Got {} requirements'.format(len(reqs)))
    for req in reqs:
        with warnings.catch_warnings():
            # h5py has an annoying FutureWarning, not our focus here
            warnings.simplefilter("ignore")
            try:
                p = pkg_resources.require(str(req))[0]
            except pkg_resources.VersionConflict:
                # print('Failed to load "{}".'.format(pkg))
                p = None
        if p is not None:
            version = str(p)
        else:
            version = 'N/A'
        msg = '{}: {} (version = {})'.format(str(req), ('ok' if p else 'FAILED'), version)
        print(args.prefix + msg)
        if not p: fails += 1
    return (fails == 0)

def verify_packages_old(args, pkgs):
    r'''Old-style package verification by iterating through a list of named modules.'''
    print('Package requirements from a stored literal')
    fails = 0
    for pkg in pkgs:
        with warnings.catch_warnings():
            # h5py has an annoying FutureWarning, not our focus here
            warnings.simplefilter("ignore")
            try:
                p = importlib.import_module(pkg)
            except:
                # print('Failed to load "{}".'.format(pkg))
                p = None
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
    
def verify_packages(args, method, pkg_list):
    r'''Return True if the pkg_list is satisfied with loadable modules.

    Switching on "method", which is either old or new.'''
    if method == 'old':
        return verify_packages_old(args, pkg_list)
    else:
        return verify_packages_new(args, pkg_list)


############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine.  Returns nonzero for trouble.'''
    if not verify_platform(args):
        return False
    # *get_package_list is a 2-tuple of a method and a requirements list
    if not verify_packages(args, *get_package_list(args)):
        return False
    if args.exosims:
        if not verify_packages(args, 'old', X_PACKAGES):
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
    parser.add_argument('-R', '--requirements', type=str, default='sdet', dest='require', metavar="[sdet|conda|FILE|URL]",
                            help='Requirements to check (literal abbreviations "sdet" or "conda-v2", or requirements.txt-format FILE, or URL of same, beginning with http[s]://).')
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

    # tempfile slot in case of -R http://...
    args.tempfile = ''

    ok = main(args)

    # delete the tempfile we might have needed for -R
    if args.tempfile and os.path.isfile(args.tempfile):
        try:
            os.remove(args.tempfile)
        except FileNotFoundError:
            pass

    print('Overall: {}.'.format('OK' if ok else 'FAILED'))
    sys.exit(0 if ok else 1)


