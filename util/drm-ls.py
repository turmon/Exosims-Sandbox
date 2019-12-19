#!/usr/bin/env python
#
# Pickles containing DRM results are summarized to the terminal in a format
# similar to "ls".
#
# Usage:
#   drm-ls.py [-lqrc] [-d N] [FILE_OR_DIRECTORY ...]
# Simplest usage:
#   drm-ls.py 
#
# where:
#  -l gives long-format output (extra columns)
#  -q gives short, summary output (no per-DRM output)
#  -r performs recursive descent (otherwise, named files/dirs are examined)
#  -d N limits recursive descent to N levels
#  -c gives CSV output instead of tabular output
#
# A file, a list of files, a directory, or list thereof, can be given
# for listing.  By default, DRMs in the current working directory are
# listed.
#
# DRM files should be named with extensions .pkl, .drm, .gz, or .bz2, and
# should correspond to Python pickles.  Gzip or bzip2 compressed files
# (bearing any extension, even .pkl) will be uncompressed and examined.
#
# For more on usage, use the -h option.
# Some options may be described there but not documented here.
#

# author:
#  Michael Turmon, JPL, 2017
#
# robustness note: the signature for a pickle in MAGIC_NUMBERS may be
# sensitive to the pickle encoding used, e.g. pickle version, or
# python 2.x or 3.x.  Extra magic numbers for these other pickles
# could be added.


import argparse
import sys
import os
import os.path
#import pickle
import cPickle as pickle
from collections import defaultdict
import bz2, gzip
import numpy as np
import astropy.units as u
#from astropy.time import Time


# magic numbers at the start of compressed files
MAGIC_NUMBERS = [
    ("\x1f\x8b\x08", gzip.GzipFile),
    ("\x42\x5a\x68", bz2.BZ2File),
    ("(lp", open), # "list begin" -- works for python 2.x DRMs, even empty ones
    ]
MAGIC_MAXLEN = max(len(m[0]) for m in MAGIC_NUMBERS)

# global modes
CSV_OUTPUT = False

class UnknownFileException(Exception):
    r'''Raised when encountering a file that does not appear to be a DRM.'''
    pass


############################################################
#
# Utilities
#
############################################################


def args_to_filenames(arglist, recurse, max_depth):
    # traverse root directory, and list directories as dirs and files as files
    result = []
    for arg in arglist:
        if os.path.isdir(arg):
            for root, dirs, files in os.walk(arg):
                path = root.split(os.sep)
                if (not recurse) or (len(path) >= max_depth):
                    del dirs[0:] # modify dirs in place => do not descend farther
                result.extend(drm_filter([os.path.join(root, f) for f in sorted(files)]))
        else:
            result.extend(drm_filter([arg]))
    return result


def drm_filter(files):
    r'''Filter file list to include only readable files, with DRM-like extensions.'''
    # don't peek into the file yet, it takes too long
    ok = []
    for f in files:
        if os.path.isfile(f) and os.path.splitext(f)[1] in ('.pkl', '.gz', '.bz2', '.drm'):
            if f.startswith('./'):
                ok.append(f[2:])
            else:
                ok.append(f)
    return ok


def file_accessor(fn):
    r'''Return an accessor for the filename, depending on the compression type, if any.'''
    with open(fn) as f:
        try:
            file_start = f.read(MAGIC_MAXLEN)
        except EOFError:
            file_start = ''
    for magic, accessor in MAGIC_NUMBERS:
        if file_start.startswith(magic):
            return accessor(fn)
    raise UnknownFileException("File %s not a DRM" % fn)


def load_drm(fn):
    r"""Return list of observations made, loaded from an external DRM pickle."""
    global DIAGNOSE

    def why_skip(fn, reason, diagnose):
        r'''Supply uniform message upon skipping a file.'''
        if diagnose:
            print 'Skipped %s (%s)' % (fn, reason)

    try:
        drm = pickle.load(file_accessor(fn))
    except (pickle.UnpicklingError, EOFError, UnknownFileException):
        # vanilla textfiles should route through this path
        # some pickles that are not DRMs will route through this path
        # Note: possibly, other errors should be added to the tuple above
        why_skip(fn, 'could not unpack as a DRM', DIAGNOSE)
        return None
    except:
        # an unqualified try/except can mask errors, so raise it
        why_skip(fn, 'Unpacking Error!', True)
        raise
    # we read something: examine it to verify it is a DRM.
    # we take a valid DRM to be either:
    #   (1) an empty list (could be a non-DRM, but tough luck)
    #   (2) a non-empty list containing dictionaries with
    #       detection statuses
    # non-DRMs can include any un-pickled thing, including
    # dictionaries, etc.
    if isinstance(drm, list):
        if len(drm) == 0:
            return drm # case 1 above
        if len(drm) > 0 and isinstance(drm[0], dict) and ('star_ind' in drm[0]):
            return drm # case 2 above
    why_skip(fn, 'unpacked, but not a DRM', DIAGNOSE)
    return None

 
def detail_summarize(drm):
    r'''More detailed DRM summary.'''
    # successful detections - an array entry for each observation
    dets = np.array([np.sum(np.maximum(0, obs['det_status'])) if 'det_status' in obs else 0
                     for obs in drm])
    n_det = np.sum(dets)
    # successful characterizations - an array entry for each observation
    chars = np.array([np.sum(np.maximum(0, obs['char_status'])) if 'char_status' in obs else 0
                      for obs in drm ])
    n_char = np.sum(chars)
    # final mass
    masses = [obs['scMass'] for obs in drm if 'scMass' in obs]
    mass = masses[-1].value if masses else -99
    # final arrival time
    arr_times = [obs['arrival_time'] for obs in drm if 'arrival_time' in obs]
    arr_time = strip_units(arr_times[-1]) if arr_times else -99
    # stars we visited, in order
    all_stars = np.array([obs['star_ind'] for obs in drm])
    all_star_det = set(all_stars[dets > 0]) # de-duplicate
    n_star_det = len(all_star_det)
    # package and return
    result = dict(Nobs=len(drm),Ndet=n_det,Nchar=n_char,Nstar_det=n_star_det,Fmass=mass, ArrTime=arr_time)
    return result


def key_print_order(s):
    default_key_order = ['Nobs', 'Ndet', 'Nchar', 'Nstar_det', 'Fmass','ArrTime']
    ordered_keys = []
    # take keys in defined order first
    for k in default_key_order:
        if k in s:
            ordered_keys.append(k)
    # next, add in anything else
    for k in sorted(s):
        if k not in ordered_keys:
            ordered_keys.append(k)
    return ordered_keys

def strip_units(x):
    r'''Remove astropy units from x, if present.'''
    try:
        return x.value
    except AttributeError:
        return x

def summarize_print(fn, s, verbosity, fmt='%.0f', first_call=[]):
    # fmt0: file line; delim1: inter-number delimiter
    if CSV_OUTPUT:
        fmt0 = '{0},'.format
        delim1 = ','
    else:
        fmt0 = '{0:<31}\t'.format
        delim1 = '\t'
    # print header, if desired - uses mutable argument to detect first call
    if len(first_call) is 0:
        first_call.append('called') # header is done
        line = delim1.join(['%s' % (k,) for k in key_print_order(s)])
        if verbosity:
            print fmt0('DRM') + line
    if verbosity > 1:
        line = delim1.join([fmt % (strip_units(s[k]),) for k in key_print_order(s)])
    elif verbosity == 1:
        line = fmt % (s['Nobs'],)
    print fmt0(fn) + line


def summarize(fn, drm, accum, verbosity):
    r'''Summarize each DRM, add up the total valid DRMs'''
    # number of observations
    if verbosity > 1:
        # inspect the drm
        s = detail_summarize(drm)
    else:
        # don't look within the drm
        s = dict(Nobs=len(drm))
    if verbosity:
        summarize_print(fn, s, verbosity)
    # update accumulator
    for k in s:
        accum[k] += s[k]
        
    
def summ_accum(Ndrm, accum, verbosity):
    r'''Summarize the accumulated counts.'''
    #print 'TOTAL NDRM =', accum['n_drm']
    #print 'TOTAL NOBS =', accum['n_obs']
    verbosity_used = max(verbosity, 1)
    if Ndrm > 0: 
        summarize_print('*TOTAL', accum, verbosity_used, first_call=['called'])
    if Ndrm > 0:
        mean = {k: accum[k]/(1.0*Ndrm) for k in accum}
        summarize_print('*MEAN', mean, verbosity_used, fmt='%.2f', first_call=['called'])
    # this is not handled per-DRM, because it makes no sense there
    if not CSV_OUTPUT:
        print '%d DRMs examined' % Ndrm
    

############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine: load and process DRMs.'''

    global CSV_OUTPUT
    CSV_OUTPUT = args.csv_output
    # this DIAGNOSE is intended to be for error-reporting to diagnose
    # issues with the code or the files, but not for ordinary use
    global DIAGNOSE
    DIAGNOSE = args.DIAGNOSE

    # Get filenames of DRMs
    fns = args_to_filenames(args.drm, args.recurse, args.max_depth)

    # Summarize each DRM
    Ndrm = 0
    accum = defaultdict(int)
    for fn in fns:
        drm = load_drm(fn)
        if drm is None:
            continue # not an actual DRM file, eg a textfile
        summarize(fn, drm, accum, args.verbose)
        Ndrm += 1

    # Produce cumulative summary across all DRMs
    summ_accum(Ndrm, accum, args.verbose)


    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Summarize EXOSIMS DRM(s).",
                                     epilog='')
    parser.add_argument('drm', metavar='DRM', nargs='*', default='.',
                            help='drm file, or directory thereof')
    parser.add_argument('-d', '--depth', help='depth limit to recursive descent',
                      dest='max_depth', type=int, default=1000)
    parser.add_argument('-l', '--long', help='long-format listing',
                      dest='verbose', action='count', default=1)
    parser.add_argument('-D', '--diagnose', help='diagnostic (debugging) output',
                      dest='DIAGNOSE', action='count', default=0)
    parser.add_argument('-q', '--quiet', help='cumulative summary only (quiet)',
                      dest='verbose', action='store_const', const=0)
    parser.add_argument('-c', '--csv', help='CSV output', default=False,
                      dest='csv_output', action='store_true')
    parser.add_argument('-r', '--recursive', help='recursively list', default=False,
                      dest='recurse', action='store_true')
    args = parser.parse_args()
    
    main(args)
    sys.exit(0)



