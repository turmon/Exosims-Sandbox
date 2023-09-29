#!/usr/bin/env python
r'''
select_runs.py: select runs from a list

usage:
  select_runs.py [ -q ] [ -n N ] [ -m MODE ] SCENARIO

where:
  SCENARIO is a location (for now, directory) containing run information

and:
  -m MODE tells what mode
     by default this is "drm" (output the DRM filenames within the SCENARIO)
  -n N tells how many runs to select.  Default is 10.  Non-negative, or T for all runs.
     It is not an error if N is larger than the number of runs.
  -q means to exit quietly if SCENARIO is not present

Typical usage:
  util/select_runs.py -n 5 sims/HabEx_4m_TSDD

turmon aug 2023
'''

from operator import itemgetter
import argparse
import sys
import hashlib
import os
import glob

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

########################################
###
###  Utility Functions
###
########################################

def load_runs(scenario):
    # get run names
    try:
        runs = glob.glob(os.path.join(scenario, 'drm', '*.pkl'))
    except IOError:
        print(f"Error: IO error finding DRMs in `{scenario}'.  Quitting.", file=sys.stderr)
        raise
    # make a key for each row
    new_info = []
    for run in runs:
        # make the hash from the basename only
        # this way, different scenarios will have matching run orders
        run_base = os.path.basename(run)
        hasher = hashlib.md5()
        hasher.update(run_base.encode())
        new_info.append({'run': run, '_key': hasher.hexdigest()})
    return new_info


def select(table, n, sortkey):
    table_s = sorted(table, key=itemgetter(sortkey))
    return table_s[:n]


def main(args):
    # load the run-list
    table = load_runs(args.scenario)
    # prevent special cases
    if args.n == 0 or len(table) == 0: return
    # handle the "n = -1" flag
    n_select = args.n if args.n >= 0 else len(table)

    # select some rows
    table_s = select(table, n_select, '_key')

    # output run names
    for row in table_s:
        print(row['run'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Select Runs from Scenario.  Print the selected files to stdout.",
                                     epilog='')
    parser.add_argument('-m', '--mode', metavar='MODE', default='drm', help='mode: "drm" only, for now')
    parser.add_argument('scenario', metavar='SCENARIO', help='scenario name, including sims/ part')
    parser.add_argument('-n', help='number of runs to output, T for all, default = %(default)s',
                      type=str, dest='n_orig', default='10', metavar='N')
    parser.add_argument('-q', help='quiet', action='store_true', dest='quiet', default=False)
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose', default=0)
    args = parser.parse_args()
    
    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # n_orig = 'T' => n = -1 (as a flag for "everything")
    # n_orig converts to int => ensure >= 0
    # else, throw error
    if args.n_orig == 'T':
        args.n = -1
    else:
        try:
            args.n = int(args.n_orig)
        except ValueError:
            print('%s: Error: n must be an integer or "T".' % args.progname, file=sys.stderr)
            sys.exit(1)
        if args.n < 0:
            print('%s: Error: n must be nonnegative.' % args.progname, file=sys.stderr)
            sys.exit(1)
        
    if args.mode not in ('drm', 'DRM'):
        print("%s: Error: MODE must be `drm'." % args.progname, file=sys.stderr)
        sys.exit(1)
    dir_for_mode = 'drm'

    if not os.path.isdir(args.scenario):
        if args.quiet: sys.exit(0)
        print(f"{args.progname}: Error: Supplied scenario directory `{args.scenario}' is not readable.",
                  file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(d_try := os.path.join(args.scenario, dir_for_mode)):
        if args.quiet: sys.exit(0)
        print(f"{args.progname}: Error: Implied DRM directory `{d_try}' is not readable.",
                  file=sys.stderr)
        sys.exit(1)

    main(args)
    sys.exit(0)
