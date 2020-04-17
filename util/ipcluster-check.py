#!/usr/bin/env python
#
# Check an ipython parallel controller/engine.
#
# Usage:
#   ...
# Simplest usage:
#   ipcluster-check.py
#
#
# For more on usage, use the -h option.
# Some options may be described there but not documented here.
#

# author:
#  Michael Turmon, JPL, 2017
#

from __future__ import print_function
import ipyparallel as ipp
import argparse
import sys

def verify(fn, profile):
    if fn:
        arglist = dict(url_file=fn)
    elif profile:
        arglist = dict(profile=profile)
    else:
        arglist = {}

    c = ipp.Client(**arglist)
    ids = c.ids
    # a trivial computation
    result = list(c[:]).apply_sync(lambda : "Hello")
    # how many came back ok
    ok = sum([1 for s in result if s == "Hello"])
    message = 'Connected to %d engines, and jobs returned from %d.  %s.' % \
      (len(ids), ok, 'All OK' if len(ids) == ok else 'Error')
    return message


############################################################
#
# Main routine
#
############################################################

def main(args):
    r'''Main routine.  Returns nonzero for trouble.'''
    message = verify(args.file, args.profile)
    print(message)
    return 0 if 'OK' in message else 1

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Check ipcluster.",
                                     epilog='')
    parser.add_argument('--profile', type=str, default=None,
                            help='Profile name (as used in ipcluster) for ipython engines.')
    parser.add_argument('--file', type=str, default=None,
                            help='ipcontroller-client.json file (as used in ipcluster) for ipython engines.')
    args = parser.parse_args()
    
    if args.file and '.json' not in args.file:
        raise ValueError('If given, file must be an ipcontroller json file')

    status = main(args)
    sys.exit(status)


