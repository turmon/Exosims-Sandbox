#!/usr/bin/env python
r'''
select_ensembles.py: select from a list of simulation-ensembles

usage:
  select_ensembles.py [ -q ] [ -n N ] [ -k key ] [ -o outkey] MODE CSVFILE

where:
  MODE is "top", "bottom", or "mix"
  CSVFILE is a CSV file containing yields and experiment names
     (it is as generated by reduce-drms.py)
and:
  -n N tells how many to select.  Default is 10.  Non-negative.
  -k key tells which key to use for "top", "bottom", or "mix".  Default is first key.
  -o outkey tells which key to print
  -q means to exit quietly if CSVFILE is not present

The interpretation of MODE is:
  top (bottom): select the rows containing the top (resp. bottom) N values of KEY, 
    with values compared numerically if possible, and alphabetically otherwise.
  mix: select an arbitrary mix of N rows, with the rows selected according to the
    md5 hash of the value of KEY.  This provides stability across invocations if 
    rows are added to the CSVFILE between runs, provided the KEY is unchanged.

In addition, the pseudokey _index (for -k or -o) gives the ordinal number
of the CSV row.  Also, for MODE of mix, when -k key is used, the pseudokey 
_key is available as output, and yields the md5 hash of "key".

Typical usage:
  util/select_ensembles.py -o experiment -k chars_earth_unique top sims/HabEx_4m_TSDDold_DD_TF17_promo6_20190123.exp/reduce-yield-plus.csv

turmon jan 2019
'''

from __future__ import division
from __future__ import print_function
from operator import itemgetter
import argparse
import sys
import hashlib
import os
import csv

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

########################################
###
###  Utility Functions
###
########################################

def load_csv(infile):
    # load the CSV
    info = []
    try:
        with open(infile, 'r') as csvfile:
            reader = csv.DictReader(csvfile, delimiter=',')
            for row in reader:
                info.append(row)
    except IOError:
        sys.stderr.write("Error: could not open `%s'.  Quitting." % infile)
        raise
    # (no conversion possible)
    if len(info) == 0: return info
    # find conversion (just string-to-number where possible)
    # (we just inspect the first row)
    c_func = {}
    for k in info[0].keys():
        try:
            float(info[0][k])
            c_func[k] = float
        except ValueError:
            c_func[k] = str
    # do the conversion
    new_info = []
    for index,row in enumerate(info):
        new_row = {}
        for k in row:
            new_row[k] = c_func[k](row[k])
        new_row['_index'] = index
        new_info.append(new_row)
    return new_info

def decorate(table, key):
    r'''Adds hashed value to table, in-place.'''
    for row in table:
        # ensure a string, then convert to bytes
        # note, repr(float) is the same between Py2/Py3, but str(float) is not
        val = repr(row[key]).encode()
        h = hashlib.md5()
        h.update(val)
        # hash as a hex string, so it can be output with -o if desired
        row['_' + key] = h.hexdigest()

def select(table, n, sortkey, reverse):
    table_s = sorted(table, key=itemgetter(sortkey), reverse=reverse)
    return table_s[:n]

def main(args):
    # load CSV
    table = load_csv(args.infile)
    # prevent special cases
    if args.n == 0 or len(table) == 0: return
    # assign sort key
    if args.key:
        args.sortkey = args.key
    else:
        args.sortkey = list(table[0].keys())[0]
    if args.sortkey not in table[0]:
        print("%s: Error: supplied sort key not in CSV." % args.progname)
        sys.exit(1)
    # case analysis to handle selection mode
    if args.mode == 'mix':
        # decorate the table with a new key, _ + sortkey, which has the md5 hash of
        # the sortkey.  We will sort on the new key to get a "random" but stable
        # selection of rows.  (By stable, we mean that adding a new row in a subsequent
        # simulation run will not change the selection radically.)
        decorate(table, args.sortkey)
        args.sortkey = '_' + args.sortkey
        reverse = False # sort order does not matter in this case
    elif args.mode == 'bottom':
        reverse = False
    elif args.mode == 'top':
        reverse = True
    else:
        print("%s: Error: bad mode." % args.progname)
        sys.exit(1)
    # select some rows
    table_s = select(table, args.n, args.sortkey, reverse)
    # assign output key -- error check now because decorate() added a key
    if args.outkey:
        outkey = args.outkey
    else:
        outkey = list(table[0].keys())[0]
    if outkey not in table_s[0]:
        print("%s: Error: requested output key not in CSV." % args.progname)
        sys.exit(1)
    # output the requested key
    for row in table_s:
        print(row[outkey])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Select Ensembles from Experiment.  Print the selected keys to stdout.",
                                     epilog='')
    parser.add_argument('mode', metavar='MODE', help='mode: "top" or "mix"')
    parser.add_argument('infile', metavar='FILE', help='mode name')
    parser.add_argument('-k', '--key', type=str, default='',
                            help='Key for selection.')
    parser.add_argument('-o', type=str, default='',
                            dest='outkey', help='Key for output.')
    parser.add_argument('-n', help=' Number of ensemble keys output, default = %(default)d',
                      type=int, dest='n', default=10)
    parser.add_argument('-q', help='quiet', action='store_true', dest='quiet',
                            default=False)
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    args = parser.parse_args()
    
    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    if args.n < 0:
        print('%s: Error: n must be nonnegative.' % args.progname)
        sys.exit(1)

    if args.mode not in ('mix', 'bottom', 'top'):
        print("%s: Error: MODE must be `mix' 'bottom' or `top'." % args.progname)
        sys.exit(1)

    if not os.access(args.infile, os.R_OK):
        if not args.quiet:
            print("%s: Error: Supplied csv file `%s' is not readable." % (args.progname, args.infile))
        sys.exit(0 if args.quiet else 1)

    main(args)
    sys.exit(0)
