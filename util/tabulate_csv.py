#!/usr/bin/env python
r'''
tabulate_csv.py: tabulate CSV file or files according to standard processes

usage:
  `tabulate_csv.py [ -q ] [ -p ] [ -o outfile ] MODE CSVFILE`

where:

* MODE is a comma-separated list of modes, or "all"
    (at present, "funnel" is the only mode)
* CSVFILE is a CSV file (or file pattern) with input fields for the table

* -o outfile gives an output file or output file pattern (stdout default)
* -p means plain text output versus the default markdown
* -q means to continue quietly if CSVFILE is not present

The MODE is intended to be one or several keywords indicating the kind
of tabulation to be done.  There are really two cases:

- single-mode operation.  The CSVFILE can be given explicitly, and so can
    the optional outfile.
- multi-mode operation.  Here the CSVFILE should be a pattern with two
    %s's: the first for the relevant MODE, the second for the "csv" suffix.
    And the outfile, if given, should also contain two %s placeholders: the first
    for the mode and the second for the suffix (md or txt).

The setup is intended to support many tabulations at once.  Now, there is just one:
   `reduce-funnel.csv -> table-funnel.md`

Typical usage:
``` shell
  util/tabulate_csv.py funnel sims/HabEx_4m_foo/reduce-funnel.csv
  util/tabulate_csv.py -q -o sims/HabEx_4m_foo/tbl/table-%s.%s all sims/HabEx_4m_foo/reduce-%s.%s
```

turmon jan 2019
'''

from __future__ import division
from __future__ import print_function
import argparse
import sys
import os
import csv

# global verbosity mode, also usable for debugging print's
VERBOSITY = 0

########################################
###
###  Utility Functions
###
########################################

def ensure_permissions(fn):
    r'''Ensure correct permissions on the named data file.  
    We use rw-rw-r-- = 664 (octal), to allow group-write.'''
    try:
        os.chmod(fn, 0o664)
    except OSError:
        pass # e.g., don't own the file

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


########################################
###
###  Dump funnel-type tables
###
########################################

targets = ('star', 'allplan', 'hzone', 'earth')
table_row_spec = [
    ['Status', 'Star', 'Planet', 'Hab. Zone', 'Earth'],
    ['  ', '[count]', '[count]', '[count]', '[count]'],
    [], # end of table header
    ['Promoted'        ] + ['%%s_%s'            % t for t in targets],
    ['Attempt: Cume'   ] + ['%%s_tries_%s_cume' % t for t in targets],
    ['Attempt: Unique' ] + ['%%s_tries_%s_uniq' % t for t in targets],
    ['Observed: Cume'  ] + ['%%s_chars_%s_cume' % t for t in targets],
    ['Observed: Unique'] + ['%%s_chars_%s_uniq' % t for t in targets],
    ]

def dump_funnel_worker(fp, plain, table, key, name):
    # delimiters at head and foot
    (head, foot) = ('=== ', '\n') if plain else ('# ', '\n')
    # delimiters within rows
    (delim0, delim1) = ('', '') if plain else ('| ', ' | ')
    # field width
    width = 16 if plain else 24
    pm_sign = '+/-' if plain else '&plusmn;'
    # number of columns
    ncol = len(table_row_spec[0])
    # make the file
    fp.write(head)
    fp.write('Summary: %s\n' % name)
    for nr, row in enumerate(table_row_spec):
        # end of table head
        if len(row) == 0:
            if not plain:
                fp.write(delim0)
                fp.write(delim1.join(['-'*width for r in range(ncol)]))
                fp.write(delim0[::-1] + '\n')
            continue
        # row-begin
        fp.write(delim0)
        # column 1
        fp.write('%-*s' % (width, row[0]))
        fp.write(delim1)
        for cell in row[1:]:
            if '%' in cell:
                key_mean = 'funnel_' + (cell % key) + '_mean'
                key_err  = 'funnel_' + (cell % key) + '_sem'
                out = '%.2f %s %.2f' % (table[key_mean], pm_sign, table[key_err])
            else:
                out = cell
            fp.write('%*s' % (width, out))
            fp.write(delim0[::-1]) # reversed row-begin
        fp.write('\n')
    fp.write(foot)

def dump_funnel(fp, args, table):
    dump_funnel_worker(fp, args.plain, table[0], 'deep', 'Deep Dive')
    dump_funnel_worker(fp, args.plain, table[0], 'promo', 'Promotion')

DUMP_DISPATCH_TABLE = {
    'funnel':  dump_funnel,
    # (insert other table-types as needed)
    }

########################################
###
###  Main routine
###
########################################

def main(args):
    # load CSV
    for mode in args.modes:
        # input file based on mode
        infile = args.infile % (mode, 'csv') if '%s' in args.infile else args.infile
        if not os.access(infile, os.R_OK):
            if not args.quiet:
                sys.stderr.write("%s: Error: Supplied csv file `%s' is not readable.\n" % (args.progname, infile))
                sys.exit(1) # stop now
            else:
                continue # if quiet, just move on - not an error
        # get the data
        table = load_csv(infile)
        # output file based on mode
        if not args.outfile:
            output_to_file = False
            out_fp = sys.stdout
        else:
            output_to_file = True
            suffix = 'txt' if args.plain else 'md'
            outfile = args.outfile % (mode, suffix) if '%s' in args.outfile else args.outfile
            try:
                out_fp = open(outfile, 'w')
            except:
                sys.stderr.write("%s: Error: Indicated output file `%s' is not writable.\n" % (args.progname, outfile))
                sys.exit(1)
            print('\tTabulation to %s' % outfile) # only if not to stdout, obviously
        # dump the data
        DUMP_DISPATCH_TABLE[mode](out_fp, args, table)
        if output_to_file:
            out_fp.close()
            ensure_permissions(outfile)
    # indicate success if using template
    # this is a sentinel file indicating the whole set was made
    if '%s' in args.outfile:
        outfile = args.outfile % ('status', 'txt')
        with open(outfile, 'w') as fp:
            fp.write('Tables written by %s\n' % os.getenv('USER'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tabulate CSV file into text or markdown.",
                                     epilog='')
    parser.add_argument('mode', metavar='MODE',
                            help='comma-separated from <%s>, or just "all"' % '|'.join(DUMP_DISPATCH_TABLE.keys()))
    parser.add_argument('infile', metavar='FILE', help='input file name (or template)')
    parser.add_argument('-o', '--outfile', metavar='FILE', help='output file name (or template)',
                            default='')
    parser.add_argument('-p', '--plain', action='store_true', dest='plain',
                            help='Plain-text output.')
    parser.add_argument('-q', action='store_true', dest='quiet',
                            default=False, help='Exit quietly if file not present.')
    parser.add_argument('-v', help='verbosity', action='count', dest='verbose',
                            default=0)
    args = parser.parse_args()
    
    # set umask in hopes that files will be group-writable
    os.umask(0o002)

    VERBOSITY = args.verbose

    # program name, for convenience
    args.progname = os.path.basename(sys.argv[0])
    
    # set args.modes from args.mode
    #   mode = comma-separated modes, or * for all
    if args.mode == 'all':
        args.modes = DUMP_DISPATCH_TABLE.keys()
    else:
        args.modes = args.mode.split(',')

    # now check these modes
    for mode in args.modes:
        if mode not in DUMP_DISPATCH_TABLE:
            sys.stderr.write("%s: Error: MODE (%s) must be legal (%s).\n" % (
                args.progname, mode, ','.join(DUMP_DISPATCH_TABLE.keys())))
            sys.exit(1)

    # if multi-modes, ensure template is in the given output file
    if len(args.modes) > 0 and len(args.outfile) > 0 and '%s' not in args.outfile:
        sys.stderr.write("%s: Need a template (%s) in the output file if it is given with multiple modes.\n" % (
            args.progname, ))
        sys.exit(1)

    # ensure enclosing dir exists
    if len(args.outfile) > 0 and '%s' in args.outfile:
        directory = os.path.dirname(args.outfile % ('dummy', 'txt'))
        if not os.path.exists(directory):
            os.makedirs(directory)

    main(args)
    sys.exit(0)
