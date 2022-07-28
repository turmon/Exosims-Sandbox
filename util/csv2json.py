#!/usr/bin/env python
#
# csv2json.py: Convert table in CSV format to a list of JSON objects
#
# Usage:
#     csv2json-v1.py test.csv test.json
#   (or)
#     csv2json-v1.py test.csv > anything.json
#
#
# Given a CSV file, with a single header line, like this:
# _____________________
# Letter,Frequency,Percentage
# A,24373121,8.1
# B,4762938,1.6
# C,8982417,3.0
# _____________________
# converts it to a JSON file, which is a list of objects, one per
# row, like this:
# _____________________
#   [
#       {
#           "RowNum": 1,
#           "Letter": "A",
#           "Frequency": 24373121,
#           "Percentage": 8.1
#       },
#       {
#           "RowNum": 2,
#           "Letter": "B",
#           "Frequency": 4762938,
#           "Percentage": 1.6
#       },
#       {
#           "RowNum": 3,
#           "Letter": "C",
#           "Frequency": 8982417,
#           "Percentage": 3
#       }
#   ]
# _____________________
# * Note that the RowNum field will be added for convenience.
# * Not tested with complex/quoted strings within the input CSV
#
# turmon 2022-july


import os
import sys
import csv
import json
import argparse


def remap(row, n):
    r'''Map values in a row from the CSV-reader into a new row,
    converting some types if possible.'''
    # ensures this is the first field in the output
    row_new = {'RowNum': n+1}
    for key, val in row.items():
        try:
            val_new = float(val)
            # up-convert to int, if possible
            if int(val_new) == val_new:
                val_new = int(val_new)
        except:
            # it was a general string ... leave it
            val_new = val
        row_new[key] = val_new
    return row_new

def make_json(fn_csv, fn_json):
    r'''Convert a named CSV file to a named JSON file.
    Over-writes the second file.'''
	
    # populate a list with dictionaries, one per CSV row
    data = []
    with open(fn_csv, encoding='utf-8') as fp:
        # sniff for a header
        if not csv.Sniffer().has_header(fp.read(4096)):
            sys.stderr.write('Warning: Could not find header in %s.\n' % fn_csv)
            sys.stderr.write('Proceeding anyway.\n')
        fp.seek(0)
		# append all rows to data
        csvReader = csv.DictReader(fp)
        for n, row in enumerate(csvReader):
            data.append(remap(row, n))
    # write json to the file
    if fn_json == '-':
        fp = sys.stdout
    else:
        fp = open(fn_json, 'w', encoding='utf-8')
    fp.write(json.dumps(data, indent=4))
    fp.close()
		
# Decide the two file paths according to your
# computer system
csvFilePath = r'Names.csv'
jsonFilePath = r'Names.json'

def main(args):
    # Call the make_json function
    make_json(args.file, args.outfile)
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert CSV-with-header to a row-list in JSON.')
    parser.add_argument('file', help='CSV input filename', metavar='FILE')
    parser.add_argument('outfile', nargs='?', default='-', help='JSON output filename (omit for stdout)', metavar='OUTFILE')
    args = parser.parse_args()
    if not os.path.isfile(args.file):
        sys.exit('Fatal: Input file "%s" not readable.' % args.file)
    sys.exit(main(args))

