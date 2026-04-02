#!/usr/bin/env python3
"""
Script to concatenate yields from the first 24 lines from sims/ across multiple CSV files,
into a single output file.

Synopsis:
   util/ad-hoc-yield-cat.py [-N nrow] SCENARIO [...]

Take initial lines from "reduce-yield-time.csv" within the given list of SCENARIOs, and make a
new yield file from them.

- Some yield columns are deleted. By default, all _red_ columns.
- The time index (first 2 columns) is updated on the second and subsequent files
- A "source" column is added saying where the yield came from

Typical usage:
  $ util/ad-hoc-yield-cat.py sims/SAG23_AAS2026a.fam/Vis_60m.fam/Vis_60m_20251229 sims/SAG23_AAS2026a.fam/NUV35m.fam/UV_35mSS_OMNI sims/SAG23_AAS2026a.fam/Vis_60m.fam/Vis_60m_20251229 
  ad-hoc-yield-cat.py: Vis_60m_20251229: Took 24 of 59 rows
  ad-hoc-yield-cat.py: UV_35mSS_OMNI: Took 24 of 59 rows
  ad-hoc-yield-cat.py: Vis_60m_20251229: Took 24 of 59 rows
  ad-hoc-yield-cat.py: Success: Processed 3 directories.
  ad-hoc-yield-cat.py: Output: Analyses/2026-01-AAS/Vis_60m_20251229-UV_35mSS_OMNI-Vis_60m_20251229/reduce-yield-time.csv

Claude.ai 2025-12-31
"""

import argparse
import sys
import os
from datetime import datetime
from pathlib import Path
import pandas as pd
import numpy as np

def generate_output_dirname(args):
    """
    Generate output directory name by concatenating basenames of all input directories.
    
    Args:
        args: Parsed command line arguments
    
    Returns:
        str: Generated output directory name
    """
    basenames = []
    for directory in args.directories:
        # Get the basename of each directory path
        basename = os.path.basename(directory.rstrip('/'))
        # Handle case where directory is given as '.' or empty
        if not basename or basename == '.':
            basename = os.path.basename(os.path.abspath(directory))
        basenames.append(basename)
    
    # Concatenate all basenames
    outdir = '-'.join(basenames)
    
    # Store in args for use by other functions
    args.outdir = outdir
    

def concatenate(args):
    """
    Read CSV files from directories and concatenate them into a combined dataframe.
    
    Args:
        args: Parsed command line arguments
    
    Returns:
        pd.DataFrame: Combined dataframe with source column
    """
    # Define constants
    csv_filename = "reduce-yield-time.csv"
    
    # List to store dataframes
    dataframes = []
    
    # Process each directory
    prior_endpoint = None
    for Ndir, (directory, Nrow) in enumerate(zip(args.directories, args.N_list)):
        dir_base = os.path.basename(directory.rstrip('/'))
        dir_path = Path(directory)
        csv_path = dir_path / csv_filename
        
        # Check if file exists and is readable
        if not csv_path.exists():
            print(f"{args.script_name}: Error: File {csv_path} does not exist", file=sys.stderr)
            sys.exit(1)
        if not csv_path.is_file():
            print(f"{args.script_name}: Error: {csv_path} is not a file", file=sys.stderr)
            sys.exit(1)
        if not os.access(csv_path, os.R_OK):
            print(f"{args.script_name}: Error: File {csv_path} is not readable", file=sys.stderr)
            sys.exit(1)
        
        # Load the CSV file
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            print(f"{args.script_name}: Error: Failed to read {csv_path}: {e}", file=sys.stderr)
            sys.exit(1)
        
        # Check if file has at least Nrow lines
        if len(df) < Nrow:
            print(f"{args.script_name}: Error: File {csv_path} has only {len(df)} lines, needs at least {Nrow}", 
                  file=sys.stderr)
            sys.exit(1)
        
        # Take first Nrow lines
        df_subset = df.head(Nrow).copy()
        print(f"{args.script_name}: {dir_base}: Took {Nrow} of {len(df)} rows")
        
        # exclude some columns, if desired
        dropped = []
        for x_1 in args.exclude:
            dropped += [col for col in df_subset.columns if x_1 in col]
        df_subset.drop(columns=dropped, inplace=True)

        # Apply column modifications unless --slice
        if not args.slice:
            for col in df_subset.columns:
                # modify the time index
                if col.endswith("time_lo") or col.endswith("time_hi"):
                    # start the offset with zero, otherwise (prior + 1 month)
                    prior = 0.0 if prior_endpoint is None else (prior_endpoint[col] + 30.5)
                    # Add the prior value to the whole column
                    df_subset[col] = df_subset[col] + prior
        
        # Add source column with dataframe number Ndir
        df_subset['source'] = Ndir + 1
        
        dataframes.append(df_subset)

        # record the last entry
        prior_endpoint = df_subset.iloc[-1]
    
    # Concatenate all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    return combined_df


def dump(args, combined_df):
    """
    Write the combined dataframe and README file to the output directory.
    
    Args:
        args: Parsed command line arguments (including args.outdir)
        combined_df: Combined dataframe to write
    """
    # Define constants
    csv_filename = "reduce-yield-time.csv"
    output_dir = Path("Analyses/2026-01-AAS") / args.outdir
    output_csv = output_dir / csv_filename
    readme_file = output_dir / "README.txt"
    
    # Create output directory if it doesn't exist
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        # Set directory permissions (group read/write)
        os.chmod(output_dir, 0o775)
    except Exception as e:
        print(f"{args.script_name}: Error: Failed to create directory {output_dir}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write the combined dataframe to CSV
    try:
        combined_df.to_csv(output_csv, index=False)
        # Set file permissions (group read/write)
        os.chmod(output_csv, 0o664)
    except Exception as e:
        print(f"{args.script_name}: Error: Failed to write output CSV {output_csv}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Create README.txt
    try:
        with open(readme_file, 'w') as f:
            f.write("This file was auto-generated\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Command: {' '.join(sys.argv)}\n")
            f.write(f"Source sims:\n")
            for directory in args.directories:
                f.write(f"\t{directory}\n")
        # Set file permissions (group read/write)
        os.chmod(readme_file, 0o664)
    except Exception as e:
        print(f"{args.script_name}: Error: Failed to write README file {readme_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Success message
    print(f"{args.script_name}: Success: Processed {len(args.directories)} directories.")
    print(f"{args.script_name}: Output: {output_csv}")


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Concatenate heads of reduce-yield-time.csv in multiple sims'
    )
    parser.add_argument(
        'directories',
        nargs='+',
        help='sim/ directories containing reduce-yield-time.csv'
    )
    # "N" is a string
    parser.add_argument("-N", dest='N_given', default="24", help="Line counts from each file")
    parser.add_argument(
        '--slice',
        dest='slice',
        default=False,
        action='store_true',
        help='Slice out, but do not cumulate, columns (for testing)'
    )
    
    args = parser.parse_args()

    # compose the list-of-lines in "N_given" -> list-of-integers in "N_list"
    N_list = args.N_given.split(',')
    if len(N_list) == 1:
        N_list = N_list * len(args.directories)
    args.N_list = [int(x) for x in N_list]
    assert len(args.N_list) == len(args.directories), "Requested lengths must match input arguments"

    # columns to exclude ... if the column name contains this string,
    # it will be omitted from the output
    args.exclude = ['_red_']

    # Store the script basename for message prefixing
    args.script_name = os.path.basename(sys.argv[0])
    
    # Generate output directory name
    generate_output_dirname(args)
    
    # Read and concatenate the dataframes
    combined_df = concatenate(args)
    
    # Write the output files
    dump(args, combined_df)


if __name__ == "__main__":
    main()

