#!/usr/bin/env python
"""
Driver script to generate all DRM plots

This script orchestrates the execution of all DRM plotting functions in a
table-driven manner, loading CSV files and calling each plotting function
with the appropriate arguments.
"""

import sys
import os
import time
import argparse
import pandas as pd
import importlib
from pathlib import Path
# this import must work: fail fast if it doesn't
import plot_drm_gallery

# Program name for error messages
PROGNAME = os.path.basename(sys.argv[0])

# Plot registry - table-driven configuration
# Each entry specifies:
#   - name: short identifier for the plot
#   - module: Python module name to import
#   - function: function name to call within the module
#   - csv_files: list of CSV file identifiers (used with src_tmpl)
#   - enabled: whether to run this plot
#   - mode: additional mode dict to merge with overall mode (optional)
PLOT_REGISTRY = [
    {
        'name': 'yield_times',
        'module': 'plot_drm_yield_times',
        'function': 'plot_drm_yield_times',
        'csv_files': ['yield-time'],
        'enabled': True,
        'mode': {},  # Additional mode settings for this plot
    },
    {
        'name': 'fuel_used',
        'module': 'plot_drm_fuel_used',
        'function': 'plot_drm_fuel_used',
        'csv_files': ['times'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'events',
        'module': 'plot_drm_events',
        'function': 'plot_drm_events',
        'csv_files': ['events'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'event_counts',
        'module': 'plot_drm_event_counts',
        'function': 'plot_drm_event_counts',
        'csv_files': ['event-counts', 'earth-char-count'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'visit_times',
        'module': 'plot_drm_visit_times',
        'function': 'plot_drm_visit_times',
        'csv_files': ['visit-time'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'time_used',
        'module': 'plot_drm_time_used',
        'function': 'plot_drm_time_used',
        'csv_files': ['times'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'promote',
        'module': 'plot_drm_promote',
        'function': 'plot_drm_promote',
        'csv_files': ['promote', 'promote-hist'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'star_targets',
        'module': 'plot_drm_star_targets',
        'function': 'plot_drm_star_targets',
        'csv_files': ['star-target'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'earth_chars',
        'module': 'plot_drm_earth_chars',
        'function': 'plot_drm_earth_chars',
        'csv_files': ['earth-char-list'],
        'enabled': True,
        'mode': {},
    },
    {
        'name': 'radlum',
        'module': 'plot_drm_radlum',
        'function': 'plot_drm_radlum',
        'csv_files': ['radlum', 'earth'],
        'enabled': True,
        'mode': {},
    },
]


def signal_plot_end(args):
    '''Signal success to a sentinel file.

    The content of the file is unimportant, but its existence tells
    the build system that the plots are complete, or sufficiently
    complete.
    '''
    fn_output = args.dest_tmpl % ('info', 'txt')
    with open(fn_output, 'w') as fp:
        print(f"Graphics written by {os.getenv('USER')}", file=fp);
        print(f"Ensemble of N={args.reduce_info.get('ensemble_size', 'N/A')} runs", file=fp)
        print(f"Runs were reduced on {args.reduce_info.get('runtime', 'N/A')}",
              file=fp, flush=True);
    # ensure group-writable
    os.chmod(fn_output, 0o664)


def load_csv_files(src_tmpl, csv_files, plot_name):
    """Load CSV files, returning list of DataFrames or None if any missing."""
    dataframes = []
    for csv_file in csv_files:
        filepath = src_tmpl % (csv_file, 'csv')
        if not os.path.exists(filepath):
            print(f"{PROGNAME}: Warning: CSV file '{filepath}' not found, skipping {plot_name} plot",
                  file=sys.stderr)
            return None
        try:
            dataframes.append(pd.read_csv(filepath))
        except Exception as e:
            print(f"{PROGNAME}: Warning: Could not read '{filepath}': {e}, skipping {plot_name} plot",
                  file=sys.stderr)
            return None
    return dataframes


def run_plot(plot_config, reduce_info, src_tmpl, dest_tmpl, overall_mode):
    """
    Run a single plot function

    Parameters
    ----------
    plot_config : dict
        Configuration dictionary for the plot
    reduce_info : dict
        Metadata dict from reduce-info.csv, passed to each plotter
    src_tmpl : str
        Template string for input CSV files
    dest_tmpl : str
        Template string for output graphics files
    overall_mode : dict
        Universal mode settings, passed to each plotter

    Returns
    -------
    tuple of (bool, list)
        (True, files_written) if plot succeeded, (False, []) otherwise
    """
    name = plot_config['name']
    module_name = plot_config['module']
    function_name = plot_config['function']
    csv_files = plot_config['csv_files']
    plot_mode = plot_config.get('mode', {})
    
    # this function's verbosity verbosity = overall verbosity
    verbose = overall_mode['verbose']

    # Merge global mode with plot-specific mode
    # Plot-specific mode takes precedence
    merged_mode = {**overall_mode, **plot_mode}
    
    # Skip if disabled
    if not plot_config.get('enabled', True):
        if verbose > 1:
            print(f"Skipping {name} (disabled in registry)")
        return True, []

    # Load required CSV files
    plot_data = load_csv_files(src_tmpl, csv_files, name)
    if plot_data is None:
        return False, []

    try:
        # Import the module
        if verbose > 1:
            print(f"Running {name}...")

        module = importlib.import_module(f"plot_drm_gallery.{module_name}")
        plot_function = getattr(module, function_name)

        # Call the plot function
        files_written = plot_function(reduce_info, plot_data, dest_tmpl, merged_mode)
        if files_written is None:
            files_written = []

        if verbose > 1:
            print(f"  {name} completed successfully")

        return True, files_written

    except ImportError as e:
        print(f"{PROGNAME}: Fatal: Could not import {module_name}: {e}",
              file=sys.stderr)
        return False, []
    except AttributeError as e:
        print(f"{PROGNAME}: Fatal: Function {function_name} not found in {module_name}: {e}",
              file=sys.stderr)
        return False, []
    except Exception as e:
        print(f"{PROGNAME}: Fatal: Error running {name}: {e}",
              file=sys.stderr)
        return False, []


def main():
    """
    Main driver function
    """
    parser = argparse.ArgumentParser(
        description='Generate all DRM plots from CSV data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python driver.py "data/%s.%s" "output/det-%s.%s"
    
The first argument is the source template with two %%s placeholders that will be
filled with CSV file identifiers (e.g., "info", "yield_time") and "csv".

The second argument is the destination template with two %%s placeholders for
the plot name and file extension.

Optional arguments:
    --mode_op OP        Set the global mode operation string (default: "*")
    --only PLOT         Run only the specified plot (by name)
    --skip PLOT         Skip the specified plot (by name), can be repeated
    --list              List all available plots and exit
    --verbose           Print verbose progress messages
        """
    )
    
    parser.add_argument('src_tmpl', type=str,
                       help='Source template string (e.g., "data/%%s.%%s")')
    parser.add_argument('dest_tmpl', type=str,
                       help='Destination template string (e.g., "output/det-%%s.%%s")')
    parser.add_argument('--mode_op', type=str, default='*',
                       help='Global mode operation string (default: "*")')
    parser.add_argument('--only', type=str, metavar='PLOT',
                       help='Run only the specified plot (by name)')
    parser.add_argument('--skip', type=str, action='append', default=[], metavar='PLOT',
                       help='Skip the specified plot (by name), can be repeated')
    parser.add_argument('--list', action='store_true',
                       help='List all available plots and exit')
    parser.add_argument('--verbose', '-v', action='count', default=1, 
                       help='Verbosity')
    parser.add_argument('--quiet', '-q', action='store_true', help='Minimal verbosity')
    
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    start_time = time.perf_counter()
    
    if args.quiet: args.verbose = 0

    # List plots if requested
    if args.list:
        print("Available plots:")
        for plot in PLOT_REGISTRY:
            status = "enabled" if plot.get('enabled', True) else "disabled"
            csv_list = ", ".join(plot['csv_files'])
            print(f"  {plot['name']:20s} ({status:8s}) - CSVs: {csv_list}")
        return 0
    
    # Create global mode dictionary
    overall_mode = {'op': args.mode_op, 'verbose': args.verbose}
    
    # Determine which plots to run
    plots_to_run = []
    
    if args.only:
        # Run only the specified plot
        found = False
        for plot in PLOT_REGISTRY:
            if plot['name'] == args.only:
                plots_to_run.append(plot)
                found = True
                break
        if not found:
            print(f"{PROGNAME}: Fatal: Plot '{args.only}' not found in registry",
                  file=sys.stderr)
            return 1
    else:
        # Run all plots except those in skip list
        for plot in PLOT_REGISTRY:
            if plot['name'] not in args.skip:
                plots_to_run.append(plot)
            elif args.verbose > 1:
                print(f"Skipping {plot['name']} (--skip)")
    
    # get basic information into args.reduce_info
    fn_info = args.src_tmpl % ('info', 'csv')
    df_info = pd.read_csv(fn_info)
    args.reduce_info = df_info.iloc[0].to_dict()
    del df_info
    
    # ensure the directory
    dir_path = os.path.dirname(args.src_tmpl % ('dummy', 'txt'))
    try:
        os.makedirs(dir_path, exist_ok=True)
    except FileExistsError:
        # runs if the path exists, but is a file
        print(f"{args.progname}: Fatal: A file exists already at '{dir_path}'.", file=sys.stedrr)
        raise

    # Make the plots
    if args.verbose > 1:
        print(f"Running {len(plots_to_run)} plots...")
        print(f"Source template: {args.src_tmpl}")
        print(f"Destination template: {args.dest_tmpl}")
        print(f"Overall mode: {overall_mode}")
        print()
    
    success_count = 0
    fail_count = 0
    skip_count = 0
    all_records = []

    for plot in plots_to_run:
        result, files_written = run_plot(plot, args.reduce_info, args.src_tmpl, args.dest_tmpl, overall_mode)
        if result:
            success_count += 1
        elif result is False:
            # False means it ran but failed
            fail_count += 1
        else:
            # None or other means it was skipped
            skip_count += 1
        # Collect records: associate each file with its routine and csv_files
        routine = plot['function']
        csv_files = plot['csv_files']
        for gfx_file in files_written:
            all_records.append((gfx_file, routine, csv_files))

    # Write the tabulation CSV
    if all_records:
        fn_csv = args.dest_tmpl % ('plot-list', 'csv')
        with open(fn_csv, 'w') as fp:
            fp.write('graphics_file,routine,csv_files\n')
            for gfx_file, routine, csv_files in all_records:
                fp.write(f'{gfx_file},{routine},{";".join(csv_files)}\n')

    # Summary
    if args.verbose > 0 or fail_count > 0:
        print(f"{args.progname}: Plot Sets:  {success_count} ok, {fail_count} failed, {skip_count} skipped")
        print(f"{args.progname}: Plot Files: {len(all_records)}")
    run_time = time.perf_counter() - start_time
    # unconditionally print 
    print(f"{args.progname}: Done. (Elapsed time: {run_time:.2f} s)")
    
    # FIXME: need a better criterion than "perfect", see also just below
    if fail_count == 0:
        signal_plot_end(args)

    # FIXME: this exit code is inconsistent with the signaling above
    # Return non-zero if any plots failed
    return 1 if fail_count > 0 else 0


if __name__ == '__main__':
    # ensure group-writable from the start
    os.umask(0o002)
    sys.exit(main())
