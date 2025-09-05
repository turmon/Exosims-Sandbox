#!/usr/bin/env python
"""
Script to generate multiple JSON files from a template with different parameter combinations
and create an index file tracking all generated files.
"""

import json
import sys
import os

###########################################################

# set umask in hopes that created files will be group-writable
os.umask(0o002)

# Define the lists of parameters to sweep
n_det_removes = [1, 2, 3, 4, 5, 6, 7]
n_det_mins = [1, 2, 3, 4, 5]

# experiment-specific place for cache files
cachedir = 'sims/test.fam/Sweep_det_202508.exp/Caches'

# staging directory for trial-and-error
STAGE = 'Staging'

if not os.path.exists(STAGE):
    os.makedirs(STAGE, 0o775)


###########################################################

def main():
    # Read the template JSON file
    try:
        with open('exosims_template.json', 'r') as f:
            template_dict = json.load(f)
    except FileNotFoundError:
        print("Error: exosims_template.json not found")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in exosims_template.json - {e}")
        sys.exit(1)
    
    # accumulate file information
    file_index = []
    
    # Nested loops to iterate over all combinations
    for n_det_min in n_det_mins:
        for n_det_remove in n_det_removes:
            # Create a copy of the template dictionary
            current_dict = template_dict.copy()
            
            # update with swpet values
            current_dict['n_det_min'] = n_det_min
            current_dict['n_det_remove'] = n_det_remove

            # plug in new cache file
            current_dict['cachedir'] = cachedir
            
            # templated script name
            script_name = f"s_det_min{n_det_min}_rm{n_det_remove}.json"
            
            # Write the modified dictionary to a new JSON file
            with open(STAGE + '/' + script_name, 'w') as f:
                json.dump(current_dict, f, indent=2)
            
            # Create index entry for this file
            index_entry = {
                'script_name': script_name,
                'n_det_min': n_det_min,
                'n_det_remove': n_det_remove,
                'run_name': script_name[:-5]  # Remove .json suffix
            }
            
            # Add to the index list
            file_index.append(index_entry)
            
            print(f"Created: {script_name}")
    
    # Write the index file
    with open(STAGE + '/' + 's_index.json', 'w') as f:
        json.dump(file_index, f, indent=2)
    
    print(f"\nGenerated {len(file_index)} files")
    print("Index file created: s_index.json")
    print(f"* All results in {STAGE}/")

if __name__ == "__main__":
    main()

