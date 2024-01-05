# script (no arguments) that:
#   -- loads the s_index file for an Experiment, as a pandas table
#      (one row per Ensemble)
#   -- iterates over the runs therein, and loads the reduce-info.csv for each
#   -- adds the reduction info as more columns in the run list
#   -- writes the pandas file as a .pkl
#
# ...it could have just been a csv, but I wanted to allow more complex outputs if desired

import pandas as pd
from pathlib import Path
import pickle

#exp_name = 'aas_run_051823.exp'
#exp_name = 'aas_run_052523.exp'
exp_name = 'aas_run_v3.1.4_053123.exp'

t_index = pd.read_json(open(Path('Scripts') / exp_name / 's_index.json'))
t_index.set_index('run_name', inplace=True)

t_infos = []
for run_name in t_index.index:
    info_name = Path('sims') / exp_name / run_name / 'reduce-info.csv'
    t1 = pd.read_csv(open(info_name))
    # add this distinguishing column
    t1['run_name'] = run_name
    t_infos.append(t1)

# stack list-of-tables => single table
t_info_as_read = pd.concat(t_infos).set_index('run_name')

# combine into: [t_index t_info_as_read]
t_info = pd.concat([t_index, t_info_as_read], axis=1)

print(f'Got {t_info.shape[0]} rows of {t_info.shape[1]} columns.')

t_out = Path(exp_name).stem + '-t_info.pkl'
with open(t_out, 'wb') as fp:
    pickle.dump(t_info, fp)

print(f'Dumped to "{t_out}"')
