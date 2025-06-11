#
# demo of a few plots using drm_tabulate as a module

# %matplotlib tk

import glob

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from_drm = False
if from_drm:
    import drm_tabulate

    in_spec = 'util/demo_plot_tab_args.json'
    dt_args = drm_tabulate.parse_arglist(['-f', in_spec])
    # in_pattern = 'sims/H6C_CODulzE_A_IFSp7_AAS2023_20230102.exp/s_H6C_base_750nm/drm/1[12]*.pkl'
    in_pattern = 'sims/H6C_CODulzE_A_IFSp7_AAS2023_20230102.exp/s_H6C_base_750nm/drm/*.pkl'
    in_files = glob.glob(in_pattern)
    df = pd.DataFrame(drm_tabulate.main(dt_args, in_files))
else:
    ensemble_csv = 'drm_info.csv'
    spc_csv = 'spc.csv'
    df = pd.read_csv(ensemble_csv)
    spc = pd.read_csv(spc_csv)
    

# number of ensemble members
Nens = len(pd.unique(df.seed))

# chars vs. time
bins_missiontime = np.arange(0, 360, 30)
plt.hist(df.arrival_time, bins=bins_missiontime)
plt.title('Characterizations vs. Elapsed Time')
plt.xlabel('Mission Time [d]')
plt.ylabel('Characterizations')

# completeness vs. time
arrival_binned = np.digitize(x=df.arrival_time, bins=bins_missiontime)
ax = sns.swarmplot(data=df, y='comp0', x=arrival_binned, hue='char_time', orient='v', palette='cool')

norm = plt.Normalize(0, max(df.char_time))
sm = plt.cm.ScalarMappable(cmap='cool', norm=norm)
sm.set_array([])
ax.legend_.remove()
ax.figure.colorbar(sm, label='Integration Time [d]')

plt.title('Target Completeness vs. Elapsed Time\n(Point color: Integration time [d])')
plt.xlabel('Mission Time [month]')
plt.ylabel('Completeness')


# #chars histogram
plt.clf()
bins_chars = np.arange(0, 20, 1)
plt.hist(df.seed.value_counts(), bins=bins_chars)
plt.title('Histogram of Number of Characterizations')
plt.xlabel('Characterizations')
plt.ylabel('Probability')

# distance/luminosity
sns.scatterplot(data=df, x='dist', y='luminosity', hue='comp0')
plt.yscale('log')
plt.title('Completeness vs. Host Star Properties')
plt.xlabel('Distance [pc]')
plt.ylabel('Luminosity [Lsun]')

