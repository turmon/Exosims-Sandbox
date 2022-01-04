# early script which appears to extract some fields about some targets from 
# the drm/spc, and export to CSV.

from __future__ import print_function
import pickle
import numpy as np
import csv
import glob

targets = ["HIP 37279","HIP 97649","HIP 67927","HIP 2021","HIP 107556",
           "HIP 46853","HIP 22449","HIP 95501","HIP 102422","HIP 86974"]


def full_summary(spk_dir, drm_dir, hip_list, outfile):
    header = ["Target", "RA", "dec", "distance", "Lstar", "Vmag", "avg_num_visits", 
              "comp", "avg_det_times", "avg_char_times", "avg_total_det_times","avg_total_char_times", 
              "avg_total_all_det_times", "avg_total_all_char_times"]
    if spk_dir[-1] != '/':
        spks = glob.glob(spk_dir + "/*.spc")
    else:
        spks = glob.glob(spk_dir + "*.spc")
    if drm_dir[-1] != '/':
        drms = glob.glob(drm_dir + "/*.pkl")
    else:
        drms = glob.glob(drm_dir + "*.pkl")
    spks = sorted(spks)
    drms = sorted(drms)
    summaries = [header]
    summ = {}
    for i, spk in enumerate(spks):
        if spk.split("/")[-1].split(".")[0] == drms[i].split("/")[-1].split(".")[0]:
            lines = summarize_results(spk, drms[i], hip_list, outfile)
            for line in lines[1:]:
                if summ.get(line[0]) is None:
                    summ[line[0]] = [line]
                else:
                    summ[line[0]] = summ.get(line[0]) + [line]
        else:
            print("DIR mismatch: {} is not {}".format(spk, drms[i]))
    for hip_key, lines in summ.items():
        lines = np.array(lines)
        line = np.append(np.append([hip_key], lines[0,1:6]), [np.mean(lines[:,6]), lines[0,7]])
        line = np.append(line, [np.mean(np.concatenate(lines[:,8])), np.mean(np.concatenate(lines[:,9])), np.mean(lines[:,10])])
        line = np.append(line, [np.mean(lines[:,11]), np.sum(lines[:,12]), np.sum(lines[:,13])])
        summaries.append(line.tolist())
    return summ, summaries


def summarize_results(spk_fpath, drm_fpath, hip_list, outfile):
    spk = pickle.load(open(spk_fpath,'r'))
    drm = pickle.load(open(drm_fpath,'r'))
    header = ["Target", "RA", "dec", "distance", "Lstar", "Vmag", "num_visits", 
              "comp", "det_times", "char_times", "total_det_times","total_char_times", 
              "total_all_det_times", "total_all_char_times"]
    lines = [header]
    for t in hip_list:
        spk_t_ind = np.where((np.array(spk['Name']) == t))[0]
        if np.any(spk_t_ind):
            coords = spk['coords'][spk_t_ind]
            (RA, dec, distance) = (coords.ra.value[0], coords.dec.value[0], coords.distance.value[0])
            Lstar = spk['L'][spk_t_ind][0]
            sInd = spk['sInds'][spk_t_ind][0]
            line = []
            all_det_times = []
            all_char_times = []
            det_times = []
            char_times = []
            comp = None
            num_visits = 0
            for d in drm:
                if 'det_time' in d.keys():
                    all_det_times.append(d['det_time'].value)
                    if str(d['star_ind']) == str(sInd):
                        det_times.append(d['det_time'].value)
                        num_visits += 1
                        if comp is None:
                            try:
                                comp = d['det_comp']
                            except:
                                pass
                if 'char_time' in d.keys():
                    all_char_times.append(d['char_time'].value)
                    if str(d['star_ind']) == str(sInd):
                        char_times.append(d['char_time'].value)
            total_det_times = sum(det_times)
            total_char_times = sum(char_times)
            total_all_det_times = sum(all_det_times)
            total_all_char_times = sum(all_char_times)
            line = [t, RA, dec, distance, Lstar, "N/A", num_visits, 
                    comp, det_times, char_times, total_det_times, 
                    total_char_times, total_all_det_times, total_all_char_times]
            lines.append(line)
    with open(outfile, 'w') as csvfile:
        c = csv.writer(csvfile)
        c.writerows(lines)
    return lines

