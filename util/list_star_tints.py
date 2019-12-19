import pickle
import numpy as np
import csv
import os
import glob
import pdb

targets = ["HIP 37279","HIP 97649","HIP 67927","HIP 2021","HIP 107556",
           "HIP 46853","HIP 22449","HIP 95501","HIP 102422","HIP 86974"]

# make own tint output directory in current dir
# averages should be in file

def full_summary(datadir, hip_list):
    header = ["Target", "RA", "dec", "distance", "Lstar", "Vmag", "avg_num_visits", 
              "avg_comp", "avg_det_times", "avg_char_times", "avg_total_star_det_times","avg_total_star_char_times", 
              "avg_total_ens_det_times", "avg_total_ens_char_times"]

    spks = glob.glob(os.path.join(datadir, "spc", "*.spc"))
    drms = glob.glob(os.path.join(datadir, "drm", "*.pkl"))

    spks = sorted(spks)
    drms = sorted(drms)
    summaries = [header]
    summ = {}

    for i, spk in enumerate(spks):
        if spk.split("/")[-1].split(".")[0] == drms[i].split("/")[-1].split(".")[0]:

            if not os.path.exists(os.path.join(datadir, "tint_out")):
                os.mkdir(os.path.join(datadir, "tint_out"))

            outfile = os.path.join(datadir, "tint_out", spk.split("/")[-1].split(".")[0] + ".csv")
            lines = summarize_results(spk, drms[i], hip_list, outfile)
            for line in lines[1:]:
                if summ.get(line[0]) is None:
                    summ[line[0]] = [line]
                else:
                    summ[line[0]] = summ.get(line[0]) + [line]
        else:
            print("DIR mismatch: {} is not {}".format(spk, drms[i]))

    hips = summ.keys()
    hips = sorted(hips)
    for hip_key in hips:
        # pdb.set_trace()
        lines = summ[hip_key]
        target = hip_key
        ra, dec, dist = lines[0][1:4]
        l = lines[0][4]
        vmag = lines[0][5]
        avg_num_vis = np.mean([line[6] for line in lines[:]])
        avg_comp = np.mean([line[7] for line in lines[:]])
        avg_det_times = np.mean(np.concatenate([line[8] for line in lines[:]]))
        avg_char_times = np.mean(np.concatenate([line[9] for line in lines[:]]))
        avg_total_star_det_times = np.mean([line[10] for line in lines[:]])
        avg_total_star_char_times = np.mean([line[11] for line in lines[:]])
        avg_total_ens_det_times = np.mean([line[12] for line in lines[:]])
        avg_total_ens_char_times = np.mean([line[13] for line in lines[:]])

        line = [target, ra, dec, dist, l, vmag, avg_num_vis, 
                avg_comp, avg_det_times, avg_char_times, avg_total_star_det_times, 
                avg_total_star_char_times, avg_total_ens_det_times, avg_total_ens_char_times]
        summaries.append(line)

    outfile = os.path.join(datadir, "tint_out", "tint_summary.csv")
    with open(outfile, 'w') as csvfile:
        c = csv.writer(csvfile)
        c.writerows(summaries)

    return summ, summaries


def summarize_results(spk_fpath, drm_fpath, hip_list, outfile):
    spk = pickle.load(open(spk_fpath,'r'))
    drm = pickle.load(open(drm_fpath,'r'))
    header = ["Target", "RA", "dec", "distance", "Lstar", "Vmag", "num_visits", 
              "comp", "det_times", "char_times", "total_star_det_times","total_star_char_times", 
              "total_ens_det_times", "total_ens_char_times"]
    lines = [header]
    print(spk_fpath)
    hip_list = sorted(hip_list)
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
            comp = 0
            num_visits = 0
            for d in drm:
                if 'det_time' in d.keys():
                    all_det_times.append(d['det_time'].value)
                    if str(d['star_ind']) == str(sInd):
                        det_times.append(d['det_time'].value)
                        num_visits += 1
                        if comp == 0:
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

