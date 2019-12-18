import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import astropy.constants as const
import astropy.units as u

import EXOSIMS,os.path
import EXOSIMS.MissionSim as msim
import EXOSIMS.StarCatalog.EXOCAT1 as exc

import sys

sys.path.append('/Users/wdula/Documents/exosims/')

scriptfile = '/Users/wdula/Documents/exosims/HabEx_4m_TSDD_top200DD_52m_5100_20180723b/HabEx_4m_TSDD_top200DD_52m_5100_20180723b.json'
sim = EXOSIMS.MissionSim.MissionSim(scriptfile)

drms = glob.glob("drm/*")
spcs = glob.glob("spc/*")

Ms = sim.TargetList.MsTrue   # star masses
num_detections_hist = []
num_visits_4_zero_det_stars = []
num_planet_detections_all = [] 
num_plan_gte_3det = 0   # number of planets with >= 3 visits
num_plan_gt_half_T = 0  # number of planets with >=3 visits that occur over half their period
num_plan_in_hz = 0      # total number of planets with >3 det in the habitable zone
total_planets = 0       # total number of planets over all sims
total_stars = 0         # total number of stars
sInds = []

for i,drm_f in enumerate(drms):
    num_planet_detections = np.array([])
    drm = pickle.load(open(drm_f,'r'))
    seed_num = os.path.basename(drm_f).split('.')[0]
    spc = pickle.load(open(glob.glob("spc/{}*".format(seed_num))[0],'r'))
    star_visits = {}
    plan2time = {}
    total_planets += spc['nPlans']
    total_stars += len(spc['sInds'])
    for d in drm:
        if "det_status" in d.keys():
            star_visits[d['star_ind']] = (star_visits.get(d['star_ind']) or []) + [d, ]
            num_planet_detections = np.concatenate([num_planet_detections, np.array(d['plan_inds'])[np.where(d['det_status']==1)[0]]])
    for skey, item in star_visits.items():
        if len(item) >= 3:
            det_num = 0
            for s_drm in item:
                det_num += np.sum(s_drm['det_status'][s_drm['det_status']==1])
                plan_inds = np.array(s_drm['plan_inds'])[np.where(s_drm['det_status']==1)[0]]
                for plan in plan_inds:
                    plan2time[plan] = (plan2time.get(plan) or []) + [s_drm['arrival_time']]
            T = np.sqrt((spc['a'][s_drm['plan_inds']]**3 * 4 * np.pi)/(const.G * (Ms[skey] + spc['Mp'][s_drm['plan_inds']]))).to('d')
            for p_i, plan in enumerate(s_drm['plan_inds']):
                if plan2time.get(plan) != None:
                    if (plan2time.get(plan)[-1] - plan2time.get(plan)[0]) > T[p_i].value:
                        num_plan_gt_half_T += 1
                        if .95*u.AU < spc['a'][plan] < 1.67*u.AU:
                            num_plan_in_hz += 1
            num_detections_hist.append(det_num)
            if det_num == 0:
                num_visits_4_zero_det_stars.append(len(item))
                if len(item) == 6:
                    sInds.append(item[0]['star_ind'])
    plan_gte_3det = np.where(np.unique(num_planet_detections, return_counts=True)[1] >= 3)[0]
    num_plan_gte_3det += len(plan_gte_3det)
    num_planet_detections_all.append(num_planet_detections)

avg_num_plan_gte_3det = num_plan_gte_3det / len(drms)
avg_num_plan_gt_half_T = num_plan_gt_half_T / len(drms)
avg_num_plan_in_hz = num_plan_in_hz / len(drms)

eta = total_planets/float(total_stars)

plt.subplot(2, 1, 1)
plt.hist(num_detections_hist, np.arange(0,50))
plt.xlabel("Number of Total Detections")
plt.ylabel("Number of Stars")
plt.title("Number of Detections for Stars with >= 3 Visits")


plt.subplot(2, 1, 2)
plt.hist(num_visits_4_zero_det_stars, np.arange(3, 10))
plt.xlabel("Number of Total Visits")
plt.ylabel("Number of Stars")
plt.title("Number of Visits for Stars with >= 3 Visits and No Detecitons")

