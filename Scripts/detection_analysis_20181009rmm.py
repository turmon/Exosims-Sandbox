
import matplotlib.pyplot as plt
import numpy as np
import pickle
import glob
import astropy.constants as const
import astropy.units as u
import sys

sys.path.append('/proj/exep/rhonda/Sandbox/HabEx/EXOSIMS')

import EXOSIMS,os.path
import EXOSIMS.MissionSim as msim
import EXOSIMS.StarCatalog.EXOCAT1 as exc

import pdb

def produce_detection_analysis(drms, spcs, scriptfile):
    sim = EXOSIMS.MissionSim.MissionSim(scriptfile,nopar=True)

    Ms = sim.TargetList.MsTrue   # star masses
    num_detections_hist = []
    Rp_bins = np.array([0.5, 1.0, 1.75, 3.5, 6.0, 14.3])
    plan_rad_hist_half_T = []
    plan_rad_hist_in_hz = []


    num_visits_4_zero_det_stars = []
    num_planet_detections_all = []

    num_plan_gte_3det = 0   # number of planets with >= 3 visits
    num_plan_gt_half_T = 0  # number of planets with >=3 visits that occur over half their period
    num_plan_in_hz = 0      # total number of planets with >3 det in the habitable zone
    num_plan_earthlike = 0

    num_stars_gte_3det = 0
    num_stars_gt_half_T = 0
    num_stars_in_hz = 0

    total_planets = 0       # total number of planets over all sims
    total_stars = 0         # total number of stars
    sInds = []
    t_det_earthlike = []

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
            num_stars_gt_half_T_flag = False
            num_stars_in_hz_flag = False
            if len(item) >= 3:
                num_stars_gte_3det += 1
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
                            num_stars_gt_half_T_flag = True
                            num_plan_gt_half_T += 1
                            plan_rad_hist_half_T.append(spc['Rp'][plan].value)
                            if .95*u.AU < spc['a'][plan] < 1.67*u.AU:
                                num_stars_in_hz_flag = True
                                num_plan_in_hz += 1
                                plan_rad_hist_in_hz.append(spc['Rp'][plan].value)
                                if .8*(spc['a'][plan]**-.5).value < spc['Rp'][plan].value < 1.4:
                                    num_plan_earthlike += 1
                                    for s_drm in item:
                                        t_det_earthlike.append(s_drm['det_time'].value)

                if num_stars_gt_half_T_flag == True:
                    num_stars_gt_half_T += 1
                if num_stars_in_hz_flag == True:
                    num_stars_in_hz += 1

                num_detections_hist.append(det_num)

                if det_num == 0:
                    num_visits_4_zero_det_stars.append(len(item))
                    if len(item) == 6:
                        sInds.append(item[0]['star_ind'])

        plan_gte_3det = np.where(np.unique(num_planet_detections, return_counts=True)[1] >= 3)[0]
        num_plan_gte_3det += len(plan_gte_3det)
        num_planet_detections_all.append(num_planet_detections)

    avg_num_plan_gte_3det = num_plan_gte_3det / float(len(drms))
    avg_num_plan_gt_half_T = num_plan_gt_half_T / float(len(drms))
    avg_num_plan_in_hz = num_plan_in_hz / float(len(drms))
    avg_num_plan_earthlike = num_plan_earthlike / float(len(drms))

    avg_num_stars_gte_3det = num_stars_gte_3det / float(len(drms))
    avg_num_stars_gt_half_T = num_stars_gt_half_T / float(len(drms))
    avg_num_stars_in_hz = num_stars_in_hz / float(len(drms))

    eta = total_planets/float(total_stars)

    avg_t_det_earthlike = sum(t_det_earthlike)/len(t_det_earthlike)

    print('\n===================\n')
    print('Avg Num Planets >= 3 Detections: {}').format(avg_num_plan_gte_3det)
    print('Avg Num Planets Detected Over a Half-Period: {}').format(avg_num_plan_gt_half_T)
    print('Avg Num Planets in The Habitable Zone: {}').format(avg_num_plan_in_hz)
    print('Avg Num Earthlike Planets: {}').format(avg_num_plan_earthlike)
    print('Avg Detection Time for Earthlike Planets: {}'.format(avg_t_det_earthlike))
    print('\n===================\n')
    print('Avg Num Stars >= 3 Detections: {}').format(avg_num_stars_gte_3det)
    print('Avg Num Stars with Detections Over a Half-Period: {}').format(avg_num_stars_gt_half_T)
    print('Avg Num Stars with Planets in The Habitable Zone: {}').format(avg_num_stars_in_hz)
    print('\n===================\n')
    print('eta: {}').format(eta)

    plt.subplot(4, 1, 1)
    plt.hist(num_detections_hist, np.arange(0,50))
    plt.xlabel("Number of Total Detections")
    plt.ylabel("Number of Stars")
    plt.title("Number of Detections for Stars with >= 3 Visits, {} DRMs".format(len(drms)))

    plt.subplot(4, 1, 2)
    plt.hist(num_visits_4_zero_det_stars, np.arange(3, 10))
    plt.xlabel("Number of Total Visits")
    plt.ylabel("Number of Stars")
    plt.title("Number of Visits for Stars with >= 3 Visits and No Detections")

    plt.subplot(4, 1, 3)
    plt.hist(plan_rad_hist_half_T, Rp_bins)
    plt.xlabel("Planet Radii")
    plt.ylabel("Number of Planets")
    plt.title("Radius of Planets Detected over Half their Periods")

    plt.subplot(4, 1, 4)
    plt.hist(plan_rad_hist_in_hz, Rp_bins)
    plt.xlabel("Planet Radii")
    plt.ylabel("Number of Planets")
    plt.title("Radius of Planets Detected over Half their Periods in the Habitable Zone")
    plt.gcf().set_size_inches(10.5, 18.5)
    plt.savefig('det_analysis_plot.png', dpi=100)


if __name__ == "__main__":

    sys.path.append('/proj/exep/rhonda/Sandbox/HabEx/sims/')
    sys.path.append('/proj/exep/rhonda/Sandbox/HabEx/Scripts/')

    drms = glob.glob("drm/*")
    spcs = glob.glob("spc/*")
    scriptfile = '/proj/exep/rhonda/Sandbox/HabEx/Scripts/HabEx_4m_TSDD_Soto_maxvisit10_5384_20181017.json'

    produce_detection_analysis(drms, spcs, scriptfile)    


