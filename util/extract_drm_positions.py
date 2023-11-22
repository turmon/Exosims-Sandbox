#!/usr/bin/env python
#
# extract_drm_positions: Extract targets, observatory, and planet positions from DRM
#
# For usage, use the -h option.  Some options may be described there but not here.
#
# This script produces a CSV of target and observatory positions, with one entry for
# every observation in the give DRM.
# Note: this 2023 version is updated for micro-thruster analysis (POC: Cameron Haag)
#
# Algorithm and notes:
# - A pickle summarizing a DRM is loaded, and the corresponding MissionSim 
# object is instantiated from the supplied script.  The observations are 
# stepped through, the sim is moved to the observation time, and properties 
# (like keepout and observatory orbit) are queried.  
# - The DRM file (.pkl extension, name supplied via --drm) is the one produced
# by the standard simulation setup.  The seed is taken from the name of the
# DRM file, and the sim object is instantiated by giving it that seed.
#
# Usage:
#   extract_drm_positions.py [-p FILE] [-x SCRIPT] SCRIPT.json SEED.pkl
#
# where the TWO required arguments are the script, and the DRM, and optionally:
#  -p FILE -- a filename template (containing a %s) for the output .csv file; the
#             %s hold the seed.  By default, ./%s-position.csv is used.
#  -x SCRIPT -- the given SCRIPT filename is loaded on top of the argument script;
#               if the given SCRIPT name begins with !, it is treated as a
#               json literal rather than a filename
#
# Actual usage example for reference:
#   This command adds the Local/ subdirectory here to PYTHONPATH, and then
#   runs the command.  The Local/ directory must be added because this particular script
#   inherits from a SurveyEnsemble (IPClusterEnsembleJPL2) which is in Local/ here.
#   The Exosims instance should come from your current Python venv.
#     S=HabEx_4m_TSDD_top200DD_52m_5184_20181015_obs
#     PYTHONPATH=Local util/extract_drm_positions.py -p sims/$S/pos/%s-position.csv Scripts/$S.json sims/$S/drm/1351679.pkl
#
# Note: there is a sh driver for this script, and if using the Sandbox, it's easier
# to just use the driver, which enforces the filename conventions.
# See: extract_drm_positions_driver.sh

# history:
#  turmon oct 2018: created, from keepout_path_graphics
#  turmon oct 2023: updated for python3 and current venv practices
#  turmon oct 2023: placed the legacy extract_drm_positions.py in a _2018.py file

from __future__ import print_function
import sys
import argparse
import os
import os.path
import csv
import json
import pickle
from time import time
import EXOSIMS
import EXOSIMS.MissionSim
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord

# unpickling python2/numpy pickles within python3 requires this
PICKLE_ARGS = {} if sys.version_info.major < 3 else {'encoding': 'latin1'}

############################################################
#
# Utility Functions
#
############################################################

def extract_time(tm):
    r'''Unbox a time, if it is boxed, for compatibility with multiple DRMs.'''
    try:
        return tm.value
    except AttributeError:
        return tm


# Container class for loading canned observational data from the DRM
class ObserveInfo(object):
    r"""Observations made by a mission, in order, as loaded from an external DRM pickle."""

    def load_from_spc_file(self, loc, spc):
        r'''Load the DRM, and star/planet info from a "spc" file given as an argument.
        This drm + spc file transfer is compatible the Exosims ipyparallel output.'''
        # these will fail noisily if there is no file present
        print('Loading DRM from', loc)
        self.drm = pickle.load(open(loc, 'rb'), **PICKLE_ARGS)
        # (for now, this is always False)
        if spc:
            # spc file contains a dict with many fields - save them all
            self.spc = pickle.load(open(spc, 'rb'), **PICKLE_ARGS)
            # and then extract the ones we need
            self.snm = self.spc['Name']
            # self.sco = self.spc['coords'] # we have the script, so coords are not needed
            # planet luminosity, scaled for distance, one entry per planet, in SolarLuminosity/AU^2
            self.L_planet = (self.spc['L'][self.spc['plan2star']]/(self.spc['a'].to('au')**2)).value
            # planet radius, in EarthRad
            self.Rp_planet = self.spc['Rp'].value

    def __init__(self, loc, spc):
        # load DRM and Star-Planet info
        self.load_from_spc_file(loc, spc)
        # number of observations
        N = len(self.drm)
        # stars we visited, in DRM order
        self.drm_s_ind = [r['star_ind'] for r in self.drm]
        # visit times
        self.times = np.array([extract_time(r['arrival_time']) for r in self.drm])
        # summaries

    def summary(self):
        s = ('DRM of %d observations, %d stars observed'
                 %
                 (len(self.drm), len(set(self.drm_s_ind))))
        return s


############################################################
#
# Main routine
#
############################################################


def dictify_coords(tag, bodies, r_obs, vr_obs, vd_obs, r_targ, r_body, v_body):
    r'''Turns a group of coordinate sets into entries in a dict, which is returned.

    This is here to get the clutter out of the main routine.'''
    coords = ('x', 'y', 'z')
    s_unit = 'pc' # for stars
    b_unit = 'au' # for solar system bodies and observatory
    v_unit = u.au/u.year
    d = {}
    d.update({'obs%s_%s'   %(tag,c):r_obs      [0,i].to(b_unit).value for i,c in enumerate(coords)})
    d.update({'star%s_%s'  %(tag,c):r_targ     [0,i].to(s_unit).value for i,c in enumerate(coords)})
    d.update({'%s%s_%s'    %(b,tag,c):r_body[b][0,i].to(b_unit).value for i,c in enumerate(coords) for b in bodies})
    d.update({'vrobs%s_%s'  %(tag,c):vr_obs    [0,i].to(v_unit).value for i,c in enumerate(coords)})
    d.update({'vdobs%s_%s'  %(tag,c):vd_obs    [0,i].to(v_unit).value for i,c in enumerate(coords)})
    d.update({'v%s%s_%s'   %(b,tag,c):v_body[b][0,i].to(v_unit).value for i,c in enumerate(coords) for b in ('earth',)})
    return d
                         
def get_char_time(obs):
    r'''Utility function, gets char time from a drm observation.

    The tricky issue here is that some characterizations (with starshade)
    have an embedded char_time key, and others (coronagraph + filters) have
    a list of char_info's with a series of char_times, which must be
    added up since they occur serially.'''
    if 'char_time' in obs:
        # starshade char - has a unitary char_time key
        return obs['char_time']
    else:
        # coronagraph-only char - sum integration windows
        return sum(cslice['char_time'] for cslice in obs['char_info'])

def obs_to_obs_time(obs):
    r'''Return the time taken by the given observation (det or char).

    Returns an astropy Quantity in days.'''
    if 'det_time' in obs:
        rv = obs['det_time']
    elif 'char_time' in obs or 'char_info' in obs:
        rv = get_char_time(obs)
    return rv

def angle_between(left, center, right):
    r'''Angle theta in between three points, left-center-right.'''
    v_left  = left  - center
    v_right = right - center
    u_left  = v_left  / np.linalg.norm(v_left)
    u_right = v_right / np.linalg.norm(v_right)
    cos_theta = np.dot(u_right, u_left)
    # clip to [-1,+1] to suppress numerical issues
    theta = np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))
    return theta

def fd_obs_velocity(Obs, tm):
    # finite-difference velocity for the observatory, based on position
    # for Obs:
    #    -- huge movement in 1d => dt << 1d
    #    -- typical tiniest obs-time ~ 0.01 day => dt << 0.01 day
    dt = 0.0001 * u.d
    time_1 = tm - 0.5 * dt
    time_2 = tm + 0.5 * dt
    r_obs_1 = Obs.orbit(time_1)
    r_obs_2 = Obs.orbit(time_2)
    # place in AU/year
    v_body = (r_obs_2 - r_obs_1) / dt.to(u.year)
    return v_body

def body_position(Obs, bodies, time, velocity=True):
    r'''Body positions vector in heliocentric equatorial frame.

    For related usage, see Observatory.keepout().'''
    # map our body names to Exosims body names
    body_map = dict(sun='Sun',
                    earth='Earth',
                    moon='Moon')
    r_body = {b: Obs.solarSystem_body_position(time, body_map[b]).to('AU')
                  for b in bodies}
    if velocity:
        # finite-difference velocity
        # for Earth:
        #    -- huge movement in 1d => dt << 1d
        #    -- typical tiniest obs-time ~ 0.01 day => dt << 0.01 day
        dt = 0.0001 * u.d
        time_1 = time - 0.5 * dt
        time_2 = time + 0.5 * dt
        r_body_1 = {b: Obs.solarSystem_body_position(time_1, body_map[b]).to('AU')
                    for b in bodies}
        r_body_2 = {b: Obs.solarSystem_body_position(time_2, body_map[b]).to('AU')
                    for b in bodies}
        # place in AU/year
        v_body = {b: (r_body_2[b] - r_body_1[b])/dt.to(u.year) for b in bodies}
    else:
        v_body = None
    return r_body, v_body

def extract_positions(OI, args, xspecs):
    r'''Extract target and observatory positions across the DRM in OI.'''

    # load script, allowing for extra specs
    sim = EXOSIMS.MissionSim.MissionSim(args.script, **xspecs)
    # shorthand for objects
    OS = sim.OpticalSystem
    Obs = sim.Observatory
    TL = sim.TargetList

    ########################################
    # 2023-10: This is no longer being used -- commented out
    #
    # set star coordinates from TL
    #   heliocentric lon/lat
    #   NB: this 
    #OI.xpos = TL.coords.barycentrictrueecliptic.lon.value
    #OI.ypos = TL.coords.barycentrictrueecliptic.lat.value
    #
    # # use (equatorial) ra/dec coordinates vs. (ecliptic) lat/lon coordinates
    # ra_dec_coord = False
    # if ra_dec_coord:
    #     # attribute name for SkyCoord's
    #     xattr_name, yattr_name = 'ra', 'dec'
    # else:
    #     # attribute name for SkyCoord's
    #     xattr_name, yattr_name = 'lon', 'lat'
    ########################################
    
    ########################################
    # Output file
    if not os.path.isdir(os.path.dirname(args.outfile)):
        os.makedirs(os.path.dirname(args.outfile), 0o775)
    try:
        csv_position = open(args.outfile, 'w')
    except IOError:
        print("Error: could not write position-file '%s'" % args.outfile)
        raise
    # determine nouns to track, fields to output, etc.
    coords = ('x', 'y', 'z')
    # solar system bodies to locate
    # bodies = ('sun', 'earth', 'moon')
    bodies = ('sun', 'earth')
    # nouns = all things to report positions/velocities on
    nouns = ('obs', 'vrobs', 'vdobs', 'star') + bodies + ('vearth', )
    # removed 'star_vmag'
    fieldnames = ['timeA_iso', 'timeA_mjd', 'time0_mjd', 'timeH_mjd', 'time1_mjd',
                      'dt_obs', 'obs_mode',
                      'star_name', 'theta_obsA', 'theta_obs0', 'theta_obsH', 'theta_obs1', ]
    # currently using the "Arrival" + "IntTimeStart" + "Halfway" times for coords
    times_used = ['A', '0', 'H', '1']
    for t in times_used:
        fieldnames.extend([f'{x}{t}_{c}' for x in nouns for c in coords])
    # open the CSV writer (ignores extra fields in dict)
    writer = csv.DictWriter(csv_position, fieldnames=fieldnames, extrasaction='ignore')
    writer.writeheader()

    ########################################
    # Set up time-by-time loop
    
    # Choose the first observingMode that is a detection mode
    detMode = list(filter(lambda mode: mode['detectionMode'], OS.observingModes))[0]
    # we are not handling the multi-mode char case now!
    charMode = list(filter(lambda mode: mode['instName'].startswith('spec'), OS.observingModes))[0]

    # these ohTimes are Quantity's in days
    det_ohTime  = detMode[ 'syst']['ohTime']
    char_ohTime = charMode['syst']['ohTime']

    # find a total overhead [qty days]
    det_ohTotal  =  det_ohTime + Obs.settlingTime
    char_ohTotal = char_ohTime + Obs.settlingTime

    print(f"{args.progname}: Modes and overheads:")
    print(f"  Detection: { detMode['instName']} with { detMode['systName']} setup+overhead { det_ohTotal:.3f}")
    print(f"  Char-tion: {charMode['instName']} with {charMode['systName']} setup+overhead {char_ohTotal:.3f}")

    # for progress indication
    system_t0 = time()
    system_dt = 10.0 # seconds
    
    ########################################
    # Loop over times
    
    print('%s: Stepping through observations.' % args.progname)
    Ntime = len(OI.drm)
    for i, obs in enumerate(OI.drm):
        # print progress if needed
        if time() > system_t0 + system_dt:
            sys.stdout.write('  [%3d/%3d]\n' % (i+1, Ntime))
            system_t0 = time()

        # observing mode
        obs_mode = 'Unknown'
        if 'det_info' in obs or 'det_time' in obs:
            obs_mode = 'det'
            obs_ohTotal = det_ohTotal
        if 'char_mode' in obs or 'char_info' in obs:
            obs_mode = 'char'
            obs_ohTotal = char_ohTotal

        # times and time strings
        # time model for one obs is: (arrival) (time for overhead) (time for integration)
        # thus:
        #   timeA = arrival
        #   time0 = start of integration
        #   timeH = halfway point for integration
        #   time1 = end of integration
        arrival_time = extract_time(obs['arrival_time']) # time-from-start
        obs_time = obs_to_obs_time(obs) # time spent observing
        # observing endpoints
        timeA = sim.TimeKeeping.missionStart + arrival_time*u.day
        time0 = timeA + obs_ohTotal
        timeH = time0 + obs_time*0.5
        time1 = time0 + obs_time
        # suppress ms on string output
        if False:
            time0.out_subfmt = 'date_hms'
            time1.out_subfmt = 'date_hms'
        # chop millisec so MS Excel parses a date
        timeA_iso = timeA.iso[:-4]
        time0_iso = time0.iso[:-4]
        timeH_iso = timeH.iso[:-4]
        time1_iso = time1.iso[:-4]
        timeA_mjd = '%.3f' % timeA.mjd
        time0_mjd = '%.3f' % time0.mjd
        timeH_mjd = '%.3f' % timeH.mjd
        time1_mjd = '%.3f' % time1.mjd
        
        sInd = obs['star_ind']

        # keepout calculation yields geometry (koangles are *not* the object angles)
        #kogood0, r_body0, r_targ0, culprit0, koangles0 = Obs.keepout(TL, sInd, time0, detMode, returnExtra=True)
        #kogoodH, r_bodyH, r_targH, culpritH, koanglesH = Obs.keepout(TL, sInd, timeH, detMode, returnExtra=True)
        #kogood1, r_body1, r_targ1, culprit1, koangles1 = Obs.keepout(TL, sInd, time1, detMode, returnExtra=True)
        
        r_bodyA, v_bodyA = body_position(Obs, bodies, timeA)
        r_body0, v_body0 = body_position(Obs, bodies, time0)
        r_bodyH, v_bodyH = body_position(Obs, bodies, timeH)
        r_body1, v_body1 = body_position(Obs, bodies, time1)

        # target star positions vector in heliocentric equatorial frame
        r_targA = TL.starprop(sInd, timeA)
        r_targ0 = TL.starprop(sInd, time0)
        r_targH = TL.starprop(sInd, timeH)
        r_targ1 = TL.starprop(sInd, time1)

        # observatory position in HE coordinates
        r_obsA = Obs.orbit(timeA)
        r_obs0 = Obs.orbit(time0)
        r_obsH = Obs.orbit(timeH)
        r_obs1 = Obs.orbit(time1)

        # observatory velocity in rotating frame [AU/yr]
        vr_obsA = Obs.haloVelocity(timeA)
        vr_obs0 = Obs.haloVelocity(time0)
        vr_obsH = Obs.haloVelocity(timeH)
        vr_obs1 = Obs.haloVelocity(time1)

        # observatory velocity in heliocentric equatorial frame [AU/yr]
        vd_obsA = fd_obs_velocity(Obs, timeA)
        vd_obs0 = fd_obs_velocity(Obs, time0)
        vd_obsH = fd_obs_velocity(Obs, timeH)
        vd_obs1 = fd_obs_velocity(Obs, time1)

        # Calculate desired (ra,dec or lon,lat) coordinates of visible targets, kept-out targets, and bright bodies
        # Recall that r_targ (the star positions) and r_body (the solar-system
        # body positions) are both in HE (heliocentric equatorial) coordinates.
        #targH = SkyCoord(r_targH[:,0],   r_targH[:,1],   r_targH[:,2],   representation='cartesian')
        #bodyH = SkyCoord(r_bodyH[:,0,0], r_bodyH[:,0,1], r_bodyH[:,0,2], representation='cartesian')
        # extract galactic lat/lon
        #targH_conv = targH.barycentrictrueecliptic
        #bodyH_conv = bodyH.barycentrictrueecliptic

        # find angle, sun-obs-target (pass as 1x3 tuples)
        # NB: sun = index #0
        theta_obsA = angle_between(r_bodyA['sun'][0,:], r_obsA[0,:], r_targA[0,:])
        theta_obs0 = angle_between(r_body0['sun'][0,:], r_obs0[0,:], r_targ0[0,:])
        theta_obsH = angle_between(r_bodyH['sun'][0,:], r_obsH[0,:], r_targH[0,:])
        theta_obs1 = angle_between(r_body1['sun'][0,:], r_obs1[0,:], r_targ1[0,:])

        # convert to dict for ease-of-output
        info = {}
        info.update(dictify_coords('A', bodies, r_obsA, vr_obsA, vd_obsA, r_targA, r_bodyA, v_bodyA))
        info.update(dictify_coords('0', bodies, r_obs0, vr_obs0, vd_obs0, r_targ0, r_body0, v_body0))
        info.update(dictify_coords('H', bodies, r_obsH, vr_obsH, vd_obsH, r_targH, r_bodyH, v_bodyH))
        info.update(dictify_coords('1', bodies, r_obs1, vr_obs1, vd_obs1, r_targ1, r_body1, v_body1))
        # previously present
        #   star_vmag=TL.Vmag[sInd],
        info.update(
            star_name=TL.Name[sInd],
            timeA_iso=timeA_iso,
            timeA_mjd=timeA_mjd,
            time0_mjd=time0_mjd,
            timeH_mjd=timeH_mjd,
            time1_mjd=time1_mjd,
            dt_obs=obs_time.value,
            obs_mode=obs_mode,
            theta_obsA=theta_obsA.to('deg').value,
            theta_obs0=theta_obs0.to('deg').value,
            theta_obsH=theta_obsH.to('deg').value,
            theta_obs1=theta_obs1.to('deg').value,
            )
        writer.writerow(info)

    # (end loop over observations)
    print('%s: Done.' % args.progname)
    csv_position.close()

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract EXOSIMS observing geometry.")
    parser.add_argument('script', metavar='SCRIPT', help='json script')
    parser.add_argument('drm', help='file where pickle of observation sequence lives',
                      metavar='DRM')
    #parser.add_argument('--spc', help='name of star-planet-characteristics (spc) file',
    #                  dest='spc', metavar='FILE', default='')
    parser.add_argument('-x', '--xspecs', help='extra json spec-file to add on top of SCRIPT',
                      dest='xspecs', metavar='FILE', default='')
    parser.add_argument('-p', '--position', help='file template for position info (contains %s)',
                      dest='out_position', metavar='FILE', default='')
    
    # set umask in hopes that files/dirs will be group-writable
    os.umask(0o002)
    # process arguments
    args = parser.parse_args()
    args.progname = os.path.basename(sys.argv[0])
    args.spc = None # unused

    # position csv file
    drm_seed = os.path.splitext(os.path.basename(args.drm))[0]
    if args.out_position:
        if '%s' not in args.out_position:
            raise ValueError('Output file template must contain "%s" for seed.')
        args.outfile = args.out_position % drm_seed
    else:
        args.outfile = './%s-position.csv' % drm_seed

    # load extra specs, if any
    # argument is a filename, or if argument begins with "!", a literal.
    if args.xspecs and not args.xspecs.startswith('!'):
        assert os.path.isfile(args.xspecs), "%s is not a file." % args.xspecs
        try:
            script = open(args.xspecs).read()
            xspecs = json.loads(script)
        except ValueError as err:
            print("Error: Input xspec file `%s' improperly formatted." % (args.xspecs,))
            print("Error: JSON error was: ", err)
            # re-raise here to suppress the rest of the backtrace.
            # it is only confusing details about the bowels of json.loads()
            raise ValueError(err)
        except:
            print("Error: %s", (sys.exc_info()[0],))
            raise
    elif args.xspecs and args.xspecs.startswith('!'):
        xspecs = json.loads(args.xspecs[1:])
    else:
        xspecs = {}

    print('%s: setting standalone mode with seed = %s' % (args.progname, drm_seed))
    # ensure simulation seed is reset
    xspecs['seed'] = int(drm_seed)
    # put into stand-alone mode (no ipyparallel)
    xspecs['ensemble_mode'] = 'standalone'

    # load info about observations (DRM) from given location
    OI = ObserveInfo(args.drm, args.spc)
    print('%s: %s' % (args.progname, OI.summary()))

    extract_positions(OI, args, xspecs)

