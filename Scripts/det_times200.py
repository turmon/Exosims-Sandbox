import os.path,EXOSIMS,EXOSIMS.MissionSim
import numpy as np
import astropy.units as u

scriptfile = '/proj/exep/rhonda/Sandbox/HabEx/Scripts/HabEx_4m_TS_dmag25p5_20180321.json'
sim = EXOSIMS.MissionSim.MissionSim(scriptfile,nopar=True)
SU=sim.SimulatedUniverse
TK=sim.TimeKeeping
TL=sim.TargetList
SS=sim.SurveySimulation
#create array of startimes
slewTimes = np.zeros(TL.nStars)*u.d
startTimes = slewTimes + TK.currentTimeAbs
sInds=SU.sInds
intTimes = np.zeros(TL.nStars)*u.d
intTimes[sInds]=sim.SurveySimulation.calc_targ_intTime(sInds,startTimes[sInds],sim.OpticalSystem.observingModes[0])
starfastInd=np.where(np.cumsum(np.sort(intTimes[sInds])).value<=0)[-1][-1]
star30Ind=np.where(np.cumsum(np.sort(intTimes[sInds])).value<=30)[-1][-1]
sortedintTimes = np.sort(intTimes[sInds])
sortedinds = np.argsort(intTimes[sInds])
TL.dist[sortedinds[starfastInd:star30Ind]]
TL.L[sortedinds[starfastInd:star30Ind]]
SS.dMagint[sortedinds[starfastInd:star30Ind]]
SS.WAint[sortedinds[starfastInd:star30Ind]]
TL.comp0[sortedinds[starfastInd:star30Ind]]

fZ = sim.ZodiacalLight.fZ(sim.Observatory, TL, sInds, startTimes, sim.OpticalSystem.observingModes[2])
fEZ = sim.ZodiacalLight.fEZ0
charTimes = sim.OpticalSystem.calc_intTime(TL, sInds, fZ, fEZ, SS.dMagint, SS.WAint, sim.OpticalSystem.observingModes[2])
C_p, C_b, C_sp = sim.OpticalSystem.Cp_Cb_Csp(TL, sInds, fZ, fEZ, SS.dMagint, SS.WAint, sim.OpticalSystem.observingModes[2])
SNR=10
char_intTime = np.true_divide(SNR**2*C_b, (C_p**2)) #no C_sp for occulter
