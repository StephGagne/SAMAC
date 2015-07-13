from pylab import *
import numpy as np
import copy
import scipy.stats.stats as st

import samac

#########################################################################        
############################   avsizedist   #############################
#########################################################################        
def avsizedist(CL,prof='belowcloud',scan=0,inst='PCASP',filler=0):
    """ Returning method. This method returns the average size distribution for a given instrument and scan type and number. The returned average size distribution is an array with the mean concentration of particles in the first row and the size of the particles in the second. 
    Example: R=CloudObj.avsizedist(prof='abovecloud',scan=1,inst='PCASP') would return in R the average size distribution for the above cloud #1 for the PCASP.
    The defaults are prof='belowcloud',scan=0,inst='PCASP', filler=0
    filler: filling with a given time step (useful for 2d type data) =1 ON, =0 OFF. Time step can be adjusted in the call of the method filler
    Special profile "belowvert" can also be used, this we return the average size distribution below cloud, where a vertical profile was measured."""

    cancel=0
    instn=[i for i,x in enumerate(CL.sd) if x["Distname"].upper()==inst.upper()]
    if len(instn)<=0: cancel=1; raise ValueError("This instrument name has not been found (inst=%s)." % inst)
    if filler==1: sd=samac.filler(CL,CL.times["cloud"][0],inst=inst)
    else: sd=copy.deepcopy(CL.sd[instn[0]])
    if type(sd["time"])==float: cancel=1
    elif len(sd["time"])==0: cancel=1
    if cancel==1: avsd=[]; return avsd
    if prof=='belowvert':
        if len(sd["time"])<=1:      # interpolation cannot take place if there is not at least 2 points to interpolate from
            print("[avsizedist] Only 1 point available: special profile belowvert cannot be used.")
            avsd=[]; return avsd
        ctime=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctime=ctime[0]
        basictime=CL.data[ctime]
        if sum(basictime==sd["time"])==len(basictime): pass      # the instrument share the timeline of the basic data from which belowvert is calculated
        else:       # the instrument does not share the basic timeline and must be interpolated
            if sum(sd["data"])>0: M=np.min(sd["data"][sd["data"]>0.0])    # minimum non-zero concentration 
            else: M=0
            f=interpolate.RectBivariateSpline(sd["bins"],sd["time"],sd["data"])
            sd["data"]=f(sd["bins"],basictime);  sd["time"]=basictime;
            sd["data"][sd["data"]<0.01*M]=0.0        # all interpolated data smaller than 1% of the minimum.
        R=belowprof(CL)
        try:
            if np.size(R)>0:
                Tix=R[scan]['areaindex']
                concs=sd["data"][:,nonzero(Tix)[0]]
            else: print("[avsizedist] There is no data below cloud under vertical profile #%d." % scan)
        except: print("[avsizedist] The instrument number (%d) may be wrong, or no data corresponds to the profile." % instn[0])            
    else:
        try: 
            T=CL.times[prof][scan]
        except: print("[avsizedist] The profile %s [%d] doesn't exist." % (prof,scan))
        try:
            concs=sd["data"][:,nonzero((sd["time"]>=T[0])*(sd["time"]<=T[1]))[0]]
        except: print("[avsizedist] The instrument number (%d) may be wrong, or no data corresponds to the profile." % instn[0])
    if "concs" in locals():
        if len(sd["time"])==1: avsd=[concs,sd["bins"]]
        else: avsd=np.array([st.nanmean(concs,axis=1), sd["bins"]])
        return avsd
    else: 
        print("[avsizedist] Failed to return the averaged size distribution.")
        avsd=[]; return avsd
    

