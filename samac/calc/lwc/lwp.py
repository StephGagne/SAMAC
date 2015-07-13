
from pylab import *
import numpy as np
from scipy import interpolate

#########################################################################        
#################################  lwp  #################################
#########################################################################                
def lwp(CL):
    """ This property calculates the liquid water path from all vertical scans. One LWP per vertical scan will be returned in an array. 
        The LWP for vertical scan #0 can be accessed with CloudObj.lwp[0][0], #1 with CloudObj.lwp[1][0]. """
    # lwc nan and wrong number only are cleared. Negative and small numbers are kept ==> the integral should cancel the noise out.
    LWP=[]
    ct=[i for i,x in enumerate(CL.dttl) if x=='time'][0]
    Calt=[i for i,x in enumerate(CL.dttl) if x=='altitude']
    if len(Calt)==1: alt=CL.data[Calt[0]];
    else: print("[lwp] No or multiple altitudes found."); return LWP
    t=CL.data[ct]
    Clwc=[i for i,x in enumerate(CL.dttl) if x=='LWC']
    if len(Clwc)>1: print("[lwp] Multiple LWC found."); return LWP
    elif len(Clwc)==1: lwc=CL.data[Clwc[0]]; ta=t
    elif len(Clwc)==0:
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
        if len(posx)==1: 
            lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            # adapting for too short data for interpolation
            if len(lwctime)<2: lwc=np.ones((2,))*NaN; lwctime=np.array([CL.times["cloud"][0][0],CL.times["cloud"][0][1]]);
            # adapting the time vector to a common time vector
            ta1=np.max([t[0],lwctime[0]]); ta2=np.min([t[-1],lwctime[-1]]);
            ta=t[nonzero((t>=ta1)*(t<=ta2))[0]]
            alt=alt[nonzero((t>=ta1)*(t<=ta2))[0]]
            f=interpolate.interp1d(lwctime,lwc,kind='linear')
            if 'ma' not in str(type(lwc)).lower(): lwc=np.ma.array(lwc,mask=False)
            fma=interpolate.interp1d(lwctime,lwc.mask,kind='linear')    # interpolating the mask
            lwc=np.ma.array(f(ta),mask=fma(ta))
        else: print("[lwp] No LWC was found in the basic data or multiple were found in the extra data."); return LWP
    for k,j in enumerate(CL.times["verticloud"]):
        alti=alt[np.nonzero((ta>=j[0])*(ta<=j[1]))]; dalt=np.diff(alti);
        lwci=lwc[np.nonzero((ta>=j[0])*(ta<=j[1]))]; lwci=lwci[1:];
        if sum(lwci)<0: LWP.append(nan); print("[lwp]: Too much noise, not enough signal in LWC. No LWP calculated for scan %d." %(k))
        else: LWP.append(abs(sum(lwci*dalt)))
    LWP=np.reshape(np.array(LWP),(-1,1))
    return LWP

