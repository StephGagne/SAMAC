################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####


from pylab import *
import numpy as np

from scipy import interpolate

import samac

#########################################################################        
###############################  avCDNC  ################################
#########################################################################

def avCDNC(CL,abovesize=nan,uppersize=nan,inst="FSSP96",scan=0):
    """ This method returns the average cloud droplet concentration in a vertical scan weighted according to the LWC.
    The default options are avCDNC(abovesize=nan,uppersize=nan,inst="FSSP96",scan=0). This method handles the extradata module for lwc.
    "abovesize" is the size above which the concentration is calculated, "uppersize" is the upper size limit. By default the entire available spectrum is integrated. 
    "inst" is the name of the instrument that should be used. Generally an FSSP100.
    "scan" is the scan number, if there are more than one vertical scan.
    The units are the same as for the instrument's concentration (as long as the altitude is in meters)."""
    if len(CL.times["verticloud"])<=scan:
        print("[avCDNC] No vertical scan for this cloud.")
        return nan
    else:
        ts=CL.times["verticloud"][scan]       # times of the vertical scan
        pos=[i for i,x in enumerate(CL.sd) if inst.lower() in x["Distname"].lower()]; 
        if len(pos)<=0: cancel=1; raise ValueError("This instrument name has not been found (inst=%s)." % inst)
        else: pos=pos[0]
        calt=[i for i,x in enumerate(CL.dttl) if x == 'altitude']; calt=calt[0]
        ctim=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctim=ctim[0]
        alt=CL.data[calt]
        tim=CL.data[ctim]
        clwc=[i for i,x in enumerate(CL.dttl) if x == 'LWC']; 
        if len(clwc)>1: print("[avCDNC] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=ts     # LWC found in the basic data
        elif len(clwc)==1: lwc=CL.data[clwc[0]]; lwctime=tim; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(CL.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[avCDNC] No LWC (or multiple) was found in the basic data or the extra data."); return nan
        
        
        timtc=CL.sd[pos]["time"]       # time corresponding to the concentration tc
        # find shortest common time stamp
        ta1=np.max([lwctime[0],tim[0],timtc[0],ts[0]]); ta2=np.min([lwctime[-1],tim[-1],timtc[-1],ts[1]]); 
        taix=[nonzero((timtc>=ta1)*(timtc<=ta2))[0]]; 
        usedtime=timtc[taix]
        if len(usedtime)==0: WA=nan    # weighted average is not calculable (no data for the vertical scan)
        else:
            tc=samac.dNdlogDp2N(CL.sd[pos],abovesize,uppersize)[taix] # total concentration between above- and upper-size
            test=tim.__array__()==timtc
            if len(tim)==len(timtc): 
                if (tim.__array__()==timtc).all()==True: alt=alt[taix]
                else:
                    print("[avCDNC] nonequal times (although same length): altitude was interpolated.")
                    # interpolating on size distribution data.
                    falt=interpolate.interp1d(tim,alt,kind='linear')   
                    alt=falt(usedtime); 
            else:
                print("[avCDNC] nonequal times: altitude was interpolated.")
                # interpolating on size distribution data.
                falt=interpolate.interp1d(tim,alt,kind='linear')   
                alt=falt(usedtime); 
            if len(lwctime)==len(timtc):
                if (lwctime.__array__()==timtc).all()==True: lwc=lwc[taix]
            else: flwc=interpolate.interp1d(lwctime,lwc,kind='linear'); lwc=flwc(usedtime)
            # integrating
            dalt=np.diff(alt);
            lwc=lwc[1:];
            tc=tc[1:];
            WA=abs(sum(tc*lwc*dalt))/(samac.lwp(CL)[scan][0])        # if inst conc is cm-3     (cm-3 lwc*meters / lwc*meters = cm-3)
        return WA
        

