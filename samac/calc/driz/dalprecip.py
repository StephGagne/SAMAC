
from pylab import *
import numpy as np

from scipy import interpolate

import samac


#########################################################################        
#############################  totalprecip  #############################
#########################################################################

def totalprecip(CL,abovesize=100,scan=0,filler=0):
    """ This method returns the integrated (with respect to altitude) amount of precipitation drops (number) in a vertical profile (in #/m2 OR in # m L-1, need to edit code to change) above the size given in options. A 2d instrument is used.
    By default, this in done for the first scan (#0), other scans can be accessed using options.
    CloudObj.totalprecip(abovesize=100, scan=0) size in microns, is the default."""
    #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
    if len(CL.times["verticloud"])<=scan:
        print("[totalprecip] No vertical scan for this cloud.")
        return nan
    else:
        ts=CL.times["verticloud"][scan]       # times of the vertical scan
        # columns & positions
        calt=[i for i,x in enumerate(CL.dttl) if x == 'altitude']; calt=calt[0]
        ctim=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctim=ctim[0]
        pos2d=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]; 
        if len(pos2d)==0: pos2d=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()];
        try: pos2d=pos2d[0]
        except: print("[totalprecip] No precipitation size distribution found.")
        # loading the relevant distribution
        if filler==1:
            fill2d=samac.filler(CL,ts=ts)
        elif filler==0: 
            try:
                fill2d=dict()
                tix=np.nonzero((CL.sd[pos2d]["time"]>=ts[0])*(CL.sd[pos2d]["time"]<=ts[1]))[0]
                if len(tix)==0: return 0.00
                fill2d["data"]=CL.sd[pos2d]["data"][:,tix]; fill2d["time"]=CL.sd[pos2d]["time"][tix];
                fill2d["bins"]=CL.sd[pos2d]["bins"]; fill2d["total"]=CL.sd[pos2d]["total"][tix];
            except: print("[totalprecip] Data not found"); return nan
        else:
            print("[totalprecip] filler must be 1 or 0.")
            return nan
        # controlling for empty distributions
        if type(fill2d["time"])==float: Int=nan     # data does not exist
        elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
        
        else:   # proceeding!
            prec=samac.dNdlogDp2N(fill2d,abovesize,nan)
            alt=CL.data[calt]
            tim=CL.data[ctim]
            
            # interpolating altitude on 2d data.
            f=interpolate.interp1d(tim,alt,kind='linear')   
            alt2=f(fill2d["time"])
            # integrating
            dalt=np.diff(alt2);
            prec=prec[1:];
            Int=abs(sum(prec*dalt))    # in litres-1 * m => 1L = 0.001m3 ==> 1000 #/m2 = 1 #m/L
            Int=Int*1000    # now in #/m2
        return Int       


#########################################################################        
###########################  totalvolprecip  ############################
#########################################################################

def totalvolprecip(CL,abovesize=100,scan=0,filler=0):
    """ This method returns the integrated (with respect to altitude) amount of precipitation drops (volume) in a vertical profile (in meters OR in um3 m L-1, need to edit code to change) above the size given in options. A 2d instrument is used. Note that it is assume that the bins sizes are diameters and not radii.
    scan: default scan=0. Access to different vertical scans in a cloud.
    filler: 1 if ON and 0 if OFF, fills missed time steps with zeros, meant to fill 2d data. Default: 0, OFF. filler works with default 2dc instrument and time steps of 10 secs. Feel free to modify as needed.
    CloudObj.totalvolprecip(abovesize=100, scan=0) size in microns, is the default."""
    #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
    if len(CL.times["verticloud"])<=scan:
        print("[totalvolprecip] No vertical scan for this cloud.")
        return nan
    else:
        ts=CL.times["verticloud"][scan]       # times of the vertical scan
        # columns & positions
        calt=[i for i,x in enumerate(CL.dttl) if x == 'altitude']; calt=calt[0]
        ctim=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctim=ctim[0]
        pos2d=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]; 
        if len(pos2d)==0: pos2d=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()];
        try: pos2d=pos2d[0]
        except: print("[totalvolprecip] No precipitation size distribution found.")
        # loading the relevant distribution
        if filler==1:
            fill2d=samac.filler(Cl,ts=ts)
        elif filler==0: 
            try:
                fill2d=dict()
                tix=np.nonzero((CL.sd[pos2d]["time"]>=ts[0])*(CL.sd[pos2d]["time"]<=ts[1]))[0]
                if len(tix)==0: return 0.00
                fill2d["data"]=CL.sd[pos2d]["data"][:,tix]; fill2d["time"]=CL.sd[pos2d]["time"][tix];
                fill2d["bins"]=CL.sd[pos2d]["bins"]; fill2d["total"]=CL.sd[pos2d]["total"][tix];
            except: print("[totalvolprecip] Data not found"); return nan
        else:
            print("[totalvolprecip] filler must be 1 or 0.")
            return nan
        # controlling for empty distributions
        if type(fill2d["time"])==float: Int=nan     # data does not exist
        elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
        else:       # proceeding!
            ix=fill2d["bins"]>=abovesize
            prec=np.sum(samac.dNdlogDp2dN(fill2d)[ix].transpose()*(pi*(fill2d["bins"][ix])**3)/6,1)        # in #/L * um3, assuming the bins are diameters.
            alt=CL.data[calt]
            tim=CL.data[ctim]
            # interpolating altitude on 2d data.
            f=interpolate.interp1d(tim,alt,kind='linear')   
            alt2=f(fill2d["time"])
            # integrating
            dalt=np.diff(alt2);
            prec=prec[1:];
            Int=abs(sum(prec*dalt))    # in litres-1 * um3 * m => 1L = 0.001m3; 1um=10^-6m ==> 1e-15 m = 1 um3/L * m
            Int=Int*1e-15    # now in meters (integrated (w/ resp. altitude) volume of water per volume of air)
        return Int  

#########################################################################        
#############################  precipconc  ##############################
#########################################################################

def precipconc(CL,abovesize=100,scan=0,minlwc=0.01,filler=0):
    """ This method returns the integrated (with respect to altitude) amount of precipitation drops (number) in a vertical profile in-cloud (lwc>minlwc (default=0.01 g/m3)) for particle greater than "abovesize" and divide by the depth of the cloud. A 2d probe is used. Units are number per litre (if the 2d used has units per litre).
    By default, this in done for the first scan (#0), other scans can be accessed using options.
    filler: 1 if ON and 0 if OFF, fills missed time steps with zeros, meant to fill 2d data. Default: 0, OFF. filler works with default 2dc instrument and time steps of 10 secs. Feel free to modify as needed.
    CloudObj.precipconc(abovesize=100, scan=0, minlwc=0.01) size in microns, lwc in g/m3 is the default."""
    #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
    if len(CL.times["verticloud"])<=scan:
        print("[precipconc] No vertical scan for this cloud.")
        return nan
    else:
        ts=CL.times["verticloud"][scan]       # times of the vertical scan
        # columns & positions
        calt=[i for i,x in enumerate(CL.dttl) if x == 'altitude']; calt=calt[0]
        ctim=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctim=ctim[0]
        alt=CL.data[calt]; tim=CL.data[ctim];
        clwc=[i for i,x in enumerate(CL.dttl) if x == 'LWC']; 
        if len(clwc)>1: print("[precipconc] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=ts     # LWC found in the basic data
        elif len(clwc)==1: lwc=CL.data[clwc[0]]; lwctime=tim; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(CL.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[precipconc] No LWC (or multiple) was found in the basic data or the extra data."); return nan
        # loading the distribution
        if filler==1:
            fill2d=samac.filler(CL,ts=ts)
        elif filler==0: 
            try:
                pos2d=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]; 
                if len(pos2d)==0: pos2d=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()];
                try: pos2d=pos2d[0]
                except: print("[precipconc] No precipitation size distribution found.")
                fill2d=dict()
                tix=np.nonzero((CL.sd[pos2d]["time"]>=ts[0])*(CL.sd[pos2d]["time"]<=ts[1]))[0]
                if len(tix)==0: return 0.00
                fill2d["data"]=CL.sd[pos2d]["data"][:,tix]; fill2d["time"]=CL.sd[pos2d]["time"][tix];
                fill2d["bins"]=CL.sd[pos2d]["bins"]; fill2d["total"]=CL.sd[pos2d]["total"][tix];
            except: print("[precipconc] Data not found"); return nan
        else:
            print("[precipconc] filler must be 1 or 0.")
            return nan
        # controlling for empty distributions
        if type(fill2d["time"])==float: Int=nan     # data does not exist
        elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
        else:       # proceeding!
            prec=samac.dNdlogDp2N(fill2d,abovesize,nan)
            # interpolating lwc on basic data
            tc1=np.max([lwctime[0],tim[0]]); tc2=np.min([lwctime[-1],tim[-1]]);
            tcix=[nonzero((tim>=tc1)*(tim<=tc2))[0]]
            flwc=interpolate.interp1d(lwctime,lwc,kind='linear')
            lwc=flwc(tim[tcix])
            cloudtop=np.max(alt[tcix][np.nonzero((tim[tcix]>=ts[0])*(tim[tcix]<=ts[1])*(lwc>minlwc))])
            cloudbase=np.min(alt[tcix][np.nonzero((tim[tcix]>=ts[0])*(tim[tcix]<=ts[1])*(lwc>minlwc))])
            # interpolating altitude on 2d data.
            f=interpolate.interp1d(tim,alt,kind='linear')
            alt2=f(fill2d["time"]);
            # integrating
            alti=alt2[(alt2>cloudbase)*(alt2<cloudtop)]; 
            dalt=np.diff(alti);
            dalt=hstack((dalt,dalt[-1]))
            preci=prec[(alt2>cloudbase)*(alt2<cloudtop)]; 
            Int=abs(sum(preci*dalt))    # in l-1 * m => m l-1
            Depth=cloudtop-cloudbase; 
            #print "Depth used:", Depth, "Difference in Depth:", Depth-(mean(CL.props["height"][scan][2:4])-mean(CL.props["height"][scan][0:2]))
            Int=Int/Depth   # in units l-1
        return Int


