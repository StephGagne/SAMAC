################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####

from pylab import *
import numpy as np
import samac

#########################################################################        
#############################   overview   ##############################
#########################################################################
def overview(CL,interact=0,Rtime=None):
    """ This method plots the cloud object's Altitude, Latitude, Longitude, Particle concentration, Droplet concentration, Liquid Water Content and Precipitation drop concentration as a function of time for the whole cloud period. Extradata module is fully handled.
        Use CloudObj.overview() for a non-interactive plot;
        Use CloudObj.overview(1) or CloudObj.overview(interact=1) for an interactive plot where you can magnify, induce offsets, and zoom.
        Use CloudObj.overview(Rtime=1) to display Hour:Minutes instead of fraction of day on the x-axis.
        The color bar at the bottom of the upper plot matches the color coding of mapcloud.
        The color bars at the top of the bottom plot correspond to the different legs as follows:
            Magenta: verticloud
            Yellow:  abovecloud
            Cyan:    belowcloud
            Black:   horicloud
        Note that they can be overlaying each other."""
    from plotcloud import plotcloud2
    
    # time
    pos=[i for i,x in enumerate(CL.dttl) if x == 'time']
    if len(pos)==1: tim=CL.data[pos[0]];
    else: print("[overview] No time (or multiple) was found. Forget about this plot! It won't happen unless you have (a) time")
    # altitude
    pos=[i for i,x in enumerate(CL.dttl) if x == 'altitude']
    if len(pos)==1: alt=CL.data[pos[0]]
    else: print("[overview] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude"); alt=NaN;
    # lat
    pos=[i for i,x in enumerate(CL.dttl) if x == 'latitude']
    if len(pos)==1: lat=CL.data[pos[0]]
    else: print("[overview] Latitude not found"); lat=NaN;
    # lon
    pos=[i for i,x in enumerate(CL.dttl) if x == 'longitude']
    if len(pos)==1: lon=CL.data[pos[0]]
    else: print("[overview] Longitude not found"); lon=NaN
    # LWC
    pos=[i for i,x in enumerate(CL.dttl) if x.upper() == 'LWC']
    if len(pos)>1: print("[overview] Multiple LWC found."); lwc=np.ones(len(tim))*NaN; lwctime=tim; lwcunits=''     # LWC found in the basic data
    elif len(pos)==1: lwc=CL.data[pos[0]]; lwctime=tim; lwcunits=CL.dunit[pos[0]]
    else:       # looking for LWC in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
        if len(posx)==1: 
            lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
            lwcunits=CL.extraunit[posx[0][0]][posx[0][1]]
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[overview] No LWC (or multiple) was/were found in the extra data."); lwc=np.ones(len(tim))*NaN; lwctime=tim; lwcunits=''
    ### Size Distributions ###
    TakenSDs=list()     # names list for size distributions, passed to set labels in plot.
    # Aerosol particles 
    As=[x for x in CL.sd if 'A' in x["sdtype"].upper()] 
    if len(As)==0:
        for i in range(len(CL.dttl)): print("%d- %s (%s)" % (i, CL.dttl[i], CL.dunit[i]))
        Q=raw_input("No aerosol sdtype found. Is there a column in the basic data corresponding to the particle number concentration? [n]/#  ")
        try: 
            Q=int(Q); ptc=CL.data[Q]; ptctime=tim; TakenSDs.append([CL.dttl[Q],CL.dunit[Q]])
        except: 
            for i,t in enumerate(CL.extrattl):
                print("%d- %s" %(i,t)); 
            Q=raw_input("Is the total aerosol concentration data in one of these lists? [n]/#  ")
            try:
                Q=int(Q); 
                for j,tt in enumerate(t):
                    print("%d- %s" %(j,tt)); 
                Qj=raw_input("Choose the position of the aerosol number concentration data. [n]/#  ")
                try: 
                    Qj=int(Qj); 
                    ptc=CL.extradata[Q][Qj]; TakenSDs.append([CL.extrattl[Q][Qj],CL.extraunit[Q][Qj]]);
                    k=[k for k,x in enumerate(CL.extrattl[Q]) if x.lower() == 'time'][0]; ptctime=CL.extradata[Q][k]; 
                except: ptc=np.ones(len(tim))*NaN; ptctime=tim; TakenSDs.append(['','']);
            except: ptc=np.ones(len(tim))*NaN; ptctime=tim; TakenSDs.append(['','']);
    elif len(As)==1: ptc=As[0]["total"]; ptctime=As[0]["time"]; TakenSDs.append([As[0]["Distname"],As[0]["units"]])
    elif len(As)>1: 
        As1=[x for x in CL.sd if 'A1' in x["sdtype"].upper()]  # looking for primary distribution
        if len(As1)==1: ptc=As1[0]["total"]; ptctime=As1[0]["time"]; TakenSDs.append([As1[0]["Distname"],As1[0]["units"]])
        else: 
            print("Aerosol distributions:")
            for i,x in enumerate(As): print("%d- %s (sdtype=%s)" % (i,x["Distname"],x["sdtype"]))
            Q=raw_input("Which is your primary Aerosol distribution (consider editing sdtype if this is permanent)? [n]/#  ")
            try: 
                Q=int(Q); ptc=As[Q]["total"]; ptctime=As[Q]["time"]; TakenSDs.append([As[Q]["Distname"],As[Q]["units"]])
            except: ptc=np.ones(len(tim))*NaN; ptctime=tim; print("[overview] Aerosol distribution not found.");  TakenSDs.append(['',''])
    else: print("[overview] impossible aerosols")
    # Cloud droplets 
    Cs=[x for x in CL.sd if 'C' in x["sdtype"].upper()] 
    if len(Cs)==0:
        for i in range(len(CL.dttl)): print("%d- %s (%s)" % (i, CL.dttl[i], CL.dunit[i]))
        Q=raw_input("No cloud droplet sdtype found. Is there a column in the basic data corresponding to the particle number concentration? [n]/#  ")
        try: 
            Q=int(Q); cdnc=CL.data[Q]; cdnctime=tim;  TakenSDs.append([CL.dttl[Q],CL.dunit[Q]])
        except: 
            for i,t in enumerate(CL.extrattl):
                print("%d- %s" %(i,t)); 
            Q=raw_input("Is the total cloud droplet number concentration data in one of these lists? [n]/#  ")
            try:
                Q=int(Q); 
                for j,tt in enumerate(t):
                    print("%d- %s" %(j,tt)); 
                Qj=raw_input("Choose the position of the cloud droplet number concentration data. [n]/#  ")
                try: 
                    Qj=int(Qj); 
                    cdnc=CL.extradata[Q][Qj]; TakenSDs.append([CL.extrattl[Q][Qj],CL.extraunit[Q][Qj]]);
                    k=[k for k,x in enumerate(CL.extrattl[Q]) if x.lower() == 'time'][0]; cdnctime=CL.extradata[Q][k]; 
                except: cdnc=np.ones(len(tim))*NaN; cdnctime=tim; TakenSDs.append(['','']);
            except: cdnc=np.ones(len(tim))*NaN; cdnctime=tim; TakenSDs.append(['','']);
    elif len(Cs)==1: cdnc=Cs[0]["total"]; cdnctime=Cs[0]["time"]; TakenSDs.append([Cs[0]["Distname"],Cs[0]["units"]])
    elif len(Cs)>1: 
        Cs1=[x for x in CL.sd if 'C1' in x["sdtype"].upper()]  # looking for primary distribution
        if len(Cs1)==1: cdnc=Cs1[0]["total"]; cdnctime=Cs1[0]["time"]; TakenSDs.append([Cs1[0]["Distname"],Cs1[0]["units"]])
        else: 
            print("Cloud droplet distributions:")
            for i,x in enumerate(Cs): print("%d- %s (sdtype=%s)" % (i,x["Distname"],x["sdtype"]))
            Q=raw_input("Which is your primary Cloud droplet distribution (consider editing sdtype if this is permanent)? [n]/#  ")
            try: 
                Q=int(Q); cdnc=Cs[Q]["total"]; cdnctime=Cs[Q]["time"]; TakenSDs.append([Cs[Q]["Distname"],Cs[Q]["units"]])
            except: cdnc=np.ones(len(tim))*NaN; cdnctime=tim; print("[overview] Cloud droplet distribution not found."); TakenSDs.append(['',''])
    else: print("[overview] impossible cloud droplets")
    # Precipitation
    pr=dict()
    p=[x for x in CL.sd if 'P' in x["sdtype"].upper()]
    if len(p)==0:
        # >100 microns
        for i in range(len(CL.dttl)): print("%d- %s (%s)" % (i, CL.dttl[i], CL.dunit[i]))
        Q1=raw_input("No precipitation sdtype found. Is there a column in the basic data corresponding to the precipitation number concentration >100 microns? [n]/#  ")
        try: 
            Q1=int(Q1); pr["100"]=CL.data[Q1]; pr["time100"]=tim; TakenSDs.append([CL.dttl[Q1],CL.dunit[Q1]]);
        except: 
            for i,t in enumerate(CL.extrattl):
                print("%d- %s" %(i,t)); 
            if len(CL.extradata)>0: Q1=raw_input("Is the >100 data in one of these lists? [n]/#  ")
            try:
                Q1=int(Q1); 
                for j,tt in enumerate(CL.extrattl[Q1]):
                    print("%d- %s" %(j,tt)); 
                Q1j=raw_input("Choose the position of the >100 micron data. [n]/#  ")
                try: 
                    Q1j=int(Q1j); 
                    pr["100"]=CL.extradata[Q1][Q1j]; TakenSDs.append([CL.extrattl[Q1][Q1j],CL.extraunit[Q1][Q1j]]);
                    k=[k for k,x in enumerate(CL.extrattl[Q1]) if x.lower() == 'time'][0]; pr["time100"]=CL.extradata[Q1][k]; 
                except: pr["100"]=np.ones(len(tim))*NaN; pr["time100"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
            except: pr["100"]=np.ones(len(tim))*NaN; pr["time100"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
        # >200 microns
        for i in range(len(CL.dttl)): print("%d- %s (%s)" % (i, CL.dttl[i], CL.dunit[i]))
        Q2=raw_input("No precipitation sdtype found. Is there a column in the basic data corresponding to the precipitation number concentration >200 microns? [n]/#  ")
        try: 
            Q2=int(Q2); pr["200"]=CL.data[Q2]; pr["time200"]=tim; TakenSDs.append([CL.dttl[Q2],CL.dunit[Q2]]);
        except: 
            for i,t in enumerate(CL.extrattl):
                print("%d- %s" %(i,t)); 
            if len(CL.extradata)>0: Q1=raw_input("Is the >200 data in one of these lists? [n]/#  ")
            try:
                Q1=int(Q1); 
                for j,tt in enumerate(CL.extrattl[Q1]):
                    print("%d- %s" %(j,tt)); 
                Q1j=raw_input("Choose the position of the >200 micron data. [n]/#  ")
                try: 
                    Q1j=int(Q1j); 
                    pr["200"]=CL.extradata[Q1][Q1j]; TakenSDs.append([CL.extrattl[Q1][Q1j],CL.extraunit[Q1][Q1j]]);
                    k=[k for k,x in enumerate(CL.extrattl[Q1]) if x.lower() == 'time'][0]; pr["time200"]=CL.extradata[Q1][k]; 
                except: pr["200"]=np.ones(len(tim))*NaN; pr["time200"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
            except: pr["200"]=np.ones(len(tim))*NaN; pr["time200"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
    elif len(p)==1: 
        pr["100"]=samac.dNdlogDp2N(p[0],100,nan); pr["200"]=samac.dNdlogDp2N(p[0],200,nan); pr["time100"]=p[0]["time"]; pr["time200"]=p[0]["time"];  
        TakenSDs.append([p[0]["Distname"],p[0]["units"]]); TakenSDs.append([p[0]["Distname"],p[0]["units"]]);
    elif len(p)>1: 
        p1=[x for x in CL.sd if 'P1' in x["sdtype"].upper()]    # looking for primary distribution
        if len(p1)==1: 
            pr["100"]=samac.dNdlogDp2N(p1[0],100,nan); pr["200"]=samac.dNdlogDp2N(p1[0],200,nan); pr["time100"]=p1[0]["time"]; pr["time200"]=p1[0]["time"]; 
            TakenSDs.append([p1[0]["Distname"],p1[0]["units"]]); TakenSDs.append([p1[0]["Distname"],p1[0]["units"]]);
        else: 
            print("Precipitation distributions:")
            for i,x in enumerate(p): print("%d- %s (sdtype=%s)" % (i,x["Distname"],x["sdtype"]))
            Q=raw_input("Which is your primary Precipitation distribution (consider editing sdtype if this is permanent)? [n]/#  ")
            try: 
                Q=int(Q); pr["100"]=samac.dNdlogDp2N(p[Q],100,nan); pr["200"]=samac.dNdlogDp2N(p[Q],200,nan); pr["time100"]=p[Q]["time"]; pr["time200"]=p[Q]["time"]; 
                TakenSDs.append([p[Q]["Distname"],p[Q]["units"]]); TakenSDs.append([p[Q]["Distname"],p[Q]["units"]]);
            except: 
                pr["100"]=np.ones(len(tim))*NaN; pr["200"]=np.ones(len(tim))*NaN; pr["time100"]=np.ones(len(tim))*NaN; pr["time200"]=np.ones(len(tim))*NaN;  
                TakenSDs.append(['','']); TakenSDs.append(['','']);
                print("[overview] Precipitation distribution not found.")
    else: print("[overview] impossible precipitation")
    if type(pr["time100"])==float: pr["100"]=np.ones(len(tim))*NaN; pr["200"]=np.ones(len(tim))*NaN; pr["time100"]=tim; pr["time200"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
    # adding lwc units to takenSDs: 
    TakenSDs.append(['LWC',lwcunits])
    
    Adj=[1,10,100,1000,1,1]     # Default magnifying values for [0-Altitude, 1-Ptcl conc, 2-droplet conc, 3-LWC, 4-5-2dc100,200]
    Offsts=[0,0,0,0,0,0]  # default offset values for [0-Altitude, 1-Ptcl conc, 2-droplet conc, 3-LWC, 4-5-2dc100,200]
    figi=1001            # Figure number
    Zx1=np.min(tim); Zx2=np.max(tim)      # initial time span of the figure
    Zs=np.array([Zx1,Zx2,NaN,NaN,NaN,NaN,NaN,NaN])
    if interact==1:
        plotcloud2(tim,alt,lat,lon,lwc,lwctime,ptc,ptctime,cdnc,cdnctime,pr,TakenSDs,Zs,Adj,Offsts,figi,1,proftimes=CL.times,Rtime=Rtime)
    else: plotcloud2(tim,alt,lat,lon,lwc,lwctime,ptc,ptctime,cdnc,cdnctime,pr,TakenSDs,Zs,Adj,Offsts,figi,0,proftimes=CL.times,Rtime=Rtime)        

