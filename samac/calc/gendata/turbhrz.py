################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####


from pylab import *
import numpy as np
import scipy.stats.stats as st
from scipy import interpolate

import samac

#########################################################################        
###############################  turbhrz  ###############################
#########################################################################
def turbhrz(CL,MaxCA=2,MinDist=20000,PtNo=25,Pts=5,Pts2=50):
    """ This method calculates the turbulence average over at least 20 km or more in a horizontal enough leg (the plane cannot be above an angle of 2 degrees using a running standard deviation (of 25 points) to calculate the turbulence. It is assumed that the TAS is in the basic data, along with other aircraft related data (e.g. extradata module not handled for TAS). Extradata module handled for updraft velocity.
    Options: 
        MaxCA: maximum change in the combined angle (phi+theta) in degrees. Default is 2.
        MinDist: minimum distance (in meters) an aircraft must be level enough for to calculated the turbulence (updraft). Default is 20000 or 20 km.
        PtNo: number of points from which to calculate the turbulence=stdev(updraft vel). Default is 25
        Pts: number of points taken in running statistics of angles. Default is 5.
        Pts2: number of points taken in running statistics of the combined angle. Default is 50.
    Returns a dictionary:
    R["InCloud"]: List with one item for each stretch level enough for at least MinDist km in cloud (N, first=0). Each item contains an array with the turbulence value R["InCloud"][N][0] and the distance on which it was calculated R["InCloud"][N][1].
    R["BlwCloud"]: List with one item for each stretch level enough for at least 20 km below cloud (N, first=0). Each item contains an array with the turbulence value R["InCloud"][N][0] and the distance on which it was calculated R["InCloud"][N][1].
    R["InCloudStretches"]: List with one item for each stretch level enough during at least 20 km in cloud (N, first=0). Each item contains an array with the indexes corresponding to the stretch used to calculate the corresponding turbulence value.
    R["BlwCloudStretches"]: List with one item for each stretch level enough during at least 20 km below cloud (N, first=0). Each item contains an array with the indexes corresponding to the stretch used to calculate the corresponding turbulence value. """
    
    palt=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
    ptim=[i for i,x in enumerate(CL.dttl) if x == 'time'][0]
    ptas=[i for i,x in enumerate(CL.dttl) if 'TAS' in x.upper()][0]
    pudv=[i for i,x in enumerate(CL.dttl) if x == 'udvel']
    if len(pudv)==1: pudv=pudv[0]; udv=CL.data[pudv][0:-1]
    elif len(pudv)>1: print("[turbhrz] multiple updraft velocities found in the basic data.")
    else:
        posx=[] 
        for i,ttl in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'udvel']    # check all titles matching with updraft velocity
        if len(posx)==1: 
            udv=CL.extradata[posx[0][0]][posx[0][1]]    # loading the updraft velocity data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            udvt=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            # padding so the interpolation works
            if udvt[0]>CL.data[ptim][0]: udv=hstack((nan,udv)); udvt=hstack((CL.data[ptim][0],udvt));
            if udvt[-1]<CL.data[ptim][-1]: udv=hstack((udv,nan)); udvt=hstack((udvt,CL.data[ptim][-1]));
            fudv=interpolate.interp1d(udvt,udv,kind='linear')
            udv=fudv(CL.data[ptim])[0:-1]
        else: print("[turbhrz] No updraft velocity (or multiple) found in the basic or the extra data."); return dict()
    Alt=CL.data[palt][0:-1]
    A=samac.angles(CL)
    theta=samac.runstats(A[1],Pts)[0]  # change of altitude
    phi=diff(A[2])
    phi=np.ma.masked_where((phi>350)+(isnan(phi)), phi)
    phi=np.ma.masked_where((phi<-350+(isnan(phi))), phi)
    phi=hstack([NaN,samac.runstats(phi,Pts)[0]])     # change in direction
    phi[(abs(phi)>=89)]=0   # accounts for jump from 0 to 360 degrees. After smoothing, the limit cannot be too high. 
    CAstats=samac.runstats(abs(theta)+abs(phi),Pts2)      # combined angle
    CA=CAstats[0]
#        udv=CL.data[pudv][0:-1]        # updraft
    Sp=CL.data[ptas][0:-1]
    dt=diff(CL.data[ptim]*24*60*60)
    D=cumsum(Sp*dt)
    # Initializing indexes for below cloud and in-cloud measurements.
    bc=ones(np.shape(CL.data[ptim]))*False
    ic=ones(np.shape(CL.data[ptim]))*False
    for j in range(len(CL.times["belowcloud"])):
        bctmp=(CL.data[ptim]>=CL.times["belowcloud"][j][0])*(CL.data[ptim]<=CL.times["belowcloud"][j][1])
        try: bc=bc+bctmp;
        except: bc=bctmp;
    for j in range(len(CL.times["horicloud"])):
        ictmp=(CL.data[ptim]>=CL.times["horicloud"][j][0])*(CL.data[ptim]<=CL.times["horicloud"][j][1])
        try: ic=ic+ictmp;
        except: ic=ictmp;
    # Below Cloud #
    Stret=list()
    i=0
    while i<len(CA):
        ix=list()
        if (CA[i]<=MaxCA and bc[i]==True):
            while (i<len(CA) and CA[i]<=MaxCA and bc[i]==True):
                ix.append(i)
                i=i+1
            Stret.append(ix)
        else: i=i+1
    BCStretches=list()
    for n,leg in enumerate(Stret):
        Dist=D[leg]; 
        Dist=Dist[(Dist.mask==False)]
        if (Dist[-1]-Dist[0]) >= MinDist:     # if the stretch is longer than the minimum distance.
            BCStretches.append(np.array(leg))
        del Dist
    # In-Cloud #
    Stret=list()
    i=0
    while i<len(CA):
        ix=list()
        if (CA[i]<=MaxCA and ic[i]==True):
            while (i<len(CA) and CA[i]<=MaxCA and ic[i]==True):
                ix.append(i)
                i=i+1
            Stret.append(ix)
        else: i=i+1
    ICStretches=list()
    for n,leg in enumerate(Stret):
        Dist=D[leg]; Dist=Dist[(Dist.mask==False)]
        if (Dist[-1]-Dist[0]) >= MinDist:     # if the strecht is longer than the minimum distance.
            ICStretches.append(np.array(leg))       
        del Dist
    # Calculating the average turbulence.
    ICturb=list(); BCturb=list()
    for leg in ICStretches:
        Dist=D[leg]; Dist=Dist[(Dist.mask==False)]
        if (sum(leg)>PtNo and sum(~isnan(udv[leg]))>PtNo):
            turb=st.nanmean(samac.runstats(udv[leg],PtNo)[1])
            #print("In-Cloud %d: %.4f km. w' = %0.4f m/s" %(cn,(Dist[-1]-Dist[0])/1000,turb))
            ICturb.append(np.array([turb,(Dist[-1]-Dist[0])]))
            del Dist
    for leg in BCStretches:
        Dist=D[leg]; Dist=Dist[(Dist.mask==False)]
        if (sum(leg)>PtNo and sum(~isnan(udv[leg]))>PtNo):
            turb=st.nanmean(samac.runstats(udv[leg],PtNo)[1])
            #print("Below-Cloud %d: %.4f km. w' = %0.4f m/s" %(cn,(Dist[-1]-Dist[0])/1000,turb))
            BCturb.append(np.array([turb,(Dist[-1]-Dist[0])]))
        else: BCturb.append(np.array([nan,(Dist[-1]-Dist[0])]))
        del Dist
            
    R=dict()
    R["InCloud"]=ICturb; R["BlwCloud"]=BCturb; R["InCloudStretches"]=ICStretches; R["BlwCloudStretches"]=BCStretches
    return R        

