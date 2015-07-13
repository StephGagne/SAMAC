################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####

from pylab import *
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats.stats as st
from scipy import interpolate
import datetime as dt

import samac

#########################################################################        
##############################   vprof   ################################
#########################################################################
def vprof(CL,interact=None,axtype1='log',axtype2='log'):
    """ This method plots the Temperature, Theta-Q, effective diameter, LWC, RH, Particle and Droplet and Precipitation drop concentrations, in the cloud as a function of Altitude. One plot for each vertical scan is produced.
        Use CloudObj.vprof(0), CloudObj.vprof() or CloudObj.vprof(interact=0)for a plot NOT in interactive mode;   
        Use CloudObj.vprof(1) or CloudObj.vprof(interact=1)for a plot in interactive mode.
        axtype1: (default='log') 'log' if the bottom axis of the left hand plot should be logarithmic, 'lin' if it should be linear.
        axtype2: (default='log') 'log' if the bottom axis of the right hand plot should be logarithmic, 'lin' if it should be linear. 
        Fully handles the extradata module."""

    
    if interact==None or interact==0: interact=0
    elif interact==1: plt.ion()
    
    # time
    post=[i for i,x in enumerate(CL.dttl) if x == 'time']
    try: 
        talt=CL.data[post]
        talt=talt.reshape(max(np.shape(talt)),)
    except:
        print("[vprof] No time (or multiple) was found. Forget about this plot! It won't happen unless you have a time")
    # altitude
    pos=[i for i,x in enumerate(CL.dttl) if x == 'altitude']
    try: 
        Alt=CL.data[pos]
    except:
        print("[vprof] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude in the main data")
    # temperature
    pos=[i for i,x in enumerate(CL.dttl) if x == 'temperature']
    if len(pos)>1:  print("[vprof] Multiple temperatures found."); T=np.ones((2,))*NaN    # temperature found in the basic data
    elif len(pos)==1: T=CL.data[pos[0]]; Ttime=talt; 
    else:       # looking for temperature in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
        if len(posx)==1: 
            T=CL.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            Ttime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[vprof] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN;
    # theta-Q
    pos=[i for i,x in enumerate(CL.dttl) if x == 'theta-q']
    if len(pos)>1: print("[vprof] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN     # Theta-Q found in the basic data
    elif len(pos)==1: 
        ThQ=CL.data[pos[0]]; ThQtime=talt; 
        if 'k' in CL.dunit[pos[0]].lower(): ThQ=ThQ-273.15        # converting to celsius
    else:       # looking for Theta-Q in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'theta-q']    # check all titles matching with Theta-Q
        if len(posx)==1: 
            ThQ=CL.extradata[posx[0][0]][posx[0][1]]    # loading the Theta-Q data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            ThQtime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            if 'k' in CL.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else: print("[vprof] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN;
    # LWC
    pos=[i for i,x in enumerate(CL.dttl) if x == 'LWC']
    if len(pos)>1: print("[vprof] Multiple LWC found."); lwc=np.ones((2,))*NaN     # LWC found in the basic data
    elif len(pos)==1: lwc=CL.data[pos[0]]; lwctime=talt; 
    else:       # looking for LWC in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
        if len(posx)==1: 
            lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[vprof] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN;
    # RH
    pos=[i for i,x in enumerate(CL.dttl) if x.upper() == 'RH']
    if len(pos)>1: print("[vprof] Multiple RH found."); rh=np.ones((2,))*NaN     # RH found in the basic data
    elif len(pos)==1: rh=CL.data[pos[0]]; rhtime=talt; 
    else:       # looking for RH in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'rh']    # check all titles matching with RH
        if len(posx)==1: 
            rh=CL.extradata[posx[0][0]][posx[0][1]]    # loading the RH data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            rhtime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[vprof] No RH (or multiple) was found in the basic data or the extra data."); rh=np.ones((2,))*NaN;

       
    ### Interpolating the altitude onto the data ###
    f=interpolate.interp1d(talt,np.ma.filled(Alt,nan),kind='linear')   # f=interp(x=time of altitude, y=alt, linear interpolation)
    
    A=list(); C=list(); D=list(); 
    for sd in CL.sd:
        if 'A' in sd["sdtype"].upper():
            A.append(sd)
        elif 'C' in sd["sdtype"].upper():
            C.append(sd)
        elif 'P' in sd["sdtype"].upper():
            D.append(sd)
        else: 
            qa=0
            while qa==0:
                Q=raw_input("To which category does this distribution belong? %s \nA (aerosols), C (cloud droplets), or P (precipitation). U (undetermined) if this distribution does not belong to any of these categories.:  " % sd["Distname"])
                if Q.upper()=='A': A.append(sd); qa=1
                elif Q.upper()=='C': C.append(sd); qa=1
                elif Q.upper()=='P': D.append(sd); qa=1
                elif Q.upper()=='U': qa=1
                else: print("You must chose A,C,P or U.")
                

    for v,ts in enumerate(CL.times['verticloud']):
        [ts1,ts2]=ts        # start and end time of the profile
        # adapting T, ThQ, lwc, rh
        try: 
            T1=T[(Ttime<=ts2)*(Ttime>=ts1)].reshape(-1,); Ttime1=Ttime[(Ttime<=ts2)*(Ttime>=ts1)].reshape(-1,)
            if len(Ttime1)<2: T1=np.ones((2,))*NaN; Ttime1=np.array([ts1,ts2]);
        except: T1=T; Ttime1=np.array([ts1,ts2])
        try: 
            ThQ1=ThQ[(ThQtime<=ts2)*(ThQtime>=ts1)].reshape(-1,); ThQtime1=ThQtime[(ThQtime<=ts2)*(ThQtime>=ts1)].reshape(-1,)
            if len(ThQtime1)<2: ThQ1=np.ones((2,))*NaN; ThQtime1=np.array([ts1,ts2]);
        except: ThQ1=ThQ; ThQtime1=np.array([ts1,ts2])
        try: 
            lwc1=lwc[(lwctime<=ts2)*(lwctime>=ts1)].reshape(-1,); lwctime1=lwctime[(lwctime<=ts2)*(lwctime>=ts1)].reshape(-1,)
            if len(lwctime1)<2: lwc1=np.ones((2,))*NaN; lwctime1=np.array([ts1,ts2]);
        except: lwc1=lwc; lwctime1=np.array([ts1,ts2]); 
        try: 
            rh1=rh[(rhtime<=ts2)*(rhtime>=ts1)].reshape(-1,); rhtime1=rhtime[(rhtime<=ts2)*(rhtime>=ts1)].reshape(-1,)
            if len(rhtime1)<2: rh1=np.ones((2,))*NaN; rhtime1=np.array([ts1,ts2]);
        except: rh1=rh; rhtime1=np.array([ts1,ts2]); 
        
        AA=list();  CC=list(); DD100=list(); DD200=list(); DDtotSize=list(); DDtot=list();
        for sd in A:        # Aerosols
            try: 
                newalt=f(sd["time"][(sd["time"]<=ts2)*(sd["time"]>=ts1)])      # new alt=f(new time)
                AA.append(np.ma.vstack((newalt,sd["total"][(sd["time"]<=ts2)*(sd["time"]>=ts1)],samac.effrad(CL,inst=sd["Distname"],bindist='lin')[(sd["time"]<=ts2)*(sd["time"]>=ts1)])))
            except: AA.append(np.ones((3,1))*NaN); print("[vprof] %s was not successfully added to the plot." % sd["Distname"])
        for sd in C:        # Cloud droplets
            try: 
                newalt=f(sd["time"][(sd["time"]<=ts2)*(sd["time"]>=ts1)])      # new alt=f(new time)
                CC.append(np.ma.vstack((newalt,sd["total"][(sd["time"]<=ts2)*(sd["time"]>=ts1)],samac.effrad(CL,inst=sd["Distname"],bindist='lin')[(sd["time"]<=ts2)*(sd["time"]>=ts1)])))
            except: CC.append(np.ones((3,1))*NaN); print("[vprof] %s was not successfully added to the plot." % sd["Distname"])
        for sd in D:        # Drizzle drops
            try: 
                ix=(sd["time"]<=ts2)*(sd["time"]>=ts1)
                if sum(ix)>0:
                    newalt=f(sd["time"][ix])
                    if sum(~isnan(sd["data"]))==0 and sum(~isnan(sd["total"]))>0:
                        DDtot.append(np.ma.vstack((newalt,sd["total"][ix])))
                        DDtotSize.append(raw_input("What is the size above which the total was made?  "))
                    elif sum(~isnan(sd["data"]))==0 and sum(~isnan(sd["total"]))==0: print("%s data full of NaNs" % sd["Distname"]); crash  # --> going into except
                    else: 
                        DD100.append(np.ma.vstack((newalt,samac.dNdlogDp2N(sd,100,nan)[ix])))
                        DD200.append(np.ma.vstack((newalt,samac.dNdlogDp2N(sd,200,nan)[ix])))
                else: print("[vprof] %s data not defined for this vertical profile (%d)" % (sd["Distname"],50+v)); crash     # --> going into except
            except: 
                DD100.append(np.ones((2,1))*NaN); DD200.append(np.ones((2,1))*NaN); DDtot.append(np.ones((2,1))*NaN); DDtotSize.append('None')
                print("[vprof] %s was not successfully added to the plot (%d)." % (sd["Distname"],50+v))

        LnStlz=['-','--',':','-.','None']
        Mrkrz=['*','o','D','s','v','p','+','x','h']
        phandlez=list();
        
        ### plotting ###
        fig=plt.figure(50+v, figsize=(9.7,6.5))
        fig.subplots_adjust(right=0.70,left=0.1,top=0.83,wspace=0.3)        # Way to make space for the legend.
        ax1=fig.add_subplot(121)
        phandle,=ax1.plot(lwc1,f(lwctime1)[0],'g.-',label='LWC'); phandlez.append(phandle)
        for ii,p in enumerate(CC):
            phandle,=ax1.plot(p[1],p[0],linestyle=LnStlz[ii], marker='None', color='k',label=C[ii]["Distname"]); phandlez.append(phandle)
        for ii,p in enumerate(AA):
            phandle,=ax1.plot(p[1],p[0],linestyle=LnStlz[ii], marker='None', color='r',label=A[ii]["Distname"]); phandlez.append(phandle)
        if axtype1=='lin': ax1.set_xscale('linear')
        else:
            try: ax1.set_xscale('log')
            except: plot(1,1,'w.'); ax1.set_xscale('log')
        ax1.set_ylabel('Altitude (m)')
        try:
            if len(A+C)==len([x["units"] for x in A+C if x["units"]==A[0]["units"]]):     # list of units for aerosol and cloud droplet instruments combined
                ax1.set_xlabel('Aero & CDNC Conc. (%s), LWC (g/m$^{3}$)' %(A[0]["units"]))
            else:
                ax1.set_xlabel('Aero & CDNC Conc. (units vary), LWC (g/m$^{3}$)')
        except: pass
        ax2=ax1.twiny()
        if sum(~isnan(rh1))>=2: phandle,=ax2.plot(rh1,f(rhtime1)[0],'b--',label='RH'); phandlez.append(phandle)
        else: plot(50,st.nanmean(Alt[0]),'w.');
        ax2.set_xlabel("RH (%)")
        [Zx1, Zx2, Zy1, Zy2]=ax2.axis()

        ax3=fig.add_subplot(122)
        for ii,p in enumerate(DD100):
            phandle,=ax3.plot(p[1],p[0],color='c',linestyle='None',marker=Mrkrz[ii],label='Drizzle ' + D[ii]["Distname"]+ ' >100um');             
            phandlez.append(phandle)
        for ii,p in enumerate(DD200):
            phandle,=ax3.plot(p[1],p[0],color='m',linestyle='None',marker=Mrkrz[ii],label='Drizzle ' + D[ii]["Distname"]+ ' >200um'); 
            phandlez.append(phandle)
        for ii,p in enumerate(DDtot):
            phandle,=ax3.plot(p[1],p[0],color='c',linestyle='None',marker=Mrkrz[ii],label='Drizzle ' + D[ii]["Distname"] + ' (tot>' + DDtotSize[ii] + ')'); 
            phandlez.append(phandle)
        if axtype2=='lin': ax3.set_xscale('linear')
        else:
            try: ax3.set_xscale('log')
            except: plot(1,1,'w.'); ax3.set_xscale('log')
        try:
            if len(D)==len([x["units"] for x in D if x["units"]==D[0]["units"]]):     # list of units for precipitation instruments
                ax3.set_xlabel('Drizzle Conc. (%s)' %(D[0]["units"]))
            else:
                ax3.set_xlabel('Drizzle Conc. (units vary)')
        except: pass
        ax4=ax3.twiny()
        phandle,=ax4.plot(T1,f(Ttime1)[0],'b', label="Temperature"); phandlez.append(phandle)
        phandle,=ax4.plot(ThQ1,f(ThQtime1)[0],'m-', label="ThetaQ"); phandlez.append(phandle)
        for ii,p in enumerate(CC):
            phandle,=ax4.plot(p[2],p[0],linestyle=LnStlz[ii], marker='None', color='k',label='Deff (' + C[ii]["Distname"] + ')'); phandlez.append(phandle)
        for ii,p in enumerate(AA):
            phandle,=ax4.plot(p[2],p[0],linestyle=LnStlz[ii], marker='None', color='r',label='Deff (' + A[ii]["Distname"] + ')'); phandlez.append(phandle)
        ax4.set_xlabel("Temperature (C), ThetaQ (C), De (um)"); 
        ax3.set_ylim(Zy1,Zy2)
        

        fig.text(0.35, 0.92, "%s @ %s - %s \n(vert. scan #%d)" % (CL.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=CL.data[post][0][0])).strftime('%H:%M'), CL.desc["humanplace"], v),horizontalalignment='center')
        lines = phandlez
        ax3.legend(lines, [l.get_label() for l in lines],loc=(1.01,0))
        
    if interact==0:
        if np.shape(CL.times['verticloud'])[0]>=1: 
            print("[vprof] The sequence may be blocked until you close all figures.")
            plt.show()
        else: print("[vprof] No vertical scan was measured in this cloud.")


