from pylab import *
import numpy as np
import matplotlib.pyplot as plt

import scipy.stats.stats as st
from scipy import interpolate
import datetime as dt



#########################################################################        
#############################   wholeprof   #############################
#########################################################################
def wholeprof(CL,interact=None,axtype1='log',axtype2='log'):
    """ This method plots the Temperature, LWC, Particle and Droplet and Precipitation drop concentrations (and more), in the cloud as a function of Altitude. All the available data, regardless of profiles are plotted.
        Use CloudObj.wholeprof(0), CloudObj.wholeprof() or CloudObj.wholeprof(interact=0)for a plot NOT in interactive mode;   
        Use CloudObj.wholeprof(1) or CloudObj.wholeprof(interact=1)for a plot in interactive mode.
        axtype1: (default='log') 'log' if the bottom axis of the left hand plot should be logarithmic, 'lin' if it should be linear.
        axtype2: (default='log') 'log' if the bottom axis of the right hand plot should be logarithmic, 'lin' if it should be linear. """
    
    if  interact==0: pass
    elif interact==None or interact==1: plt.ion()
    ts1=CL.times["cloud"][0][0]; ts2=CL.times["cloud"][0][1]
    
    # time
    post=[i for i,x in enumerate(CL.dttl) if x == 'time']
    try: 
        talt=CL.data[post]
        talt=talt.reshape(max(np.shape(talt)),)
    except:
        print("[wholeprof] No time (or multiple) was found. Forget about this plot! It won't happen unless you have (a) time")
    # altitude
    pos=[i for i,x in enumerate(CL.dttl) if x == 'altitude']
    try: 
        Alt=CL.data[pos]
    except:
        print("[wholeprof] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude in the main data")
    # temperature
    pos=[i for i,x in enumerate(CL.dttl) if x == 'temperature']
    if len(pos)>1:  print("[wholeprof] Multiple temperatures found."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN    # temperature found in the basic data
    elif len(pos)==1: T=CL.data[pos[0]]; Ttime=talt; 
    else:       # looking for temperature in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
        if len(posx)==1: 
            T=CL.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            Ttime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            T=T[(Ttime>=ts1)*(Ttime<=ts2)].reshape(-1,); Ttime=Ttime[(Ttime>=ts1)*(Ttime<=ts2)].reshape(-1,); 
        else: print("[wholeprof] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN
    # theta-Q
    pos=[i for i,x in enumerate(CL.dttl) if x == 'theta-q']
    if len(pos)>1: print("[wholeprof] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN     # Theta-Q found in the basic data
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
            ThQ=ThQ[(ThQtime>=ts1)*(ThQtime<=ts2)].reshape(-1,); ThQtime=ThQtime[(ThQtime>=ts1)*(ThQtime<=ts2)].reshape(-1,); 
            if 'k' in CL.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else: print("[wholeprof] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN
    # LWC
    pos=[i for i,x in enumerate(CL.dttl) if x == 'LWC']
    if len(pos)>1: print("[wholeprof] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN     # LWC found in the basic data
    elif len(pos)==1: lwc=CL.data[pos[0]]; lwctime=talt; 
    else:       # looking for LWC in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
        if len(posx)==1: 
            lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            lwc=lwc[(lwctime>=ts1)*(lwctime<=ts2)].reshape(-1,); lwctime=lwctime[(lwctime>=ts1)*(lwctime<=ts2)].reshape(-1,); 
        else: print("[wholeprof] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
    # RH
    pos=[i for i,x in enumerate(CL.dttl) if x.upper() == 'RH']
    if len(pos)>1: print("[wholeprof] Multiple RH found."); rh=np.ones((2,))*NaN; rhtime=np.ones((2,))*NaN     # RH found in the basic data
    elif len(pos)==1: rh=CL.data[pos[0]]; rhtime=talt; 
    else:       # looking for RH in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'rh']    # check all titles matching with RH
        if len(posx)==1: 
            rh=CL.extradata[posx[0][0]][posx[0][1]]    # loading the RH data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            rhtime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            rh=rh[(rhtime>=ts1)*(rhtime<=ts2)].reshape(-1,); rhtime=rhtime[(rhtime>=ts1)*(rhtime<=ts2)].reshape(-1,); 
        else: print("[wholeprof] No RH (or multiple) was found in the basic data or the extra data."); rh=np.ones((2,))*NaN; rhtime=np.ones((2,))*NaN
        
    # adapting T, ThQ, lwc, rh (if empty or too short)
    if len(Ttime)<2: T=np.ones((2,))*NaN; Ttime=np.array([ts1,ts2]);
    if len(ThQtime)<2: ThQ=np.ones((2,))*NaN; ThQtime=np.array([ts1,ts2]);
    if len(lwctime)<2: lwc=np.ones((2,))*NaN; lwctime=np.array([ts1,ts2]);
    if len(rhtime)<2: rh=np.ones((2,))*NaN; rhtime=np.array([ts1,ts2]);
            
        
    ### Interpolating the altitude onto the data ###
    f=interpolate.interp1d(talt,np.ma.filled(Alt,nan),kind='linear')   # f=interp(x=time of altitude, y=alt, linear interpolation)
    
    # classifying the size distributions according to the type of particles they measure
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

    AA=list();  CC=list(); DD100=list(); DD200=list(); DDtotSize=list(); DDtot=list();
    for sd in A:        # Aerosols
        try: 
            ix=nonzero((sd["time"]>=talt[0])*(sd["time"]<=talt[-1]))[0]
            sdtime=sd["time"][ix]
            newalt=f(sdtime)
            AA.append(np.ma.vstack((newalt,sd["total"][ix])))
        except: AA.append(np.ones((2,1))*NaN); print("[wholeprof] %s was not successfully added to the plot." % sd["Distname"])
    for sd in C:        # Cloud droplets
        try: 
            ix=nonzero((sd["time"]>=talt[0])*(sd["time"]<=talt[-1]))[0]
            sdtime=sd["time"][ix]
            newalt=f(sdtime)
            CC.append(np.ma.vstack((newalt,sd["total"][ix])))
        except: CC.append(np.ones((2,1))*NaN); print("[wholeprof] %s was not successfully added to the plot." % sd["Distname"])
    for sd in D:        # Drizzle drops
        try: 
            ix=nonzero((sd["time"]>=talt[0])*(sd["time"]<=talt[-1]))[0]
            sdtime=sd["time"][ix]
            newalt=f(sdtime)
            if sum(~isnan(sd["data"]))==0 and sum(~isnan(sd["total"]))>0:
                DDtot.append(np.ma.vstack((newalt,sd["total"][ix])))
                DDtotSize.append(raw_input("What is the size above which the total was made?  "))
            elif sum(~isnan(sd["data"]))==0 and sum(~isnan(sd["total"]))==0: crash
            else: 
                DD100.append(np.ma.vstack((newalt,dNdlogDp2N(sd,100,nan)[ix])))
                DD200.append(np.ma.vstack((newalt,dNdlogDp2N(sd,200,nan)[ix])))
        except: 
            DD100.append(np.ones((2,1))*NaN); DD200.append(np.ones((2,1))*NaN); DDtot.append(np.ones((2,1))*NaN); DDtotSize.append('None')
            print("[wholeprof] %s was not successfully added to the plot." % sd["Distname"])


    # defining Line type, Marker types and predefining handles' list.
    LnStlz=['-','--',':','-.','None']
    Mrkrz=['*','o','D','s','v','p','+','x','h']
    phandlez=list();


    ### plotting ###
    fig=plt.figure(60, figsize=(9.7,6.5))
    fig.subplots_adjust(right=0.70,left=0.1,top=0.83,wspace=0.3)        # Way to make space for the legend.
    ax1=fig.add_subplot(121)
    phandle,=ax1.plot(lwc,f(lwctime)[0],'g-',label='LWC'); phandlez.append(phandle)
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
    if sum(~isnan(rh))>=2: phandle,=ax2.plot(rh,f(rhtime)[0],'b--',label='RH'); phandlez.append(phandle)
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
    phandle,=ax4.plot(T,f(Ttime)[0],'b', label="Temperature"); phandlez.append(phandle)
    phandle,=ax4.plot(ThQ,f(ThQtime)[0],'m-', label="ThetaQ"); phandlez.append(phandle)
    ax4.set_xlabel("Temperature (C), ThetaQ (C)");  
    ax3.set_ylim(Zy1,Zy2)
        

    fig.text(0.35, 0.92, "%s @ %s - %s \n(Whole cloud profile)" % (CL.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=CL.data[post][0][0])).strftime('%H:%M'), CL.desc["humanplace"]),horizontalalignment='center')
    lines = phandlez
    ax3.legend(lines, [l.get_label() for l in lines],loc=(1.01,0))
    plt.show()
            
    if interact==0:
        if np.shape(CL.times['cloud'])[0]>=1: 
            print("[wholeprof] The sequence may be blocked until you close all figures.")
            plt.show()
        else: print("[wholeprof] No graph has been generated.")


