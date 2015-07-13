################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####


from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import scipy.stats.stats as st

import samac

#########################################################################        
##############################   hrzplot   ##############################
#########################################################################
def hrzplot(CL,interact=None,Rtime=None):
    """ This method plots the LWC & Altitude, Na (aerosols) & Nd (droplets), Temperature & Theta-Q, updraft velocity and its standard deviation (turbulence) as a function of time. One plot for each horizontal scan is produced. Fully handles the extradata module.
        Use CloudObj.hrzplot(0), CloudObj.hrzplot() or CloudObj.hrzplot(interact=0) for a plot NOT in interactive mode;   
        Use CloudObj.hrzplot(1) or CloudObj.hrzplot(interact=1)for a plot in interactive mode. 
        Use CloudObd.hrzplot(1,1) or CloudObd.hrzplot(Rtime=1) to get the time in Hour:Minute format."""
    
    if interact==None or interact==0: interact=0
    elif interact==1: plt.ion()
    
    # time
    post=[i for i,x in enumerate(CL.dttl) if x == 'time']
    if len(post)==1: talt=CL.data[post[0]]
    else: print("[hrzplot] No time (or multiple) was found. Forget about this plot! It won't happen unless you have a time")
    # altitude
    pos=[i for i,x in enumerate(CL.dttl) if x == 'altitude']
    if len(pos)==1: Alt=CL.data[pos[0]]
    else: print("[hrzplot] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an a(l)titude")
    # temperature
    pos=[i for i,x in enumerate(CL.dttl) if x == 'temperature']
    if len(pos)>1:  print("[hrzplot] Multiple temperatures found."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN    # temperature found in the basic data
    elif len(pos)==1: T=CL.data[pos[0]]; Ttime=talt; 
    else:       # looking for temperature in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
        if len(posx)==1: 
            T=CL.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            Ttime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[hrzplot] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN;
    # theta-Q
    pos=[i for i,x in enumerate(CL.dttl) if x == 'theta-q']
    if len(pos)>1: print("[hrzplot] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN     # Theta-Q found in the basic data
    elif len(pos)==1: 
        ThQ=CL.data[pos[0]]; ThQtime=talt; 
        #if 'k' in CL.dunit[pos[0]].lower(): ThQ=ThQ-273.15        # converting to celsius
    else:       # looking for Theta-Q in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'theta-q']    # check all titles matching with Theta-Q
        if len(posx)==1: 
            ThQ=CL.extradata[posx[0][0]][posx[0][1]]    # loading the Theta-Q data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            ThQtime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            #if 'k' in CL.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else: print("[hrzplot] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN    
    # LWC
    pos=[i for i,x in enumerate(CL.dttl) if x == 'LWC']
    if len(pos)>1: print("[hrzplot] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN     # LWC found in the basic data
    elif len(pos)==1: lwc=CL.data[pos[0]]; lwctime=talt; 
    else:       # looking for LWC in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
        if len(posx)==1: 
            lwc=CL.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            lwctime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[hrzplot] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
    # updraft velocity
    pos=[i for i,x in enumerate(CL.dttl) if x == 'udvel']
    if len(pos)>1: print("[hrzplot] Multiple updraft velocity found."); udvel=np.ones((2,))*NaN; udveltime=np.ones((2,))*NaN     # udvel found in the basic data
    elif len(pos)==1: udvel=CL.data[pos[0]]; udveltime=talt; 
    else:       # looking for updraft velocity in the extradata
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'udvel']    # check all titles matching with LWC
        if len(posx)==1: 
            udvel=CL.extradata[posx[0][0]][posx[0][1]]    # loading the updraft velocity data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            udveltime=CL.extradata[posx[0][0]][j]     # loading associated time stamp
        else: print("[hrzplot] No LWC (or multiple) was found in the basic data or the extra data."); udvel=np.ones((2,))*NaN; udveltime=np.ones((2,))*NaN
    # turbulence calculations
    if sum(~isnan(udveltime))==0: turb=np.ones((2,))*NaN;
    else: turb=samac.runstats(udvel,25)     # returns the running average and running standard deviation of the first argument over the number of points specified in the second argument.
    try: turb[1][udvel.mask]=NaN
    except: pass
    # Cloud droplet primary
    pos=[i for i,x in enumerate(CL.sd) if 'C1' in x["sdtype"].upper()]
    if len(pos)==1: pass
    elif len(pos)==0: 
        pos=[i for i,x in enumerate(CL.sd) if 'C' in x["sdtype"].upper()]
        if len(pos)==1: pass
        elif len(pos)==0: 
            print("[hrzplot] No cloud droplet size distribution found.");
            CDrop=np.ones((2,1))*NaN; CDroptime=np.ones((2,))*NaN;  CDropName=''
        elif len(pos)>1:
            for i,p in enumerate(pos):
                print("%d - %s (%s)" %(i,CL.sd[p]["Distname"],CL.sd[p]["sdtype"]))
            try: 
                chosenp=pos[int(raw_input("Which of these distributions do you want to use for cloud droplets? (number, enter if none)  "))]
                pos=list(); pos.append(chosenp)
            except: 
                print("[hrzplot] Failed to load the distribution.")
                CDrop=np.ones((2,1))*NaN; CDroptime=np.ones((2,))*NaN;  CDropName=''
    if len(pos)==1:
        pos=pos[0]
        CDrop=CL.sd[pos]["total"]; CDroptime=CL.sd[pos]["time"]; CDropName=CL.sd[pos]["Distname"]
    else: print("[hrzplot] No distribution was found.")
    # Aerosol primary
    pos=[i for i,x in enumerate(CL.sd) if 'A1' in x["sdtype"].upper()]
    if len(pos)==1: pass
    elif len(pos)==0: 
        pos=[i for i,x in enumerate(CL.sd) if 'A' in x["sdtype"].upper()]
        if len(pos)==1: pass
        elif len(pos)==0: 
            print("[hrzplot] No aerosol size distribution found.");
            Aero=np.ones((2,1))*NaN; Aerotime=np.ones((2,))*NaN; AeroName=''
        elif len(pos)>1:
            for i,p in enumerate(pos):
                print("%d - %s (%s)" %(i,CL.sd[p]["Distname"],CL.sd[p]["sdtype"]))
            try: 
                chosenp=pos[int(raw_input("Which of these distributions do you want to use for aerosols? (number, enter if none)  "))]
                pos=list(); pos.append(chosenp)
            except: 
                print("[hrzplot] Failed to load the distribution.")
                Aero=np.ones((2,1))*NaN; Aerotime=np.ones((2,))*NaN; AeroName=''
    if len(pos)==1:
        pos=pos[0]
        Aero=CL.sd[pos]["total"]; Aerotime=CL.sd[pos]["time"]; AeroName=CL.sd[pos]["Distname"]
    else: print("[hrzplot] No distribution was found.")
    
    
    
    for v,ts in enumerate(CL.times['horicloud']):
        [ts1,ts2]=ts        # start and end time of the profile
        dayn=floor(ts1)

        ### plotting ###
        fig=plt.figure(70+v, figsize=(9.0,8.5))
        fig.subplots_adjust(right=0.63,left=0.125,top=0.90)        # Way to make space for the legend.
        
        ax11=fig.add_subplot(411)
        ix=nonzero((talt>=ts1)*(talt<=ts2))[0]
        if sum(~isnan(talt[ix]))==0:
            phandle11,=ax11.plot(ts-dayn,np.array([1,1]),'w.',label=("Altitude: m"))
        else:
            phandle11,=ax11.plot(talt[ix]-dayn,Alt[ix],'c-',label=("Altitude: %.0f m" % st.nanmean(Alt[ix])))
        ax12=ax11.twinx()
        ix=nonzero((lwctime>=ts1)*(lwctime<=ts2))[0]
        if sum(~isnan(lwctime[ix]))==0:
            phandle12,=ax12.plot(ts-dayn,np.array([1,1]),'w.',label=('LWC: g/m3'))
        else:
            phandle12,=ax12.plot(lwctime[ix]-dayn,lwc[ix],'m-',label=('LWC: %.2f g/m3' % st.nanmean(lwc[ix])))
        #ax11.set_xlabel('Time')
        ax11.set_ylabel('Altitude (m)')
        ax12.set_ylabel('LWC (g/m3)')
        lines = [phandle11,phandle12]
        ax11.legend(lines, [l.get_label() for l in lines],loc=(1.15,0))
        ### Rtime module (to show times in a hour,minute format instead of fractions of day)
        if Rtime==None: pass
        else:
            locs, labels= xticks()
            Alabels=list()
            for l in locs:
                Rhr=l*24; Rhour=int(floor(Rhr)); 
                Rmt=(Rhr-Rhour)*60; Rminute=int(floor(Rmt)); 
                Rsd=(Rmt-Rminute)*60; Rsecond=int(floor(Rsd));
                Alabels.append(dt.datetime(2000,1,1,Rhour,Rminute,Rsecond).strftime('%H:%M'))
            Alabels[0]=''; Alabels[-1]='';
            ax11.set_xticks(locs, minor=False)
            ax11.set_xticklabels(Alabels)
        ### End of Rtime module
        
        ax21=fig.add_subplot(412, sharex=ax11)
        ix=nonzero((CDroptime>=ts1)*(CDroptime<=ts2))[0]
        if sum(~isnan(CDroptime[ix]))==0:
            phandle21,=ax21.plot(ts-dayn,np.array([1,1]),'w.',label=('%s: \ncm-3' % CDropName))
        else:
            phandle21,=ax21.plot(CDroptime[ix]-dayn,CDrop[ix],'c-',label=('%s: \n%.0f cm-3' % (CDropName,st.nanmean(CDrop[ix]))))
        ax22=ax21.twinx()
        ix=nonzero((Aerotime>=ts1)*(Aerotime<=ts2))[0]
        if sum(~isnan(Aerotime[ix]))==0:
            phandle22,=ax22.plot(ts-dayn,np.array([1,1]),'w.',label=('%s: cm-3' % AeroName))
        else:
            phandle22,=ax22.plot(Aerotime[ix]-dayn,Aero[ix],'m-',label=('%s: %.0f cm-3' % (AeroName,st.nanmean(Aero[ix]))))
        #ax21.set_xlabel('Time')
        ax21.set_ylabel('Nd (cm-3)')
        ax22.set_ylabel('Na (cm-3)')
        lines = [phandle21,phandle22]
        ax21.legend(lines, [l.get_label() for l in lines],loc=(1.15,0))
        ### Rtime module (to show times in a hour,minute format instead of fractions of day)
        if Rtime==None: pass
        else:
            locs, labels= xticks()
            Alabels=list()
            for l in locs:
                Rhr=l*24; Rhour=int(floor(Rhr)); 
                Rmt=(Rhr-Rhour)*60; Rminute=int(floor(Rmt)); 
                Rsd=(Rmt-Rminute)*60; Rsecond=int(floor(Rsd));
                Alabels.append(dt.datetime(2000,1,1,Rhour,Rminute,Rsecond).strftime('%H:%M'))
            Alabels[0]=''; Alabels[-1]='';
            ax21.set_xticks(locs, minor=False)
            ax21.set_xticklabels(Alabels)
        ### End of Rtime module
        
        ax31=fig.add_subplot(413, sharex=ax11)
        ix=nonzero((ThQtime>=ts1)*(ThQtime<=ts2))[0]
        if sum(~isnan(ThQtime[ix]))==0:
            phandle31,=ax31.plot(ts-dayn,np.array([1,1]),'w.',label=('ThetaQ:K'))
        else:
            phandle31,=ax31.plot(ThQtime[ix]-dayn,ThQ[ix],'c-',label=('ThetaQ: %.2f K' % st.nanmean(ThQ[ix])))
        ax32=ax31.twinx()
        ix=nonzero((Ttime>=ts1)*(Ttime<=ts2))[0]
        if sum(~isnan(Ttime[ix]))==0:
            phandle32,=ax32.plot(ts-dayn,np.array([1,1]),'w.',label=('Temperature: \nC'))
        else:
            phandle32,=ax32.plot(Ttime[ix]-dayn,T[ix],'m-',label=('Temperature: \n%0.2f C' % st.nanmean(T[ix])))
        #ax31.set_xlabel('Time')
        ax31.set_ylabel('Theta-Q (K)')
        ax32.set_ylabel('Temp (C)')
        lines = [phandle31,phandle32]
        ax31.legend(lines, [l.get_label() for l in lines],loc=(1.15,0))
        ### Rtime module (to show times in a hour,minute format instead of fractions of day)
        if Rtime==None: pass
        else:
            locs, labels= xticks()
            Alabels=list()
            for l in locs:
                Rhr=l*24; Rhour=int(floor(Rhr)); 
                Rmt=(Rhr-Rhour)*60; Rminute=int(floor(Rmt)); 
                Rsd=(Rmt-Rminute)*60; Rsecond=int(floor(Rsd));
                Alabels.append(dt.datetime(2000,1,1,Rhour,Rminute,Rsecond).strftime('%H:%M'))
            Alabels[0]=''; Alabels[-1]='';
            ax31.set_xticks(locs, minor=False)
            ax31.set_xticklabels(Alabels)
        ### End of Rtime module
        
        ax41=fig.add_subplot(414, sharex=ax11)
        ix=nonzero((udveltime>=ts1)*(udveltime<=ts2))[0]
        if sum(~isnan(udveltime[ix]))==0:
            phandle41,=ax41.plot(ts-dayn,np.array([1,1]),'w.',label=("w\': m/s"))
        else:
            phandle41,=ax41.plot(udveltime[ix]-dayn,udvel[ix],'c-',label=("w\': %.2f m/s" % st.nanmean(udvel[ix])))
        ax42=ax41.twinx()
        if sum(~isnan(udveltime[ix]))==0:
            phandle42,=ax42.plot(ts-dayn,np.array([1,1]),'w.',label=("std(w'): m/s"))
        else:
            phandle42,=ax42.plot(udveltime[ix]-dayn,turb[1][ix],'m-',label=("std(w'): %.2f m/s" % st.nanmean(turb[1][ix])))
        ax41.set_xlabel('Time')
        ax41.set_ylabel('Uddraft velocity (m/s)')
        ax42.set_ylabel("""Turbulence (std(w'), m/s)""")
        lines = [phandle41,phandle42]
        ax41.legend(lines, [l.get_label() for l in lines],loc=(1.15,0))
        ### Rtime module (to show times in a hour,minute format instead of fractions of day)
        if Rtime==None: pass
        else:
            locs, labels= xticks()
            Alabels=list()
            for l in locs:
                Rhr=l*24; Rhour=int(floor(Rhr)); 
                Rmt=(Rhr-Rhour)*60; Rminute=int(floor(Rmt)); 
                Rsd=(Rmt-Rminute)*60; Rsecond=int(floor(Rsd));
                Alabels.append(dt.datetime(2000,1,1,Rhour,Rminute,Rsecond).strftime('%H:%M'))
            Alabels[0]=''; Alabels[-1]='';
            ax41.set_xticks(locs, minor=False)
            ax41.set_xticklabels(Alabels)
            ax41.axis()[0:2], ax42.axis()[0:2]      # for a mysterious reason, this seems to be necessary to get the axes aligned. Go figure.
        ### End of Rtime module

        fig.text(0.35, 0.92, "%s @ %s - %s \n(horiz. scan #%d)" % (CL.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=CL.data[post][0][0])).strftime('%H:%M'), CL.desc["humanplace"], v),horizontalalignment='center')
    if interact==0:
        if np.shape(CL.times['horicloud'])[0]>=1: 
            print("The sequence may be blocked until you close all figures.")
            plt.show()
        else: print("[hrzplot] No horizontal scan was measured in this cloud.")

