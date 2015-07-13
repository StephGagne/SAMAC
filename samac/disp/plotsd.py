################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####

from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import interpolate

#########################################################################        
###############################  plotsd  ################################
#########################################################################
def plotsd(CL,Rtime=None):
    """ This method is used to plot the size distrubution for a selection of instruments, and scans, as a function of time.
    The instrument and scans are chosen by the user once the method has been called (interactively).
    Call CloudObj.plotsd(Rtime=1) to get the time in Hour:Minute:Second format. """
    #MMM Make those choices in the method's call instead of prompt based?
#        plt.ion()         # turn interactive figure mode on
    for i in range(len(CL.sd)):
        print("%d - %s" % (i,CL.sd[i]["Distname"]))
    IQ=raw_input("Which instruments do you wish to plot? [all] or # separated by commas:  ")
    try: IQ=map(int,IQ.split(','))
    except: IQ=range(len(CL.sd))
    for i,x in enumerate(CL.times.keys()):
        print("%d - %s (%d scans)" % (i,x,len(CL.times[x])))
    SQ=raw_input("Which type of scan do you wish to plot? [all] or # separated by commas:  ")
    print("---------------------------")
    try: SQ=map(int,SQ.split(','))
    except: SQ=range(len(CL.times.keys()))
    k=0; plt.close('all');
    for i in IQ:
        for j in SQ:
            scanq=CL.times.keys()[j]
            for h in range(np.shape(CL.times[scanq])[0]):      # for each scan of the same type (e.g. each of 2 vertical scans)
                print("[plotsd] Now drawing: %s - %s no. %d" % (CL.sd[i]["Distname"],scanq, h))
                if type(CL.sd[i]["data"])==float: print("    [plotsd] Error due to the distribution being a NaN.")
                elif np.shape(CL.sd[i]["data"])[1]==0 or np.shape(CL.sd[i]["data"])[0]==0: print("    [plotsd] Error due to the distribution being empty.")
                else:
                    chan=CL.sd[i]["bins"]
                    Zdata=copy.copy(CL.sd[i]["data"]); Zdata[Zdata<=1e-6]=1e-6;
                    if np.shape(np.nonzero((CL.sd[i]["time"]>=CL.times[scanq][h][0])*(CL.sd[i]["time"]<=CL.times[scanq][h][1]))[0])[0]==0:
                        print("    [plotsd] No data available for this period.")
                    elif np.shape(np.nonzero((CL.sd[i]["time"]>=CL.times[scanq][h][0])*(CL.sd[i]["time"]<=CL.times[scanq][h][1]))[0])[0]<=1:
                        print("    [plotsd] Only one point available for this period.")
                    else:
                        k+=1
                        fig=plt.figure(k)
                        Zdata=Zdata[:,np.nonzero((CL.sd[i]["time"]>=CL.times[scanq][h][0])*(CL.sd[i]["time"]<=CL.times[scanq][h][1]))[0]]
                        ptime=CL.sd[i]["time"]; 
                        ptime=ptime[np.nonzero((CL.sd[i]["time"]>=CL.times[scanq][h][0])*(CL.sd[i]["time"]<=CL.times[scanq][h][1]))[0]];
                        ptime=ptime-np.floor(ptime[1]); 
                        plt.pcolor(ptime,chan,np.log10(Zdata))
                        Zmax=math.ceil(np.nanmax(np.log10(Zdata)))
                        if Zmax>0:
                            plt.clim([0,Zmax])
                            cbar=plt.colorbar(orientation='horizontal',ticks=[0,Zmax/4.,Zmax/2.,3*Zmax/4.,Zmax],format=plt.FormatStrFormatter('$10^{%2.2f}$'))
                        else:
                            plt.clim([-6,0])
                            cbar=plt.colorbar(orientation='horizontal',ticks=[-6,-4.5,-3.0,-1.5,0],format=plt.FormatStrFormatter('$10^{%2.2f}$'))
                        if '2d' in CL.sd[i]["Distname"].lower():
                            cbar.set_label('Concentration ([dN/dlogDp]Litre$^{-1}$)')
                        else: cbar.set_label('Concentration ([dN/dlogDp]cm$^{-3}$)')
                        ax=plt.gca()
                        ax.set_yscale('log')
                        plt.ylabel('Diameter/channel')
                        plt.xlabel('Time (days)')
                        plt.title("%s - %s no. %d" % (CL.sd[i]["Distname"],scanq, h))
                        if scanq=='verticloud':
                            Calt=[g for g,x in enumerate(CL.dttl) if x=='altitude']
                            ct=[g for g,x in enumerate(CL.dttl) if x=='time']
                            alt=CL.data[Calt].reshape(-1,)
                            if len(ptime)==len(alt): pass
                            else: 
                                ttmp=CL.data[ct[0]]; 
                                ttmp=ttmp-np.floor(ttmp[1])
                                f=interpolate.interp1d(ttmp,np.ma.filled(alt,nan),kind='linear')   
                                if 'ma' not in str(type(alt)).lower(): alt=np.ma.array(alt,mask=False)
                                fma=interpolate.interp1d(ttmp,alt.mask,kind='linear')    # interpolating the mask
                                ptime.reshape(max(np.shape(ptime)),)
                                ptime=ptime[(ptime>=np.min(ttmp))*(ptime<=np.max(ttmp))]
                                alt=np.ma.array(f(ptime),mask=fma(ptime))   
                                del ttmp, f
                            ax2=plt.twinx()
                            ax2.plot(ptime,alt,'k',linewidth='1.5')
                            ax2.set_ylabel('Altitude (m)')
                        ax.axis([min(ptime),max(ptime),min(chan),max(chan)])
                        ### Rtime module (to show times in a hour,minute format instead of fractions of day)
                        if Rtime==None: pass
                        else:
                            locs, labels= xticks()
                            Alabels=list()
                            for l in locs:
                                Rhr=l*24; Rhour=int(floor(Rhr)); 
                                Rmt=(Rhr-Rhour)*60; Rminute=int(floor(Rmt)); 
                                Rsd=(Rmt-Rminute)*60; Rsecond=int(floor(Rsd));
                                Alabels.append(dt.datetime(2000,1,1,Rhour,Rminute,Rsecond).strftime('%H:%M:%S'))
                            Alabels[0]=''; Alabels[-1]='';
                            ax.set_xticks(locs, minor=False)
                            ax.set_xticklabels(Alabels)
                        ### End of Rtime module                    
                        try: del ptime, chan, Zmax, Zdata;
                        except: pass
                        try: del Calt, alt, ct;
                        except: pass
    if k>=1: print("--> The sequence may be blocked until you close all figures. Sometimes."); plt.show();  
    

