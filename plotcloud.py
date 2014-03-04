
#################################################################################
##    Copyright 2013 Stephanie Gagne, Landan MacDonald                         ##
##                                                                             ##
##    This file is part of SAMAC.                                              ##
##                                                                             ##
##    SAMAC is free software: you can redistribute it and/or modify            ##
##    it under the terms of the GNU General Public License as published by     ##
##    the Free Software Foundation, either version 3 of the License, or        ##
##    (at your option) any later version.                                      ##
##                                                                             ##
##    SAMAC is distributed in the hope that it will be useful,                 ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of           ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            ##
##    GNU General Public License for more details.                             ##
##                                                                             ##
##    You should have received a copy of the GNU General Public License        ##
##    along with SAMAC.  If not, see <http://www.gnu.org/licenses/>.           ##
#################################################################################


#In this file: plotcloud2, plotcloud, maplwc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import datetime as dt


def plotcloud2(tim,alt,lat,lon,lwc,lwctime,aero,aerotime,cld,cldtime,pr,TakenSDs,cZs,Adj,Offsts,figi,Qs,proftimes=None,Rtime=None):
    """tim: times
       alt: altitude
       lat: latitude
       lon: longitude
       lwc: liquid water content
       lwctime: time associated with lwc
       aero: aerosol concentration
       aerotime: time associated with aero
       cld: cloud drop concentration
       cldtime: time associated with cld
       pr: precipitation data dictionary (pr["100"],pr["200"],pr["time100"],pr["time200"])
       TakenSDs: units to the data
       cZs: includes Zx1,Zx2,Zy11,Zy21,Zy12,Zy22 in an array
       Adj: multipliers
       Offsts: offsets
       figi: figure window number
       Qs: 0 if adjustments are needed, 1 to just display
       proftimes: contains the CloudObj.times dictionary to add a profile colour bar
       Rtime: None for ordinal time, anything else for time in HH:MM format"""
    import copy as copy
    plt.ion()         # turn interactive figure mode on
    
    # making copies
    tim=copy.copy(tim); alt=copy.copy(alt);lat=copy.copy(lat);lon=copy.copy(lon);lwc=copy.copy(lwc);lwctime=copy.copy(lwctime);aero=copy.copy(aero);aerotime=copy.copy(aerotime);cld=copy.copy(cld);cldtime=copy.copy(cldtime);pr=copy.deepcopy(pr);Zs=copy.copy(cZs);
    Zorigtime=np.floor(Zs)
    Zs=Zs-np.floor(Zs)
    Zx1=Zs[0]; Zx2=Zs[1]; Zy11=Zs[2]; Zy21=Zs[3]; Zy12=Zs[4]; Zy22=Zs[5]; Zy12L=Zs[6]; Zy22L=Zs[7];
    #[tim, alt, lwc, ptc, dropc, lat, lon] = mcols
    origZ=copy.deepcopy([Zx1, Zx2, Zy11, Zy21, Zy12, Zy22, Zy12L, Zy22L]);
                
    
    # cleaning the data
    # The plane's GPS doesn't use negative values.
    lat=np.ma.masked_outside(lat, -90, 90)
    lon=np.ma.masked_outside(lon, -180, 180)
    aero[aero<=0]=1e-2; # minimum conc.
    aero=np.ma.masked_less_equal(aero,1e-2); 
    cld[cld<=0]=1e-2;
    cld=np.ma.masked_less_equal(cld,15e-2);
    lwc[lwc<=0]=5e-4; 
    alt=alt/1000.   # now in kilometers
    pr["100"][pr["100"]<=0]=1e-2; pr["200"][pr["200"]<=0]=1e-2    # minimum concentration
    #c2d=np.ma.array(c2d)
    pr["100"]=np.ma.masked_less_equal(pr["100"],1e-2);
    pr["200"]=np.ma.masked_less_equal(pr["200"],1e-2);
    daytim=tim-np.floor(tim)
    lwctime=lwctime-np.floor(lwctime)
    aerotime=aerotime-np.floor(aerotime)
    cldtime=cldtime-np.floor(cldtime)
    pr["time100"]=pr["time100"]-np.floor(pr["time100"])
    pr["time200"]=pr["time200"]-np.floor(pr["time200"])    
        
    PGen='ask'
    while PGen=='ask':      # while we want to modify the picture (zoom, magnifying values)
        fig=plt.figure(figi, figsize=(12,6))
        fig.subplots_adjust(right=0.75,left=0.1)        # Way to make space for the legend.
        ax1=fig.add_subplot(211)    # MMM use this whenever available in pyplot: ax1=subplot2grid((2,3), (0, 0), colspan=2)
        ax2=fig.add_subplot(212)
        ############## Altitude/Ptcl conc/droplet conc ##############
        # to change to serveral independent axes: http://matplotlib.sourceforge.net/examples/pylab_examples/multiple_yaxis_with_spines.html
        plt.axes(ax1)
        phandle1=ax1.plot(daytim,alt*Adj[0]+Offsts[0],'g-',aerotime,np.log10(aero*Adj[1]+Offsts[1]),'b*',cldtime,np.log10(cld*Adj[2]+Offsts[2]),'k^',lwctime,np.log10(lwc*Adj[3]+Offsts[3]),'m-',pr["time100"],np.log10(pr["100"]*Adj[4]+Offsts[4]),'c.',pr["time200"],np.log10(pr["200"]*Adj[5]+Offsts[5]),'co')
        plt.setp(phandle1,markeredgewidth=0,markersize=5,linewidth=2)
        if np.isnan(Zy11) and np.isnan(Zy21): ax1.set_xlim(Zx1,Zx2)
        else:     ax1.axis([Zx1,Zx2,Zy11,Zy21])
        ax1.legend( ['Altitude (km)','Ptcl conc (log/%s)' %(TakenSDs[0][1]),'Dropl conc (log/%s)' %(TakenSDs[1][1]),'LWC (log/%s)' %(TakenSDs[-1][1]),'Drizzle>100 (log/%s)' %(TakenSDs[2][1]),'Drizzle>200 (log/%s)' %(TakenSDs[3][1])], loc=(1.01,0) )
        [Zx1, Zx2, Zy11, Zy21]=plt.axis()
        # time tracker bar
        ax1.scatter(daytim,Zy11*np.ones(np.shape(daytim)),c=daytim,s=40, marker='o',cmap='spectral',edgecolors='none')
        ax1.axis([Zx1,Zx2,Zy11,Zy21])
        xm1=np.mean([Zx1,Zx2]); ym1=np.mean([Zy11,Zy21]);
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
            ax1.set_xticks(locs, minor=False)
            ax1.set_xticklabels(Alabels)
        ### End of Rtime module
            

        ############## Latitude/Longitude ##############
        plt.axes(ax2)
        lns1=ax2.plot(daytim,lat,'ro', markeredgewidth=0, markersize=3.5, label='Latitude')
        plt.axes(ax2)
        if np.isnan(Zy12) and np.isnan(Zy22): ax2.set_xlim(Zx1,Zx2)
        else:     ax2.axis([Zx1,Zx2,Zy12,Zy22])
        ax22=ax2.twinx()
        lns2=plt.plot(daytim,lon,'bo', markeredgewidth=0, markersize=3.5, label='Longitude')
        ax2.grid(); 
        ax2.set_xlabel("Time"); ax2.set_ylabel("Latitude"); ax22.set_ylabel("Longitude")
        plt.axes(ax22)
        if np.isnan(Zy12L) and np.isnan(Zy22L): ax22.set_xlim(Zx1,Zx2)
        else:     ax22.axis([Zx1,Zx2,Zy12L,Zy22L])
        lns=lns1+lns2
        labs=[l.get_label() for l in lns]
        plt.legend(lns, labs, loc=(1.10,0) )
        if np.isnan(Zy12L): 
            [Zx1, Zx2, Zy12, Zy22]=ax2.axis()
            [Zy12L, Zy22L]=ax22.get_ylim()
        xm2=np.mean([Zx1,Zx2]); ym2=np.mean([Zy12,Zy22]);
        ax2.annotate("Alt M=%s; O=%s\nPtcl M=%s; O=%s\nDrop M=%s; O=%s\nLWC M=%s; O=%s\n2d100 M=%s; O=%s\n2d200 M=%s; O=%s" % (Adj[0],Offsts[0],Adj[1],Offsts[1],Adj[2],Offsts[2],Adj[3],Offsts[3],Adj[4],Offsts[4],Adj[5],Offsts[5]),(1.12, 1.10),xycoords="axes fraction", va="top", ha="left",bbox=dict(boxstyle="round, pad=0.5", fc="w"))
        # Adding the manoeuvre information
        if proftimes==None: pass
        else:
            for i in range(len(proftimes["abovecloud"])):
                ax2.plot(proftimes["abovecloud"][i,:]-np.floor(proftimes["abovecloud"][i,:]),[Zy22,Zy22],'y-',linewidth=10.0)
            for i in range(len(proftimes["belowcloud"])):
                ax2.plot(proftimes["belowcloud"][i,:]-np.floor(proftimes["belowcloud"][i,:]),[Zy22,Zy22],'c-',linewidth=10.0)
            for i in range(len(proftimes["horicloud"])):
                ax2.plot(proftimes["horicloud"][i,:]-np.floor(proftimes["horicloud"][i,:]),[Zy22,Zy22],'k-',linewidth=10.0)
            for i in range(len(proftimes["verticloud"])):
                ax2.plot(proftimes["verticloud"][i,:]-np.floor(proftimes["verticloud"][i,:]),[Zy22,Zy22],'m-',linewidth=10.0)
        # return to regular plotting
        ax2.set_xlim(Zx1,Zx2)
        ax22.set_xlim(Zx1,Zx2)
        ax2.set_ylim(Zy12,Zy22)
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
            ax2.set_xticks(locs, minor=False)
            ax2.set_xticklabels(Alabels)
            [RZx1, RZx2, RZy12, RZy22]=ax2.axis()
            ax22.set_xlim(RZx1,RZx2)
        ### End of Rtime module
        
        if Qs==1:   
            ############## Adjustments ##############
            Optie=1
            while Optie==1:     # Optie=1 as long as we don't know how to modify (or not) the plot, 0 means we are going to generate the new figure
                Pquest=raw_input("Do you wish to modify the plot? \n Zx for time zoom-in, Zy for y-axis zoom-in, Zout to return to the original size, M to change the magnifying factors, Of for offsets, and OK otherwise:\n").lower()
                # MMM also give the possibility to go back to default values if someone zoomed in too much
                if Pquest.lower()=='zx':        ####### Zoom in time #######
                    ZxQ='n'
                    print("Click twice to zoom in (time)")
                    while ZxQ=='n':
                        Lims=[]
                        Lims.append(list(plt.ginput(1,timeout=0)))
                        ax1.plot(Lims[0][0][0],ym1,'ro',ms=6)
                        ax1.axis([Zx1,Zx2,Zy11,Zy21])
                        ax2.plot(Lims[0][0][0],ym2,'ro',ms=6)
                        ax2.axis([Zx1,Zx2,Zy12,Zy22])
                        Lims.append(list(plt.ginput(1,timeout=0)))
                        ax1.plot(Lims[1][0][0],ym1,'ro',ms=6)
                        ax1.axis([Zx1,Zx2,Zy11,Zy21])
                        ax2.plot(Lims[1][0][0],ym2,'ro',ms=6)
                        ax2.axis([Zx1,Zx2,Zy12,Zy22])
                        plt.draw()
                        ZxQ=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if ZxQ=='n': print("Try again: Click twice. It ain't that hard")
                    # end of while ZxQ loop
                    if Lims[0][0][0]<=Lims[1][0][0]: Zx1=Lims[0][0][0]; Zx2=Lims[1][0][0];
                    else: Zx2=Lims[0][0][0]; Zx1=Lims[1][0][0];
                    plt.close(figi)
                    Optie=0
                elif Pquest=='zy':          ###### Zoom of the y axis ######
                    ZyQ='n'
                    SubQ=raw_input("Which subplot do you want to adjust? (1/2)  ")
                    if SubQ=='1' or SubQ=='2': pass
                    else: print("This is not a choice. The default is 2. \nClick twice on the bottom subplot to zoom in (y-axis)");    SubQ='2';
                    while ZyQ=='n':
                        Lims=[]
                        if SubQ=='1':       ### Subplot(211)
                            print("Click twice on the upper subplot to zoom in (y-axis) ")
                            Lims.append(list(plt.ginput(1,timeout=0)))
                            ax1.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                            ax1.axis([Zx1,Zx2,Zy11,Zy21])
                            Lims.append(list(plt.ginput(1,timeout=0)))
                            ax1.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                            ax1.axis([Zx1,Zx2,Zy11,Zy21])
                            plt.draw()
                        elif SubQ=='2':     ### Subplot(212)
                            LLQ=0
                            while LLQ==0:
                                LLQ=raw_input("Are you adjusting the longitude or Latitude? (Lat/Lon)  ").lower()
                                print("Click twice on the bottom subplot to zoom in (y-axis) ")
                                plt.figure(figi+1)
                                if LLQ=='lat':
                                    #plt.axes(ax2)
                                    plt.plot(daytim,lat,'ro',markeredgewidth=0, label='Latitude')
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    plt.xlabel("Time"); plt.ylabel("Latitude")
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    plt.draw()
                                elif LLQ.lower()=='lon':
                                    #plt.axes(ax22)
                                    plt.plot(daytim,lon,'bo',markeredgewidth=0,label="Longitude")
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    plt.xlabel("Time"); plt.ylabel("Longitude")
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    plt.draw()
                                else: print("This is not a choice. Try again"); LLQ=0;
                        ZyQ=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if ZyQ=='n': print("Try again: Click twice. It ain't that hard")
                        else: plt.close(figi+1)
                        # end of while ZyQ loop
                        if SubQ=='1': 
                            if Lims[0][0][1]<=Lims[1][0][1]: Zy11=Lims[0][0][1]; Zy21=Lims[1][0][1];
                            else: Zy21=Lims[0][0][1]; Zy11=Lims[1][0][1];
                        elif SubQ=='2': 
                            if LLQ=='lat':
                                if Lims[0][0][1]<=Lims[1][0][1]: Zy12=Lims[0][0][1]; Zy22=Lims[1][0][1];
                                else: Zy22=Lims[0][0][1]; Zy12=Lims[1][0][1];
                            elif LLQ=='lon':
                                if Lims[0][0][1]<=Lims[1][0][1]: Zy12L=Lims[0][0][1]; Zy22L=Lims[1][0][1];
                                else: Zy22L=Lims[0][0][1]; Zy12L=Lims[1][0][1];
                    plt.close(figi)
                    Optie=0
                elif Pquest=='zout':          ###### Zoom out all ######
                    [Zx1, Zx2, Zy11, Zy21, Zy12, Zy22, Zy12L, Zy22L]=origZ
                    plt.close(figi)
                    Optie=0
                elif Pquest=='m':           ###### Magnifying factors ######
                    print("The magnifying factors are:\n 1-altitude: %f,\n 2-Ptcl. conc.: %f,\n 3-Droplet conc.: %f\n 4-LWC: %f,\n 5-2d100: %f,\n 6-2d200: %f\n" % (Adj[0],Adj[1],Adj[2],Adj[3],Adj[4],Adj[5]))
                    ix=int(raw_input("Which factor do you want to change? (1 to 6) "))
                    if ix==1 or ix==2 or ix==3 or ix==4 or ix==5 or ix==6:
                        Adj[ix-1]=float(raw_input("You want to replace %f by:  " % Adj[ix-1]))
                        plt.close(figi)
                    else:
                        print("This is not an option. Try again!")
                    Optie=0
                    Zy11=np.min([np.min(np.log10(alt*Adj[0]+Offsts[0])),np.min(np.log10(aero*Adj[1]+Offsts[1])),np.min(np.log10(cld*Adj[2]+Offsts[2])),np.min(np.log10(lwc*Adj[3]+Offsts[3])),np.min(np.log10(pr["100"]*Adj[4]+Offsts[4])),np.min(np.log10(pr["200"]*Adj[5]+Offsts[5]))]) 
                    Zy21=np.max([np.max(alt*Adj[0]+Offsts[0]),np.max(np.log10(aero*Adj[1]+Offsts[1])),np.max(np.log10(cld*Adj[2]+Offsts[2])),np.max(np.log10(lwc*Adj[3]+Offsts[3])),np.max(np.log10(pr["100"]*Adj[4]+Offsts[4])),np.max(np.log10(pr["200"]*Adj[5]+Offsts[5]))])
                elif Pquest=='of':           ###### Offsets ######
                    print("The offsets are:\n 1-altitude: %f,\n 2-Ptcl. conc.: %f,\n 3-Droplet conc.: %f\n 4-LWC: %f,\n 5-2d100: %f,\n 6-2d200: %f\n" % (Offsts[0],Offsts[1],Offsts[2],Offsts[3],Offsts[4],Offsts[5]))
                    ix=int(raw_input("Which offset do you want to change? (1 to 6) "))
                    if ix==1 or ix==2 or ix==3 or ix==4 or ix==5 or ix==6:
                        Offsts[ix-1]=float(raw_input("You want to replace %f by:  " % Offsts[ix-1]))
                        plt.close(figi)
                    else:
                        print("This is not an option. Try again!")
                    Optie=0

                elif Pquest.lower()=='mf':        ####### add pts to map #######
                    mfq='n'
                    print("Click on a point in time to highlight it on the map")
                    while mfq=='n':
                        Lims=[]
                        Lims.append(list(fig.ginput(1,timeout=0)))
                        tim1=Lims[0][0][0]
                        tol=1.E20
                        for i in range(len(tim)):
                            time=tim[i] % 1
                            diff=abs(time-tim1)
                            if diff<tol:
                                tol=diff
                                minindex=i
                        from mpl_toolkits.basemap import Basemap
                        fig2=plt.figure(1002)
                        ax2=fig2.add_subplot(122)
                        llon=np.nanmin(lon); hlon=np.nanmax(lon)
                        llat=np.nanmin(lat); hlat=np.nanmax(lat)
                        m = Basemap(llcrnrlon=llon, llcrnrlat=llat, urcrnrlon=hlon, urcrnrlat=hlat, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=10000)
                        x, y = m(lon[minindex]*-1.,lat[minindex])
                        ax2.scatter(x,y,s=50, marker='o',color='r',edgecolors='none')
                        mfq=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if mfq=='n':
                            print("Click on a point in time to highlight it on the map")
                            pass
                        # end of while ZxQ loop
                    if Lims[0][0][0]<=Lims[1][0][0]: Zx1=Lims[0][0][0]; Zx2=Lims[1][0][0];
                    else: Zx2=Lims[0][0][0]; Zx1=Lims[1][0][0];
                    plt.close(figi)
                    Optie=0






                elif Pquest=='ok':          ###### The figure is fine ######
                    PGen='dontaskanymore'
                    Optie=0
                else:                       ###### Wrong Entry ######
                    print("This is not an option. Try again!")
                    Optie=1
            # end of while Optie loop
        else: PGen='dontaskanymore'; Optie=0;
        # end of Qs loop
    # end of the while PGen loop

    Zs=np.array([Zx1,Zx2,Zy11,Zy21,Zy12,Zy22,Zy12L,Zy22L])+Zorigtime
    return [Zs,Adj,Offsts,figi,fig]
    

def plotcloud(cmall,mcols,c2d,cZs,Adj,Offsts,figi,Qs,proftimes=None,Rtime=None):
# mall is the data
# mcols are the columns from which to take mall
# Zs includes Zx1,Zx2,Zy11,Zy21,Zy12,Zy22 in an array
# this version is used with Read_CreateClouds.py
    """cmall is the data table, 
    mcols is the relevant columns, 
    c2d is the concentration of 2d data, 
    cZs is the zooming parameters, 
    Adj is the multipliers, 
    Offsts is the offsets, 
    figi is the figure window number in which the figure will open, 
    Qs is 0 if adjustments are needed, 1 to just display, 
    proftimes contains the CloudObj.times dictionary to add a profile colour bar, 
    RTime is specified/not None if you want the time displayed in hours:minute formats."""
    import copy as copy
#    from matplotlib import rc
#    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#    rc('text', usetex=True)
    plt.ion()         # turn interactive figure mode on
    
    mall=copy.copy(cmall); Zs=copy.copy(cZs);
    Zorigtime=np.floor(Zs)
    Zs=Zs-np.floor(Zs)
    Zx1=Zs[0]; Zx2=Zs[1]; Zy11=Zs[2]; Zy21=Zs[3]; Zy12=Zs[4]; Zy22=Zs[5]; Zy12L=Zs[6]; Zy22L=Zs[7];
    [tim, alt, lwc, ptc, dropc, lat, lon] = mcols
    origZ=copy.deepcopy([Zx1, Zx2, Zy11, Zy21, Zy12, Zy22, Zy12L, Zy22L]);
                
    
    # cleaning the data
    mall=np.ma.array(mall)
    # The plane's GPS doesn't use negative values.
    mall[lat]=np.ma.masked_outside(mall[lat], -90, 90)
    mall[lon]=np.ma.masked_outside(mall[lon], -180, 180)
    mall[ptc][mall[ptc]<=0]=1e-2; 
    mall[ptc]=np.ma.masked_less_equal(mall[ptc],1e-2); 
    mall[dropc][mall[dropc]<=0]=1e-2;
    mall[dropc]=np.ma.masked_less_equal(mall[dropc],15e-2);
    c2d[1:][c2d[1:]<=0]=1e-2;  
    c2d=np.ma.array(c2d)
    c2d[1]=np.ma.masked_less_equal(c2d[1],1e-2);
    c2d[2]=np.ma.masked_less_equal(c2d[2],1e-2);
    mall[lwc][mall[lwc]<=0]=5e-4; 
    mall[alt]=mall[alt]/1000.
    malltim=mall[tim]-np.floor(mall[tim])
    c2dtim=c2d[0]-np.floor(c2d[0])
        
    PGen='ask'
    while PGen=='ask':      # while we want to modify the picture (zoom, magnifying values)
        fig=plt.figure(figi, figsize=(12,6))
        fig.subplots_adjust(right=0.75,left=0.1)        # Way to make space for the legend.
        ax1=fig.add_subplot(211)    # MMM use this whenever available in pyplot: ax1=subplot2grid((2,3), (0, 0), colspan=2)
        ax2=fig.add_subplot(212)
        ############## Altitude/Ptcl conc/droplet conc ##############
        # to change to serveral independent axes: http://matplotlib.sourceforge.net/examples/pylab_examples/multiple_yaxis_with_spines.html
        plt.axes(ax1)
        phandle1=ax1.plot(malltim,mall[alt]*Adj[0]+Offsts[0],'g-',malltim,np.log10(mall[ptc]*Adj[1]+Offsts[1]),'b*',malltim,np.log10(mall[dropc]*Adj[2]+Offsts[2]),'k^',malltim,np.log10(mall[lwc]*Adj[3]+Offsts[3]),'m-',c2dtim,np.log10(c2d[1]*Adj[4]+Offsts[4]),'c.',c2dtim,np.log10(c2d[2]*Adj[5]+Offsts[5]),'co')
        plt.setp(phandle1,markeredgewidth=0,markersize=5,linewidth=2)
        if np.isnan(Zy11) and np.isnan(Zy21): ax1.set_xlim(Zx1,Zx2)
        else:     ax1.axis([Zx1,Zx2,Zy11,Zy21])
        ax1.legend( ['Altitude (km)','Ptcl conc (log/cm-3)','Dropl conc (log/cm-3)','LWC (log)','Drizzle>100 (log/l-1)','Drizzle>200 (log/l-1)'], loc=(1.01,0) )
        [Zx1, Zx2, Zy11, Zy21]=plt.axis()
        ax1.scatter(malltim,Zy11*np.ones(np.shape(malltim)),c=malltim,s=40, marker='o',cmap='spectral',edgecolors='none')
        ax1.axis([Zx1,Zx2,Zy11,Zy21])
        xm1=np.mean([Zx1,Zx2]); ym1=np.mean([Zy11,Zy21]);
        #ax1.annotate("Alt M=%s; O=%s\nPtcl M=%s; O=%s\nDrop M=%s; O=%s\nLWC M=%s; O=%s\n2d100 M=%s; O=%s\n2d200 M=%s; O=%s" % (Adj[0],Offsts[0],Adj[1],Offsts[1],Adj[2],Offsts[2],Adj[3],Offsts[3],Adj[4],Offsts[4],Adj[5],Offsts[5]),(1.017, 1.00),xycoords="axes fraction", va="top", ha="left",bbox=dict(boxstyle="round, pad=0.5", fc="w"))
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
            ax1.set_xticks(locs, minor=False)
            ax1.set_xticklabels(Alabels)
        ### End of Rtime module
            

        ############## Latitude/Longitude ##############
        plt.axes(ax2)
        lns1=ax2.plot(malltim,mall[lat],'ro', markeredgewidth=0, markersize=3.5, label='Latitude')
        plt.axes(ax2)
        if np.isnan(Zy12) and np.isnan(Zy22): ax2.set_xlim(Zx1,Zx2)
        else:     ax2.axis([Zx1,Zx2,Zy12,Zy22])
        ax22=ax2.twinx()
        lns2=plt.plot(malltim,mall[lon],'bo', markeredgewidth=0, markersize=3.5, label='Longitude')
        ax2.grid(); 
        ax2.set_xlabel("Time"); ax2.set_ylabel("Latitude"); ax22.set_ylabel("Longitude")
        plt.axes(ax22)
        if np.isnan(Zy12L) and np.isnan(Zy22L): ax22.set_xlim(Zx1,Zx2)
        else:     ax22.axis([Zx1,Zx2,Zy12L,Zy22L])
        lns=lns1+lns2
        labs=[l.get_label() for l in lns]
        plt.legend(lns, labs, loc=(1.10,0) )
        if np.isnan(Zy12L): 
            [Zx1, Zx2, Zy12, Zy22]=ax2.axis()
            [Zy12L, Zy22L]=ax22.get_ylim()
        xm2=np.mean([Zx1,Zx2]); ym2=np.mean([Zy12,Zy22]);
        ax2.annotate("Alt M=%s; O=%s\nPtcl M=%s; O=%s\nDrop M=%s; O=%s\nLWC M=%s; O=%s\n2d100 M=%s; O=%s\n2d200 M=%s; O=%s" % (Adj[0],Offsts[0],Adj[1],Offsts[1],Adj[2],Offsts[2],Adj[3],Offsts[3],Adj[4],Offsts[4],Adj[5],Offsts[5]),(1.12, 1.10),xycoords="axes fraction", va="top", ha="left",bbox=dict(boxstyle="round, pad=0.5", fc="w"))
        # Adding the profile information #
        if proftimes==None: pass
        else:
            for i in range(len(proftimes["abovecloud"])):
                ax2.plot(proftimes["abovecloud"][i,:]-np.floor(proftimes["abovecloud"][i,:]),[Zy22,Zy22],'y-',linewidth=10.0)
            for i in range(len(proftimes["belowcloud"])):
                ax2.plot(proftimes["belowcloud"][i,:]-np.floor(proftimes["belowcloud"][i,:]),[Zy22,Zy22],'c-',linewidth=10.0)
            for i in range(len(proftimes["horicloud"])):
                ax2.plot(proftimes["horicloud"][i,:]-np.floor(proftimes["horicloud"][i,:]),[Zy22,Zy22],'k-',linewidth=10.0)
            for i in range(len(proftimes["verticloud"])):
                ax2.plot(proftimes["verticloud"][i,:]-np.floor(proftimes["verticloud"][i,:]),[Zy22,Zy22],'m-',linewidth=10.0)
        # return to regular plotting
        ax2.set_xlim(Zx1,Zx2)
        ax22.set_xlim(Zx1,Zx2)
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
            ax2.set_xticks(locs, minor=False)
            ax2.set_xticklabels(Alabels)
            [RZx1, RZx2, RZy12, RZy22]=ax2.axis()
            ax22.set_xlim(RZx1,RZx2)
        ### End of Rtime module
        
        if Qs==1:   
            ############## Adjustments ##############
            Optie=1
            while Optie==1:     # Optie=1 as long as we don't know how to modify (or not) the plot, 0 means we are going to generate the new figure
                Pquest=raw_input("Do you wish to modify the plot? \n Zx for time zoom-in, Zy for y-axis zoom-in, Zout to return to the original size, M to change the magnifying factors, Of for offsets, and OK otherwise:\n").lower()
                # MMM also give the possibility to go back to default values if someone zoomed in too much
                if Pquest.lower()=='zx':        ####### Zoom in time #######
                    ZxQ='n'
                    print("Click twice to zoom in (time)")
                    while ZxQ=='n':
                        Lims=[]
                        Lims.append(list(plt.ginput(1,timeout=0)))
                        ax1.plot(Lims[0][0][0],ym1,'ro',ms=6)
                        ax1.axis([Zx1,Zx2,Zy11,Zy21])
                        ax2.plot(Lims[0][0][0],ym2,'ro',ms=6)
                        ax2.axis([Zx1,Zx2,Zy12,Zy22])
                        Lims.append(list(plt.ginput(1,timeout=0)))
                        ax1.plot(Lims[1][0][0],ym1,'ro',ms=6)
                        ax1.axis([Zx1,Zx2,Zy11,Zy21])
                        ax2.plot(Lims[1][0][0],ym2,'ro',ms=6)
                        ax2.axis([Zx1,Zx2,Zy12,Zy22])
                        plt.draw()
                        ZxQ=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if ZxQ=='n': print("Try again: Click twice. It ain't that hard")
                    # end of while ZxQ loop
                    if Lims[0][0][0]<=Lims[1][0][0]: Zx1=Lims[0][0][0]; Zx2=Lims[1][0][0];
                    else: Zx2=Lims[0][0][0]; Zx1=Lims[1][0][0];
                    plt.close(figi)
                    Optie=0
                elif Pquest=='zy':          ###### Zoom of the y axis ######
                    ZyQ='n'
                    SubQ=raw_input("Which subplot do you want to adjust? (1/2)  ")
                    if SubQ=='1' or SubQ=='2': pass
                    else: print("This is not a choice. The default is 2. \nClick twice on the bottom subplot to zoom in (y-axis)");    SubQ='2';
                    while ZyQ=='n':
                        Lims=[]
                        if SubQ=='1':       ### Subplot(211)
                            print("Click twice on the upper subplot to zoom in (y-axis) ")
                            Lims.append(list(plt.ginput(1,timeout=0)))
                            ax1.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                            ax1.axis([Zx1,Zx2,Zy11,Zy21])
                            Lims.append(list(plt.ginput(1,timeout=0)))
                            ax1.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                            ax1.axis([Zx1,Zx2,Zy11,Zy21])
                            plt.draw()
                        elif SubQ=='2':     ### Subplot(212)
                            LLQ=0
                            while LLQ==0:
                                LLQ=raw_input("Are you adjusting the longitude or Latitude? (Lat/Lon)  ").lower()
                                print("Click twice on the bottom subplot to zoom in (y-axis) ")
                                plt.figure(figi+1)
                                if LLQ=='lat':
                                    #plt.axes(ax2)
                                    plt.plot(malltim,mall[lat],'ro',markeredgewidth=0, label='Latitude')
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    plt.xlabel("Time"); plt.ylabel("Latitude")
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12,Zy22])
                                    plt.draw()
                                elif LLQ.lower()=='lon':
                                    #plt.axes(ax22)
                                    plt.plot(malltim,mall[lon],'bo',markeredgewidth=0,label="Longitude")
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    plt.xlabel("Time"); plt.ylabel("Longitude")
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[0][0][0],Lims[0][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    Lims.append(list(plt.ginput(1,timeout=0)))
                                    plt.plot(Lims[1][0][0],Lims[1][0][1],'ro',ms=6)
                                    plt.axis([Zx1,Zx2,Zy12L,Zy22L])
                                    plt.draw()
                                else: print("This is not a choice. Try again"); LLQ=0;
                        ZyQ=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if ZyQ=='n': print("Try again: Click twice. It ain't that hard")
                        else: plt.close(figi+1)
                        # end of while ZyQ loop
                        if SubQ=='1': 
                            if Lims[0][0][1]<=Lims[1][0][1]: Zy11=Lims[0][0][1]; Zy21=Lims[1][0][1];
                            else: Zy21=Lims[0][0][1]; Zy11=Lims[1][0][1];
                        elif SubQ=='2': 
                            if LLQ=='lat':
                                if Lims[0][0][1]<=Lims[1][0][1]: Zy12=Lims[0][0][1]; Zy22=Lims[1][0][1];
                                else: Zy22=Lims[0][0][1]; Zy12=Lims[1][0][1];
                            elif LLQ=='lon':
                                if Lims[0][0][1]<=Lims[1][0][1]: Zy12L=Lims[0][0][1]; Zy22L=Lims[1][0][1];
                                else: Zy22L=Lims[0][0][1]; Zy12L=Lims[1][0][1];
                    plt.close(figi)
                    Optie=0
                elif Pquest=='zout':          ###### Zoom out all ######
                    [Zx1, Zx2, Zy11, Zy21, Zy12, Zy22, Zy12L, Zy22L]=origZ
                    plt.close(figi)
                    Optie=0
                elif Pquest=='m':           ###### Magnifying factors ######
                    print("The magnifying factors are:\n 1-altitude: %f,\n 2-Ptcl. conc.: %f,\n 3-Droplet conc.: %f\n 4-LWC: %f,\n 5-2d100: %f,\n 6-2d200: %f\n" % (Adj[0],Adj[1],Adj[2],Adj[3],Adj[4],Adj[5]))
                    ix=int(raw_input("Which factor do you want to change? (1 to 6) "))
                    if ix==1 or ix==2 or ix==3 or ix==4 or ix==5 or ix==6:
                        Adj[ix-1]=float(raw_input("You want to replace %f by:  " % Adj[ix-1]))
                        plt.close(figi)
                    else:
                        print("This is not an option. Try again!")
                    Optie=0
                    Zy11=np.min([np.min(np.log10(mall[alt]*Adj[0]+Offsts[0])),np.min(np.log10(mall[ptc]*Adj[1]+Offsts[1])),np.min(np.log10(mall[dropc]*Adj[2]+Offsts[2])),np.min(np.log10(mall[lwc]*Adj[3]+Offsts[3])),np.min(np.log10(c2d[1]*Adj[4]+Offsts[4])),np.min(np.log10(c2d[2]*Adj[5]+Offsts[5]))]) 
                    Zy21=np.max([np.max(mall[alt]*Adj[0]+Offsts[0]),np.max(np.log10(mall[ptc]*Adj[1]+Offsts[1])),np.max(np.log10(mall[dropc]*Adj[2]+Offsts[2])),np.max(np.log10(mall[lwc]*Adj[3]+Offsts[3])),np.max(np.log10(c2d[1]*Adj[4]+Offsts[4])),np.max(np.log10(c2d[2]*Adj[5]+Offsts[5]))])
                elif Pquest=='of':           ###### Offsets ######
                    print("The offsets are:\n 1-altitude: %f,\n 2-Ptcl. conc.: %f,\n 3-Droplet conc.: %f\n 4-LWC: %f,\n 5-2d100: %f,\n 6-2d200: %f\n" % (Offsts[0],Offsts[1],Offsts[2],Offsts[3],Offsts[4],Offsts[5]))
                    ix=int(raw_input("Which offset do you want to change? (1 to 6) "))
                    if ix==1 or ix==2 or ix==3 or ix==4 or ix==5 or ix==6:
                        Offsts[ix-1]=float(raw_input("You want to replace %f by:  " % Offsts[ix-1]))
                        plt.close(figi)
                    else:
                        print("This is not an option. Try again!")
                    Optie=0

                elif Pquest.lower()=='mf':        ####### add pts to map #######
                    mfq='n'
                    print("Click on a point in time to highlight it on the map")
                    while mfq=='n':
                        Lims=[]
                        Lims.append(list(fig.ginput(1,timeout=0)))
                        tim1=Lims[0][0][0]
                        tol=1.E20
                        for i in range(len(mall[tim])):
                            time=mall[tim][i] % 1
                            diff=abs(time-tim1)
                            if diff<tol:
                                tol=diff
                                minindex=i
                        from mpl_toolkits.basemap import Basemap
                        fig2=plt.figure(1002)
                        ax2=fig2.add_subplot(122)
                        llon=np.nanmin(mall[lon]*-1.); hlon=np.nanmax(mall[lon]*-1.)
                        llat=np.nanmin(mall[lat]); hlat=np.nanmax(mall[lat])
                        m = Basemap(llcrnrlon=llon, llcrnrlat=llat, urcrnrlon=hlon, urcrnrlat=hlat, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=10000)
                        x, y = m(mall[lon][minindex]*-1.,mall[lat][minindex])
                        ax2.scatter(x,y,s=50, marker='o',color='r',edgecolors='none')
                        mfq=raw_input("Are you satisfied with your selection? ([y]/n)   ")
                        if mfq=='n':
                            print("Click on a point in time to highlight it on the map")
                            pass
                        # end of while ZxQ loop
                    if Lims[0][0][0]<=Lims[1][0][0]: Zx1=Lims[0][0][0]; Zx2=Lims[1][0][0];
                    else: Zx2=Lims[0][0][0]; Zx1=Lims[1][0][0];
                    plt.close(figi)
                    Optie=0






                elif Pquest=='ok':          ###### The figure is fine ######
                    PGen='dontaskanymore'
                    Optie=0
                else:                       ###### Wrong Entry ######
                    print("This is not an option. Try again!")
                    Optie=1
            # end of while Optie loop
        else: PGen='dontaskanymore'; Optie=0;
        # end of Qs loop
    # end of the while PGen loop

    Zs=np.array([Zx1,Zx2,Zy11,Zy21,Zy12,Zy22,Zy12L,Zy22L])+Zorigtime
    return [Zs,Adj,Offsts,figi,fig]
    

def maplwc(tim,lwc,lat,lon,imode,mf=0):
    import copy
    from math import isnan
    """ This function maps an aircraft trajectory. """
    # this plots the cloud's lwc on a map using basemap.
    # all inputs must be the same length (except for imode, mf)
    # imode=0 => interactive mode off, imode=1 => interactive mode on
    # Stephanie Gagne, Dal, March 2012
    from mpl_toolkits.basemap import Basemap
    lat=np.ma.masked_outside(lat, -90, 90)
    lon=np.ma.masked_outside(lon, -180, 180)
    if imode==1: plt.ion()
    elif imode==0: plt.ioff()
        
        # regional subplot
    fig2=plt.figure(1002)
    ax1=fig2.add_subplot(121)
    try:
        llon=np.floor(np.nanmin(lon)); hlon=np.ceil(np.nanmax(lon))
        llat=np.floor(np.nanmin(lat)); hlat=np.ceil(np.nanmax(lat))
    except:
        llon=np.floor(lon.min()); hlon=np.ceil(lon.max())
        llat=np.floor(lat.min()); hlat=np.ceil(lat.max())

    m = Basemap(llcrnrlon=llon-2, llcrnrlat=llat-2, urcrnrlon=hlon+2, urcrnrlat=hlat+2, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=10000)
    x, y = m(lon, lat)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.2)
    m.drawmapboundary()
    m.drawcountries()
    m.drawmapboundary(fill_color='lightcyan')
    m.fillcontinents(color='tan',lake_color='lightcyan')
    ax1.plot(x,y,'ro',markeredgewidth=0)
    m.drawparallels(np.arange(llat-2,hlat+2,1),labels=[1,0,0,0])
    m.drawmeridians(np.arange(llon-2,hlon+2,1),labels=[0,0,0,1])
    ax1.set_title('Situation')

    # zoomed on cloud subplot
    if mf==0:
        if len(np.shape(lwc))==2:
            lwc=lwc.data[0]
        else:
            lwc=lwc.data
    ax2=fig2.add_subplot(122)
    try:
        llon=np.nanmin(lon); hlon=np.nanmax(lon)
        llat=np.nanmin(lat); hlat=np.nanmax(lat)
    except: 
        if type(lon)==np.ma.core.MaskedArray: llon=np.min(lon); hlon=np.max(lon)
        if type(lon)==np.ma.core.MaskedArray: llat=np.min(lat); hlat=np.max(lat)
    m = Basemap(llcrnrlon=llon, llcrnrlat=llat, urcrnrlon=hlon, urcrnrlat=hlat, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=10000)
    x, y = m(lon, lat)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.2)
    m.drawmapboundary()
    m.drawcountries()
    m.drawmapboundary(fill_color='lightcyan')
    m.fillcontinents(color='tan',lake_color='lightcyan',zorder=0)
    if hlat-llat>0.5: m.drawparallels(np.arange(np.floor(llat),np.ceil(hlat),0.5),labels=[0,1,0,0])
    elif hlat-llat<=0.5: m.drawparallels(np.arange(np.floor(llat),np.ceil(hlat),0.05),labels=[0,1,0,0])
    if hlon-llon>1.: m.drawmeridians(np.arange(np.floor(llon),np.ceil(hlon),1),labels=[0,0,0,1])
    elif hlon-llon<=1.: m.drawmeridians(np.arange(np.floor(llon),np.ceil(hlon),0.1),labels=[0,0,0,1])
    for i in range(len(lwc)):
        if isnan(lwc[i]):
            lwc[i]=0.
    ax2.scatter(x,y,c=tim,s=70*lwc+1.5, marker='o',cmap='spectral',edgecolors='none')
    ax2.set_title('Flight map: \nsize=lwc \ncolor=time (black to white)')
    
    if imode==0: plt.show()





