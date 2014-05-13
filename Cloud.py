#! /usr/bin/python

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


## LIST OF METHODS AND PROPERTIES

# methods:
# __init__(self,data,dttl,dunit,times,sd)
# __del__(self)
# overview(self,interact=0, Rtime=None)
# mapcloud(self,interact=1)
# add2basic(self,newdata,newttl,newunit)
# describe(self)
# plotsd(self,Rtime=0)
# addsizedistxl(self)
# addsizedist(self)
# vprof(self,interact=None, axtype1='log', axtype2='log')
# wholeprof(self,interact=None, axtype1='log', axtype2='log')
# defheight(self)
# defBGheight(self)
# lwcflag(self,base=bg)
# adlwc(self,base=bg)
# timechange(self)
# MaskData(self)
# avsizedist(self, prof='belowcloud', scan=0, inst='PCASP', filler=0)
# plotavsd(self,prof='belowcloud',num=0,inst='PCASP', filler=0)
# plotallavsd(self, prof='belowcloud', num=0, interact=None, exc='2d')
# belowprof(self)
# hrzplot(self,interact=None, Rtime=None)
# precipincloudTF(self,abovesize=100)
# totalprecip(self,abovesize=100,scan=0, filler=0)
# totalvolprecip(self,abovesize=100,scan=0, filler=0)
# precipconc(self,abovesize=100,scan=0,minlwc=0.01, filler=0)
# avCDNC(self,abovesize=nan,uppersize=nan,inst="FSSP96",scan=0)
# vpinfo(self,param,base='bg')
# effrad(self,inst, bindist='lin')
# path3d(self)
# filler(self,ts,timestep=10.0,inst='2dc')
# turbhrz(self,MaxCA=2,MinDist=20000,PtNo=25,Pts=5,Pts2=50)
# writeouthdf5(fName=None,fdir=None)

# properties:
# lwp
# thickness
# tempinfo
# presinfo
# precipblwTF
# precipblwvpTF
# angles

from pylab import *
import numpy as np
import copy as copy
import datetime as dt
import matplotlib.pyplot as plt
import math
import xlrd
import os
import scipy.stats.stats as st
from CloudDataSorter import runstats
from CloudDataSorter import dNdlogDp2N
from scipy import interpolate
np.seterr(divide='ignore',invalid='ignore')        # ignoring warning of divisions by zero, or invalid values (such as NaNs)


### SAMAC version ###
__version__="0.9.2" #
### SAMAC version ###

class Cloud(dict):
    numclouds=0
    def __init__(self,data,dttl,dunit,times,sd):
        # inputs
        self.data=copy.copy(data)
        self.dttl=copy.copy(dttl)
        self.dunit=copy.copy(dunit)
        self.times=copy.deepcopy(times)
        self.sd=copy.deepcopy(sd)
        self.extradata=list(); self.extrattl=list(); self.extraunit=list()

        
        ###### others ######
        Cloud.numclouds+=1      # number of cloud objects. NB: this only works if the python session stays open.
        self.extradatanum=0     # initializing the number of data that will be added afterwards (not during initilization)
        self.descedit=0         # initializing the number of time the cloud description was edited
        self.props=dict()       # diverse physical and chemical properties of the cloud
        
        # orientation of the data (in rows or columns)
        if len(self.dttl)>=2:
            if np.shape(self.data)[0]==len(self.dttl):
                self.orient='row'
            elif np.shape(self.data)[1]==len(self.dttl):
                self.orient='col'
            else: 
                print("Odd shaped data. Orientation has not been defined.")
                self.orient='unknown'
        else: 
            print("[Cloud initialization] The data should have at least 2 variables, one of which should be time. Orientation has not been defined.")
            self.orient='unknown'
        
            
    def __del__(self):
        Cloud.numclouds-=1


###########
# methods #
###########


    #########################################################################        
    #############################   overview   ##############################
    #########################################################################
    def overview(self,interact=0,Rtime=None):
        """ This method plots the cloud object's Altitude, Latitude, Longitude, Particle concentration, Droplet concentration, Liquid Water Content and Precipitation drop concentration as a function of time for the whole cloud period. Extradata module is fully handled.
            Use CloudObj.overview() for a non-interactive plot;
            Use CloudObj.overview(1) or CloudObj.overview(interact=1) for an interactive plot where you can magnify, induce offsets, and zoom.
            Use CloudObj.overview(Rtime=1) to display Hour:Minutes instead of fraction of day on the x-axis."""
        from plotcloud import plotcloud2
        
        # time
        pos=[i for i,x in enumerate(self.dttl) if x == 'time']
        if len(pos)==1: tim=self.data[pos[0]];
        else: print("[overview] No time (or multiple) was found. Forget about this plot! It won't happen unless you have (a) time")
        # altitude
        pos=[i for i,x in enumerate(self.dttl) if x == 'altitude']
        if len(pos)==1: alt=self.data[pos[0]]
        else: print("[overview] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude"); alt=NaN;
        # lat
        pos=[i for i,x in enumerate(self.dttl) if x == 'latitude']
        if len(pos)==1: lat=self.data[pos[0]]
        else: print("[overview] Latitude not found"); lat=NaN;
        # lon
        pos=[i for i,x in enumerate(self.dttl) if x == 'longitude']
        if len(pos)==1: lon=self.data[pos[0]]
        else: print("[overview] Longitude not found"); lon=NaN
        # LWC
        pos=[i for i,x in enumerate(self.dttl) if x.upper() == 'LWC']
        if len(pos)>1: print("[overview] Multiple LWC found."); lwc=np.ones(len(tim))*NaN; lwctime=tim; lwcunits=''     # LWC found in the basic data
        elif len(pos)==1: lwc=self.data[pos[0]]; lwctime=tim; lwcunits=self.dunit[pos[0]]
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                lwcunits=self.extraunit[posx[0][0]][posx[0][1]]
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[overview] No LWC (or multiple) was/were found in the extra data."); lwc=np.ones(len(tim))*NaN; lwctime=tim; lwcunits=''
        ### Size Distributions ###
        TakenSDs=list()     # names list for size distributions, passed to set labels in plot.
        # Aerosol particles 
        As=[x for x in self.sd if 'A' in x["sdtype"].upper()] 
        if len(As)==0:
            for i in range(len(self.dttl)): print("%d- %s (%s)" % (i, self.dttl[i], self.dunit[i]))
            Q=raw_input("No aerosol sdtype found. Is there a column in the basic data corresponding to the particle number concentration? [n]/#  ")
            try: 
                Q=int(Q); ptc=self.data[Q]; ptctime=tim; TakenSDs.append([self.dttl[Q],self.dunit[Q]])
            except: 
                for i,t in enumerate(self.extrattl):
                    print("%d- %s" %(i,t)); 
                Q=raw_input("Is the total aerosol concentration data in one of these lists? [n]/#  ")
                try:
                    Q=int(Q); 
                    for j,tt in enumerate(t):
                        print("%d- %s" %(j,tt)); 
                    Qj=raw_input("Choose the position of the aerosol number concentration data. [n]/#  ")
                    try: 
                        Qj=int(Qj); 
                        ptc=self.extradata[Q][Qj]; TakenSDs.append([self.extrattl[Q][Qj],self.extraunit[Q][Qj]]);
                        k=[k for k,x in enumerate(self.extrattl[Q]) if x.lower() == 'time'][0]; ptctime=self.extradata[Q][k]; 
                    except: ptc=np.ones(len(tim))*NaN; ptctime=tim; TakenSDs.append(['','']);
                except: ptc=np.ones(len(tim))*NaN; ptctime=tim; TakenSDs.append(['','']);
        elif len(As)==1: ptc=As[0]["total"]; ptctime=As[0]["time"]; TakenSDs.append([As[0]["Distname"],As[0]["units"]])
        elif len(As)>1: 
            As1=[x for x in self.sd if 'A1' in x["sdtype"].upper()]  # looking for primary distribution
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
        Cs=[x for x in self.sd if 'C' in x["sdtype"].upper()] 
        if len(Cs)==0:
            for i in range(len(self.dttl)): print("%d- %s (%s)" % (i, self.dttl[i], self.dunit[i]))
            Q=raw_input("No cloud droplet sdtype found. Is there a column in the basic data corresponding to the particle number concentration? [n]/#  ")
            try: 
                Q=int(Q); cdnc=self.data[Q]; cdnctime=tim;  TakenSDs.append([self.dttl[Q],self.dunit[Q]])
            except: 
                for i,t in enumerate(self.extrattl):
                    print("%d- %s" %(i,t)); 
                Q=raw_input("Is the total cloud droplet number concentration data in one of these lists? [n]/#  ")
                try:
                    Q=int(Q); 
                    for j,tt in enumerate(t):
                        print("%d- %s" %(j,tt)); 
                    Qj=raw_input("Choose the position of the cloud droplet number concentration data. [n]/#  ")
                    try: 
                        Qj=int(Qj); 
                        cdnc=self.extradata[Q][Qj]; TakenSDs.append([self.extrattl[Q][Qj],self.extraunit[Q][Qj]]);
                        k=[k for k,x in enumerate(self.extrattl[Q]) if x.lower() == 'time'][0]; cdnctime=self.extradata[Q][k]; 
                    except: cdnc=np.ones(len(tim))*NaN; cdnctime=tim; TakenSDs.append(['','']);
                except: cdnc=np.ones(len(tim))*NaN; cdnctime=tim; TakenSDs.append(['','']);
        elif len(Cs)==1: cdnc=Cs[0]["total"]; cdnctime=Cs[0]["time"]; TakenSDs.append([Cs[0]["Distname"],Cs[0]["units"]])
        elif len(Cs)>1: 
            Cs1=[x for x in self.sd if 'C1' in x["sdtype"].upper()]  # looking for primary distribution
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
        p=[x for x in self.sd if 'P' in x["sdtype"].upper()]
        if len(p)==0:
            # >100 microns
            for i in range(len(self.dttl)): print("%d- %s (%s)" % (i, self.dttl[i], self.dunit[i]))
            Q1=raw_input("No precipitation sdtype found. Is there a column in the basic data corresponding to the precipitation number concentration >100 microns? [n]/#  ")
            try: 
                Q1=int(Q1); pr["100"]=self.data[Q1]; pr["time100"]=tim; TakenSDs.append([self.dttl[Q1],self.dunit[Q1]]);
            except: 
                for i,t in enumerate(self.extrattl):
                    print("%d- %s" %(i,t)); 
                if len(self.extradata)>0: Q1=raw_input("Is the >100 data in one of these lists? [n]/#  ")
                try:
                    Q1=int(Q1); 
                    for j,tt in enumerate(self.extrattl[Q1]):
                        print("%d- %s" %(j,tt)); 
                    Q1j=raw_input("Choose the position of the >100 micron data. [n]/#  ")
                    try: 
                        Q1j=int(Q1j); 
                        pr["100"]=self.extradata[Q1][Q1j]; TakenSDs.append([self.extrattl[Q1][Q1j],self.extraunit[Q1][Q1j]]);
                        k=[k for k,x in enumerate(self.extrattl[Q1]) if x.lower() == 'time'][0]; pr["time100"]=self.extradata[Q1][k]; 
                    except: pr["100"]=np.ones(len(tim))*NaN; pr["time100"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
                except: pr["100"]=np.ones(len(tim))*NaN; pr["time100"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
            # >200 microns
            for i in range(len(self.dttl)): print("%d- %s (%s)" % (i, self.dttl[i], self.dunit[i]))
            Q2=raw_input("No precipitation sdtype found. Is there a column in the basic data corresponding to the precipitation number concentration >200 microns? [n]/#  ")
            try: 
                Q2=int(Q2); pr["200"]=self.data[Q2]; pr["time200"]=tim; TakenSDs.append([self.dttl[Q2],self.dunit[Q2]]);
            except: 
                for i,t in enumerate(self.extrattl):
                    print("%d- %s" %(i,t)); 
                if len(self.extradata)>0: Q1=raw_input("Is the >200 data in one of these lists? [n]/#  ")
                try:
                    Q1=int(Q1); 
                    for j,tt in enumerate(self.extrattl[Q1]):
                        print("%d- %s" %(j,tt)); 
                    Q1j=raw_input("Choose the position of the >200 micron data. [n]/#  ")
                    try: 
                        Q1j=int(Q1j); 
                        pr["200"]=self.extradata[Q1][Q1j]; TakenSDs.append([self.extrattl[Q1][Q1j],self.extraunit[Q1][Q1j]]);
                        k=[k for k,x in enumerate(self.extrattl[Q1]) if x.lower() == 'time'][0]; pr["time200"]=self.extradata[Q1][k]; 
                    except: pr["200"]=np.ones(len(tim))*NaN; pr["time200"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
                except: pr["200"]=np.ones(len(tim))*NaN; pr["time200"]=tim; TakenSDs.append(['','']); TakenSDs.append(['','']);
        elif len(p)==1: 
            pr["100"]=dNdlogDp2N(p[0],100,nan); pr["200"]=dNdlogDp2N(p[0],200,nan); pr["time100"]=p[0]["time"]; pr["time200"]=p[0]["time"];  
            TakenSDs.append([p[0]["Distname"],p[0]["units"]]); TakenSDs.append([p[0]["Distname"],p[0]["units"]]);
        elif len(p)>1: 
            p1=[x for x in self.sd if 'P1' in x["sdtype"].upper()]    # looking for primary distribution
            if len(p1)==1: 
                pr["100"]=dNdlogDp2N(p1[0],100,nan); pr["200"]=dNdlogDp2N(p1[0],200,nan); pr["time100"]=p1[0]["time"]; pr["time200"]=p1[0]["time"]; 
                TakenSDs.append([p1[0]["Distname"],p1[0]["units"]]); TakenSDs.append([p1[0]["Distname"],p1[0]["units"]]);
            else: 
                print("Precipitation distributions:")
                for i,x in enumerate(p): print("%d- %s (sdtype=%s)" % (i,x["Distname"],x["sdtype"]))
                Q=raw_input("Which is your primary Precipitation distribution (consider editing sdtype if this is permanent)? [n]/#  ")
                try: 
                    Q=int(Q); pr["100"]=dNdlogDp2N(p[Q],100,nan); pr["200"]=dNdlogDp2N(p[Q],200,nan); pr["time100"]=p[Q]["time"]; pr["time200"]=p[Q]["time"]; 
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
            plotcloud2(tim,alt,lat,lon,lwc,lwctime,ptc,ptctime,cdnc,cdnctime,pr,TakenSDs,Zs,Adj,Offsts,figi,1,proftimes=self.times,Rtime=Rtime)
        else: plotcloud2(tim,alt,lat,lon,lwc,lwctime,ptc,ptctime,cdnc,cdnctime,pr,TakenSDs,Zs,Adj,Offsts,figi,0,proftimes=self.times,Rtime=Rtime)        
        

    #########################################################################        
    #############################   mapcloud   ##############################
    #########################################################################
    def mapcloud(self,interact=1):
        """ This method plots the trajectory of the aircraft during the whole cloud period over a map of the region. Time is colour-coded and the thickness of the line reflects the LWC.
            Use CloudObj.mapcloud() or CloudObj.mapcloud(1) to continue using ipython (interactive mode);
            Use CloudObj.mapcloud(0) if you would rather close the figure before continuing to use ipython (the plot is only plotted after the whole code is read).
            Extradata module is handled for LWC.  """
        from plotcloud import maplwc
        
        # time
        pos=[i for i,x in enumerate(self.dttl) if x == 'time']
        if len(pos)==1: tim=self.data[pos[0]]
        elif len(pos)>1: print("[mapcloud] Multiple time entries were found. Go tidy up your data.")
        elif len(pos)==0: print("[mapcloud] No time was found in the basic data. Go tidy up your data.") 
        else: print("[mapcloud] Unknown Error loading time.")
       # LWC
        pos=[i for i,x in enumerate(self.dttl) if x == 'LWC']
        if len(pos)>0:
            try: 
                lwc=self.data[pos]
            except: 
                print("[mapcloud] Multiple LWC were found in the basic data. Go tidy up your data.")
                lwc=NaN
        else:
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.upper() == 'LWC']    # check all titles matching with lwc
            if len(posx)==1: 
                lwc1=self.extradata[posx[0][0]][posx[0][1]]    # loading the lwc data
                j=[j for j,x in enumerate(self.extrattl[i]) if x == 'time'][0]
                t1=self.extradata[posx[0][0]][j]     # associated time stamp
                # interpolation on the general time stamp
                lwc=nan*ones([len(tim)])    # preallocating
                ix=((tim>=t1[0])*(tim<=t1[-1]))    # index of interpolating times within the limit of the original time stamp
                f=interpolate.interp1d(t1,lwc1,'linear')    
                lwc[ix]=f(tim[ix])
                lwc=np.ma.masked_where(isnan(lwc), lwc)
            else: 
                print("[mapcloud] No LWC (or multiple) was found in the basic data, and not in the extra data either.")
                lwc=NaN
        # lat
        pos=[i for i,x in enumerate(self.dttl) if x == 'latitude']
        try: 
            lat=self.data[pos]
        except: 
            print("[mapcloud] Latitude not found in basic data")
            lat=NaN
        # lon
        pos=[i for i,x in enumerate(self.dttl) if x == 'longitude']
        try: 
            lon=self.data[pos]
        except: 
            print("[mapcloud] Longitude not found in basic data")
            lon=NaN
        
        if interact==1: maplwc(tim,lwc,lat,lon,1)     # plotting program
        else: maplwc(tim,lwc,lat,lon,0)
      

        
    #########################################################################        
    ############################   add2basic   ##############################
    #########################################################################
    def add2basic(self,newdata,newttl,newunit):
        """ This method is used to add data to after the cloud's initialization.
            If the new data has the same time vector as the original basic aircraft data, it will be appended to CloudObj.data, dttl, dunit. 
            If the new data has a different time vector as the original basic aircraft data, it will be added in the extradata module.
            newdata: 2d numpy array (same format as CloudObj.data) containing the data to be added. If the data goes to extradata, it is important to add the time stamp with it.
            newttl: list of titles corresponding to newdata (same format as CloudObj.dttl). If the data goes to the BASIC data, all data with the title "time" will be removed to avoid conflict with the main time data in the basic data. The time stamp is considered to be the same and there is no point having two copies. If you wish, however, to keep the extra time stamp in the basic data, the title can be named anything else than 'time' (like 'time1') and will be kept. 
            If the new data goes to EXTRADATA the time stamp should be passed in newdata and the title must be exactly 'time'.
            newunit: list of units corresponding to newdata (same format as CloudObj.dunit)."""
        # copying
        newdata=copy.copy(newdata); newttl=copy.deepcopy(newttl); newunit=copy.deepcopy(newunit)
        newdata=np.array(newdata)
        if len(np.shape(newdata))==1: newdata.reshape(1,-1) # if 1d array, convert to 2d array
        elif len(np.shape(newdata))==2: pass
        else: print("[add2basic] newdata must be a 2d array.")
        if (len(newttl)==len(newunit) and len(newttl)==np.shape(newdata)[0]) or (len(newttl)==len(newunit) and len(newttl)==np.shape(newdata)[1]): pass
        else: 
            print("[add2basic] newttl and newunit must be the same length as newdata. There must be a title and unit for each new data."); 
            return "[add2basic] Error: inputs non conform"
        
        # newdata includes the time (if different from the basic data) and the new data (in rows or columns)
        # newttl are the titles of the new data (same format as dttl)
        print("[add2basic] shape newdata: %s, shape basic data: %s" % (np.shape(newdata), np.shape(self.data)))
        
        if self.orient=='row':
            tlen=np.shape(self.data)[1]     #time stamp length of basic data
            # adapting newdata to row orientation
            if np.shape(newdata)[0]==len(newttl): pass
            elif np.shape(newdata)[1]==len(newttl): newdata=newdata.transpose()
            dlen=np.shape(newdata)[0]     #data length (number of different data)
            if np.shape(newdata)[1]==tlen:    # same time stamp length
                print("[add2basic] The new data has the same length as the basic data time stamp.")
                X=raw_input("Do you agree that the new data shares the basic data's time stamp, and agree to have the new data appended to the basic data ( and its own time stamp deleted)? Y/[N]  ")
                if X.lower()=='y':
                    ix=list(arange(0,dlen))
                    post=[i for i,x in enumerate(newttl) if x=='time']
                    for i in range(len(post)):
                        j=post.pop(-1);
                        garbage=newttl.pop(j); garbage=newunit.pop(j);      # removing title and units related to time
                        garbage=ix.pop(j)
                    self.data=np.vstack((self.data,newdata[ix,:]));
                    for count,title in enumerate(newttl):
                        self.dttl.append(title)
                        self.dunit.append(newunit[count])
                    AP='B'      # appended to basic
                else: AP='X'    # otherwise needs to be appended to extradata
            else: AP='X'
            if AP=='X':
                self.extradatanum+=1
                self.extradata.append(newdata)
                self.extrattl.append(newttl)
                self.extraunit.append(newunit)
        elif self.orient=='col':    
            tlen=np.shape(self.data)[0]     #time stamp length
            # adapting newdata to column orientation
            if np.shape(newdata)[1]==len(newttl): pass
            elif np.shape(newdata)[0]==len(newttl): newdata=newdata.transpose()
            dlen=np.shape(newdata)[1]     #data length (number of different data)
            if np.shape(newdata)[0]==tlen:    # same time stamp length
                print("[add2basic] The new data has the same length as the basic data time stamp.")
                X=raw_input("Do you agree that the new data shares the basic data's time stamp, and agree to have the new data appended to the basic data ( and its own time stamp deleted)? Y/[N]  ")
                if X.lower()=='y':
                    ix=list(arange(0,dlen))
                    post=[i for i,x in enumerate(newttl) if x=='time']
                    for i in range(len(post)):
                        j=post.pop(-1);
                        garbage=newttl.pop(j); garbage=newunit.pop(j);      # removing title and units related to time
                        garbage=ix.pop(j)
                    self.data=np.vstack((self.data,newdata[:,ix]));
                    for count,title in enumerate(newttl):
                        self.dttl.append(title)
                        self.dunit.append(newunit[count])
                    AP='B'      # appended to basic
                else: AP='X'    # otherwise needs to be appended to extradata
            else: AP='X'
            if AP=='X':
                self.extradatanum+=1
                self.extradata.append(newdata.transpose())
                self.extrattl.append(newttl)
                self.extraunit.append(newunit)        
        elif self.orient=='unknown': print("[add2basic] The orientation is unknown, nothing was done.")
        else: print("[add2basic] CloudObj.orient is not 'row','col' or 'unknown'. This is strange and unusual - check __init__.")
        if AP=='B': print("[add2basic] Data was appended to the basic data at CloudObj.data.")
        if AP=='X': print("[add2basic] Data was added to the extradata at CloudObj.extradata.")



    #########################################################################        
    ##############################  describe  ###############################
    #########################################################################
    def describe(self):
        """ This method is used to enter the description of a cloud if it is used for the first time, otherwise to change previous entries. Entries can also be modified by simply accessing CloudObj.desc (python dictionary)"""
        # call this method to enter the description of the cloud.
        if self.descedit==0:        # there is no data in the description yet
            self.descedit+=1
            self.desc=dict()
            self.desc["campaign"]=raw_input("Enter the name of the campaign:  ")
            self.desc["date"]=(dt.datetime(1,1,1)+dt.timedelta(days=self.times['cloud'][0][0]-1)).strftime('%Y%m%d')
            self.desc["time"]=(dt.datetime(1,1,1)+dt.timedelta(days=self.times['cloud'][0][0])).strftime('%H:%M')
            self.desc["cloudtype"]=raw_input("What kind of cloud was it? (ex. Cumulus)  ")
            self.desc["environment"]=raw_input("What kind of climate is this cloud in? (ex. Maritime/Continental/Arctic)  ")
            self.desc["flight#"]=raw_input("Flight number (enter NA if not known).  ")
            self.desc["humanplace"]=raw_input("Name the place (meaningful to humans) where this cloud occured.\nIf you would like to see a map, type 'map':  ")
            if self.desc["humanplace"].lower()=='map':
                self.mapcloud(1)
                self.desc["humanplace"]=raw_input("Name the place (meaningful to humans) where this cloud was sampled.  ")
                plt.close(1002)
            self.desc["warmcold"]=raw_input("Is this a warm, cold, or mixed cloud?  ")
            self.desc["landwater"]=raw_input("Is this cloud over land, water or both?  ")
            self.desc["layered"]="no"       # MMM Not layered by default. If the cloud is layered the user should copy the cloud and write 1 for the lowest layer, 2 for the one above, and so on.
            
            
        else:                       # there is already some data in the description
            for k in self.desc.keys():
                print("%s: %s" % (k, self.desc[k]))
            RR=raw_input("What would you like to change? key name / all (erases all entry and asks again) / add (add new entry):  ")
            if RR=='all':
                self.descedit=0
                self.describe()
            if RR=='add':
                newkey=raw_input("What is the name of the new entry's key (will be accessed by c.desc['keyname'])?  ")
                self.desc[newkey]=raw_input("what is your entry for c.desc['%s']?  " % newkey)
            else:
                if RR=='humanplace': 
                    self.desc["humanplace"]=raw_input("Name the place (meaningful to humans) where this cloud occured.\nIf you would like to see a map, type 'map':  ")
                    if self.desc["humanplace"].lower()=='map':
                        self.mapcloud(1)
                        self.desc["humanplace"]=raw_input("Name the place (meaningful to humans) where this cloud occured.  ")
                        plt.close(1002)
                elif self.desc.has_key(RR):
                    print("Old entry is %s" % self.desc[RR])
                    self.desc[RR]=raw_input("New entry for %s:  " % RR)
                else: print("Typo? Could not find the mentionned key.")



    #########################################################################        
    ###############################  plotsd  ################################
    #########################################################################
    def plotsd(self,Rtime=None):
        """ This method is used to plot the size distrubution for a selection of instruments, and scans, as a function of time.
        The instrument and scans are chosen by the user once the method has been called (interactively).
        Call CloudObj.plotsd(Rtime=1) to get the time in Hour:Minute:Second format. """
        #MMM Make those choices in the method's call instead of prompt based?
#        plt.ion()         # turn interactive figure mode on
        for i in range(len(self.sd)):
            print("%d - %s" % (i,self.sd[i]["Distname"]))
        IQ=raw_input("Which instruments do you wish to plot? [all] or # separated by commas:  ")
        try: IQ=map(int,IQ.split(','))
        except: IQ=range(len(self.sd))
        for i,x in enumerate(self.times.keys()):
            print("%d - %s (%d scans)" % (i,x,len(self.times[x])))
        SQ=raw_input("Which type of scan do you wish to plot? [all] or # separated by commas:  ")
        print("---------------------------")
        try: SQ=map(int,SQ.split(','))
        except: SQ=range(len(self.times.keys()))
        k=0; plt.close('all');
        for i in IQ:
            for j in SQ:
                scanq=self.times.keys()[j]
                for h in range(np.shape(self.times[scanq])[0]):      # for each scan of the same type (e.g. each of 2 vertical scans)
                    print("[plotsd] Now drawing: %s - %s no. %d" % (self.sd[i]["Distname"],scanq, h))
                    if type(self.sd[i]["data"])==float: print("    [plotsd] Error due to the distribution being a NaN.")
                    elif np.shape(self.sd[i]["data"])[1]==0 or np.shape(self.sd[i]["data"])[0]==0: print("    [plotsd] Error due to the distribution being empty.")
                    else:
                        chan=self.sd[i]["bins"]
                        Zdata=copy.copy(self.sd[i]["data"]); Zdata[Zdata<=1e-6]=1e-6;
                        if np.shape(np.nonzero((self.sd[i]["time"]>=self.times[scanq][h][0])*(self.sd[i]["time"]<=self.times[scanq][h][1]))[0])[0]==0:
                            print("    [plotsd] No data available for this period.")
                        elif np.shape(np.nonzero((self.sd[i]["time"]>=self.times[scanq][h][0])*(self.sd[i]["time"]<=self.times[scanq][h][1]))[0])[0]<=1:
                            print("    [plotsd] Only one point available for this period.")
                        else:
                            k+=1
                            fig=plt.figure(k)
                            Zdata=Zdata[:,np.nonzero((self.sd[i]["time"]>=self.times[scanq][h][0])*(self.sd[i]["time"]<=self.times[scanq][h][1]))[0]]
                            ptime=self.sd[i]["time"]; 
                            ptime=ptime[np.nonzero((self.sd[i]["time"]>=self.times[scanq][h][0])*(self.sd[i]["time"]<=self.times[scanq][h][1]))[0]];
                            ptime=ptime-np.floor(ptime[1]); 
                            plt.pcolor(ptime,chan,np.log10(Zdata))
                            Zmax=math.ceil(np.max(np.log10(Zdata)))
                            if Zmax>0:
                                plt.clim([0,Zmax])
                                cbar=plt.colorbar(orientation='horizontal',ticks=[0,Zmax/4.,Zmax/2.,3*Zmax/4.,Zmax],format=plt.FormatStrFormatter('$10^{%2.2f}$'))
                            else:
                                plt.clim([-6,0])
                                cbar=plt.colorbar(orientation='horizontal',ticks=[-6,-4.5,-3.0,-1.5,0],format=plt.FormatStrFormatter('$10^{%2.2f}$'))
                            if '2d' in self.sd[i]["Distname"].lower():
                                cbar.set_label('Concentration ([dN/dlogDp]Litre$^{-1}$)')
                            else: cbar.set_label('Concentration ([dN/dlogDp]cm$^{-3}$)')
                            ax=plt.gca()
                            ax.set_yscale('log')
                            plt.ylabel('Diameter/channel')
                            plt.xlabel('Time (days)')
                            plt.title("%s - %s no. %d" % (self.sd[i]["Distname"],scanq, h))
                            if scanq=='verticloud':
                                Calt=[g for g,x in enumerate(self.dttl) if x=='altitude']
                                ct=[g for g,x in enumerate(self.dttl) if x=='time']
                                alt=self.data[Calt].reshape(-1,)
                                if len(ptime)==len(alt): pass
                                else: 
                                    ttmp=self.data[ct[0]]; 
                                    ttmp=ttmp-np.floor(ttmp[1])
                                    f=interpolate.interp1d(ttmp,np.ma.filled(alt,nan),kind='linear')   
                                    ptime.reshape(max(np.shape(ptime)),)
                                    ptime=ptime[(ptime>=np.min(ttmp))*(ptime<=np.max(ttmp))]
                                    alt=f(ptime)
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
        


    #########################################################################        
    ###########################  addsizedistxl  #############################
    #########################################################################
    def addsizedistxl(self):
        """ This method is used to add a size distribution to the cloud object from an excel (.xls) file (assumed in dN/dlogDp format, conversions can be done later using codes in CloudDataSorter). 
            The added size distribution will be found appended to CloudObj.sd.
            The source file as well as the source sheet will be saved along with the distribution. """
        import Tkinter, tkFileDialog
        dirnames=[]; pathnames=[]; fnames=[];
        print("Please choose the directory:")
        pathdef=""
        while len(pathdef)==0:
            root = Tkinter.Tk()
            pathdef = tkFileDialog.askdirectory(parent=root,initialdir="/home/shared/",title='Please select a directory')
            root.withdraw()
            if len(pathdef) > 0:
                print("You chose %s" % pathdef) 
                break
            else:
                pd=raw_input("You must choose a directory. Would you like to choose again?  ")
                if pd=='y':
                    pass
                else:
                    break

        kw_main=[]
        while 1:
            kw=raw_input("Is there a search keyword (if more than one separate with commas)?  ")
            if kw=='\n':
                break
            kw=kw.lower(); kw=kw.split(',')
            for i in range(0,len(kw)):
                kw_main.append(kw[i])
            for i,x in enumerate(kw_main): kw_main[i]=str(x)
            fnames=[]
            for path, dirs, files in os.walk(pathdef):
                for i,name in enumerate(files):
                    TF=1
                    for j in range(len(kw_main)):
                        TF*=kw_main[j] in name.lower()
                    if TF==1: fnames.append(("".join([path,'/',name])))
            for i,name in enumerate(fnames): print("%d- %s" % (i,name))
            M=raw_input("Enter the file number you want to load\nOR enter (r) to refine search results\nOR enter (R) to restart search\t\t---> ")
            if M=='r':
                pass
            elif M=='R':
                fnames=[]
                kw_main=[]
            else:
                Q=int(M)
                break

        try:
            #Q=int(raw_input("Which file do you want to load (enter number) or press enter?  "))
            smpsfn=fnames[Q]
        except: print("[addsizedistxl] Failed!")
        if smpsfn[-3:].lower()=='xls':
            file = xlrd.open_workbook(smpsfn) # open excel workbook
        else: print("[addsizedistxl] Not an excel file. Use addsizedist")
        sheetnames=[]
        for i,s in enumerate(file.sheets()): 
            print("%d- %s" % (i, s.name))
            sheetnames.append(s.name)
        Q='n'
        while Q=='n':
            try: Q=int(raw_input("Which sheet do you want to load (enter number) or press enter?  "))
            except: Q='n'
        sds=file.sheet_by_name(sheetnames[Q])
        smpssheetname=sheetnames[Q]
        ttl='n'; i=0; unt='n'; dat='n'; # initiating the finding title row loop
        while ttl=='n' or unt=='n' or dat=='n':
            print("line %d- %s" % (i,sds.row_values(i,0,7)))
            Q=raw_input("Is this the title line [t], the size bin line [s], both [ts] or neither [enter]. For the first line of data: [d]. [0] to start again.  ")
            if Q=='0':
                i=0
            elif Q=='ts': ttll=i; untl=i; ttl='y'; unt='y'; i+=1
            elif Q=='t': ttll=i; i+=1; ttl='y'
            elif Q=='s': untl=i; i+=1; unt='y';
            elif Q=='d': datl=i; i+=1; dat='y'
            else: i+=1
        del ttl, i, unt, dat
        
        nupo=sds.nrows-datl;         # how long is the time vector
        sdd=[]; sddttl=[]; sddunt=[];
        for c in range(sds.ncols):
            cotype=sds.cell(datl,c).ctype
            tmpco=sds.col_values(c,datl)
            if cotype==1 or cotype==0:      # is string or empty
                colemp=1       # col is empty (or full of strings) as far as it's equal to 1
                i=0
                while colemp==1 and i<sds.nrows:
                    if sds.cell(i,c).ctype!=1 and sds.cell(i,c).ctype!=0: colemp=0;     # if there is at least one non-string or non-empty, the column is no longer empty
                    i+=1
                if colemp==0:
                    for i in range(datl,sds.nrows): 
                        if sds.cell(i,c).ctype==0 or sds.cell(i,c).ctype==1: tmpco[i-datl]=NaN
            elif cotype==2 or cotype==3: 
                colemp=0              # it's not empty, it's good.
                for i in range(datl,sds.nrows): 
                    if sds.cell(i,c).ctype==0 or sds.cell(i,c).ctype==1:
                        print "empty of string removed on line", i+1
                        tmpco[i-datl]=NaN
            else: colemp=1          # if it is none of those types, we consider the column empty.
            if colemp==0:
                sdd.append(tmpco)
                sddttl.append(sds.cell(ttll,c).value)
                sddunt.append(sds.cell(untl,c).value)
            
        # size bins
        sdtmp=dict()
        for i,x in enumerate(sddttl): print ("%d- %s (%s)" % (i, x, sddunt[i]))
        Q=raw_input("Enter the first and last column of the SIZE BINS (ex: 4,24)  ")
        try: 
            [c1, c2]=Q.split(',')
            c1=int(c1); c2=int(c2)
            sdtmp["bins"]=np.array(sddunt[c1:c2+1]);
        except: sdtmp["bins"]=NaN
        # concentrations of the size distribution
        if isnan(c1): sdtmp["data"]=np.ones(nupo)*NaN;
        else: sdtmp["data"]=np.array(sdd[c1:c2+1]);
        # time for size distribution
        Q1=raw_input("Are the times for this distribution the same times previously loaded for the cloud? [y]/n   ")
        if Q1=='n':
            Q=raw_input("Which is the time column?  ")
            try: 
                Q=int(Q)
                sdtmp["time"]=np.array(sdd[Q]);
            except: sdtmp["time"]=np.ones(nupo)*NaN;
            xltime=dt.datetime(1899,12,31,0,0,0).toordinal()-1
            if any(sdtmp["time"]>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"]<dt.datetime(2050,1,1,0,0,0).toordinal()): pass
            elif any(sdtmp["time"][0]+xltime>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"][0]+xltime<dt.datetime(2050,1,1,0,0,0).toordinal()): 
                sdtmp["time"]=sdtmp["time"]+xltime
            else: print("[addsizedistxl] Strange times...")
        else:
            sdtmp["time"]=self.data[0]
            
        # format check
        logconc=raw_input("Are the concentrations in dN [1] format, dN/dDp [2] or dN/dlogDp [3]?  ")
        if logconc=="1": 
            from CloudDataSorter import dN2dNdlogDp
            sdtmp["data"]=dN2dNdlogDp(sdtmp)
        elif logconc=="2":
            from CloudDataSorter import dNdDp2dNdlogDp
            sdtmp["data"]=dNdDp2dNdlogDp(sdtmp)
        elif logconc=="3": pass
        else: print "The concentrations have been assumed to be in dN/dlogDp format"
        
        #total concentration for size distribution
        Q=raw_input("Enter the column with the total concentration (press enter if none):  ")
        try:
            Q=int(Q)
            sdtmp["total"]=np.array(sdd[Q])
        except: 
            try: 
                sd4tot=dict(); sd4tot["bins"]=sdtmp["bins"]; sd4tot["data"]=sdtmp["data"]; sd4tot["time"]=sdtmp["time"];
                sdtmp["total"]=dNdlogDp2N(sd4tot,nan,nan)
                del sd4tot
            except: sdtmp["total"]=np.ones(nupo)*NaN;
        # MMM need to add median and mean diameters
        # Distribution name
        sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
        sdtmp["units"]=raw_input("What are the units of the size distribution? (e.g. cm-3)  ")
        sdtmp["sdtype"]=raw_input("What kind of size distribution is it? (A=aerosol, C=cloud droplet, P=precipitation)  ")
        # appending distribution to Cloud object.
        sdtmp["total"]=sdtmp["total"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["data"]=sdtmp["data"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["time"]=sdtmp["time"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["sourcefile"]=smpsfn; sdtmp["sourcesheet"]=smpssheetname

        self.sd.append(sdtmp)


        try: del dirnames, pathnames, fnames, kw, smpsfn, sheetnames, Q, sds, sdtmp, nupo, c1,c2, colemp, cotype
        except: print("[addsizedistxl] Failed to clear memory")
        
    
    #########################################################################        
    ###########################   addsizedist   ############################
    #########################################################################
    def addsizedist(self):
        """ This method is used to add a size distribution to the cloud object, from file of basic text formats (e.g. txt or dat). Excel spreadsheets are not supported by this function. Use addsizedistxl instead.
            The added size distribution will be found appended to CloudObj.sd.
            The source file will be saved along with the distribution."""
        from CloudDataSorter import todatenum
        import Tkinter, tkFileDialog
        dirnames=[]; pathnames=[]; fnames=[];
        print("Please choose a directory:")
        pathdef=""
        while len(pathdef)==0:
            root = Tkinter.Tk()
            pathdef = tkFileDialog.askdirectory(parent=root,initialdir="/home/shared/",title='Please select a directory')
            root.withdraw()
            if len(pathdef) > 0:
                print("You chose %s" % pathdef) 
                break
            else:
                pd=raw_input("You must choose a directory. Would you like to choose again?  ")
                if pd=='y':
                    pass
                else:
                    break

        kw_main=[]
        while 1:
            kw=raw_input("Is there a search keyword (if more than one separate with commas)?  ")
            if kw=='\n':
                break
            kw=kw.lower(); kw=kw.split(',')
            for i in range(0,len(kw)):
                kw_main.append(kw[i])
            for i,x in enumerate(kw_main): kw_main[i]=str(x)
            fnames=[]
            for path, dirs, files in os.walk(pathdef):
                for i,name in enumerate(files):
                    TF=1
                    for j in range(len(kw_main)):
                        TF*=kw_main[j] in name.lower()
                    if TF==1: fnames.append(("".join([path,'/',name])))
            for i,name in enumerate(fnames): print("%d- %s" % (i,name))
            M=raw_input("Enter the file number you want to load\nOR enter (r) to refine search results\nOR enter (R) to restart search\t\t---> ")
            if M=='r':
                pass
            elif M=='R':
                fnames=[]
                kw_main=[]
            else:
                Q=int(M)
                break

        try:
            #Q=int(raw_input("Which file do you want to load (enter number, 2dc is recommended) or press enter?  "))
            smpsfn=fnames[Q]
        except: print("[addsizedist] Failed!")
        if smpsfn[-3:].lower()=='xls': print("[addsizedist] This is an excel file. Please use addsizedistxl")
        else: inp = open(smpsfn,'r') # open file
        sds=[]
        for line in inp.readlines():
            ltmp=line.strip("\n"); ltmp=ltmp.split('\t');
            sds.append(ltmp)
            inp.close()       # consider using the "with open(....) as inp" instead (http://stackoverflow.com/questions/4599980/python-close-file-descriptor-question last visited Dec 2013)
        try: 
            if np.shape(sds)[0]==1 or np.shape(sds)[1]==1:   # if the delimiter was not a tab
                print("%s \n %s" % (sds[0],sds[-1]))
                delima=raw_input("What is the file's delimiter?  ")
                sds=[]; inp = open(smpsfn,'r') # open file
                for line in inp.readlines():
                    ltmp=line.strip("\n"); ltmp=ltmp.split(delima);       # trying two spaces as delimiter
                    sds.append(ltmp)
                    inp.close()
        except:
            try:
                sds=[]; inp=open(smpsfn,'r') # open file
                for line in inp.readlines():
                    ltmp=line.strip("\n"); ltmp=ltmp.split();
                    sds.append(ltmp)
                    inp.close()
            except: print("[addsizedist] Failed to find a working loading case.")

        ttl='n'; i=0; unt='n'; dat='n'; # initiating the finding title row loop
        while ttl=='n' or unt=='n' or dat=='n':
            print("line %d- %s" % (i,sds[i][0:5]))
            Q=raw_input("Is this the title line [t], the size bin line [s], both [ts] or neither [enter]. For the first line of data: [d]. [0] to start again.  ")
            if Q=='0':
                i=0
            elif Q=='ts': ttll=i; untl=i; ttl='y'; unt='y'; i+=1
            elif Q=='t': ttll=i; i+=1; ttl='y'
            elif Q=='s': untl=i; i+=1; unt='y';
            elif Q=='d': datl=i; i+=1; dat='y'
            else: i+=1
        del ttl, i, unt, dat
        sddttl=sds[ttll]
        sddunt=sds[untl]
        sdd=map(float,np.array(sds[datl]))
        for i in range(datl+1,len(sds)):
            try: sdd=vstack((sdd,map(float,np.array(sds[i]))))
            except: print i,sds[i]
        nupo=np.shape(sdd)[0]
        
        # size bins
        sdtmp=dict()
        for i,x in enumerate(sddttl): print ("%d- %s (%s)" % (i, x, sddunt[i]))
        Q=raw_input("Enter the first and last column of the SIZE BINS (ex: 4,24)  ")
        try: 
            [c1, c2]=Q.split(',')
            c1=int(c1); c2=int(c2)
            sdtmp["bins"]=np.array(map(float,sddunt[c1:c2+1]));
        except: sdtmp["bins"]=NaN
        # concentrations of the size distribution
        if isnan(c1): sdtmp["data"]=np.ones(nupo,(c2+1-c1))*NaN;
        else: sdtmp["data"]=np.array(sdd[:,c1:c2+1]);
        # time for size distribution
        Q=raw_input("Which is(are) the time column(s)?  ")
        # MMM I didn't feel like wasting a lot of time making the time retrieval more robust. Might wanna do that later.
        try: 
            Q=map(int,Q.split(','))
            if len(Q)==1: Q=Q[0]
            QQ=raw_input("Is %s an ordinal format? y/[n]  " % sdd[1,Q])
            if QQ=='y':
                sdtmp["time"]=np.array(sdd[:,Q]);
                if any(sdtmp["time"]>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"]<dt.datetime(2050,1,1,0,0,0).toordinal()): pass
                elif any(sdtmp["time"][0]+xltime>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"][0]+xltime<dt.datetime(2050,1,1,0,0,0).toordinal()): 
                    sdtmp["time"]=sdtmp["time"]+xltime
                else: print("[addsizedist] Strange ordinal...")
            else: 
                T1=raw_input("What is the date? (yyyymmdd)  ")
                T1=todatenum(dt.datetime.strptime(T1, '%Y%m%d'))
                T2f=raw_input("Enter the time format of this %f (%s) [%%Y (4 digit) %%y (2 digits) %%m %%d %%H %%M %%S (enter exact format with all -, /, spaces, etc.)]:  " % (sdd[1,Q], sddttl[Q]))
                ttmp=[]
                for i in range(nupo):
                    T2=todatenum(dt.datetime.strptime(str(int(sdd[i,Q])), T2f))
                    TT=T1+T2-np.floor(T2)
                    ttmp.append(TT)
                sdtmp["time"]=np.array(ttmp)
        except: sdtmp["time"]=np.ones(nupo)*NaN;
        xltime=dt.datetime(1899,12,31,0,0,0).toordinal()
        logconc=raw_input("Are the concentrations in dN [1] format, dN/dDp [2] or dN/dlogDp [3]?  ")
        if logconc=="1": 
            from CloudDataSorter import dN2dNdlogDp
            sdtmp["data"]=dN2dNdlogDp(sdtmp)
        elif logconc=="2":
            from CloudDataSorter import dNdDp2dNdlogDp
            sdtmp["data"]=dNdDp2dNdlogDp(sdtmp)
        elif logconc=="3": pass
        else: print "The concentrations have been assumed to be in dN/dlogDp format"
        
        # total concentration for size distribution
        Q=raw_input("Enter the column with the total concentration:  ")
        try:
            Q=int(Q)
            sdtmp["total"]=np.array(sdd[Q])
        except: 
            try: 
                sd4tot=dict(); sd4tot["bins"]=sdtmp["bins"]; sd4tot["data"]=sdtmp["data"].transpose(); sd4tot["time"]=sdtmp["time"];
                sdtmp["total"]=dNdlogDp2N(sd4tot,nan,nan)
                del sd4tot
            except: sdtmp["total"]=np.ones(nupo)*NaN;
        # MMM need to add median and mean diameters
        # Distribution name
        sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
        sdtmp["units"]=raw_input("What are the units of the size distribution? (e.g. cm-3)  ")
        sdtmp["sdtype"]=raw_input("What kind of size distribution is it? (A=aerosol, C=cloud droplet, P=precipitation)  ")
        # appending distribution to Cloud object.
        sdtmp["total"]=sdtmp["total"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        if np.shape(sdtmp["data"])[0]==len(sdtmp["time"]):
            sdtmp["data"]=sdtmp["data"].transpose()
        sdtmp["data"]=sdtmp["data"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["time"]=sdtmp["time"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["sourcefile"]=smpsfn
        self.sd.append(sdtmp)
        
    
    #########################################################################        
    ##############################   vprof   ################################
    #########################################################################
    def vprof(self,interact=None,axtype1='log',axtype2='log'):
        """ This method plots the Temperature, Theta-Q, effective diameter, LWC, RH, Particle and Droplet and Precipitation drop concentrations, in the cloud as a function of Altitude. One plot for each vertical scan is produced.
            Use CloudObj.vprof(0), CloudObj.vprof() or CloudObj.vprof(interact=0)for a plot NOT in interactive mode;   
            Use CloudObj.vprof(1) or CloudObj.vprof(interact=1)for a plot in interactive mode.
            axtype1: (default='log') 'log' if the bottom axis of the left hand plot should be logarithmic, 'lin' if it should be linear.
            axtype2: (default='log') 'log' if the bottom axis of the right hand plot should be logarithmic, 'lin' if it should be linear. 
            Fully handles the extradata module."""

        
        if interact==None or interact==0: interact=0
        elif interact==1: plt.ion()
        
        # time
        post=[i for i,x in enumerate(self.dttl) if x == 'time']
        try: 
            talt=self.data[post]
            talt=talt.reshape(max(np.shape(talt)),)
        except:
            print("[vprof] No time (or multiple) was found. Forget about this plot! It won't happen unless you have a time")
        # altitude
        pos=[i for i,x in enumerate(self.dttl) if x == 'altitude']
        try: 
            Alt=self.data[pos]
        except:
            print("[vprof] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude in the main data")
        # temperature
        pos=[i for i,x in enumerate(self.dttl) if x == 'temperature']
        if len(pos)>1:  print("[vprof] Multiple temperatures found."); T=np.ones((2,))*NaN    # temperature found in the basic data
        elif len(pos)==1: T=self.data[pos[0]]; Ttime=talt; 
        else:       # looking for temperature in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
            if len(posx)==1: 
                T=self.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                Ttime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[vprof] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN;
        # theta-Q
        pos=[i for i,x in enumerate(self.dttl) if x == 'theta-q']
        if len(pos)>1: print("[vprof] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN     # Theta-Q found in the basic data
        elif len(pos)==1: 
            ThQ=self.data[pos[0]]; ThQtime=talt; 
            if 'k' in self.dunit[pos[0]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else:       # looking for Theta-Q in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'theta-q']    # check all titles matching with Theta-Q
            if len(posx)==1: 
                ThQ=self.extradata[posx[0][0]][posx[0][1]]    # loading the Theta-Q data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                ThQtime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                if 'k' in self.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
            else: print("[vprof] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN;
        # LWC
        pos=[i for i,x in enumerate(self.dttl) if x == 'LWC']
        if len(pos)>1: print("[vprof] Multiple LWC found."); lwc=np.ones((2,))*NaN     # LWC found in the basic data
        elif len(pos)==1: lwc=self.data[pos[0]]; lwctime=talt; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[vprof] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN;
        # RH
        pos=[i for i,x in enumerate(self.dttl) if x.upper() == 'RH']
        if len(pos)>1: print("[vprof] Multiple RH found."); rh=np.ones((2,))*NaN     # RH found in the basic data
        elif len(pos)==1: rh=self.data[pos[0]]; rhtime=talt; 
        else:       # looking for RH in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'rh']    # check all titles matching with RH
            if len(posx)==1: 
                rh=self.extradata[posx[0][0]][posx[0][1]]    # loading the RH data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                rhtime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[vprof] No RH (or multiple) was found in the basic data or the extra data."); rh=np.ones((2,))*NaN;

           
        ### Interpolating the altitude onto the data ###
        f=interpolate.interp1d(talt,np.ma.filled(Alt,nan),kind='linear')   # f=interp(x=time of altitude, y=alt, linear interpolation)
        
        A=list(); C=list(); D=list(); 
        for sd in self.sd:
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
                    

        for v,ts in enumerate(self.times['verticloud']):
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
                    AA.append(np.ma.vstack((newalt,sd["total"][(sd["time"]<=ts2)*(sd["time"]>=ts1)],self.effrad(inst=sd["Distname"],bindist='lin')[(sd["time"]<=ts2)*(sd["time"]>=ts1)])))
                except: AA.append(np.ones((3,1))*NaN); print("[vprof] %s was not successfully added to the plot." % sd["Distname"])
            for sd in C:        # Cloud droplets
                try: 
                    newalt=f(sd["time"][(sd["time"]<=ts2)*(sd["time"]>=ts1)])      # new alt=f(new time)
                    CC.append(np.ma.vstack((newalt,sd["total"][(sd["time"]<=ts2)*(sd["time"]>=ts1)],self.effrad(inst=sd["Distname"],bindist='lin')[(sd["time"]<=ts2)*(sd["time"]>=ts1)])))
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
                            DD100.append(np.ma.vstack((newalt,dNdlogDp2N(sd,100,nan)[ix])))
                            DD200.append(np.ma.vstack((newalt,dNdlogDp2N(sd,200,nan)[ix])))
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
            

            fig.text(0.35, 0.92, "%s @ %s - %s \n(vert. scan #%d)" % (self.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=self.data[post][0][0])).strftime('%H:%M'), self.desc["humanplace"], v),horizontalalignment='center')
            lines = phandlez
            ax3.legend(lines, [l.get_label() for l in lines],loc=(1.01,0))
            
        if interact==0:
            if np.shape(self.times['verticloud'])[0]>=1: 
                print("[vprof] The sequence may be blocked until you close all figures.")
                plt.show()
            else: print("[vprof] No vertical scan was measured in this cloud.")

    

    #########################################################################        
    #############################   wholeprof   #############################
    #########################################################################
    def wholeprof(self,interact=None,axtype1='log',axtype2='log'):
        """ This method plots the Temperature, LWC, Particle and Droplet and Precipitation drop concentrations (and more), in the cloud as a function of Altitude. All the available data, regardless of profiles are plotted.
            Use CloudObj.wholeprof(0), CloudObj.wholeprof() or CloudObj.wholeprof(interact=0)for a plot NOT in interactive mode;   
            Use CloudObj.wholeprof(1) or CloudObj.wholeprof(interact=1)for a plot in interactive mode.
            axtype1: (default='log') 'log' if the bottom axis of the left hand plot should be logarithmic, 'lin' if it should be linear.
            axtype2: (default='log') 'log' if the bottom axis of the right hand plot should be logarithmic, 'lin' if it should be linear. """
        
        if  interact==0: pass
        elif interact==None or interact==1: plt.ion()
        ts1=self.times["cloud"][0][0]; ts2=self.times["cloud"][0][1]
        
        # time
        post=[i for i,x in enumerate(self.dttl) if x == 'time']
        try: 
            talt=self.data[post]
            talt=talt.reshape(max(np.shape(talt)),)
        except:
            print("[wholeprof] No time (or multiple) was found. Forget about this plot! It won't happen unless you have (a) time")
        # altitude
        pos=[i for i,x in enumerate(self.dttl) if x == 'altitude']
        try: 
            Alt=self.data[pos]
        except:
            print("[wholeprof] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an altitude in the main data")
        # temperature
        pos=[i for i,x in enumerate(self.dttl) if x == 'temperature']
        if len(pos)>1:  print("[wholeprof] Multiple temperatures found."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN    # temperature found in the basic data
        elif len(pos)==1: T=self.data[pos[0]]; Ttime=talt; 
        else:       # looking for temperature in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
            if len(posx)==1: 
                T=self.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                Ttime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                T=T[(Ttime>=ts1)*(Ttime<=ts2)].reshape(-1,); Ttime=Ttime[(Ttime>=ts1)*(Ttime<=ts2)].reshape(-1,); 
            else: print("[wholeprof] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN
        # theta-Q
        pos=[i for i,x in enumerate(self.dttl) if x == 'theta-q']
        if len(pos)>1: print("[wholeprof] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN     # Theta-Q found in the basic data
        elif len(pos)==1: 
            ThQ=self.data[pos[0]]; ThQtime=talt; 
            if 'k' in self.dunit[pos[0]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else:       # looking for Theta-Q in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'theta-q']    # check all titles matching with Theta-Q
            if len(posx)==1: 
                ThQ=self.extradata[posx[0][0]][posx[0][1]]    # loading the Theta-Q data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                ThQtime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                ThQ=ThQ[(ThQtime>=ts1)*(ThQtime<=ts2)].reshape(-1,); ThQtime=ThQtime[(ThQtime>=ts1)*(ThQtime<=ts2)].reshape(-1,); 
                if 'k' in self.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
            else: print("[wholeprof] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN
        # LWC
        pos=[i for i,x in enumerate(self.dttl) if x == 'LWC']
        if len(pos)>1: print("[wholeprof] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN     # LWC found in the basic data
        elif len(pos)==1: lwc=self.data[pos[0]]; lwctime=talt; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                lwc=lwc[(lwctime>=ts1)*(lwctime<=ts2)].reshape(-1,); lwctime=lwctime[(lwctime>=ts1)*(lwctime<=ts2)].reshape(-1,); 
            else: print("[wholeprof] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
        # RH
        pos=[i for i,x in enumerate(self.dttl) if x.upper() == 'RH']
        if len(pos)>1: print("[wholeprof] Multiple RH found."); rh=np.ones((2,))*NaN; rhtime=np.ones((2,))*NaN     # RH found in the basic data
        elif len(pos)==1: rh=self.data[pos[0]]; rhtime=talt; 
        else:       # looking for RH in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'rh']    # check all titles matching with RH
            if len(posx)==1: 
                rh=self.extradata[posx[0][0]][posx[0][1]]    # loading the RH data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                rhtime=self.extradata[posx[0][0]][j]     # loading associated time stamp
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
        for sd in self.sd:
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
                ix=nonzero((sd["time"]>=ts1)*(sd["time"]<=ts2))[0]
                sdtime=sd["time"][ix]
                newalt=f(sdtime)
                AA.append(np.ma.vstack((newalt,sd["total"][ix])))
            except: AA.append(np.ones((2,1))*NaN); print("[wholeprof] %s was not successfully added to the plot." % sd["Distname"])
        for sd in C:        # Cloud droplets
            try: 
                ix=nonzero((sd["time"]>=ts1)*(sd["time"]<=ts2))[0]
                sdtime=sd["time"][ix]
                newalt=f(sdtime)
                CC.append(np.ma.vstack((newalt,sd["total"][ix])))
            except: CC.append(np.ones((2,1))*NaN); print("[wholeprof] %s was not successfully added to the plot." % sd["Distname"])
        for sd in D:        # Drizzle drops
            try: 
                ix=nonzero((sd["time"]>=ts1)*(sd["time"]<=ts2))[0]
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
            

        fig.text(0.35, 0.92, "%s @ %s - %s \n(Whole cloud profile)" % (self.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=self.data[post][0][0])).strftime('%H:%M'), self.desc["humanplace"]),horizontalalignment='center')
        lines = phandlez
        ax3.legend(lines, [l.get_label() for l in lines],loc=(1.01,0))
        plt.show()
                
        if interact==0:
            if np.shape(self.times['cloud'])[0]>=1: 
                print("[wholeprof] The sequence may be blocked until you close all figures.")
                plt.show()
            else: print("[wholeprof] No graph has been generated.")

            


    #########################################################################        
    #############################   defheight   #############################
    #########################################################################
    def defheight(self):
        """ Interactive method. Sets up or clears and resets the position of the cloud. You will be asked to click 4 times per vertical profile: twice at the bottom (lower and upper guesses), and twice at the top (lower and upper guesses). This data will be stored under CloudObj.props["height"] with one such array per vertical profile. 
        If you need to add a profile, we recommend you insert a space in the list, which can then be overwritten without touching the others."""
        plt.close('all')
        import matplotlib._pylab_helpers
        figalready=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        self.vprof(1)
        fignow=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        try: 
            if len(self.props["height"])>=1:
                Q=raw_input("Lower and upper estimates of heights are already defined for %d profiles.\nWould you like to overwrite them all? y/n/a(to abort)  " % (len(self.props["height"])))
            elif len(self.props["height"])==0: Q='y'
        except: Q='y'
        if Q=='y':
            self.props["height"]=list()
            for i in range(len(figalready),len(fignow)):
                print("Click on the bottom (lower and upper estimate) and top (lower and upper estimate) of the cloud on figure %d." % (50+i))
                L=fignow[i].ginput(4,timeout=0)
                L2=np.array(L)[:,1]; L2.sort()
                self.props["height"].append(np.reshape(L2,(-1,1)))      # 4 heights: bottom 2 estimates and top 2 estimates.
                plt.close(50+i)
        if Q=='n':
            try:
                N=int(raw_input("Enter the scan number (2nd line of the figure title) you wish to modify.  "))
                f=len(figalready)+N    # figure number
                print("Click on the bottom (lower and upper estimate) and top (lower and upper estimate) of the cloud on figure %d." % (50+N))
                L=fignow[f].ginput(4,timeout=0)
                L2=np.array(L)[:,1]; L2.sort()
                self.props["height"][N]=np.reshape(L2,(-1,1))
                for i in range(len(figalready),len(fignow)):
                    plt.close(50+i)
            except: print("[defheight] Failed. Was the scan number okay?")
        else: 
            for i in range(len(figalready),len(fignow)):
                plt.close(50+i)
                
    #########################################################################        
    ############################   defBGheight   ############################
    #########################################################################
    def defBGheight(self):
        """ Interactive method. Sets up or clears and resets the position of the cloud. You will be asked to click twice per vertical profile: at the bottom and at the top. This is a best guess for cloud top and base. This data will be stored under CloudObj.props["BGheight"] with one such array per vertical profile.
        If you need to add a profile, we recommend you insert a space in the list, which can then be overwritten without touching the others."""
        plt.close('all')
        import matplotlib._pylab_helpers
        figalready=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        self.vprof(1)
        fignow=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        try: 
            if len(self.props["BGheight"])>=1:
                Q=raw_input("Best guess heights are already defined for %d profiles.\nWould you like to overwrite them all? y/n/a(to abort)  " % (len(self.props["BGheight"])))
            else: Q='y'
        except: Q='y'
        if Q=='y':
            self.props["BGheight"]=list()
            for i in range(len(figalready),len(fignow)):
                print("Click on the bottom and top of the cloud (best guess) on figure %d." % (50+i))
                L=fignow[i].ginput(2,timeout=0)
                L2=np.array(L)[:,1]; L2.sort()
                self.props["BGheight"].append(np.reshape(L2,(-1,1)))      # 2 heights: bottom and top, best guess.
                plt.close(50+i)
        if Q=='n':
            try:
                N=int(raw_input("Enter the scan number (2nd line of the figure title) you wish to modify.  "))
                f=len(figalready)+N    # figure number
                print("Click at the bottom and top of the cloud (best guess) on figure %d." % (50+N))
                L=fignow[f].ginput(2,timeout=0)
                L2=np.array(L)[:,1]; L2.sort()
                self.props["BGheight"][N]=np.reshape(L2,(-1,1))
                for i in range(len(figalready),len(fignow)):
                    plt.close(50+i)
            except: print("[defBGheight] Failed. Was the scan number okay?")
        else: 
            for i in range(len(figalready),len(fignow)):
                plt.close(50+i)                
        
            
    #########################################################################        
    ##############################   lwcflag   ##############################
    #########################################################################
    def lwcflag(self,base='bg'):
        """ Interactive method. This method asks the users the rate the LWC profiles according to their quality and features. All vertical scans are evaluated the first time, after which one can edit a particular scan. This method fully handles the extradata module.
        base: base=4point: to use base height selection between 4-point height average of the lower and upper estimates for the base height (defheight).
              base=bg: to use the best guess base height (defBGheight))"""
        plt.ion()
        calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
        ct=[i for i,x in enumerate(self.dttl) if x == 'time']; ct=ct[0]
        # Finding lwc in the basic or the extra data
        clwc=[i for i,x in enumerate(self.dttl) if x.upper() == 'LWC'];
        if len(clwc)>1: print("[lwcflag] Multiple LWC found in the basic data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
        elif len(clwc)==1: lwc=self.data[clwc[0]]; lwctime=self.data[ct]; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[lwcflag] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
        # Adiabatic lwc
        ad=self.adlwc(base=base)
        ### Interpolating the altitude onto the data ###
        f=interpolate.interp1d(self.data[ct],np.ma.filled(self.data[calt],nan),kind='linear')   # f=interp(x=time of altitude, y=alt, linear interpolation)
        if len(lwctime)<2: lwc=np.ones((2,))*NaN; lwctime=np.array([ts1,ts2]);      # adapting for if lwc is too short for interpolation.
        # plotting
        try: 
            NV=len(self.times["verticloud"])
            if NV==0: crash     # this will take us out of the try and skip to the except where flags are defined on a new cloud
            # below: testing that the number of lwcflag entries is the same as the number of vertical scans
            if len(self.props["lwcflags"])==len(self.times["verticloud"]): pass
            else: 
                QQ=raw_input("The number of lwcflag entries is not equal to the number of vertical scans! Do you want to delete the dictionary and start anew? Y/[N]  ")
                if QQ.lower()=='y':
                    del self.props["lwcflags"]
                    crash       # the previous entries are deleted, we skip to the except to define flags.
                else: print("""[lwcflag] Please abort (by pressing enter at the next question) and correct your dictionary entries in CLOUD.props["lwcflags"].""")
            figure(2)
            for j,ts in enumerate(self.times["verticloud"]):
                ix=np.nonzero((lwctime>=ts[0])*(lwctime<=ts[1]))
                subplot(1,NV,j+1)
                plot(lwc[ix],f(lwctime[ix]),'b-*',ad[j][1],ad[j][0],'r-*')
                title("vertical scan %d" % (j))
                xlabel("LWC (g/m3)")
                ylabel("Altitude (m)")
                print("scan %d flags: %s" % (j, self.props["lwcflags"][j]))
            try:
                Q=int(raw_input("Which of these scans would you like to edit (press enter to abort)?  ")); 
                if Q<0 or Q>=NV:
                    print("[lwcflag] Not a choice, the operation will be aborted.")
                else:
                    plt.close(2)
                    figure(3)
                    ts=self.times["verticloud"][Q]
                    ix=np.nonzero((lwctime>=ts[0])*(lwctime<=ts[1]))
                    plot(lwc[ix],f(lwctime[ix]),'b-*',ad[Q][1],ad[Q][0],'r-*')
                    title("Vertical scan %d" % (Q))
                    xlabel("LWC (g/m3)")
                    ylabel("Altitude (m)")
                    print("1-Good quality (young/smooth) \t2-Entrainment at the top \t3-Entrainment at the bottom \t4-Entrainment in the cloud \t5-missing data \t6-Other/Bad/Evaporated \t7-lwc>adiabatic \nEnter the flags separated with commas.  ")
                    L=np.array(map(int,raw_input("What are the vertical profile's flags belonging to scan %d.  " % (Q)).split(',')))
                    self.props["lwcflags"][Q]=L
                    plt.close(3)
            except: print("[lwcflag] The sequence was aborted.")
        except:
            self.props["lwcflags"]=list()
            if len(self.times["verticloud"])>0: print("1-Good quality (young/smooth) \t2-Entrainment at the top \t3-Entrainment at the bottom \t4-Entrainment in the cloud \t5-missing data \t6-Other/Bad/Evaporated \t7-lwc>adiabatic \nEnter the flags separated with commas.  ")
            for j,ts in enumerate(self.times["verticloud"]):
                ix=np.nonzero((lwctime>=ts[0])*(lwctime<=ts[1]))
                figure(2)
                plot(lwc[ix],f(lwctime[ix]),'b-*',ad[j][1],ad[j][0],'r-*')
                title("vertical scan %d of %d" % (j,len(self.times["verticloud"])))
                xlabel("LWC (g/m3)")
                ylabel("Altitude (m)")

                L=np.array(map(int,raw_input("What are the vertical profile's flags belonging to scan %d.  " % (j)).split(',')))
                self.props["lwcflags"].append(L) 
                plt.close(2)


    #########################################################################        
    ################################  adlwc  ################################
    #########################################################################                
    def adlwc(self,base='bg'):
        """ This method returns the adiabatic LWC based on the temperature and pressure for each vertical profile. It returns the altitude and the LWC in an array for each vertical scan. The arrays are then returned in a list. This property fully handles the extradata module.
        Example: vertical scan #1's Altitude would be R[1][0] and the Adiabatic LWC R[1][1]
        base: base=4point: to use base height selection between 4-point height average of the lower and upper estimates for the base height (defheight).
              base=bg: to use the best guess base height (defBGheight))"""
        from CloudDataSorter import AdiabLWC
        try:
            Ti=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
            A=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
            t=self.data[Ti]; # time
            alt=self.data[A]    # altitude
            # Pressure
            P=[i for i,x in enumerate(self.dttl) if x.lower() == 'pressure']
            if len(P)==1: P=P[0]; Pd=self.data[P]; Pt=t; 
            elif len(P)>1: print("[adlwc] Multiple Pressure found in the basic data."); crash;
            else:       # looking for pressure in the extradata
                posx=[] 
                for i,ttl in enumerate(self.extrattl):     # for all extra datasets available
                    posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'pressure']    # check all titles matching with pressure
                if len(posx)==1: 
                    Pd=self.extradata[posx[0][0]][posx[0][1]]    # loading the pressure data
                    j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                    Pt=self.extradata[posx[0][0]][j]     # loading associated time stamp
                else: print("[adlwc] No Pressure found in the basic or the extra data."); crash
            # Temperature         
            Te=[i for i,x in enumerate(self.dttl) if x.lower() == 'temperature']
            if len(Te)==1: Te=Te[0]; Ted=self.data[Te]; Tet=t; 
            elif len(Te)>1: print("[adlwc] Multiple temperature found in the basic data."); crash;
            else:       # looking for pressure in the extradata
                posx=[] 
                for i,ttl in enumerate(self.extrattl):     # for all extra datasets available
                    posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'temperature']    # check all titles matching with temperature
                if len(posx)==1: 
                    Ted=self.extradata[posx[0][0]][posx[0][1]]    # loading the pressure data
                    j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                    Tet=self.extradata[posx[0][0]][j]     # loading associated time stamp
                else: print("[adlwc] No temperature found in the basic or the extra data."); crash
        except: 
            print("[adlwc] One of pressure, altitude, time or temperature are missing, multiple or in the wrong format.");            
            return []
        # adapting for too short data for interpolation
        if len(Pt)<2: Pd=np.ones((2,))*NaN; Pt=np.array([self.times["cloud"][0][0],self.times["cloud"][0][1]]);
        if len(Tet)<2: Ted=np.ones((2,))*NaN; Tet=np.array([self.times["cloud"][0][0],self.times["cloud"][0][1]]);
        fP=interpolate.interp1d(Pt,Pd,kind='linear')   # f=interp(x=time of pressure, y=pressure, linear interpolation)
        fT=interpolate.interp1d(Tet,Ted,kind='linear')   # f=interp(x=time of temperature, y=temperature, linear interpolation)
        # adapting time vectors to a unique interpolable time vector.
        ta1=np.max([t[0],Pt[0],Tet[0]]); ta2=np.min([t[-1],Pt[-1],Tet[-1]]);
        ta=t[nonzero((t>=ta1)*(t<=ta2))[0]];
        alt=alt[nonzero((t>=ta1)*(t<=ta2))[0]];
        # calculations
        LWCad=[]
        dlwcdz=AdiabLWC(fT(ta),fP(ta))[1]
        for i,j in enumerate(self.times["verticloud"]):
            lwci=[]
            if base=='4point':
                cb=self.props["height"][i][0]    # lower cloud base estimate
                ct=self.props["height"][i][3]    # upper cloud top estimate
                if isnan(ct): ct=np.max(alt[np.nonzero((ta>=j[0])*(ta<=j[1]))])
            else:
                cb=self.props["BGheight"][i][0]    # best guess cloud base estimate
                ct=self.props["BGheight"][i][1]    # best guess cloud top estimate
                if isnan(ct): ct=np.max(alt[np.nonzero((ta>=j[0])*(ta<=j[1]))])
            ix=np.nonzero((ta>=j[0])*(ta<=j[1])*(alt>=cb)*(alt<=ct))
            alti=alt[ix]; dalt=np.diff(alti);
            dlwci=dlwcdz[ix];   
            dlwci=dlwci[0:-1]; alti=alti[0:-1] # cutting off the last digit
            for k in range(len(dlwci)):
                ixk=np.nonzero((alti>=cb)*(alti<=alti[k]))
                lwci.append(abs(sum(dlwci[ixk]*dalt[ixk])))
            lwci=np.array(lwci)
            LWCad.append(np.array([alti,lwci]))
        return LWCad  




    #########################################################################        
    ############################   timechange   #############################
    #########################################################################
    def timechange(self,Rtime=None):
        """ This method is used to change one of the time frames in CloudObj.times. The user will be asked to choose the scan type and scan number of which times he or she wants to modify (addition and deletion are also possible. 
        Use CloudObj(Rtime=1) to get times in Hour:Minute format."""
        plt.close('all')
        self.overview(Rtime=Rtime)
        for i,x in enumerate(self.times):
            print("%d - %s" % (i,x))
        IQ=raw_input("Which type of scan do you wish to modify?  ")
        try: IQ=int(IQ); scantype=self.times.keys()[IQ]
        except: print("You made a typo. The programmer was too lazy to make this smoother for you. Just abort and start again. Sorry!")
        txist=0     # reports if there is already a time to be modified (1) or if it needs to be added
        for i,y in enumerate(self.times[scantype]):
            print("%d - %s" % (i, y)); txist=1
        if txist==1:
            SQ=raw_input("Which of those scans would you like to modify (type A to add, D to delete):  ")
            if SQ.lower()=='a': txist=0     # adding a new scan
            elif SQ.lower()=='d': 
                try:
                    Delix=raw_input("Which of the scans would you like to delete (numbers separated with commas):  ")
                    Delix=map(int,Delix.split(','))
                    txist=3     # deletion mode
                except: print("[timechange] deletion aborted (did you list the number of the scans you wish to delete separated by commas?)")
            else:
                try: SQ=int(SQ)
                except: print("You made a typo. The programmer was too lazy to make this smoother for you. Just abort and start again. Sorry!")
        elif txist==0:
            print("You are adding a scan.")
        if txist!=3:        # not deleting
            try: Zsdaynum=np.floor(self.times["cloud"][0,0])
            except: 
                tp=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
                Zsdaynum=np.floor(self.data[tp][0])
            plt.close(1001)
            self.overview(interact=1,Rtime=Rtime)
            if txist==1: print("Please click at the start and end of the period that will replace %s #%d" % (scantype,SQ))
            elif txist==0: print("Please click at the start and end of the new %s period." % (scantype))
            Ip=ginput(2,timeout=0)
            if Rtime==None: QQ=raw_input("You are going to replace the old scan definition by or add %.5f and %.5f. Are you sure? y/[n]  " % (Ip[0][0]+Zsdaynum,Ip[1][0]+Zsdaynum))
            else: 
                Rhr1=Ip[0][0]*24; Rhour1=int(floor(Rhr1)); 
                Rmt1=(Rhr1-Rhour1)*60; Rminute1=int(floor(Rmt1)); 
                Rsd1=(Rmt1-Rminute1)*60; Rsecond1=int(floor(Rsd1));
                Rhr2=Ip[1][0]*24; Rhour2=int(floor(Rhr2)); 
                Rmt2=(Rhr2-Rhour2)*60; Rminute2=int(floor(Rmt2)); 
                Rsd2=(Rmt2-Rminute2)*60; Rsecond2=int(floor(Rsd2));
                QQ=raw_input("You are going to replace the old scan definition by or add %s and %s. Are you sure? y/[n]  " % (dt.datetime(2000,1,1,Rhour1,Rminute1,Rsecond1).strftime('%H:%M:%S'),dt.datetime(2000,1,1,Rhour2,Rminute2,Rsecond2).strftime('%H:%M:%S')))
        else:
            print("You will delete the following scans:  ")
            for i in Delix:
                if Rtime==None: print("%d - %.5f to %.5f" %(i, self.times[scantype][i][0], self.times[scantype][i][1]))
                else: 
                    Rhr1=self.times[scantype][i][0]*24; Rhour1=int(floor(Rhr1)); 
                    Rmt1=(Rhr1-Rhour1)*60; Rminute1=int(floor(Rmt1)); 
                    Rsd1=(Rmt1-Rminute1)*60; Rsecond1=int(floor(Rsd1));
                    Rhr2=self.times[scantype][i][1]*24; Rhour2=int(floor(Rhr2)); 
                    Rmt2=(Rhr2-Rhour2)*60; Rminute2=int(floor(Rmt2)); 
                    Rsd2=(Rmt2-Rminute2)*60; Rsecond2=int(floor(Rsd2));
                    print("%d - %s to %s" %(dt.datetime(2000,1,1,Rhour1,Rminute1,Rsecond1).strftime('%H:%M:%S'),dt.datetime(2000,1,1,Rhour2,Rminute2,Rsecond2).strftime('%H:%M:%S')))
            QQ=raw_input("Are you sure you want to delete them? y/[n]  ")
        if QQ=='y': 
            if txist==1: self.times[scantype][SQ]=np.reshape(np.array([Ip[0][0]+Zsdaynum,Ip[1][0]+Zsdaynum]),(-1,2)); 
            elif txist==0: self.times[scantype]=vstack((self.times[scantype],np.reshape(np.array([Ip[0][0]+Zsdaynum,Ip[1][0]+Zsdaynum]),(-1,2)))); 
            elif txist==3:      # deleting
                L=list()
                for i,ts in enumerate(self.times[scantype]):
                    if i not in Delix: L.append(ts)
                self.times[scantype]=np.reshape(L,(-1,2))
                        

    #########################################################################        
    #############################   MaskData   ##############################
    #########################################################################
    def MaskData(self):
        """This method is used to mask points that are measurement errors or that should be removed from the analysis for another reason. This method is interactive."""
        try: plt.close(30)
        except: pass
        plt.ion()         # turn interactive figure mode on
        for i in range(len(self.dttl)):
            print("%d - %s (%s)" % (i,self.dttl[i],self.dunit[i]))
        for j,XL in enumerate(self.extrattl):
            for i in range(len(XL)):
                print("X%d-%d - %s (%s)" % (j,i,XL[i],self.extraunit[j][i]))
        X1=-1000       # defining under zero for later testing
        Q=raw_input("Choose the quantity you want to inspect (type number as it appears in the list).  ")
        if Q[0].upper()=='X':
            [X1,Q]=map(int,Q.upper().strip('X').split('-'))
            if X1<0: print("[MaskData] the number cannot be negative."); crash
            try:
                self.extradata[X1]=np.ma.array(self.extradata[X1])
                dv=self.extradata[X1][Q];   # data vector
                k=[k for k,x in enumerate(self.extrattl[X1]) if x.lower() == 'time'][0]
                dvt=self.extradata[X1][k]   # time vector
            except: print("[MaskData] Problem loading data."); crash
        else: 
            try: 
                Q=int(Q)
                dv=self.data[Q]
                post=[i for i,x in enumerate(self.dttl) if x == 'time'][0] # position of time
                dvt=self.data[post]
            except: print("[MaskData] Wrong entry! This is not as expected."); crash
        
        fig=plt.figure(30)
        plot(dvt,dv,'b.')
        G=raw_input("Do you want to modify this data? [y]/n  ")
        if G.lower()=='n': plt.close(30)
        else:
            G1=raw_input("Do you want to automatically remove data more than x standard deviations from the y point running average? [y]/n  ")
            if G1=='n': G2='y'
            else: 
                try:
                    numsd=float(raw_input("Enter the amount of standard deviation (recommendation: 3):  "))
                    numra=int(raw_input("Enter the amount of point for the running average (recommendation: 25):  "))
                    M=runstats(dv,numra)
                    dv=np.ma.masked_where((dv>(M[0]+M[1]*numsd)+(isnan(dv))), dv)
                    dv=np.ma.masked_where((dv<(M[0]-M[1]*numsd)+(isnan(dv))), dv)
                    plt.close(30)
                    plt.figure(30)
                    plot(dvt,dv,'b.')    
                    G2=raw_input("Do you need to mask more data?  y/[n]  ")
                except: print("[MaskData] Operation failed")
            if G2=='y':
                #self.overview()
                Fin=0;      # Fin = are there more selections needed?
                while Fin==0:
                    Sat=0;  zq=0;     # Sat= satisfied of the selection?, zq = do we need to zoom?
                    while Sat==0:
                        plt.close(30)
                        figure(30)
                        plot(dvt,dv,'b.')
                        [Zx1, Zx2, Zy1, Zy2]=axis()
                        Q1=raw_input("Do you need to zoom in (x-axis zoom only)? [y]/n  ")
                        if Q1.lower()=='n': 
                            zq=1        # no need to zoom
                            index=nonzero((dvt>=0))[0]
                        else: zq=0
                        while zq==0:        # while the zoom is not done
                            index=nonzero((dvt>=Zx1)*(dvt<=Zx2))[0]
                            plt.close(30)
                            fig=plt.figure(30)
                            plot(dvt[index],dv[index],'b.')
                            xlim(Zx1,Zx2)
                            print("Please click twice to zoom in.")
                            Zx=ginput(2,timeout=0)
                            plot([Zx[0][0],Zx[1][0]],[Zx[0][1],Zx[1][1]],'ro')
                            xlim(Zx1,Zx2)
                            Q1=raw_input("Does this work? [y]/n  ")
                            if Q1.lower()=='n': pass
                            else: 
                                Zx1=Zx[0][0]; Zx2=Zx[1][0];
                                index=nonzero((dvt>=Zx1)*(dvt<=Zx2))[0]
                                plt.close(30)
                                fig=plt.figure(30)
                                ax=fig.add_subplot(111)
                                plot(dvt[index],dv[index],'b.')
                                xlim(Zx1,Zx2)
                                Q1=raw_input("Do you want to zoom in more? y/[n]  ")
                                if Q1.lower()=='y': zq=0
                                else: zq=1
                        print("Please define a square in which points will be masked by clicking at the upper left corner and the lower right corner.")
                        ULLR=ginput(2,timeout=0)        # click upper left, lower right
                        if ULLR[0][0]>ULLR[1][0]: x2=ULLR[0][0]; x1=ULLR[1][0];
                        else: x1=ULLR[0][0]; x2=ULLR[1][0];
                        if ULLR[0][1]>ULLR[1][1]: y2=ULLR[0][1]; y1=ULLR[1][1];
                        else: y1=ULLR[0][1]; y2=ULLR[1][1];
                        plt.close(30)
                        fig=plt.figure(30)
                        plot(dvt[index],dv[index],'b.')
                        plot([x1,x1],[y1,y2],'r-')
                        plot([x2,x2],[y1,y2],'r-')
                        plot([x1,x2],[y1,y1],'r-')
                        plot([x1,x2],[y2,y2],'r-')
                        xlim(Zx1,Zx2)
                        Q2=raw_input("Are you satisfied with the selection (no to make a new selection)? y/[n]  ")
                        if Q2.lower()=='y': Sat=1
                        if Sat==1:
                            ix=nonzero((dvt>=x1)*(dvt<=x2)*(dv>=y1)*(dv<=y2))[0]
                            dv.mask[:,ix]='True'
                            Q3=raw_input("Do you want to remove more points? [y]/n  ")
                            if Q3.lower()=='n': Fin=1
                        plt.close(30)
                    else: pass      # if there is no need to modify the data after the standard deviation treatment.
            if X1<0: self.data[Q]=dv
            else: self.extradata[X1][Q]=dv
            figure(30)
            plot(dvt,dv,'b.')


    #########################################################################        
    ############################   avsizedist   #############################
    #########################################################################        
    def avsizedist(self,prof='belowcloud',scan=0,inst='PCASP',filler=0):
        """ Returning method. This method returns the average size distribution for a given instrument and scan type and number. The returned average size distribution is an array with the mean concentration of particles in the first row and the size of the particles in the second. 
        Example: R=CloudObj.avsizedist(prof='abovecloud',scan=1,inst='PCASP') would return in R the average size distribution for the above cloud #1 for the PCASP.
        The defaults are prof='belowcloud',scan=0,inst='PCASP', filler=0
        filler: filling with a given time step (useful for 2d type data) =1 ON, =0 OFF. Time step can be adjusted in the call of the method filler
        Special profile "belowvert" can also be used, this we return the average size distribution below cloud, where a vertical profile was measured."""

        cancel=0
        instn=[i for i,x in enumerate(self.sd) if x["Distname"].upper()==inst.upper()]
        if len(instn)<=0: cancel=1; raise ValueError("This instrument name has not been found (inst=%s)." % inst)
        if filler==1: sd=self.filler(self.times["cloud"][0],inst=inst)
        else: sd=copy.deepcopy(self.sd[instn[0]])
        if type(sd["time"])==float: cancel=1
        elif len(sd["time"])==0: cancel=1
        if cancel==1: avsd=[]; return avsd
        if prof=='belowvert':
            if len(sd["time"])<=1:      # interpolation cannot take place if there is not at least 2 points to interpolate from
                print("[avsizedist] Only 1 point available: special profile belowvert cannot be used.")
                avsd=[]; return avsd
            ctime=[i for i,x in enumerate(self.dttl) if x == 'time']; ctime=ctime[0]
            basictime=self.data[ctime]
            if sum(basictime==sd["time"])==len(basictime): pass      # the instrument share the timeline of the basic data from which belowvert is calculated
            else:       # the instrument does not share the basic timeline and must be interpolated
                if sum(sd["data"])>0: M=np.min(sd["data"][sd["data"]>0.0])    # minimum non-zero concentration 
                else: M=0
                f=interpolate.RectBivariateSpline(sd["bins"],sd["time"],sd["data"])
                sd["data"]=f(sd["bins"],basictime);  sd["time"]=basictime;
                sd["data"][sd["data"]<0.01*M]=0.0        # all interpolated data smaller than 1% of the minimum.
            R=self.belowprof()
            try:
                if np.size(R)>0:
                    Tix=R[scan]['areaindex']
                    concs=sd["data"][:,nonzero(Tix)[0]]
                else: print("[avsizedist] There is no data below cloud under vertical profile #%d." % scan)
            except: print("[avsizedist] The instrument number (%d) may be wrong, or no data corresponds to the profile." % instn[0])            
        else:
            try: 
                T=self.times[prof][scan]
            except: print("[avsizedist] The profile %s [%d] doesn't exist." % (prof,scan))
            try:
                concs=sd["data"][:,nonzero((sd["time"]>=T[0])*(sd["time"]<=T[1]))[0]]
            except: print("[avsizedist] The instrument number (%d) may be wrong, or no data corresponds to the profile." % instn[0])
        if "concs" in locals():
            if len(sd["time"])==1: avsd=[concs,sd["bins"]]
            else: avsd=np.array([st.nanmean(concs,axis=1), sd["bins"]])
            return avsd
        else: 
            print("[avsizedist] Failed to return the averaged size distribution.")
            avsd=[]; return avsd
        
        
    #########################################################################        
    #############################   plotavsd   ##############################
    #########################################################################            
    def plotavsd(self,prof='belowcloud',scan=0,inst='PCASP',filler=0):
        """ This method plots the average size distribution for a given instrument and scan type and number. 
        Example: CloudOjb.plotavsd(prof='abovecloud',num=1,inst='PCASP') will plot the average size distribution for the above cloud #1 for the PCASP in the current figure window.
        The defaults are prof='belowcloud',num=0,inst='PCASP'
        The Special profile 'belowvert' can also be used"""
        asd=self.avsizedist(prof=prof,scan=scan,inst=inst,filler=filler)
        if len(asd)==0: print("[plotavsd] no data available")
        else:
            if sum(asd[0])==0: plt.plot(asd[1],asd[0],'b-*')
            else: plt.loglog(asd[1],asd[0],'b-*')
            plt.xlabel("Diameter (um)")
            unitsname=[self.sd[i]["units"] for i,x in enumerate(self.sd) if x["Distname"].upper()==inst.upper()]
            plt.ylabel('concentration [dN/dlogDp](%s)' %(unitsname[0]))
            plt.show()
        

    #########################################################################        
    ############################   plotallavsd   ############################
    #########################################################################
    def plotallavsd(self,prof='belowcloud',num=0,interact=None,exc='2d'):
        """ This method plots the average size distribution of all size distribution intruments (by default, except for 2d-type measurements) for a given scan. The maximum number of instruments that can be plotted is 8.
            prof: manoeuvre type (cloud, abovecloud, belowcloud, verticloud, horicloud and belowvert (special profile, see method belowprof). Default is 'belowcloud'.
            Use CloudObj.plotallavsd(0), CloudObj.plotallavsd() or CloudObj.plotallavsd(interact=0) for a plot NOT in interactive mode;   
            Use CloudObj.plotallavsd(1) or CloudObj.plotallavsd(interact=1)for a plot in interactive mode.  
            Use CloudObj.plotallavsd(exc='inst1,inst2') to exclude instruments inst1 and inst2 (includes 2d unless specified in the list).
            exc: exceptions not to be plotted. Exceptions should be stated in a single string, with instrument names (as found in c.sd[n]["Distname"]) separated by commas (e.g. exc='2dc,PCASP,FSSP100'). Default excepts all instruments that contain '2d' in their name.
            Use CloudObj.plotallavsd(exc=None) to exclude no instrument (includes 2d).
            The default is prof='belowcloud',num=0,interact=None,exc='2d'"""
        if interact==None or interact==0: interact=0
        elif interact==1: plt.ion()
        SD=list()
        specs=['b-*','r-o','c-v','m-d','g-s','k-^','b--o','r--*']
        instlist=[x["Distname"] for i,x in enumerate(self.sd)]
        unitlist=[x["units"] for i,x in enumerate(self.sd)]
        alllist=[[x,unitlist[i]] for i,x in enumerate(instlist)]
        if exc==None: pass
        else:
            X=exc.split(',');
            for c in range(len(X)):
                X[c]=X[c].strip()
            for c in range(len(X)):
                alllist=[x for i,x in enumerate(alllist) if X[c].lower() not in x[0].lower()]
        for i in range(len(alllist)):
            SD.append(self.avsizedist(prof=prof,scan=num,inst=alllist[i][0]))
            plt.loglog(SD[i][1],SD[i][0],specs[i])
        instunit=list()
        for i in range(len(alllist)):
            instunit.append(alllist[i][0]+' ('+alllist[i][1]+')')
        plt.legend(instunit,loc=3)
        plt.xlabel('diameter (um)')
        plt.ylabel('concentration [dN/dlogDp]')
        post=[i for i,x in enumerate(self.dttl) if x == 'time'] # position of the time
        plt.title("Average size dist. for %s no. %d on %s @ %s \n (%s)" % (prof, num, self.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=self.data[post][0][0])).strftime('%H:%M'), self.desc["humanplace"]))
        if interact==0: plt.show()


    #########################################################################        
    ############################   belowprof   ##############################
    #########################################################################            
    def belowprof(self):
        """This is a returning method. This method returns a list with one dictionary per vertical scan that includes the index of what belongs to below the vertical scan, the area index (binary), and the times. This method is used to access the special profile 'belowvert'.
        Returned list format: 
        - R[0]["index"]: index of what is below cloud and within 0.0075 degrees radius (below the vertical scan)
        - R[0]["areaindex"]: index of what is in the area (within 0.0075 degrees radius)
        - R[0]["times"]: times of what is below the vertical scan."""
        BV=list()       # BV has as many list entries as the cloud has vertical scans.
        # each list element is a dictionary including times, average aerosol distribution,  ETC.?
        ctime=[i for i,x in enumerate(self.dttl) if x == 'time']; ctime=ctime[0]
        clat=[i for i,x in enumerate(self.dttl) if x == 'latitude']; clat=clat[0]
        clon=[i for i,x in enumerate(self.dttl) if x == 'longitude']; clon=clon[0]
        calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
        res=0.0075 # resolution in degrees
        for j in range(len(self.times["belowcloud"])):
                bctimetmp=(self.data[ctime]>=self.times["belowcloud"][j][0])*(self.data[ctime]<=self.times["belowcloud"][j][1])
                try: bctime=bctime+bctimetmp; 
                except: bctime=bctimetmp; 
        if len(self.times["belowcloud"])==0: bctime=False*ones(np.shape(self.data[ctime])); print("[belowprof] No below cloud scans")
        for i in range(len(self.times["verticloud"])):      # for each vertical scan
            t=dict()
            lat=self.data[clat][nonzero((self.data[ctime]>=self.times["verticloud"][i][0])*(self.data[ctime]<=self.times["verticloud"][i][1]))[0]]
            lon=self.data[clon][nonzero((self.data[ctime]>=self.times["verticloud"][i][0])*(self.data[ctime]<=self.times["verticloud"][i][1]))[0]]
            scant=self.data[ctime][nonzero((self.data[ctime]>=self.times["verticloud"][i][0])*(self.data[ctime]<=self.times["verticloud"][i][1]))[0]]
            # Defining a number of points along the trajectory of the plane, seperated by  'res' degrees
            pt=list()
            dist=0;  # units of degrees 
            k=0     # step counter
            pt.append([scant[0],lat[0],lon[0]])
            while pt[-1][0]<scant[-1] and scant[k]<scant[-1]:
                k+=1
                dist=dist+sqrt((lat[k]-lat[k-1])**2+(lon[k]-lon[k-1])**2)
                if dist>len(pt)*res: 
                    pt.append([scant[k],lat[k],lon[k]]); 
            pt.append([scant[-1],lat[-1],lon[-1]])
            # Finding the points within a radius 'res' of at least one of these points.
            arix=(np.ones(np.shape(self.data[clat]))*0)
            for m in range(len(pt)):
                arix=arix+((self.data[clat]-pt[m][1])**2+(self.data[clon]-pt[m][2])**2<=res**2)
            # now finding the index of points within this area and below the cloud base in the same vertical or below cloud:
            Ix=nonzero( arix * ( ((self.data[ctime]>=self.times["verticloud"][i][0])*(self.data[ctime]<=self.times["verticloud"][i][1])*(self.data[calt]<self.props["height"][i][0])) + bctime ) )
            
            t["index"]=Ix[0]    # index of what is below cloud and within 0.0075 degrees radius (below the vertical scan)
            t["areaindex"]=arix     # index of what is in the area (within 0.0075 degrees radius)
            t["times"]=self.data[ctime][Ix]     # times of what is below the vertical scan
            
            BV.append(t)
        return BV

    #########################################################################        
    ##############################   hrzplot   ##############################
    #########################################################################
    def hrzplot(self,interact=None,Rtime=None):
        """ This method plots the LWC & Altitude, Na (aerosols) & Nd (droplets), Temperature & Theta-Q, updraft velocity and its standard deviation (turbulence) as a function of time. One plot for each horizontal scan is produced. Fully handles the extradata module.
            Use CloudObj.hrzplot(0), CloudObj.hrzplot() or CloudObj.hrzplot(interact=0) for a plot NOT in interactive mode;   
            Use CloudObj.hrzplot(1) or CloudObj.hrzplot(interact=1)for a plot in interactive mode. 
            Use CloudObd.hrzplot(1,1) or CloudObd.hrzplot(Rtime=1) to get the time in Hour:Minute format."""
        
        if interact==None or interact==0: interact=0
        elif interact==1: plt.ion()
        
        # time
        post=[i for i,x in enumerate(self.dttl) if x == 'time']
        if len(post)==1: talt=self.data[post[0]]
        else: print("[hrzplot] No time (or multiple) was found. Forget about this plot! It won't happen unless you have a time")
        # altitude
        pos=[i for i,x in enumerate(self.dttl) if x == 'altitude']
        if len(pos)==1: Alt=self.data[pos[0]]
        else: print("[hrzplot] No altitude (or multiple) was found. Forget about this plot! It won't happen unless you have an a(l)titude")
        # temperature
        pos=[i for i,x in enumerate(self.dttl) if x == 'temperature']
        if len(pos)>1:  print("[hrzplot] Multiple temperatures found."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN    # temperature found in the basic data
        elif len(pos)==1: T=self.data[pos[0]]; Ttime=talt; 
        else:       # looking for temperature in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'temperature']    # check all titles matching with temperature
            if len(posx)==1: 
                T=self.extradata[posx[0][0]][posx[0][1]]    # loading the temperature data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                Ttime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[hrzplot] No temperature (or multiple) was found in the basic data or the extra data."); T=np.ones((2,))*NaN; Ttime=np.ones((2,))*NaN;
        # theta-Q
        pos=[i for i,x in enumerate(self.dttl) if x == 'theta-q']
        if len(pos)>1: print("[hrzplot] Multiple Theta-Q found."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN     # Theta-Q found in the basic data
        elif len(pos)==1: 
            ThQ=self.data[pos[0]]; ThQtime=talt; 
            #if 'k' in self.dunit[pos[0]].lower(): ThQ=ThQ-273.15        # converting to celsius
        else:       # looking for Theta-Q in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'theta-q']    # check all titles matching with Theta-Q
            if len(posx)==1: 
                ThQ=self.extradata[posx[0][0]][posx[0][1]]    # loading the Theta-Q data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                ThQtime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                #if 'k' in self.extraunit[posx[0][0]][posx[0][1]].lower(): ThQ=ThQ-273.15        # converting to celsius
            else: print("[hrzplot] No Theta-Q (or multiple) was found in the basic data or the extra data."); ThQ=np.ones((2,))*NaN; ThQtime=np.ones((2,))*NaN    
        # LWC
        pos=[i for i,x in enumerate(self.dttl) if x == 'LWC']
        if len(pos)>1: print("[hrzplot] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN     # LWC found in the basic data
        elif len(pos)==1: lwc=self.data[pos[0]]; lwctime=talt; 
        else:       # looking for LWC in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[hrzplot] No LWC (or multiple) was found in the basic data or the extra data."); lwc=np.ones((2,))*NaN; lwctime=np.ones((2,))*NaN
        # updraft velocity
        pos=[i for i,x in enumerate(self.dttl) if x == 'udvel']
        if len(pos)>1: print("[hrzplot] Multiple updraft velocity found."); udvel=np.ones((2,))*NaN; udveltime=np.ones((2,))*NaN     # udvel found in the basic data
        elif len(pos)==1: udvel=self.data[pos[0]]; udveltime=talt; 
        else:       # looking for updraft velocity in the extradata
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'udvel']    # check all titles matching with LWC
            if len(posx)==1: 
                udvel=self.extradata[posx[0][0]][posx[0][1]]    # loading the updraft velocity data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                udveltime=self.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[hrzplot] No LWC (or multiple) was found in the basic data or the extra data."); udvel=np.ones((2,))*NaN; udveltime=np.ones((2,))*NaN
        # turbulence calculations
        if sum(~isnan(udveltime))==0: turb=np.ones((2,))*NaN;
        else: turb=runstats(udvel,25)     # returns the running average and running standard deviation of the first argument over the number of points specified in the second argument.
        try: turb[1][udvel.mask]=NaN
        except: pass
        # Cloud droplet primary
        pos=[i for i,x in enumerate(self.sd) if 'C1' in x["sdtype"].upper()]
        if len(pos)==1: pass
        elif len(pos)==0: 
            pos=[i for i,x in enumerate(self.sd) if 'C' in x["sdtype"].upper()]
            if len(pos)==1: pass
            elif len(pos)==0: 
                print("[hrzplot] No cloud droplet size distribution found.");
                CDrop=np.ones((2,1))*NaN; CDroptime=np.ones((2,))*NaN;  CDropName=''
            elif len(pos)>1:
                for i,p in enumerate(pos):
                    print("%d - %s (%s)" %(i,self.sd[p]["Distname"],self.sd[p]["sdtype"]))
                try: 
                    chosenp=pos[int(raw_input("Which of these distributions do you want to use for cloud droplets? (number, enter if none)  "))]
                    pos=list(); pos.append(chosenp)
                except: 
                    print("[hrzplot] Failed to load the distribution.")
                    CDrop=np.ones((2,1))*NaN; CDroptime=np.ones((2,))*NaN;  CDropName=''
        if len(pos)==1:
            pos=pos[0]
            CDrop=self.sd[pos]["total"]; CDroptime=self.sd[pos]["time"]; CDropName=self.sd[pos]["Distname"]
        else: print("[hrzplot] No distribution was found.")
        # Aerosol primary
        pos=[i for i,x in enumerate(self.sd) if 'A1' in x["sdtype"].upper()]
        if len(pos)==1: pass
        elif len(pos)==0: 
            pos=[i for i,x in enumerate(self.sd) if 'A' in x["sdtype"].upper()]
            if len(pos)==1: pass
            elif len(pos)==0: 
                print("[hrzplot] No aerosol size distribution found.");
                Aero=np.ones((2,1))*NaN; Aerotime=np.ones((2,))*NaN; AeroName=''
            elif len(pos)>1:
                for i,p in enumerate(pos):
                    print("%d - %s (%s)" %(i,self.sd[p]["Distname"],self.sd[p]["sdtype"]))
                try: 
                    chosenp=pos[int(raw_input("Which of these distributions do you want to use for aerosols? (number, enter if none)  "))]
                    pos=list(); pos.append(chosenp)
                except: 
                    print("[hrzplot] Failed to load the distribution.")
                    Aero=np.ones((2,1))*NaN; Aerotime=np.ones((2,))*NaN; AeroName=''
        if len(pos)==1:
            pos=pos[0]
            Aero=self.sd[pos]["total"]; Aerotime=self.sd[pos]["time"]; AeroName=self.sd[pos]["Distname"]
        else: print("[hrzplot] No distribution was found.")
        
        
        
        for v,ts in enumerate(self.times['horicloud']):
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

            fig.text(0.35, 0.92, "%s @ %s - %s \n(horiz. scan #%d)" % (self.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=self.data[post][0][0])).strftime('%H:%M'), self.desc["humanplace"], v),horizontalalignment='center')
        if interact==0:
            if np.shape(self.times['horicloud'])[0]>=1: 
                print("The sequence may be blocked until you close all figures.")
                plt.show()
            else: print("[hrzplot] No horizontal scan was measured in this cloud.")


    #########################################################################        
    ##########################  precipincloudTF  ############################
    #########################################################################

    def precipincloudTF(self,abovesize=100):
        """ This method tells whether precipitation in cloud (horicloud + verticloud) was observed during the entire cloud period: true or false. 
        Default value abovesize=100 for precipitation >100 microns. Noise level at 1e-2."""
        # Do not use to count the total amount of rain detected as some scans may be counted more than once.
        pos=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]        # finding P1 size dist
        if len(pos)==0: pos=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()]    # finding P size dist
        if len(pos)==1: pos=int(pos[0])
        elif len(pos)>1:
            for p in pos:
                print("%d - %s" % (p, self.sd[p]["Distname"]))
            pos=int(raw_input("What is the prefered precipitation distribution? Enter the number."))
        else: pos=nan
        if type(pos)==float:
            if isnan(pos):
                precb=nan   
        elif type(self.sd[pos]["time"])==float: 
            precb=nan; 
        elif np.shape(self.sd[pos]["time"])[0]<=1: 
            precb=nan
        else:
            Pas=dNdlogDp2N(self.sd[pos],abovesize,nan)
            #Pas=np.ma.masked_less_equal(Pas,1e-2);     # nansum doesn't seem to handle masked arrays. 
            Pas[Pas<=1e-2]=0;
            t100=0;
            for b in range(np.shape(self.times['horicloud'])[0]):
                if len(Pas[(self.sd[pos]["time"]<=self.times['horicloud'][b][1])*(self.sd[pos]["time"]>=self.times['horicloud'][b][0])])==0: pass
                else: t100=t100+np.nansum(Pas[(self.sd[pos]["time"]<=self.times['horicloud'][b][1])*(self.sd[pos]["time"]>=self.times['horicloud'][b][0])])
            for b in range(np.shape(self.times['verticloud'])[0]):
                if len(Pas[(self.sd[pos]["time"]<=self.times['verticloud'][b][1])*(self.sd[pos]["time"]>=self.times['verticloud'][b][0])])==0: pass
                else: t100=t100+np.nansum(Pas[(self.sd[pos]["time"]<=self.times['verticloud'][b][1])*(self.sd[pos]["time"]>=self.times['verticloud'][b][0])])
            if t100>0: precb=True
            else: precb=False
        return precb


    #########################################################################        
    #############################  totalprecip  #############################
    #########################################################################

    def totalprecip(self,abovesize=100,scan=0,filler=0):
        """ This method returns the integrated (with respect to altitude) amount of precipitation drops (number) in a vertical profile (in #/m2 OR in # m L-1, need to edit code to change) above the size given in options. A 2d instrument is used.
        By default, this in done for the first scan (#0), other scans can be accessed using options.
        CloudObj.totalprecip(abovesize=100, scan=0) size in microns, is the default."""
        #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
        if len(self.times["verticloud"])<=scan:
            print("[totalprecip] No vertical scan for this cloud.")
            return nan
        else:
            ts=self.times["verticloud"][scan]       # times of the vertical scan
            # columns & positions
            calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
            ctim=[i for i,x in enumerate(self.dttl) if x == 'time']; ctim=ctim[0]
            pos2d=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]; 
            if len(pos2d)==0: pos2d=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()];
            try: pos2d=pos2d[0]
            except: print("[totalprecip] No precipitation size distribution found.")
            # loading the relevant distribution
            if filler==1:
                fill2d=self.filler(ts=ts)
            elif filler==0: 
                try:
                    fill2d=dict()
                    tix=np.nonzero((self.sd[pos2d]["time"]>=ts[0])*(self.sd[pos2d]["time"]<=ts[1]))[0]
                    if len(tix)==0: return 0.00
                    fill2d["data"]=self.sd[pos2d]["data"][:,tix]; fill2d["time"]=self.sd[pos2d]["time"][tix];
                    fill2d["bins"]=self.sd[pos2d]["bins"]; fill2d["total"]=self.sd[pos2d]["total"][tix];
                except: print("[totalprecip] Data not found"); return nan
            else:
                print("[totalprecip] filler must be 1 or 0.")
                return nan
            # controlling for empty distributions
            if type(fill2d["time"])==float: Int=nan     # data does not exist
            elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
            
            else:   # proceeding!
                prec=dNdlogDp2N(fill2d,abovesize,nan)
                alt=self.data[calt]
                tim=self.data[ctim]
                
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

    def totalvolprecip(self,abovesize=100,scan=0,filler=0):
        """ This method returns the integrated (with respect to altitude) amount of precipitation drops (volume) in a vertical profile (in meters OR in um3 m L-1, need to edit code to change) above the size given in options. A 2d instrument is used. Note that it is assume that the bins sizes are diameters and not radii.
        scan: default scan=0. Access to different vertical scans in a cloud.
        filler: 1 if ON and 0 if OFF, fills missed time steps with zeros, meant to fill 2d data. Default: 0, OFF. filler works with default 2dc instrument and time steps of 10 secs. Feel free to modify as needed.
        CloudObj.totalvolprecip(abovesize=100, scan=0) size in microns, is the default."""
        #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
        from CloudDataSorter import dNdlogDp2dN
        if len(self.times["verticloud"])<=scan:
            print("[totalvolprecip] No vertical scan for this cloud.")
            return nan
        else:
            ts=self.times["verticloud"][scan]       # times of the vertical scan
            # columns & positions
            calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
            ctim=[i for i,x in enumerate(self.dttl) if x == 'time']; ctim=ctim[0]
            pos2d=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]; 
            if len(pos2d)==0: pos2d=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()];
            try: pos2d=pos2d[0]
            except: print("[totalvolprecip] No precipitation size distribution found.")
            # loading the relevant distribution
            if filler==1:
                fill2d=self.filler(ts=ts)
            elif filler==0: 
                try:
                    fill2d=dict()
                    tix=np.nonzero((self.sd[pos2d]["time"]>=ts[0])*(self.sd[pos2d]["time"]<=ts[1]))[0]
                    if len(tix)==0: return 0.00
                    fill2d["data"]=self.sd[pos2d]["data"][:,tix]; fill2d["time"]=self.sd[pos2d]["time"][tix];
                    fill2d["bins"]=self.sd[pos2d]["bins"]; fill2d["total"]=self.sd[pos2d]["total"][tix];
                except: print("[totalvolprecip] Data not found"); return nan
            else:
                print("[totalvolprecip] filler must be 1 or 0.")
                return nan
            # controlling for empty distributions
            if type(fill2d["time"])==float: Int=nan     # data does not exist
            elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
            else:       # proceeding!
                ix=fill2d["bins"]>=abovesize
                prec=np.sum(dNdlogDp2dN(fill2d)[ix].transpose()*(pi*(fill2d["bins"][ix])**3)/6,1)        # in #/L * um3, assuming the bins are diameters.
                alt=self.data[calt]
                tim=self.data[ctim]
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

    def precipconc(self,abovesize=100,scan=0,minlwc=0.01,filler=0):
        """ This method returns the integrated (with respect to altitude) amount of precipitation drops (number) in a vertical profile in-cloud (lwc>minlwc (default=0.01 g/m3)) for particle greater than "abovesize" and divide by the depth of the cloud. A 2d probe is used. Units are number per litre (if the 2d used has units per litre).
        By default, this in done for the first scan (#0), other scans can be accessed using options.
        filler: 1 if ON and 0 if OFF, fills missed time steps with zeros, meant to fill 2d data. Default: 0, OFF. filler works with default 2dc instrument and time steps of 10 secs. Feel free to modify as needed.
        CloudObj.precipconc(abovesize=100, scan=0, minlwc=0.01) size in microns, lwc in g/m3 is the default."""
        #MMM Add a check for a minimum number of points detected (time points)???? How many is enough?
        if len(self.times["verticloud"])<=scan:
            print("[precipconc] No vertical scan for this cloud.")
            return nan
        else:
            ts=self.times["verticloud"][scan]       # times of the vertical scan
            # columns & positions
            calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
            ctim=[i for i,x in enumerate(self.dttl) if x == 'time']; ctim=ctim[0]
            alt=self.data[calt]; tim=self.data[ctim];
            clwc=[i for i,x in enumerate(self.dttl) if x == 'LWC']; 
            if len(clwc)>1: print("[precipconc] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=ts     # LWC found in the basic data
            elif len(clwc)==1: lwc=self.data[clwc[0]]; lwctime=tim; 
            else:       # looking for LWC in the extradata
                posx=[] 
                for i,L in enumerate(self.extrattl):     # for all extra datasets available
                    posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
                if len(posx)==1: 
                    lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                    j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                    lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                else: print("[precipconc] No LWC (or multiple) was found in the basic data or the extra data."); return nan
            # loading the distribution
            if filler==1:
                fill2d=self.filler(ts=ts)
            elif filler==0: 
                try:
                    pos2d=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]; 
                    if len(pos2d)==0: pos2d=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()];
                    try: pos2d=pos2d[0]
                    except: print("[precipconc] No precipitation size distribution found.")
                    fill2d=dict()
                    tix=np.nonzero((self.sd[pos2d]["time"]>=ts[0])*(self.sd[pos2d]["time"]<=ts[1]))[0]
                    if len(tix)==0: return 0.00
                    fill2d["data"]=self.sd[pos2d]["data"][:,tix]; fill2d["time"]=self.sd[pos2d]["time"][tix];
                    fill2d["bins"]=self.sd[pos2d]["bins"]; fill2d["total"]=self.sd[pos2d]["total"][tix];
                except: print("[precipconc] Data not found"); return nan
            else:
                print("[precipconc] filler must be 1 or 0.")
                return nan
            # controlling for empty distributions
            if type(fill2d["time"])==float: Int=nan     # data does not exist
            elif sum(sum(fill2d["data"]))==0: Int=0.00      # data is zero
            else:       # proceeding!
                prec=dNdlogDp2N(fill2d,abovesize,nan)
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
                #print "Depth used:", Depth, "Difference in Depth:", Depth-(mean(self.props["height"][scan][2:4])-mean(self.props["height"][scan][0:2]))
                Int=Int/Depth   # in units l-1
            return Int

    #########################################################################        
    ###############################  avCDNC  ################################
    #########################################################################

    def avCDNC(self,abovesize=nan,uppersize=nan,inst="FSSP96",scan=0):
        """ This method returns the average cloud droplet concentration in a vertical scan weighted according to the LWC.
        The default options are avCDNC(abovesize=nan,uppersize=nan,inst="FSSP96",scan=0). This method handles the extradata module for lwc.
        "abovesize" is the size above which the concentration is calculated, "uppersize" is the upper size limit. By default the entire available spectrum is integrated. 
        "inst" is the name of the instrument that should be used. Generally an FSSP100.
        "scan" is the scan number, if there are more than one vertical scan.
        The units are the same as for the instrument's concentration (as long as the altitude is in meters)."""
        if len(self.times["verticloud"])<=scan:
            print("[avCDNC] No vertical scan for this cloud.")
            return nan
        else:
            ts=self.times["verticloud"][scan]       # times of the vertical scan
            pos=[i for i,x in enumerate(self.sd) if inst.lower() in x["Distname"].lower()]; 
            if len(pos)<=0: cancel=1; raise ValueError("This instrument name has not been found (inst=%s)." % inst)
            else: pos=pos[0]
            calt=[i for i,x in enumerate(self.dttl) if x == 'altitude']; calt=calt[0]
            ctim=[i for i,x in enumerate(self.dttl) if x == 'time']; ctim=ctim[0]
            alt=self.data[calt]
            tim=self.data[ctim]
            clwc=[i for i,x in enumerate(self.dttl) if x == 'LWC']; 
            if len(clwc)>1: print("[avCDNC] Multiple LWC found."); lwc=np.ones((2,))*NaN; lwctime=ts     # LWC found in the basic data
            elif len(clwc)==1: lwc=self.data[clwc[0]]; lwctime=tim; 
            else:       # looking for LWC in the extradata
                posx=[] 
                for i,L in enumerate(self.extrattl):     # for all extra datasets available
                    posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
                if len(posx)==1: 
                    lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                    j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                    lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                else: print("[avCDNC] No LWC (or multiple) was found in the basic data or the extra data."); return nan
            
            
            timtc=self.sd[pos]["time"]       # time corresponding to the concentration tc
            # find shortest common time stamp
            ta1=np.max([lwctime[0],tim[0],timtc[0],ts[0]]); ta2=np.min([lwctime[-1],tim[-1],timtc[-1],ts[1]]); 
            taix=[nonzero((timtc>=ta1)*(timtc<=ta2))[0]]; 
            usedtime=timtc[taix]
            if len(usedtime)==0: WA=nan    # weighted average is not calculable (no data for the vertical scan)
            else:
                tc=dNdlogDp2N(self.sd[pos],abovesize,uppersize)[taix] # total concentration between above- and upper-size
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
                WA=abs(sum(tc*lwc*dalt))/(self.lwp[scan][0])        # if inst conc is cm-3     (cm-3 lwc*meters / lwc*meters = cm-3)
            return WA
            
            
    ########################################################################        
    ###############################  vpinfo  ###############################
    ########################################################################
    def vpinfo(self,param,base='bg'):
        """ This method returns information on the chosen parameter from CloudObj.dttl in the cloud for all vertical scan. The averaging of the parameter is done in the particular column of the vertical scan.
            Options:
                param: string containing the title of the parameter as found in CloudObj.dttl or CloudObj.extrattl
                base: method to find the cloud base and top. Default is best guess (defBGheight) base='bg'; to use the 4-point method (defheight) base='4point'.
            Returns H:
            H["bottom"]: parameter at the cloud base
            H["top"]: parameter at the cloud top
            H["mean"]: mean parameter through the cloud
            H["median"]: median parameter through the cloud
            H["minimum"]: minimum parameter through the cloud
            H["maximum"]: maximum parameter through the cloud
            H["stdev"]: standard deviation of the parameter through the cloud
            H["delta"]: difference of parameter between the bottom and the top
            H["slope"]: delta divided by the mean thickness
            H["units"]: units of the parameter """
        if type(param)==str: pass
        else: param=str(param)
        H=dict()
        altp=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        tim=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
        T=[i for i,x in enumerate(self.dttl) if x.lower() == param.lower()]
        if len(T)==1: 
            T=T[0]
            Td=self.data[T]
            Tunits=self.dunit[T]
            alt=self.data[altp]
            ta=self.data[tim]
        elif len(T)>1: print("[vpinfo] Parameter %s was found multiple times in the basic data." %(param)); return dict()
        elif len(T)==0:
            posx=[] 
            for i,ttl in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == param.lower()]    # check all titles matching with temperature
            if len(posx)==1: 
                Td=self.extradata[posx[0][0]][posx[0][1]]    # loading the data
                Tunits=self.extraunit[posx[0][0]][posx[0][1]]
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                Tt=self.extradata[posx[0][0]][j]     # loading associated time stamp
                # adapting for too short data for interpolation
                if len(Tt)<2: Td=np.ones((2,))*NaN; Tt=np.array([self.times["cloud"][0][0],self.times["cloud"][0][1]]);
                # adapting the time vector to a common time vector
                ta1=np.max([self.data[tim][0],Tt[0]]); ta2=np.min([self.data[tim][-1],Tt[-1]]);
                ta=self.data[tim][nonzero((self.data[tim]>=ta1)*(self.data[tim]<=ta2))[0]]
                alt=self.data[altp][nonzero((self.data[tim]>=ta1)*(self.data[tim]<=ta2))[0]]
                fT=interpolate.interp1d(Tt,Td,kind='linear')
                Td=fT(ta)
            else: print("[vpinfo] No or multiple %s found in the basic or the extra data." %(param)); return dict()
        
        H["bottom"]=list(); H["top"]=list(); H["mean"]=list(); H["median"]=list(); H["stdev"]=list(); H["delta"]=list(); H["slope"]=list(); H["units"]=list(); H["minimum"]=list();  H["maximum"]=list();
        try:
            for i in range(len(self.times["verticloud"])):
                if base=='4point': cb=self.props["height"][i][1]; ct=self.props["height"][i][2];
                else: cb=self.props["BGheight"][i][0]; ct=self.props["BGheight"][i][1];
                ix=nonzero((alt>=cb)*(alt<=ct)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))[0]
                if len(ix)==0:
                    H["mean"].append(nan); H["median"].append(nan); H["stdev"].append(nan); H["minimum"].append(nan); H["maximum"].append(nan); H["top"].append(nan); H["bottom"].append(nan); H["delta"].append(nan); H["slope"].append(nan); H["units"].append(nan)
                else:
                    H["mean"].append(float(st.nanmean(Td[ix])))
                    H["median"].append(float(st.nanmedian(Td[ix])))
                    H["stdev"].append(float(st.nanstd(Td[ix])))
                    H["minimum"].append(float(np.nanmin(Td[ix])))
                    H["maximum"].append(float(np.nanmax(Td[ix])))
                    if base=='4point': 
                        if len(nonzero((alt>=ct)*(alt<=self.props["height"][i][3])*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))[0])==0: H["top"].append(nan)
                        else: H["top"].append(float(st.nanmedian(Td[nonzero((alt>=ct)*(alt<=self.props["height"][i][3])*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))])))
                        if len(nonzero((alt>=self.props["height"][i][0])*(alt<=cb)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))[0])==0: H["bottom"].append(nan)
                        else: H["bottom"].append(float(st.nanmedian(Td[nonzero((alt>=self.props["height"][i][0])*(alt<=cb)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))])))
                        H["delta"].append(H["bottom"][i]-H["top"][i])
                        H["slope"].append(H["delta"][i]/(np.mean([ct, self.props["height"][i][3]])-np.mean([self.props["height"][i][0], cb])))
                    else: 
                        R=10     # plus/minus R meters around the cloud top
                        if len(nonzero((alt>=ct-R)*(alt<=ct+R)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))[0])==0: H["top"].append(nan)
                        else: H["top"].append(float(st.nanmedian(Td[nonzero((alt>=ct-R)*(alt<=ct+R)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))])))
                        if len(nonzero((alt>=cb-R)*(alt<=cb+R)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))[0])==0: H["bottom"].append(nan)
                        else: H["bottom"].append(float(st.nanmedian(Td[nonzero((alt>=cb-R)*(alt<=cb+R)*(ta>=self.times["verticloud"][i][0])*(ta<=self.times["verticloud"][i][1]))])))
                        H["delta"].append(H["bottom"][i]-H["top"][i])
                        H["slope"].append(float(H["delta"][i]/(ct-cb)))
                    H["units"].append(Tunits)
                del ix
        except: 
            if base=='4point': print("[vpinfo] Height properties must be defined first using the defheight method.")
            else: print("[vpinfo] Height properties must be defined first using the defBGheight method.")
        return H     
        
    ########################################################################        
    ###############################  effrad  ###############################
    ########################################################################
    def effrad(self,inst,bindist='lin'):
        """ This method returns the effective radius for a given instrument for the entire cloud period. The radius is in the same units as the instrument's units (usually micrometers). Note that one might get the effective diameter if the instrument's size bins are diameters.
            example: CloudObj.effrad(inst='FSSP96',bindist='lin')
            bindist is 'lin' if the difference between bins is linearly distributed (FSSPs) and 'log' if they are logarythmically distributed (PCASP)"""
        # according to the formula in https://en.wikipedia.org/wiki/Cloud_drop_effective_radius latest access on Oct 2013.
        [pos,sd]=[[i,sd] for i,sd in enumerate(self.sd) if sd["Distname"].lower() == inst.lower()][0]
        # building dr (dradius) vector
        R=sd["bins"]; t=len(R)
        b=np.zeros([t]);
        h=zeros([t]);
        if bindist=='lin':
            for i in range(1,t):
                b[i]=(R[i-1]+R[i])/2.;
            for i in range(0,t-1):
                h[i]=(R[i+1]+R[i])/2.;
            b[0]=R[0]-(R[1]-h[0]);
            h[t-1]=R[t-1]+(b[t-1]-R[t-2]);
            dR=h-b;
        elif bindist=='log':
            for i in range(1,t):
                b[i]=10**((np.log10(R[i-1])+np.log10(R[i]))/2.);
            for i in range(0,t-1):
                h[i]=10**((np.log10(R[i+1])+np.log10(R[i]))/2.);
            b[0]=10**(np.log10(R[0])+(log10(R[1])-log10(h[1])));
            h[t-1]=10**(np.log10(R[t-1])-(log10(b[t-2])-np.log10(R[t-2])));
            dR=h-b;
        else: print("[effrad] bindist option entry is neither 'lin' or 'log'.")
        # calculating the effective radius
        ER=np.nansum((sd["bins"]**3 *dR) * sd["data"].transpose(),axis=1)/np.nansum((sd["bins"]**2 *dR) * sd["data"].transpose(),axis=1)
        return ER

    ########################################################################        
    ###############################  path3d  ###############################
    ########################################################################
    def path3d(self):
        """ This method will make a 3D plot of the flightpath of the plane during measurement.
        Colours correspond to the time tracking colours (colours may not work if using an old version of matplotlib)."""
        from mpl_toolkits.mplot3d import Axes3D
        
        plat=[i for i,x in enumerate(self.dttl) if x == 'latitude'][0]
        lat=copy.deepcopy(self.data[plat])
        plon=[i for i,x in enumerate(self.dttl) if x == 'longitude'][0]
        lon=copy.deepcopy(self.data[plon])
        palt=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        alt=copy.deepcopy(self.data[palt].data)
        pt=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
        t=copy.deepcopy(self.data[pt].data)

        #FIND THE QUADRANT
        if st.nanmean(lat)>0:
            quad_NS = 'N'
        else:
            quad_NS = 'S'

        if st.nanmean(lon)>0:
            quad_EW = 'E'
        else:
            quad_EW = 'W'


        M=runstats(alt,20)
        alt=np.ma.masked_where((alt>(M[0]+M[1]*1.)+(isnan(alt))), alt)
            
        norm = matplotlib.colors.Normalize(vmin=t[0],vmax=t[-1])

        fig = figure()
        ax = Axes3D(fig)
        majorFormatter_lon = FormatStrFormatter('%.2f '+quad_EW)
        majorFormatter_lat = FormatStrFormatter('%.2f '+quad_NS)
        try:
            if int(matplotlib.__version__[0])>0:
                ax.scatter(abs(lat),abs(lon),alt,lw=0,alpha=1,cmap='spectral',norm=norm,c=t)
                ax.view_init(28,145)
                ax.yaxis.set_major_formatter(majorFormatter_lon)
                ax.xaxis.set_major_formatter(majorFormatter_lat)
                if quad_EW == 'E':
                    ax.set_ylim(ax.get_ylim()[::-1])
                if quad_NS == 'S':
                    ax.set_xlim(ax.get_xlim()[::-1])
            else: 
                print "yup"
                ax.scatter(abs(lat),abs(lon),alt,lw=0,alpha=1,cmap='spectral',norm=norm)
                ax.view_init(28,145)
                ax.yaxis.set_major_formatter(majorFormatter_lon)
                ax.xaxis.set_major_formatter(majorFormatter_lat)
                if quad_EW == 'E':
                    ax.set_ylim(ax.get_ylim()[::-1])
                if quad_NS == 'S':
                    ax.set_xlim(ax.get_xlim()[::-1])
        except: print("[path3d] Error evaluating your version of matplotlib.")
        ax.set_xlabel('Latitude')
        ax.set_ylabel('Longitude')
        ax.set_zlabel('Altitude')
        
        plt.show()
        
    ########################################################################        
    ###############################  filler  ###############################
    ########################################################################
    def filler(self,ts,timestep=10.0,inst='2dc'):
        """This method returns a filled distribution (type dict), filled with zeros when no measurements were returned. 
        This was designed for 2d data: this instrument only reports data when it is detected (will not report zero measurements.
        ts: start and end time over which the data should be filled. Normally a couple of times from CloudObj.times['manoeuvre'][scan#].
        timestep: (default 10.0), in seconds (must be a float), is the time step with which the instrument returns data when measurements are non-zero.
        inst: (default '2dc') is a string format corresponding to the name of the distribution that needs to be filled."""
        newsd=dict()        # define the returning dictionary
        pos2d=[i for i,x in enumerate(self.sd) if x["Distname"].lower()==inst.lower()]; 
        try:                # verifying that the instrument indeed exists and can be found
            pos2d=pos2d[0]
            tim2d=self.sd[pos2d]["time"];
            prec=self.sd[pos2d]["data"];
            ptot=self.sd[pos2d]["total"];
            bins=self.sd[pos2d]["bins"]; 
        except:             # if the instrument does not exist, we go into the no file case (if statement below)
            tim2d=nan
                
        # Testing for no file or empty file
        if type(tim2d)==float:                           # measurements do not exist or not found
            if isnan(tim2d):
                print("[filler] The data were not found.")
                newsd["data"]=nan; newsd["time"]=nan; newsd["bins"]=nan; newsd["total"]=nan;
                return newsd
        else: 
            tix=np.nonzero((tim2d>=ts[0])*(tim2d<=ts[1]))[0]        # index of times in the file within ts[0] and ts[1]
            if sum(tix)==0:         # measurements exists but there is no record ==> the instrument measured nothing = zero
                tstep=timestep/(24*60*60)
                newsd["time"]=arange(ts[0],ts[1],tstep); newsd["bins"]=bins; newsd["total"]=reshape(zeros([len(newsd["time"])]),[-1]); 
                newsd["data"]=reshape(zeros([len(bins),len(newsd["time"])]),[len(bins),-1]); 
                return newsd
        
        # setting up 
        prec=prec[:,tix]; tim2d=tim2d[tix]
        t2d=np.array([]); p2d=np.array([]); pt2d=np.array([]);
        tstep=timestep/(24*60*60)
        dT=((1.20*tstep)>np.diff(tim2d));    # looking for place where the time step was >20% longer than the defined time step.
                
        # starting the new distribution
        if tim2d[0]>ts[0]:     # if the 2d record starts after the start time
            r=arange(tim2d[0]-tstep,ts[0],-tstep)
            if len(r)>0:
                t2d=r[::-1]; pt2d=zeros([len(r)]); p2d=zeros([len(r),len(bins)]);
        if len(t2d)==0:     # if t2d is still empty (e.i. we didn't fill before the first 2d time)
            t2d=tim2d[0]; pt2d=ptot[0]; p2d=reshape(prec[:,0],[1,-1]);
        else:
            t2d=np.hstack((t2d,tim2d[0])); pt2d=np.hstack((pt2d,ptot[0])); p2d=np.vstack((p2d,prec[:,0]));
        
        # continuing the new distribution
        for i,TF in enumerate(dT):   # for each time step in the 2d record
            if TF==True:    # if the step is about equal or smaller to the time step
                t2d=np.hstack((t2d,tim2d[i+1])); pt2d=np.hstack((pt2d,ptot[i+1])); p2d=np.vstack((p2d,prec[:,i+1]));
            else:       # if it needs to be patched with zeros
                r=arange(tim2d[i]+tstep,tim2d[i+1],tstep)
                if tim2d[i+1]-r[-1]<=0.2*tstep: r=r[0:-1]        # cases where the last point is too close to the next measured point (<=20% of the time step)
                t2d=np.hstack((t2d,r,tim2d[i+1])); pt2d=np.hstack((pt2d,zeros([len(r)]),ptot[i+1])); p2d=np.vstack((p2d,zeros([len(r),len(bins)]),prec[:,i+1]));
        
        # finishing the new distribution
        if tim2d[-1]<ts[1]:     # if the 2d record ends before the end of the scan time
            r=arange(tim2d[-1]+tstep,ts[1],tstep)
            if len(r)>0:
                t2d=np.hstack((t2d,r)); pt2d=np.hstack((pt2d,zeros([len(r)]))); p2d=np.vstack((p2d,zeros([len(r),len(bins)])));
        
        # returning
        newsd["data"]=p2d.transpose(); newsd["time"]=t2d; newsd["bins"]=bins; newsd["total"]=pt2d;
        return newsd
        
    #########################################################################        
    ###############################  turbhrz  ###############################
    #########################################################################
    def turbhrz(self,MaxCA=2,MinDist=20000,PtNo=25,Pts=5,Pts2=50):
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
        
        palt=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        ptim=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
        ptas=[i for i,x in enumerate(self.dttl) if 'TAS' in x.upper()][0]
        pudv=[i for i,x in enumerate(self.dttl) if x == 'udvel']
        if len(pudv)==1: pudv=pudv[0]; udv=self.data[pudv][0:-1]
        elif len(pudv)>1: print("[turbhrz] multiple updraft velocities found in the basic data.")
        else:
            posx=[] 
            for i,ttl in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'udvel']    # check all titles matching with updraft velocity
            if len(posx)==1: 
                udv=self.extradata[posx[0][0]][posx[0][1]]    # loading the updraft velocity data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                udvt=self.extradata[posx[0][0]][j]     # loading associated time stamp
                # padding so the interpolation works
                if udvt[0]>self.data[ptim][0]: udv=hstack((nan,udv)); udvt=hstack((self.data[ptim][0],udvt));
                if udvt[-1]<self.data[ptim][-1]: udv=hstack((udv,nan)); udvt=hstack((udvt,self.data[ptim][-1]));
                fudv=interpolate.interp1d(udvt,udv,kind='linear')
                udv=fudv(self.data[ptim])[0:-1]
            else: print("[turbhrz] No updraft velocity (or multiple) found in the basic or the extra data."); return dict()
        Alt=self.data[palt][0:-1]
        A=self.angles
        theta=runstats(A[1],Pts)[0]  # change of altitude
        phi=diff(A[2])
        phi=np.ma.masked_where((phi>350)+(isnan(phi)), phi)
        phi=np.ma.masked_where((phi<-350+(isnan(phi))), phi)
        phi=hstack([NaN,runstats(phi,Pts)[0]])     # change in direction
        phi[(abs(phi)>=89)]=0   # accounts for jump from 0 to 360 degrees. After smoothing, the limit cannot be too high. 
        CAstats=runstats(abs(theta)+abs(phi),Pts2)      # combined angle
        CA=CAstats[0]
#        udv=self.data[pudv][0:-1]        # updraft
        Sp=self.data[ptas][0:-1]
        dt=diff(self.data[ptim]*24*60*60)
        D=cumsum(Sp*dt)
        # Initializing indexes for below cloud and in-cloud measurements.
        bc=ones(np.shape(self.data[ptim]))*False
        ic=ones(np.shape(self.data[ptim]))*False
        for j in range(len(self.times["belowcloud"])):
            bctmp=(self.data[ptim]>=self.times["belowcloud"][j][0])*(self.data[ptim]<=self.times["belowcloud"][j][1])
            try: bc=bc+bctmp;
            except: bc=bctmp;
        for j in range(len(self.times["horicloud"])):
            ictmp=(self.data[ptim]>=self.times["horicloud"][j][0])*(self.data[ptim]<=self.times["horicloud"][j][1])
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
                turb=st.nanmean(runstats(udv[leg],PtNo)[1])
                #print("In-Cloud %d: %.4f km. w' = %0.4f m/s" %(cn,(Dist[-1]-Dist[0])/1000,turb))
                ICturb.append(np.array([turb,(Dist[-1]-Dist[0])]))
                del Dist
        for leg in BCStretches:
            Dist=D[leg]; Dist=Dist[(Dist.mask==False)]
            if (sum(leg)>PtNo and sum(~isnan(udv[leg]))>PtNo):
                turb=st.nanmean(runstats(udv[leg],PtNo)[1])
                #print("Below-Cloud %d: %.4f km. w' = %0.4f m/s" %(cn,(Dist[-1]-Dist[0])/1000,turb))
                BCturb.append(np.array([turb,(Dist[-1]-Dist[0])]))
            else: BCturb.append(np.array([nan,(Dist[-1]-Dist[0])]))
            del Dist
                
        R=dict()
        R["InCloud"]=ICturb; R["BlwCloud"]=BCturb; R["InCloudStretches"]=ICStretches; R["BlwCloudStretches"]=BCStretches
        return R        
        


    #########################################################################        
    ###########################   writeouthdf5   ############################
    #########################################################################
    def writeouthdf5(self,fName=None, fdir=None):
        """This method writes out the content of the cloud in a string representation. This method can be used if one wants to export the cloud's data to a different programming language or make a language-independent backup of their cloud-object's data. Uses Python's repr function.
        fName: Name of the file where the cloud will be saved. Default saves with name SAMAC_Cloud_date_time.
        fdir: directory in which the file will be saved.
        WARNING: this method is designed to handle the current official structure of the object, if the structure has been customized by users, these custom structures will not be saved. This method will need to be modified to accomodate the new structure members.
        WARNING: any masked data will lost their masks. The feature to save masks is in planning for future versions."""
        #MMM writeout the masks of the arrays also!
        import Tkinter, tkFileDialog
        import h5py
        
        if fName==None: 
            if "layered" in self.desc:
                if self.desc["layered"]=='no': layno=''
                else: layno='_'+str(self.desc["layered"])
            else: layno=''
            fName="SAMAC_Cloud_" + self.desc["date"] + "_" + self.desc["time"] + layno      # name of the file where the cloud will be saved
        if fdir==None: fdir=os.getcwd()+'/'
        try:        # verifying if the file already exists
            f=open(fdir+fName+'.hdf5','r'); f.close();
            X=True
        except: X=False;
        Abort=False     # not aborting by default
        if X:
            print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir))
            while X==True:
                Q=raw_input("Do you want to change the file name (N), change directory (D), both (ND) or overwrite (O)?  ")
                if 'n' in Q.lower():
                    fName=str(raw_input("Enter the new file name:  "))
                    try: 
                        f=open(fName+'.hdf5','r'); f.close();
                        X=True
                        print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir))
                    except: X=False
                if 'd' in Q.lower():
                    root = Tkinter.Tk()
                    print("Please choose the directory in which the file should be saved:")
                    pathdef=""
                    while len(pathdef)==0:
                        pathdef = tkFileDialog.askdirectory(parent=root,title='Please select a directory')
                        root.withdraw()
                        if len(pathdef) > 0:
                            print("You chose %s" % pathdef) 
                            fdir=pathdef+'/'
                            try:        # verifying if the file already exists
                                f=open(fdir+fName+'.hdf5','r'); f.close();
                                print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir));
                                X=True
                            except: X=False;
                            break
                        else:
                            pd=raw_input("You must choose a directory. Would you like to choose again?  ")
                            if pd=='y': pass
                            else: break
                if Q.lower()=='o':
                    X=False
                if Q.lower()=='o' or 'n' in Q.lower() or 'd' in Q.lower(): pass
                else: X=False; Abort=True
        if X: print("You are caught in a bug, try fixing the code or try something else. Sorry for the inconvenience.")
        else: 
            if Abort: print("Aborted.")
            else:
                # saving data in file
                f=h5py.File(fdir+fName+'.hdf5','w')
                data=f.create_dataset("data",data=self.data)
                desc=f.create_group("desc"); 
                for key in self.desc: 
                    desc.attrs[key]=np.string_(self.desc[key])
                descedit=f.create_dataset("descedit",data=self.descedit)
                dttl=f.create_dataset("dttl",data=np.string_(self.dttl))
                dunit=f.create_dataset("dunit",data=np.string_(self.dunit))
                for i,actdata in enumerate(self.extradata):
                    extradata=f.create_dataset("extradata"+str(i),data=actdata)
                xdnum=f.create_dataset("extradatanum",data=self.extradatanum)
                for i, actttl in enumerate(self.extrattl):
                    xttl=f.create_dataset("extrattl"+str(i),data=np.string_(actttl))
                for i, actunit in enumerate(self.extraunit):
                    xunit=f.create_dataset("extraunit"+str(i),data=np.string_(actunit))
                orient=f.create_dataset("orient",data=np.string_(self.orient))
                props=f.create_group("props");
                for key in self.props:
                    if isinstance(self.props[key],np.ndarray):
                        dset=props.create_dataset(key,data=self.props[key])
                    elif isinstance(self.props[key],list):
                        for i,L in enumerate(self.props[key]):
                            if isinstance(L,np.ndarray): dset=props.create_dataset(key+str(i),data=L)
                            elif isinstance(L,basestring): props.attrs[key+str(i)]=np.string_(L)
                            else: print("C.props[%s] is a list C.props[%s][%s] is neither an array nor a string" % (key,i))
                    elif isinstance(self.props[key],basestring):
                        props.attrs[key]=np.string_(self.props[key])
                    else: print("C.props[%s] not saved: type %s not handled" %(key,type(self.props[key])))
                for i,sd in enumerate(self.sd):
                    szdst=f.create_group("sd"+str(i))
                    for key in sd:
                        if isinstance(sd[key],np.ndarray):
                            dset=szdst.create_dataset(key,data=sd[key])
                        elif isinstance(sd[key],basestring):
                            szdst.attrs[key]=np.string_(sd[key])
                        else: print("C.sd[%s][%s] is of type %s, this type is not expected and not handled by this method." %(i,key,type(sd[key])))
                times=f.create_group("times")
                for key in self.times:
                    dset=times.create_dataset(key,data=self.times[key])
                f.close()
                print("File saved at %s." % (fdir+fName+'.hdf5'))
                
        # Tips to open and read the hdf5 file in python:
        #    g=h5py.File('filename.hdf5','r')
        # For arrays:
        #    data=g["data"][:]      The whole array is returned in variable "data"
        # For dictionaries containing arrays:
        #    t=g["times"]
        #    t["abovecloud"][:]
        # For dictionaries containing strings:
        #    desc=g["desc"]
        #    desc2=dict()
        #    for key in desc.attrs:
        #        desc2[key]=desc.attrs[key]
        # For lists, now np.strings:
        #    dttl=g["dttl"]
        #    L=str(dttl[...]).replace(" u'","").replace(" '","").replace("'","").strip("[").strip("]").split(",")

            

##############
# Properties #
##############

    #########################################################################        
    #################################  lwp  #################################
    #########################################################################                
    @property
    def lwp(self):
        """ This property calculates the liquid water path from all vertical scans. One LWP per vertical scan will be returned in an array. 
            The LWP for vertical scan #0 can be accessed with CloudObj.lwp[0][0], #1 with CloudObj.lwp[1][0]. """
        # lwc nan and wrong number only are cleared. Negative and small numbers are kept ==> the integral should cancel the noise out.
        LWP=[]
        ct=[i for i,x in enumerate(self.dttl) if x=='time'][0]
        Calt=[i for i,x in enumerate(self.dttl) if x=='altitude']
        if len(Calt)==1: alt=self.data[Calt[0]];
        else: print("[lwp] No or multiple altitudes found."); return LWP
        t=self.data[ct]
        Clwc=[i for i,x in enumerate(self.dttl) if x=='LWC']
        if len(Clwc)>1: print("[lwp] Multiple LWC found."); return LWP
        elif len(Clwc)==1: lwc=self.data[Clwc[0]]; ta=t
        elif len(Clwc)==0:
            posx=[] 
            for i,L in enumerate(self.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(L) if x.lower() == 'lwc']    # check all titles matching with LWC
            if len(posx)==1: 
                lwc=self.extradata[posx[0][0]][posx[0][1]]    # loading the LWC data
                j=[j for j,x in enumerate(self.extrattl[i]) if x.lower() == 'time'][0]
                lwctime=self.extradata[posx[0][0]][j]     # loading associated time stamp
                # adapting for too short data for interpolation
                if len(lwctime)<2: lwc=np.ones((2,))*NaN; lwctime=np.array([self.times["cloud"][0][0],self.times["cloud"][0][1]]);
                # adapting the time vector to a common time vector
                ta1=np.max([t[0],lwctime[0]]); ta2=np.min([t[-1],lwctime[-1]]);
                ta=t[nonzero((t>=ta1)*(t<=ta2))[0]]
                alt=alt[nonzero((t>=ta1)*(t<=ta2))[0]]
                f=interpolate.interp1d(lwctime,lwc,kind='linear')
                lwc=f(ta)
            else: print("[lwp] No LWC was found in the basic data or multiple were found in the extra data."); return LWP
        for j in self.times["verticloud"]:
            alti=alt[np.nonzero((ta>=j[0])*(ta<=j[1]))]; dalt=np.diff(alti);
            lwci=lwc[np.nonzero((ta>=j[0])*(ta<=j[1]))]; lwci=lwci[1:];
            LWP.append(abs(sum(lwci*dalt)))
        LWP=np.reshape(np.array(LWP),(-1,1))
        return LWP

    #########################################################################        
    #############################  thickness  ###############################
    #########################################################################
    @property
    def thickness(self):
        """ This property calculates the max and min cloud thickness (in meters) based on the height properties. The method defheight must be run at least once on the object for this property to work.
        The property can be accessed as CloudObj.thickness["min"] or CloudObj.thickness (returns a dict) """
        H=dict()
        try:
            H["max"]=list(); H["min"]=list();
            for i in range(len(self.props["height"])):
                H["max"].append(self.props["height"][i][3]-self.props["height"][i][0])
                H["min"].append(self.props["height"][i][2]-self.props["height"][i][1])
        except: print("[thickness] Height properties must be defined first using the defheight method.")
        return H

    #########################################################################        
    ###############################  tempinfo  ##############################
    #########################################################################
    @property
    def tempinfo(self):
        """ This property returns information on the temperature in the cloud (all given in Celsius). Note that the temperature is averaged over the entire cloud time at the altitude required (bottom, top or in-cloud) - not the case using CloudObj.vpinfo(temperature).
            CloudObj.tempinfo["bottom"]: temperature at the cloud base
            CloudObj.tempinfo["top"]: temperature at the cloud top
            CloudObj.tempinfo["mean"]: mean temperature through the cloud (in cloud)
            CloudObj.tempinfo["median"]: median temperature through the cloud (in cloud)
            CloudObj.tempinfo["stdev"]: standard deviation of the temperature through the cloud (in cloud)
            CloudObj.tempinfo["delta"]: difference of temperature between the bottom and the top
            CloudObj.tempinfo["slope"]: delta divided by the mean thickness
            The property can be accessed as e.g. CloudObj.tempinfo["bottom"] or CloudObj.tempinfo (dictionary) """
        H=dict()
        H["bottom"]=list(); H["top"]=list(); H["mean"]=list(); H["median"]=list(); H["stdev"]=list(); H["delta"]=list(); H["slope"]=list(); 
        alt=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        T=[i for i,x in enumerate(self.dttl) if x == 'temperature'][0]
        try:
            for i in range(len(self.props["height"])):
                ix=nonzero((self.data[alt]>=self.props["height"][i][1])*(self.data[alt]<=self.props["height"][i][2]))
                H["bottom"].append(float(st.nanmedian(self.data[T][nonzero((self.data[alt]>=self.props["height"][i][0])*(self.data[alt]<=self.props["height"][i][1]))])))
                H["top"].append(float(st.nanmedian(self.data[T][nonzero((self.data[alt]>=self.props["height"][i][2])*(self.data[alt]<=self.props["height"][i][3]))])))
                H["mean"].append(float(st.nanmean(self.data[T][ix])))
                H["median"].append(float(st.nanmedian(self.data[T][ix])))
                H["stdev"].append(float(st.nanstd(self.data[T][ix])))
                H["delta"].append(H["bottom"][i]-H["top"][i])
                H["slope"].append(H["delta"][i]/(np.mean([self.props["height"][i][2], self.props["height"][i][3]])-np.mean([self.props["height"][i][0], self.props["height"][i][1]])))     # Celsius/meter
                del ix
        except: print("[tempinfo] Height properties must be defined first using the defheight method.")
        return H
        
    ########################################################################        
    #############################  presinfo  ###############################
    ########################################################################
    @property
    def presinfo(self):
        """ This property returns information on the pressure (in mb) in the cloud. Note that the pressure is averaged over the entire cloud time at the altitude required (bottom, top or in-cloud). - not the case using CloudObj.vpinfo(pressure)
            CloudObj.presinfo["bottom"]: pressure at the cloud base
            CloudObj.presinfo["top"]: pressure at the cloud top
            CloudObj.presinfo["mean"]: mean pressure through the cloud
            CloudObj.presinfo["median"]: median pressure through the cloud
            CloudObj.presinfo["stdev"]: standard deviation of the pressure through the cloud
            CloudObj.presinfo["delta"]: difference of pressure between the bottom and the top
            CloudObj.presinfo["slope"]: delta divided by the mean thickness """
        H=dict()
        H["bottom"]=list(); H["top"]=list(); H["mean"]=list(); H["median"]=list(); H["stdev"]=list(); H["delta"]=list(); H["slope"]=list(); 
        alt=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        T=[i for i,x in enumerate(self.dttl) if x == 'Pressure'][0]
        try:
            for i in range(len(self.props["height"])):
                ix=nonzero((self.data[alt]>=self.props["height"][i][1])*(self.data[alt]<=self.props["height"][i][2]))
                H["bottom"].append(float(st.nanmedian(self.data[T][nonzero((self.data[alt]>=self.props["height"][i][0])*(self.data[alt]<=self.props["height"][i][1]))])))
                H["top"].append(float(st.nanmedian(self.data[T][nonzero((self.data[alt]>=self.props["height"][i][2])*(self.data[alt]<=self.props["height"][i][3]))])))
                H["mean"].append(float(st.nanmean(self.data[T][ix])))
                H["median"].append(float(st.nanmedian(self.data[T][ix])))
                H["stdev"].append(float(st.nanstd(self.data[T][ix])))
                H["delta"].append(H["bottom"][i]-H["top"][i])
                H["slope"].append(H["delta"][i]/(np.mean([self.props["height"][i][2], self.props["height"][i][3]])-np.mean([self.props["height"][i][0], self.props["height"][i][1]])))     # Celsius/meter
                del ix
        except: print("[presinfo] Height properties must be defined first using the defheight method.")
        return H        

    #########################################################################        
    ############################  precipblwTF  ##############################
    #########################################################################
    @property
    def precipblwTF(self):
        """ This property tells whether precipitation below cloud was observed during the entire cloud period: true or false for >100 & >200 microns.
            The dictionary entries can be accessed through CloudObj.precipblwTF["100"] or CloudObj.precipblwTF["200"]. """
        pos=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]       
        if len(pos)==0: pos=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()]
        if len(pos)==0: raise ValueError("[precipblwTF] No precipitation distribution found.")
        elif len(pos)==1: pos=int(pos[0])
        elif len(pos)>1:
            for p in pos:
                print("%d - %s" % (p, self.sd[p]["Distname"]))
            pos=int(raw_input("What is the prefered precipitation distribution? Enter the number."))
        else: pos=nan   
        if type(self.sd[pos]["time"])==float or pos==nan: 
            precb=dict(); precb["100"]=list(); precb["200"]=list();
            precb["100"].append(nan); precb["200"].append(nan);
        elif np.shape(self.sd[pos]["time"])[0]<=1: 
            precb=dict(); precb["100"]=list(); precb["200"]=list();
            precb["100"].append(nan); precb["200"].append(nan);
        else:
            P100=dNdlogDp2N(self.sd[pos],100,nan)
            #P100=np.ma.masked_less_equal(P100,1e-2);
            P100[P100<=1e-2]=0
            P200=dNdlogDp2N(self.sd[pos],200,nan)
            P200[P200<=1e-2]=0
            #P200=np.ma.masked_less_equal(P200,1e-2);
            t100=0; t200=0; precb=dict();
            for b in range(np.shape(self.times['belowcloud'])[0]):
                if len(P100[(self.sd[pos]["time"]<=self.times['belowcloud'][b][1])*(self.sd[pos]["time"]>=self.times['belowcloud'][b][0])])==0: pass
                else:
                    t100=t100+np.nansum(P100[(self.sd[pos]["time"]<=self.times['belowcloud'][b][1])*(self.sd[pos]["time"]>=self.times['belowcloud'][b][0])])
                    t200=t200+np.nansum(P200[(self.sd[pos]["time"]<=self.times['belowcloud'][b][1])*(self.sd[pos]["time"]>=self.times['belowcloud'][b][0])])    
            if t100>0: precb["100"]=True
            else: precb["100"]=False
            if t200>0: precb["200"]=True
            else: precb["200"]=False
        return precb

    #########################################################################        
    ###########################  precipblwvpTF  #############################
    #########################################################################
    @property
    def precipblwvpTF(self):
        """ This property tells whether precipitation below cloud was observed for a given vertical scan and below that scan: true or false for >100 & >200 microns.
            It is important to note that the plane not detecting rain when it passed below the vertical scan doesn't mean that the shaft has not precipitated or will not precipitate.
            The dictionary entries can be accessed through CloudObj.precipblwvpTF["100"] or CloudObj.precipblwvpTF["200"]. """
        pos=[i for i,x in enumerate(self.sd) if 'p1' in x["sdtype"].lower()]  
        if len(pos)==0: pos=[i for i,x in enumerate(self.sd) if 'p' in x["sdtype"].lower()]  
        if len(pos)==1: pos=int(pos[0])
        else: pos=nan; print("[precipblwvpTF] Multiple precipitation distribution found, no priority is set.")
        if type(self.sd[pos]["time"])==float or pos==nan: 
            precb=dict(); precb["100"]=list(); precb["200"]=list();
            for r in range(len(self.times["verticloud"])):
                precb["100"].append(nan); precb["200"].append(nan);
        elif np.shape(self.sd[pos]["time"])[0]<=1: 
            precb=dict(); precb["100"]=list(); precb["200"]=list();
            for r in range(len(self.times["verticloud"])):
                precb["100"].append(nan); precb["200"].append(nan);
        else:
            precb=dict();
            precb["100"]=list(); precb["200"]=list();
            P100=dNdlogDp2N(self.sd[pos],100,nan)
            #P100=np.ma.masked_less_equal(P100,1e-2);   # nansum doesn't handle masked arrays.
            P100[P100<=1e-2]=0
            P200=dNdlogDp2N(self.sd[pos],200,nan)
            P200[P200<=1e-2]=0
            #P200=np.ma.masked_less_equal(P200,1e-2);
            for k in range(len(self.times["verticloud"])):
                tsix=list()
                G=self.belowprof()[k]["index"]
                Gt=self.belowprof()[k]["times"]
                H=np.nonzero(diff(G)>5)[0]
                if np.shape(H)[0]==0: tsix.append([Gt[0],Gt[-1]])
                else:
                    for c in range(len(H)+1):
                        if c==0: tsix.append([Gt[0],Gt[H[c]]])
                        elif c<=len(H)-1: tsix.append([Gt[H[c-1]+1],Gt[H[c]]])
                        elif c==len(H): tsix.append([Gt[H[c-1]+1],Gt[-1]])
                t100=0; t200=0;
                for b in range(len(tsix)):
                    if len(P100[(self.sd[pos]["time"]<=tsix[b][1])*(self.sd[pos]["time"]>=tsix[b][0])])==0: pass
                    else:
                        t100=t100+np.nansum(P100[(self.sd[pos]["time"]<=tsix[b][1])*(self.sd[pos]["time"]>=tsix[b][0])])
                        t200=t200+np.nansum(P200[(self.sd[pos]["time"]<=tsix[b][1])*(self.sd[pos]["time"]>=tsix[b][0])])
                if t100>0: precb["100"].append(True)
                else: precb["100"].append(False)
                if t200>0: precb["200"].append(True)
                else: precb["200"].append(False)
        return precb


    #########################################################################        
    ###############################  angles  ################################
    #########################################################################
    @property
    def angles(self):
        """ This property calculates the vertical angle of the aircraft (or other real or hypothetical aerial vehicle) based on TAS and Altitude change.
        It also calculates the horizontal angle of the aircraft based on the longitude and latitude. The TAS must be in the basic data (extradata module not handled).
        Returns a list with: time, vertical angle, horizontal angle (with one point less than the time in CLOUD.data)."""
        
        # vertical angle calculations
        palt=[i for i,x in enumerate(self.dttl) if x == 'altitude'][0]
        ptim=[i for i,x in enumerate(self.dttl) if x == 'time'][0]
        ptas=[i for i,x in enumerate(self.dttl) if 'TAS' in x.upper()]
        if len(ptas)==1: ptas=ptas[0]
        elif len(ptas)==0:
            print("[angles] True Air Speed (TAS) not found.")
            for i,title in enumerate(self.dttl):
                print("%d - %s (%s)" % (i,title,self.dunit[i]))
            try: ptas=int(raw_input("Enter the number of the TAS (press return if none correspond):  "))
            except: print("[angles] no TAS available: failure of method angles.")
        elif len(ptas)>=2:
            print("[angles] Multiple possible True Air Speed (TAS) were found.")
            for i in range(len(ptas)):
                print("%d - %s (%s)" % (i,self.dttl[ptas[i]],self.dunit[ptas[i]]))
            p=int(raw_input("Enter the number of the TAS (press return if none correspond):  "))
            try: ptas=ptas[0]
            except: print("[angles] no TAS available: failure of method.")
        else: print("[angles] Very strange indeed.")
        if type(ptas)==int: pass
        else: return "Error: clould not find TAS position!"
        dalt=np.diff(runstats(self.data[palt],5)[0])
        dt=np.diff(self.data[ptim])*24*60*60        # time in seconds
        tas=self.data[ptas]; tas=tas[0:-1]      # removing the last point so it has the same length as the differentiated vectors
        tas=runstats(tas,5)[0]      # running average of tas
        # the angle zero is when the aircraft is parallel to the ground
        VA=360*np.arcsin(dalt/(tas*dt))/(2*pi)      # plane angle, in degrees (compared to horizontal)
        
        # horizontal angle calculations
        plat=[i for i,x in enumerate(self.dttl) if x == 'latitude'][0]
        plon=[i for i,x in enumerate(self.dttl) if x == 'longitude'][0]
        # Let's define East as the angle zero.
        dlat=np.diff(runstats(self.data[plat],5)[0]); 
        dlon=np.diff(runstats(self.data[plon],5)[0]); 
        # assuming NA quadrant  MMM - this needs to be fixed (although right now the data has only positive latitudes/longitudes)
        # assuming parallel enough longitude lines
        HA=360*np.arctan(abs(dlat/dlon))/(2*pi)     # calculating everything as a first quadrant
        HA[((dlat>=0)*(dlon>0))]=180-HA[((dlat>=0)*(dlon>0))]   # 2nd quadrant
        HA[((dlat<0)*(dlon>0))]=180+HA[((dlat<0)*(dlon>0))]     # 3rd quadrant
        HA[((dlat<0)*(dlon<0))]=360-HA[((dlat<0)*(dlon<0))]   # 4th quadrant
        HA[((dlat<0)*(dlon==0))]=270
        HA[((dlat>0)*(dlon==0))]=90
        
        return [self.data[ptim][0:-1], VA, HA]
        
        

# mainifization!
if __name__ == 'main':
    a = Cloud({})

