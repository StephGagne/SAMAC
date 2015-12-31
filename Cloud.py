#! /usr/bin/python

#################################################################################
##    Copyright 2013-2015 Stephanie Gagne, 2013-2014 Landan MacDonald          ##
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
# add2basic(self,newdata,newttl,newunit)
# describe(self)
# addsd(self)
# defheight(self)
# defBGheight(self)
# lwcflag(self,base=bg)
# timechange(self)
# MaskData(self)


from pylab import *     # syntax not recommended consider changing for import pylab (and check if all names/functions still work)
import copy as copy
import os

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import math
import xlrd
import scipy.stats.stats as st
from scipy import interpolate

import samac

np.seterr(divide='ignore',invalid='ignore')        # ignoring warning of divisions by zero, or invalid values (such as NaNs)


### SAMAC version ###
__version__="1.0.0" #
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
                samac.mapcloud(self,1)
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
                        samac.mapcloud(self,1)
                        self.desc["humanplace"]=raw_input("Name the place (meaningful to humans) where this cloud occured.  ")
                        plt.close(1002)
                elif self.desc.has_key(RR):
                    print("Old entry is %s" % self.desc[RR])
                    self.desc[RR]=raw_input("New entry for %s:  " % RR)
                else: print("Typo? Could not find the mentionned key.")


    #########################################################################        
    ###############################  addsd  #################################
    #########################################################################
    def addsd(self):
        """ This method is used to add a size distribution to the cloud object from an excel (.xls) file (assumed in dN/dlogDp format, conversions can be done later using codes in CloudDataSorter). 
            The added size distribution will be found appended to CloudObj.sd.
            The source file as well as the source sheet will be saved along with the distribution. """
        import Tkinter, tkFileDialog
        # choosing directory and file
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
        except: print("[addsd] Failed!")
        
        # checking format
        if smpsfn[-3:].lower()=='xls' or smpsfn[-4:].lower()=='xlsx': xlf=1
        elif smpsfn[-3:].lower()=='ods': raise IOError("Open document spreadsheets (ods) are not accepted. Convert to xls/xlsx or to a text file.")
        else: xlf=0
        
        # transforming excel time zero into an ordinal format
        xltime=dt.datetime(1899,12,31,0,0,0).toordinal()-1
        
        ##### Loading File #####
        
        ### Excel File Loop ###
        if xlf==1:
            import xlrd
            file = xlrd.open_workbook(smpsfn) # open excel workbook
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
            
        ### Non-Excel File Loop ### 
        else:
            inp = open(smpsfn,'r') # open file
            sds=[]
            for line in inp.readlines():
                ltmp=line.strip("\n"); ltmp=ltmp.split('\t');
                sds.append(ltmp)
                inp.close()       # consider using the "with open(....) as inp" instead (see 
                                  # http://stackoverflow.com/questions/4599980/python-close-file-descriptor-question last visited Dec 2013)
            try: 
                if np.shape(sds)[0]==1 or np.shape(sds)[1]==1:   # if the delimiter was not a tab
                    print("%s \n %s" % (sds[0],sds[-1]))
                    delima=raw_input("What is the file's delimiter?  ")
                    sds=[]; inp = open(smpsfn,'r') # open file
                    for line in inp.readlines():
                        ltmp=line.strip("\n"); ltmp=ltmp.split(delima);     # consider quicker loading functions such as np.loadtxt
                        sds.append(ltmp)
                        inp.close()
            except:
                try:
                    sds=[]; inp=open(smpsfn,'r') # open file
                    for line in inp.readlines():
                        ltmp=line.strip("\n"); ltmp=ltmp.split();
                        sds.append(ltmp)
                        inp.close()
                except: print("[addsd] Failed to find a working loading file.")

        ##### Common Code: finding titles, units, data #####
        ttl='n'; i=0; unt='n'; dat='n'; # initiating the finding title row loop
        while ttl=='n' or unt=='n' or dat=='n':
            if xlf==1: print("line %d- %s" % (i,sds.row_values(i,0,7)))
            else: print("line %d- %s" % (i,sds[i][0:5]))
            Q=raw_input("Is this the title line [t], the size bin line [s], both [ts] or neither [enter]. For the first line of data: [d]. [0] to start again.  ")
            if Q=='0':
                i=0
            elif Q=='ts': ttll=i; untl=i; ttl='y'; unt='y'; i+=1
            elif Q=='t': ttll=i; i+=1; ttl='y'
            elif Q=='s': untl=i; i+=1; unt='y';
            elif Q=='d': datl=i; i+=1; dat='y'
            else: i+=1
        del ttl, i, unt, dat

        ##### Fishing for the Data #####
        
        ### Excel Format ###
        if xlf==1:
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
                else: print("[addsd] Strange times...")
            else:
                sdtmp["time"]=self.data[0]

        ### Non-Excel Format ###
        else:
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
            # MMM I didn't feel like wasting a lot of time making the time retrieval more robust. Might want to do that later.
            try: 
                Q=map(int,Q.split(','))
                if len(Q)==1: Q=Q[0]
                QQ=raw_input("Is %s an ordinal format? y/[n]  " % sdd[1,Q])
                if QQ=='y':
                    sdtmp["time"]=np.array(sdd[:,Q]);
                    if any(sdtmp["time"]>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"]<dt.datetime(2050,1,1,0,0,0).toordinal()): pass
                    elif any(sdtmp["time"][0]+xltime>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"][0]+xltime<dt.datetime(2050,1,1,0,0,0).toordinal()): 
                        sdtmp["time"]=sdtmp["time"]+xltime
                    else: print("[addsd] Strange ordinal...")
                else: 
                    T1=raw_input("What is the date? (yyyymmdd)  ")
                    T1=samac.todatenum(dt.datetime.strptime(T1, '%Y%m%d'))
                    T2f=raw_input("Enter the time format of this %f (%s) [%%Y (4 digit) %%y (2 digits) %%m %%d %%H %%M %%S (enter exact format with all -, /, spaces, etc.)]:  " % (sdd[1,Q], sddttl[Q]))
                    ttmp=[]
                    for i in range(nupo):
                        T2=samac.todatenum(dt.datetime.strptime(str(int(sdd[i,Q])), T2f))
                        TT=T1+T2-np.floor(T2)
                        ttmp.append(TT)
                    sdtmp["time"]=np.array(ttmp)
            except: sdtmp["time"]=np.ones(nupo)*NaN;
        

        ##### Common to all formats: treating & placing the data #####

        # format check
        logconc=raw_input("Are the concentrations in dN [1] format, dN/dDp [2] or dN/dlogDp [3]?  ")
        if logconc=="1": 
            sdtmp["data"]=samac.dN2dNdlogDp(sdtmp)
        elif logconc=="2":
            sdtmp["data"]=samac.dNdDp2dNdlogDp(sdtmp)
        elif logconc=="3": pass
        else: print "[addsd] The concentrations have been assumed to be in dN/dlogDp format"
        
        #total concentration for size distribution
        Q=raw_input("Enter the column with the total concentration (press enter if none):  ")
        try:
            Q=int(Q)
            if len(sdd[Q])==nupo: sdtmp["total"]=np.array(sdd[Q])       # this if statement handles different versions of numpy
            else: sdtmp["total"]=np.array(sdd[:,Q])
        except: 
            try: 
                sd4tot=dict(); sd4tot["bins"]=sdtmp["bins"]; sd4tot["data"]=sdtmp["data"]; sd4tot["time"]=sdtmp["time"];
                sdtmp["total"]=samac.dNdlogDp2N(sd4tot,nan,nan)
                del sd4tot
            except: sdtmp["total"]=np.ones(nupo)*NaN;

        # MMM need to add median and mean diameters

        # Distribution name
        sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
        sdtmp["units"]=raw_input("What are the units of the size distribution? (e.g. cm-3)  ")
        sdtmp["sdtype"]=raw_input("What kind of size distribution is it? (A=aerosol, C=cloud droplet, P=precipitation)  ")
        # appending distribution to Cloud object.
        sdtmp["total"]=sdtmp["total"][nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        if np.shape(sdtmp["data"])[0]==len(sdtmp["time"]):
            sdtmp["data"]=sdtmp["data"].transpose()
        sdtmp["data"]=sdtmp["data"][:,nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["time"]=sdtmp["time"][nonzero((sdtmp["time"]>=self.times["cloud"][0][0])*(sdtmp["time"]<=self.times["cloud"][0][1]))[0]]
        sdtmp["sourcefile"]=smpsfn; 
        
        # excel file: last-minute additions
        if xlf==1: 
            sdtmp["sourcesheet"]=smpssheetname
            try: del dirnames, pathnames, fnames, kw, smpsfn, sheetnames, Q, sds, sdd, nupo, c1,c2, colemp, cotype
            except: print("[addsd] Failed to clear memory")
        
        # Adding to the cloud object
        self.sd.append(sdtmp)

    

    #########################################################################        
    #############################   defheight   #############################
    #########################################################################
    def defheight(self):
        """ Interactive method. Sets up or clears and resets the position of the cloud. You will be asked to click 4 times per vertical profile: twice at the bottom (lower and upper guesses), and twice at the top (lower and upper guesses). This data will be stored under CloudObj.props["height"] with one such array per vertical profile. 
        If you need to add a profile, we recommend you insert a space in the list, which can then be overwritten without touching the others."""
        plt.close('all')
        import matplotlib._pylab_helpers
        figalready=[manager.canvas.figure for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        samac.vprof(self,1)
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
        samac.vprof(self,1)
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
              base=bg: to use the best guess base height (defBGheight))
        flags: 1-Good quality (young/smooth) \t2-Entrainment at the top \t3-Entrainment at the bottom \t4-Entrainment in the cloud \t5-missing data \t6-Other/Bad/Evaporated/Incomplete Profile \t7-lwc>adiabatic"""
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
        ad=samac.adlwc(self,base=base)
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
                    print("1-Good quality (young/smooth) \t2-Entrainment at the top \t3-Entrainment at the bottom \t4-Entrainment in the cloud \t5-missing data \t6-Other/Bad/Evaporated/Incomplete Profile \t7-lwc>adiabatic \nEnter the flags separated with commas.  ")
                    L=np.array(map(int,raw_input("What are the vertical profile's flags belonging to scan %d.  " % (Q)).split(',')))
                    self.props["lwcflags"][Q]=L
                    plt.close(3)
            except: print("[lwcflag] The sequence was aborted.")
        except:
            self.props["lwcflags"]=list()
            if len(self.times["verticloud"])>0: print("1-Good quality (young/smooth) \t2-Entrainment at the top \t3-Entrainment at the bottom \t4-Entrainment in the cloud \t5-missing data \t6-Other/Bad/Evaporated/Incomplete Profile \t7-lwc>adiabatic \nEnter the flags separated with commas.  ")
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
    ############################   timechange   #############################
    #########################################################################
    def timechange(self,Rtime=None):
        """ This method is used to change one of the time frames in CloudObj.times. The user will be asked to choose the scan type and scan number of which times he or she wants to modify (addition and deletion are also possible. 
        Use CloudObj(Rtime=1) to get times in Hour:Minute format."""
        plt.close('all')
        samac.overview(self,Rtime=Rtime)
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
            samac.overview(self,interact=1,Rtime=Rtime)
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
            if X1<0: raise ValueError("[MaskData] the number cannot be negative."); 
            try:
                self.extradata[X1]=np.ma.array(self.extradata[X1])
                dv=self.extradata[X1][Q];   # data vector
                k=[k for k,x in enumerate(self.extrattl[X1]) if x.lower() == 'time'][0]
                dvt=self.extradata[X1][k]   # time vector
            except: raise RuntimeError("[MaskData] Problem loading data. Was your entry correct?"); 
        else: 
            try: 
                Q=int(Q)
                if Q<0: print("[MaskData] No negative numbers accepted."); crash
                dv=np.ma.array(self.data[Q])
                post=[i for i,x in enumerate(self.dttl) if x == 'time'][0] # position of time
                dvt=self.data[post]
            except: raise RuntimeError("[MaskData] Wrong entry or unexpected behaviour. Was your entry the right integer? Is there an entry called 'time'?");
        
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
                    M=samac.runstats(dv,numra)
                    dv=np.ma.masked_where((dv>(M[0]+M[1]*numsd)+(isnan(dv))), dv)
                    dv=np.ma.masked_where((dv<(M[0]-M[1]*numsd)+(isnan(dv))), dv)
                    plt.close(30)
                    plt.figure(30)
                    plot(dvt,dv,'b.')    
                    G2=raw_input("Do you need to mask more data?  y/[n]  ")
                except: print("[MaskData] Operation failed")
            if G2=='y':
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
                            dv.mask[ix]='True'  
                            Q3=raw_input("Do you want to remove more points? [y]/n  ")
                            if Q3.lower()=='n': Fin=1
                        plt.close(30)
                    else: pass      # if there is no need to modify the data after the standard deviation treatment.
            if X1<0: self.data[Q]=dv
            else: self.extradata[X1][Q]=dv
            figure(30)
            plot(dvt,dv,'b.')


       

# mainifization!
if __name__ == 'main':
    a = Cloud({})

