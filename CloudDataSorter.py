

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



# In this file:
# CloudDataSorter
# todatenum
# dNdlogDp2N
# dNdlogDp2dN
# dN2dNdlogDp
# runstats
# AdiabLWC


from pylab import *
import datetime as dt
import numpy as np
import copy as copy    
import scipy.stats.stats as st
import os


  
#MMM improvements could include changing titles that are identical to standard titles.  

def CloudDataSorter(data_,dttl_,dunit_):
# data is the data, dttl are titles, dunit are the units
# this function aims at taking non-standard data and giving it a standard shape.


    from CloudDataSorter import todatenum
    import sys
    sys.path.append("/home/shared/Convair_Maritimes/Calibrations_Keys/")
    from InstrumentSpecs import InstrumentSpecs
    
   
    # copying
    data=copy.copy(data_); dttl=copy.copy(dttl_); dunit=copy.copy(dunit_);
    # orientation: having each value on a line (default with xlread and make code shorter)
    if np.shape(data)[0]==len(dttl):pass
    elif np.shape(data)[1]==len(dttl): data=data.transpose()
    # number of points of each value (number of times)
    nupo=np.shape(data)[1]
    # converting into a list to be able to use pop
    data=data.tolist()

    
    ####### Time #######
    Jttl=" ".join(dttl).lower(); Junit=" ".join(dunit).lower();
    TQ=0
    if ("ordinal" in Jttl) or ("ordinal" in Junit) or ("datetime" in Jttl) or ("datetime" in Junit):
        for i in range(len(dttl)):
            if TQ==0:
                if (("ordinal" in dttl[i]) or ("ordinal" in dunit[i]) or ("datetime" in dttl[i]) or ("datetime" in dunit[i])):
                    Sdata=np.array(data.pop(i)); TQ=1;   # this has to be the first thing we do to Sdata
                    trash=dttl.pop(i); trash=dunit.pop(i);
    else:
        if ("time" in Jttl) or ("time" in Junit):
            for i in range(len(dttl)):
                if TQ==0:
                    if (("time" in dttl[i]) or ("time" in dunit[i])):
                        Q=raw_input("Does this %s look like a ordinal time format to you? y/[n]  " % data[i][10])
                        if Q=='y':
                            Sdata=np.array(data.pop(i)); TQ=1   # this has to be the first thing we do to Sdata
                            trash=dttl.pop(i); trash=dunit.pop(i);
        if TQ==0:
            for i in range(len(dttl)):
                    print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            Q=raw_input("Are any of these describing ordinal time formats? [n]/number:  ")
            try: 
                Sdata=np.array(data.pop(int(Q)))   # this has to be the first thing we do to Sdata
                Sttl=dttl.pop(int(Q)); Sunit=dunit.pop(int(Q));
            except: 
                Q=raw_input("Is there a 'date (# of days)' and a 'time of day (# of days)' column?  [n]/date#,time#  ")
                try: 
                    Ptime=map(int,Q.split(','))
                    Sdata=np.array(data[Ptime[0]])+np.array(data[Ptime[1]])
                    Ptime.sort(); Ptime.reverse();
                    trash=data.pop(Ptime[1]); trash=data.pop(Ptime[0]);
                    trash=dttl.pop(Ptime[1]); trash=dunit.pop(Ptime[1]);
                    trash=dttl.pop(Ptime[0]); trash=dunit.pop(Ptime[0]);
                    del trash, Ptime, Q
                except: 
                    Q=raw_input("Are there year, month, day, hour, min, sec columns?  [n]/year#,month#,day#,hour#,min#,sec#  ")
                    try: 
                        Ptime=map(int,Q.split(','))
                        Y=data[Ptime[0]]; M=data[Ptime[1]]; D=data[Ptime[2]]; H=data[Ptime[3]]; Mi=data[Ptime[4]]; S=data[Ptime[5]];
                        Ptime.sort(); Ptime.reverse();
                        for col in Ptime:
                            trash=data.pop(col); trash=dttl.pop(col); trash=dunit.pop(col); 
                        tmp=[]
                        for i in range(len(Y)):
                            tmp.append(todatenum(dt.datetime(int(Y[i]),int(M[i]),int(D[i]),int(H[i]),int(Mi[i]),int(S[i]))))
                        Sdata=np.array(tmp)
                        del trash, tmp, Ptime, Y, M, D, H, Mi, S
                    except:
                        Q=raw_input("Are there another combination of dates and time?  [n]/number of ALL and ONLY necessary columns separated by commas ','  ")
                        try: 
                            Ptime=map(int,Q.split(','))
                            datetype=[];
                            tmp=[]
                            Datecol=0; Timecol=0;
                            for i in range(len(Ptime)):
                                DaTi=raw_input("What type of date is this: %s (%s)? 'datenum' (ordinal), 'timenum' (fraction of day) or %%Y (4 digit) %%y (2 digits) %%m %%d %%H %%M %%S (enter exact format with all -, /, spaces, etc.):  " % (int(data[Ptime[i]][10]), dttl[Ptime[i]]))
                                if DaTi.lower()=='datenum': Datecol=Ptime.pop(i)
                                elif DaTi.lower()=='timenum': Timecol=Ptime.pop(i)
                                else: datetype.append(DaTi)
                            dateform=" ".join(datetype);
                            if len(Ptime)>0:
                                req=0
                                for i in range(nupo):
                                    timedata=[]
                                    for j in Ptime:
                                        try: timedata.append(str(int(data[j][i])))
                                        except: req=1; print("Incongruity in the file at time #%d! Please correct file." % i)
                                    datedata=" ".join(timedata)
                                    if 'y' in dateform.lower():     # date is in dateform
                                        if Timecol==0:      # time and date are in dateform
                                            try: tmp.append(todatenum(dt.datetime.strptime(datedata, dateform)))
                                            except: pass
                                        else:       # date is in dateform, time comes as a timenum
                                            try: tmp.append(todatenum(dt.datetime.strptime(datedata, dateform))+data[Timecol][i])
                                            except: pass
                                    elif Datecol!=0:        # date is an ordinal
                                        if Timecol==0:       # date is an ordinal but time is in dateform
                                            t0=dt.datetime.strptime(datedata, dateform); 
                                            t1=dt.time(t0.hour,t0.minute,t0.second); 
                                            tmp.append(data[Datecol][i]+todatenum(t1))
                                if req==1: trash=raw_input("Problems in the file. Cancel this operation and correct it.")
                            else:               # date is an ordinal and time is a timenum (len(Ptime)=0)
                                print("datenum and timenum should have been declared in a previous question. Tsk tsk tsk.")
                                tmp.append(data[Datecol][i]+data[Timecol][i])
                            Sdata=np.array(tmp)
                            del tmp, Ptime, datetype, Datecol, Timecol, DaTi, dateform
                        except: 
                            trash=raw_input("Seriously? I give up! Get better data! (press enter to continue)")
                            Sdata=np.ones(nupo)*NaN
    xltime=dt.datetime(1899,12,31,0,0,0).toordinal()
    if Sdata[0]>dt.datetime(1950,1,1,0,0,0).toordinal() and Sdata[0]<dt.datetime(2050,1,1,0,0,0).toordinal(): pass
    elif Sdata[0]+xltime>dt.datetime(1950,1,1,0,0,0).toordinal() and Sdata[0]+xltime<dt.datetime(2050,1,1,0,0,0).toordinal(): Sdata=Sdata+xltime-1
    else: print("Strange times...")
    Sttl=["time"]; Sunit=["ordinal"];
    for i in range(len(dttl)):
        print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
    Q=raw_input("Cleaning time columns: Enter all column numbers involving time: they are going to be deleted:  [enter]/column # separated by commas: ")
    try:
        DelCol=map(int,Q.split(","))
        DelCol.sort(); DelCol.reverse();
        for col in DelCol:
            trash=data.pop(col); trash=dttl.pop(col); trash=dunit.pop(col);
        del trash
    except: pass
    del Jttl, Junit
       
    
    ####### size distributions #######
    for i in range(len(dttl)): print("%d- %s (%s)" % (i, dttl[i],dunit[i]))
    while 1:
        sdnum_i=raw_input("How many different size distributions do you have?  ")      # number of size distributions
        try:
            sdnum=int(sdnum_i)
            break
        except:
            print("That wasn't a number, dear. Please try again.\n")
            pass
    sd=[];
    keys_list=[]; 
    for j in range(sdnum):
        fortrash=[];    # setting the trashbin empty
        for i in range(len(dttl)): print("%d- %s (%s)" % (i, dttl[i],dunit[i]))
        sdtmp=dict()
        print("Size distribution #%d!" % (j+1))
        sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
        ctim=raw_input("Enter the column number corresponding to the TIME if different from the default:  ")
        try: ctim=int(ctim);
        except: ctim=NaN; 
        ctot=raw_input("Enter the column number corresponding to the TOTAL PARTICLE concentration if any:  ")
        try: ctot=int(ctot)
        except: ctot=NaN
        logconc=raw_input("Are the concentrations in dN [1] format or dN/dlogDp [2]?  ")
        cmdp=raw_input("Enter the column number corresponding to the MEAN VOLUME DIAMETER of particles if any:  ");
        try: cmdp=int(cmdp)
        except: cmdp=NaN
        cmedp=raw_input("Enter the column number corresponding to the MEDIAN VOLUME DIAMETER of particles if any:  ")
        try: cmedp=int(cmedp)
        except: cmedp=NaN
        Q=raw_input("Enter the first and last column of the SIZE DISTribution (ex: 4,24)  ")
        try: 
            [c1, c2]=Q.split(',')
            c1=int(c1); c2=int(c2)
            szbin=dunit[c1:c2+1]
        except: 
            c1=NaN; c2=NaN; szbin=NaN;
        # size bins:
        Ins=InstrumentSpecs(sdtmp["Distname"])
        if raw_input("%s \nAre these the sizes? [y]/n  " % szbin)=='n': 
            Q='n'
            if raw_input("%s \nAre these the sizes? [y]/n  " % dttl[c1:c2+1])=='n': Q='n'
            else: szbin=dttl[c1:c2+1]; Q='y'
        else: Q='y'
        while Q=='n':
            if sdtmp["Distname"].lower()=='fssp300': szbin=np.array(map(float,Ins['range 1'].split(','))); Q='y'
            elif sdtmp["Distname"].lower()=='pcasp': szbin=np.array(map(float,Ins['range 1'].split(','))); Q='y'
            elif sdtmp["Distname"].lower()=='fssp124': 
                pos=[i for i,x in enumerate(dttl) if (x.lower() == 'r124' or x.lower() == 'r24')]
                try:
                    if len(pos)==1: pos=int(pos[0]); 
                    R=("range %d" % data[pos][1]);
                    szbin=np.array(map(float,Ins[R].split(','))); Q='y'
                except: 
                    print("No size dist was found in the instrument specification file for instrument %s." % sdtmp["Distname"])
                    for i in range(len(Ins)):
                        keys_list(Ins.keys()[i][-1])
                    for i in keys_list:
                        print("range "+str(i)+" = "+str(Ins['range '+str(i)]))
                    qr=raw_input("Are any of the above ranges the one you're looking for? [n]/#  ")
                    try:
                        R=("range "+str(qr))
                        szbin=np.array(map(float,Ins[R].split(',')))
                        Q='y'
                    except:
                        pass   
            elif sdtmp["Distname"].lower()=='fssp96': 
                pos=[i for i,x in enumerate(dttl) if (x.lower() == 'r96' or x.lower() == 'r6')]
                try:
                    if len(pos)==1: pos=int(pos[0]); 
                    R=("range %d" % data[pos][1]); 
                    szbin=np.array(map(float,Ins[R].split(','))); Q='y'
                except: print("No size dist was found in the instrument specification file for instrument %s." % sdtmp["Distname"])
            elif sdtmp["Distname"].lower()=='fssp100': 
                pos=[i for i,x in enumerate(dttl) if (x.lower() == 'r100' or x.lower() == 'r10')]
                try:
                    if len(pos)==1: pos=int(pos[0]); 
                    R=("range %d" % data[pos][1]);
                    szbin=np.array(map(float,Ins[R].split(','))); Q='y'
                except: 
                    print("No size dist was found in the instrument specification file for instrument %s." % sdtmp["Distname"])
                    for i in range(len(Ins)):
                        keys_list.append(Ins.keys()[i][-1])
                    for i in keys_list:
                        print("range "+str(i)+" = "+str(Ins['range '+str(i)]))
                    qr=raw_input("Are any of the above ranges the one you're looking for? [n]/#  ")
                    try:
                        R=("range "+str(qr))
                        szbin=np.array(map(float,Ins[R].split(',')))
                        Q='y'
                    except:
                        pass  
            if Q=='n':
                if raw_input("Do you want to assign automatic channel numbers?  y/[n]  ")=='y': szbin=range(1,(c2-c1+2)); Q='y'
                else: 
                    try: szbin=np.array(map(float,raw_input("Enter the sizes separated by commas on one line here: ").split(',')))
                    except: szbin=[]
                    if len(szbin)==c2+1-c1: Q='y'
                    else: print("You need to have %d channels. Try again." % (c2+1-c1)); Q='n'
        sdtmp["bins"]=np.array(szbin); 
        # concentrations of the size distribution
        if isnan(c1): sdtmp["data"]=np.ones(nupo)*NaN;
        else: 
            sdtmp["data"]=np.array(data[c1:c2+1]);
            sdtmp["data"][sdtmp["data"]<0]=NaN
        if logconc=="1": 
            sdtmp["data"]=dN2dNdlogDp(sdtmp); 
            sdtmp["data"][sdtmp["data"]<0]=NaN
            logconc="2" # now it is in dndlog format.
        elif logconc=="2": pass
        else: print "The concentrations have been assumed to be in dN/dlogDp format"
        # time for size distribution
        if isnan(ctim): sdtmp["time"]=Sdata;
        else: sdtmp["time"]=np.array(data[ctim]);
        # total concentration for size distribution
        if isnan(ctot) and logconc=="2": 
            try: 
                sd4tot=dict(); sd4tot["bins"]=copy.deepcopy(sdtmp["bins"]); sd4tot["data"]=copy.deepcopy(sdtmp["data"]); sd4tot["time"]=copy.deepcopy(sdtmp["time"]);
                from CloudDataSorter import dNdlogDp2N
                sdtmp["total"]=dNdlogDp2N(sd4tot,nan,nan)
                del sd4tot
            except: sdtmp["total"]=np.ones(nupo)*NaN;
        elif isnan(ctot) and logconc=="1":
            try: 
                sdtmp["total"]=np.nansum(sdtmp["data"],axis=0)
            except: sdtmp["total"]=np.ones(nupo)*NaN;      
        else: 
            try:
                sdtmp["total"]=np.array(data[ctot])
            except: sdtmp["total"]=np.ones(nupo)*NaN;
        # mean diameter of particles
        if isnan(cmdp):   # assume dN/dlogDp (should be at this stage)
            try: 
                ttmp=copy.deepcopy(sdtmp["total"])
                ttmp[ttmp<=0]=NaN
                sdtmp["mdp"]=np.nansum(dNdlogDp2dN(sdtmp).transpose()*sdtmp["bins"],axis=1)/ttmp
                sdtmp["mdp"][sdtmp["mdp"]<=0]=NaN
            except: sdtmp["mdp"]=np.ones(nupo)*NaN   
        else: 
            sdtmp["mdp"]=np.array(data[cmdp])         
            sdtmp["mdp"][sdtmp["mdp"]<=0]=NaN
#        # median diameter of particles
        if isnan(cmedp):
            try: 
                dN=dNdlogDp2dN(sdtmp)
                ttmp=np.nansum(dN,axis=0)
                ttmp[ttmp<=0]=NaN
                sdtmp["medp"]=np.ones(np.shape(sdtmp["total"]))*NaN
                for t in range(0,len(sdtmp["time"])):
                    if isnan(ttmp[t]): 
                        sdtmp["medp"][t]=NaN
                    else:
                        ntot=0; bin=-1
                        while ntot<ttmp[t]/2.0 and bin<=len(sdtmp["bins"])+1:
                            bin+=1; 
                            ntot=ntot+dN[bin,t]; 
                        sdtmp["medp"][t]=(sdtmp["bins"][bin]);
                sdtmp["medp"][sdtmp["medp"]<=0]=NaN
            except: sdtmp["medp"]=np.ones(nupo)*NaN   
        else: 
            sdtmp["medp"]=np.array(data[cmedp])      
            sdtmp["medp"][sdtmp["medp"]<=0]=NaN
        sdtmp["units"]=raw_input("What are the units of the size distribution? (e.g. cm-3)  ")
        sdtmp["sdtype"]=raw_input("What kind of size distribution is it? (A=aerosol, C=cloud droplet, P=precipitation)  ")
        sdtmp["sourcefile"]='basic'
        
        sd.append(sdtmp)
        # trashing the size distribution
        fortrash.append(ctot); fortrash.append(ctim); fortrash.append(cmdp); fortrash.append(cmedp);
        
        try: fortrash.append(pos); 
        except: pass
        if not isnan(c1): 
            for i in range(c1,c2+1):
                fortrash.append(i)
        try: 
            g=[i for i in fortrash if not isnan(i)]
            g.sort(); g.reverse();
            try: g.remove([]);
            except: pass
            for col in g:
                trash=data.pop(col); trash=dttl.pop(col); trash=dunit.pop(col);
            del trash
        except: print("This size distribution is empty!")
        # end of single size dist loop
    # access through for example: sd[0]["total"]

    ####### 2D data #######
    from CloudDataSorter import addsizedist2
    print("Importing 2D data (precipitation-sized droplets)")
    try: 
        sdtmp=addsizedist2()
        sd.append(sdtmp)
    except: print "Something went wrong with the 2D data. You can try adding it later."
    
    ####### Altitude #######
    Jttl=" ".join(dttl).lower(); Junit=" ".join(dunit).lower();
    if "alt" in Jttl:
        Q=0     # the right altitude has not yet been found
        for i in range(len(dttl)):
            if Q==0:
                if ("alt" in dttl[i].lower()) and (dunit[i].strip(" ")).lower()=='m':
                    Stmp=np.array(data.pop(i)); Q=1     # the right altitude has been found
                    trash=dttl.pop(i); trash=dunit.pop(i);
        if Q==0:
            i=0        # initiating
            while Q==0 and i<=(len(dttl)-1):
                if "alt" in dttl[i].lower() and Q==0:
                    if (dunit[i].strip(" ")).lower()=='km': cv=1000
                    elif (dunit[i].strip(" ")).lower()=='f' or (dunit[i].strip(" ")).lower()=='ft': cv=0.3048
                    elif (dunit[i].strip(" ")).lower()=='kf' or (dunit[i].strip(" ")).lower()=='kft': cv=1000*0.3048
                    else: cv=float(raw_input("Enter the factor of conversion for units %s:  " % dunit[i].strip(" ")))
                    Stmp=np.array(data.pop(i))*cv; Q=1;    # feet/kfeet converted in m (1f=0.3048 m)
                    trash=dttl.pop(i); trash=dunit.pop(i);
                    del cv
                i+=1
        del trash, Q, i
    else: 
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        i=raw_input("No altitude was found. Can you find it? [n]/number ")
        try:
            i=int(i)
            cv=raw_input("Enter the conversion factor to meters if necessary m=1, km=1000, f=0.3048, kf=304.8, etc.:  ")
            #Sdata=np.vstack((Sdata,np.array(data.pop(i))*cv)); 
            Stmp=np.array(data.pop(i))*cv; 
            trash=dttl.pop(i); trash=dunit.pop(i);
        except: 
            print("No altitude found.")    
            #Sdata=np.vstack((Sdata,np.ones(nupo)*NaN))
            Stmp=np.ones(nupo)*NaN
    Sttl.append("altitude"); Sunit.append("m");
    del Jttl, Junit
    M=runstats(Stmp,10)
    Stmp[Stmp==88888]=NaN
    Stmp[Stmp==-9999]=NaN
    Stmp=np.ma.masked_where((Stmp>(M[0]+M[1]*1.5)+(isnan(Stmp))), Stmp)
    Sdata=np.ma.vstack((Sdata,Stmp))
    

 
    ####### Latitude & Longitude #######
    # Latitude
    P=[i for i,x in enumerate(dttl) if 'lat' in x.lower()]
    if len(P)==1:
        P=int(P[0])
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
        trash=dttl.pop(P); trash=dunit.pop(P);
    else:
        if len(P)>1:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered latitude:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); trash=dunit.pop(i);
            except: P=[]
        if len(P)==0:
            for i in range(len(dttl)):
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("No latitude was found. Can you see it? [n]/number  ")
            try:
                P=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
                trash=dttl.pop(P); trash=dunit.pop(P);
            except: 
                print("No latitude found.")
                Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("latitude"); Sunit.append("deg");
    Sdata[-1][Sdata[-1]==-88.88888]=NaN
    Sdata[-1][Sdata[-1]>888.]=NaN
    Sdata[-1][Sdata[-1]==0.]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    # Longitude
    P=[i for i,x in enumerate(dttl) if 'lon' in x.lower()]
    if len(P)==1:
        P=int(P[0])
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
        trash=dttl.pop(P); trash=dunit.pop(P);
    else:
        if len(P)>1:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered longitude:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); trash=dunit.pop(i);
            except: P=[]
        if len(P)==0:
            for i in range(len(dttl)):
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("No longitude was found. Can you see it? [n]/number  ")
            try:
                P=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
                trash=dttl.pop(P); trash=dunit.pop(P);
            except: 
                print("No longitude found.")
                Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("longitude"); Sunit.append("deg");
    Sdata[-1][Sdata[-1]==-88.88888]=NaN
    Sdata[-1][Sdata[-1]>888.]=NaN
    Sdata[-1][Sdata[-1]==0.]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
        
    ####### Potential temperatures #######
    # Theta-D
    P=[i for i,x in enumerate(dttl) if ('the' in x.lower() and 'd' in x.lower() and 'k' in dunit[i].lower())]
    if len(P)==1:
        try:
            i=int(P[0]); 
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i))))
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    elif len(P)>1:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the pefered Theta-D:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No Theta-D was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); thunit=dunit.pop(P);
        except: 
            print("No Theta-D found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
            thunit='k'
    if 'k' in thunit: pass
    elif 'c' in thunit: Sdata[-1]=Sdata[-1]+273.15
    Sttl.append("theta-d"); Sunit.append("kelvin");
    Sdata[-1][Sdata[-1]==888.8]=NaN
    # Theta-Q
    P=[i for i,x in enumerate(dttl) if ('the' in x.lower() and 'q' in x.lower() and 'k' in dunit[i].lower())]
    if len(P)==1:
        try:
            i=int(P[0]); 
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i))))
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    elif len(P)>1:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the pefered Theta-Q:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No Theta-Q was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); thunit=dunit.pop(P);
        except: 
            print("No Theta-Q found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
            thunit='k'
    if 'k' in thunit: pass
    elif 'c' in thunit: Sdata[-1]=Sdata[-1]+273.15
    Sttl.append("theta-q"); Sunit.append("kelvin");
    Sdata[-1][Sdata[-1]==888.8]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    # Theta-E
    P=[i for i,x in enumerate(dttl) if ('the' in x.lower() and 'e' in x.lower() and 'k' in dunit[i].lower())]
    if len(P)==1:
        try:
            i=int(P[0]); 
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i))))
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    elif len(P)>1:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the pefered Theta-E:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); thunit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No Theta-E was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); thunit=dunit.pop(P);
        except: 
            print("No Theta-E found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
            thunit='k'
    if 'k' in thunit: pass
    elif 'c' in thunit: Sdata[-1]=Sdata[-1]+273.15
    Sttl.append("theta-e"); Sunit.append("kelvin");
    Sdata[-1][Sdata[-1]==888.8]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    ####### Temperature #######
    P=[i for i,x in enumerate(dttl) if ('t' in x.lower() and ('c' in dunit[i].lower() or 'k' in dunit[i].lower()))]
    if len(P)>0:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered temperature:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); trash=dunit.pop(i);
            except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No temperature was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); trash=dunit.pop(P);
        except: 
            print("No temperature found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("temperature"); Sunit.append("celsius");
    Sdata[-1][Sdata[-1]==-888.8]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
    ####### Relative Humidity #######
    P=[i for i,x in enumerate(dttl) if ('r' in x.lower() and 'h' in x.lower() and "%" in dunit[i].lower())]
    if len(P)==1:
        P=int(P[0])
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
        trash=dttl.pop(P); trash=dunit.pop(P);
    else:
        if len(P)>1:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered relative humidity:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); trash=dunit.pop(i);
            except: P=[]
        if len(P)==0:
            for i in range(len(dttl)):
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("No RH was found. Can you see it? [n]/number  ")
            try:
                P=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
                trash=dttl.pop(P); trash=dunit.pop(P);
            except: 
                print("No RH found.")
                Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("RH"); Sunit.append("%");
    Sdata[-1][Sdata[-1]==888.8]=NaN
    Sdata[-1][Sdata[-1]>500]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
    ####### Pressure #######
    P=[i for i,x in enumerate(dttl) if ('p' in x.lower() and ('b' in dunit[i].lower() or 'atm' in dunit[i].lower() or 'pa' in dunit[i].lower()))]
    if len(P)==1:
        P=int(P[0])
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
        trash=dttl.pop(P); Punit=dunit.pop(P);
    else:
        if len(P)>1:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered pressure:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); Punit=dunit.pop(i);
            except: P=[]
        if len(P)==0:
            for i in range(len(dttl)):
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("No pressure was found. Can you see it? [n]/number  ")
            try:
                P=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
                trash=dttl.pop(P); Punit=dunit.pop(P);
            except: 
                print("No pressure found.")
                Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
                Punit="empty"
    Punit=Punit.strip(" ")
    Punit=Punit.lower()
    if Punit=='mb' or Punit=='mbar' or Punit=="empty": cv=1
    elif Punit=='b' or Punit=='bar': cv=1000
    elif Punit=='atm': cv=1013.25
    elif Punit=='pa': cv=0.01
    elif Punit=='kpa': cv=10
    Sdata[-1]=Sdata[-1]*cv
    Sttl.append("Pressure"); Sunit.append("mb");
    Sdata[-1][Sdata[-1]<0]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])

    ####### liquid water content #######
    P=[i for i,x in enumerate(dttl) if ('lw' in x.lower() or 'kszr' in x.lower() or 'klwc' in x.lower() and 'g/m3' in dunit[i].lower())]
    if len(P)==1:
        P=int(P[0])
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
        trash=dttl.pop(P); trash=dunit.pop(P);
    else:
        if len(P)>1:
            for i in P:
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("Choose which of these is the prefered liquid water content:  number/[n]  ")
            try:
                i=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
                trash=dttl.pop(i); trash=dunit.pop(i);
            except: P=[]
        if len(P)==0:
            for i in range(len(dttl)):
                print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
            P=raw_input("No LWC was found. Can you see it? [n]/number  ")
            try:
                P=int(P)
                Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
                trash=dttl.pop(P); trash=dunit.pop(P);
            except: 
                print("No LWC found.")
                Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("LWC"); Sunit.append("g/m3");
    Sdata[-1][Sdata[-1]==-8.888]=NaN
    Sdata[-1][Sdata[-1]==8.888]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])

    
    ####### surface temperature #######
    P=[i for i,x in enumerate(dttl) if ('t' in x.lower() and ('c' in dunit[i].lower() or 'k' in dunit[i].lower()))]
    if len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered surface temperature:  number/[n]  ")
        try:
            i=int(P)
            if 'k' in dunit[i].lower(): cv=273.15
            else: cv=0
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i))-cv)); 
            trash=dttl.pop(i); trash=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No surface temperature was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            if 'k' in dunit[i].lower(): cv=273.15
            else: cv=0
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P))-cv)); 
            trash=dttl.pop(P); trash=dunit.pop(P);
        except: 
            print("No surface temperature found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("surface temp"); Sunit.append("celsius");
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
#    ####### median volume diameter #######
#    P=[i for i,x in enumerate(dttl) if ('mevd' in x.lower() and 'm' in dunit[i].lower())]
#    if len(P)==1:
#        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P[0])))); 
#        trash=dttl.pop(P[0]); Punit=dunit.pop(P[0]);
#    elif len(P)>0:
#        for i in P:
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("Choose which of these is the prefered median volume diameter:  number/[n]  ")
#        try:
#            i=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
#            trash=dttl.pop(i); Punit=dunit.pop(i);
#        except: P=[]
#    if len(P)==0:
#        for i in range(len(dttl)):
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("No median volume diameter was found. Can you see it? [n]/number  ")
#        try:
#            P=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
#            trash=dttl.pop(P); Punit=dunit.pop(P);
#        except: 
#            print("No median volume diameter found.")
#            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
#            Punit='um'
#    Punit=Punit.strip(" ")
#    Punit=Punit.lower()
#    if Punit=='um': cv=1
#    else: cv=float(raw_input("Units are %s. Enter the correction factor to get um:  " % Punit))
#    Sdata[-1]=Sdata[-1]*cv
#    Sttl.append("medvd"); Sunit.append("um");
#    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
#    ####### median volume diameter #######
#    P=[i for i,x in enumerate(dttl) if ('mevd' in x.lower() and 'm' in dunit[i].lower())]
#    if len(P)==1:
#        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P[0])))); 
#        trash=dttl.pop(P[0]); Punit=dunit.pop(P[0]);
#    elif len(P)>0:
#        for i in P:
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("Choose which of these is the prefered median volume diameter:  number/[n]  ")
#        try:
#            i=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
#            trash=dttl.pop(i); Punit=dunit.pop(i);
#        except: P=[]
#    if len(P)==0:
#        for i in range(len(dttl)):
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("No median volume diameter was found. Can you see it? [n]/number  ")
#        try:
#            P=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
#            trash=dttl.pop(P); Punit=dunit.pop(P);
#        except: 
#            print("No median volume diameter found.")
#            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
#            Punit='um'
#    Punit=Punit.strip(" ")
#    Punit=Punit.lower()
#    if Punit=='um': cv=1
#    else: cv=float(raw_input("Units are %s. Enter the correction factor to get um:  " % Punit))
#    Sdata[-1]=Sdata[-1]*cv
#    Sttl.append("medvd"); Sunit.append("um");
#    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
    ####### Mixing ratio air-water #######
    P=[i for i,x in enumerate(dttl) if ('mixra' in x.lower() and 'g/kg' in dunit[i].lower())]
    if len(P)==1:
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P[0])))); 
        trash=dttl.pop(P[0]); trash=dunit.pop(P[0]);
    elif len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered mixing ratio air-water:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); trash=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No mixing ratio air-water was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); trash=dunit.pop(P);
        except: 
            print("No mixing ratio air-water found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("mixrair"); Sunit.append("g/kg");
    Sdata[-1][Sdata[-1]==88.888]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    
#    ####### Mixing ratio of liquid #######
#    P=[i for i,x in enumerate(dttl) if ('mixrl' in x.lower() and 'g/kg' in dunit[i].lower())]
#    if len(P)==1:
#        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P[0])))); 
#        trash=dttl.pop(P[0]); trash=dunit.pop(P[0]);
#    elif len(P)>0:
#        for i in P:
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("Choose which of these is the prefered mixing ratio of liquid:  number/[n]  ")
#        try:
#            i=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
#            trash=dttl.pop(i); trash=dunit.pop(i);
#        except: P=[]
#    if len(P)==0:
#        for i in range(len(dttl)):
#            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
#        P=raw_input("No mixing ratio of liquid was found. Can you see it? [n]/number  ")
#        try:
#            P=int(P)
#            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
#            trash=dttl.pop(P); trash=dunit.pop(P);
#        except: 
#            print("No mixing ratio of liquid found.")
#            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
#    Sttl.append("mixrliquid"); Sunit.append("g/kg"); 
#    Sdata[-1][Sdata[-1]==88.888]=NaN   
#    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    

    ####### Total mixing ratio #######
    P=[i for i,x in enumerate(dttl) if ('mixrt' in x.lower() and 'g/kg' in dunit[i].lower())]
    if len(P)==1:
        Sdata=np.ma.vstack((Sdata,np.array(data.pop(P[0])))); 
        trash=dttl.pop(P[0]); trash=dunit.pop(P[0]);
    elif len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered total mixing ratio:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); trash=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No total mixing ratio was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); trash=dunit.pop(P);
        except: 
            print("No total mixing ratio found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("mixrtotal"); Sunit.append("g/kg"); 
    Sdata[-1][Sdata[-1]==88.888]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    

    ####### Dew Point #######
    P=[i for i,x in enumerate(dttl) if ('d' in x.lower() and ('c' in dunit[i].lower() or 'k' in dunit[i].lower()))]
    if len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered dew point:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); Punit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No dew point was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); Punit=dunit.pop(P);
        except: 
            print("No dew point found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Punit=Punit.strip(" ")
    Punit=Punit.lower()
    if 'c' in Punit: cv=0
    elif 'k' in Punit: cv=273.15
    Sdata[-1]=Sdata[-1]+cv
    Sttl.append("dewpt"); Sunit.append("celsius");
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    del Punit


    ####### Wind speed #######
    P=[i for i,x in enumerate(dttl) if (('w' in x.lower() and 's' in x.lower()) and ('m' in dunit[i].lower() or 'k' in dunit[i].lower()))]
    if len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered wind speed:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); Punit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No wind speed was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); Punit=dunit.pop(P);
        except: 
            print("No wind speed found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
            Punit='empty'
    Punit=Punit.strip(" ")
    Punit=Punit.lower()
    if Punit=='m/s' or Punit=='m/sec' or Punit=='ms-1' or Punit=='m s-1' or Punit=='msec-1' or Punit=='m sec-1': cv=1
    elif Punit=='km/s' or Punit=='km/sec' or Punit=='kms-1' or Punit=='km s-1' or Punit=='kmsec-1' or Punit=='km sec-1': cv=1000
    elif Punit=='knots' or Punit=='kt' or ('k' in Punit and 't' in Punit): cv=0.514444444
    elif Punit=='km/h' or Punit=='kmh-1' or Punit=='km h-1': cv=0.277777778
    elif Punit=='empty': cv=0
    else: cv=float(raw_input("Units are %s. Enter the correction factor to get m/s:  " % Punit))
    Sdata[-1]=Sdata[-1]*cv
    Sttl.append("windspeed"); Sunit.append("m/s");
    Sdata[-1][Sdata[-1]==888.8]=NaN
    Sdata[-1][Sdata[-1]==-99.99]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    del Punit
    
    ####### Wind direction #######
    P=[i for i,x in enumerate(dttl) if ('w' in x.lower() and 'd' in x.lower())]
    if len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered wind direction:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); trash=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No wind direction was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); trash=dunit.pop(P);
        except: 
            print("No wind direction found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
    Sttl.append("winddir"); Sunit.append("deg");
    Sdata[-1][abs(Sdata[-1])>360]=NaN
    Sdata[-1][Sdata[-1]==-99.99]=NaN
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    
    ####### updraft velocity (gust) #######
    P=[i for i,x in enumerate(dttl) if ('w' in x.lower() and ('m' in dunit[i].lower() or 'k' in dunit[i].lower()) and ('s' in dunit[i].lower() or 'h' in dunit[i].lower()))]
    if len(P)>0:
        for i in P:
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("Choose which of these is the prefered updraft velocity:  number/[n]  ")
        try:
            i=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(i)))); 
            trash=dttl.pop(i); Punit=dunit.pop(i);
        except: P=[]
    if len(P)==0:
        for i in range(len(dttl)):
            print("%d- %s (%s)" % (i, dttl[i], dunit[i]))
        P=raw_input("No updraft velocity was found. Can you see it? [n]/number  ")
        try:
            P=int(P)
            Sdata=np.ma.vstack((Sdata,np.array(data.pop(P)))); 
            trash=dttl.pop(P); Punit=dunit.pop(P);
        except: 
            print("No updraft velocity found.")
            Sdata=np.ma.vstack((Sdata,np.ones(nupo)*NaN))
            Punit="empty"
    Sttl.append("udvel"); Sunit.append("m/s");
    Punit=Punit.strip(" ")
    Punit=Punit.lower()
    if Punit=='m/s' or Punit=='m/sec' or Punit=='ms-1' or Punit=='m s-1' or Punit=='msec-1' or Punit=='m sec-1' or Punit=="empty": cv=1
    elif Punit=='km/s' or Punit=='km/sec' or Punit=='kms-1' or Punit=='km s-1' or Punit=='kmsec-1' or Punit=='km sec-1': cv=1000
    elif Punit=='knots' or Punit=='kt' or ('k' in Punit and 't' in Punit): cv=0.514444444
    elif Punit=='km/h' or Punit=='kmh-1' or Punit=='km h-1': cv=0.277777778
    else: cv=float(raw_input("Units are %s. Enter the correction factor to get m/s:  " % Punit))
    Sdata[-1][Sdata[-1]==8888.8]=NaN
    Sdata[-1][Sdata[-1]==888.8]=NaN
    Sdata[-1][Sdata[-1]==9999.9]=NaN
    Sdata[-1]=Sdata[-1]*cv
    Sdata[-1]=np.ma.masked_where(isnan(Sdata[-1]),Sdata[-1])
    del Punit


    
    #################################################################
    ###  Add leftover original data at the end of the sorted data ###
    #################################################################
    Sdata=np.ma.vstack((Sdata,np.array(data)))
    for i in range(len(dttl)):
        Sttl.append(dttl[i])
        Sunit.append(dunit[i])
    return [Sdata, Sttl, Sunit, sd]     # return sorted data with legend


def addsizedist2():
# loads a size dist in a text file (not excel)
# returns a size distribution
    from CloudDataSorter import todatenum
    import Tkinter, tkFileDialog
    dirnames=[]; pathnames=[]; fnames=[];
    print("Please choose the directory using the following dialogue")
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
    except: 
        print("Failed! Maybe there are no files corresponding to your request?")
        smpsfn='no file'
    if smpsfn=='no file':
        sdtmp=dict()
        sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
        sdtmp["total"]=NaN; sdtmp["time"]=NaN; sdtmp["sourcefile"]='No file found'; sdtmp["bins"]=NaN; sdtmp["data"]=NaN; sdtmp["mdp"]=NaN; sdtmp["medp"]=NaN;
    else:
        if smpsfn[-3:]=='xls': print("This is an excel file. Please use addsizedist1")
        else: inp = open(smpsfn,'r') # open file
        sds=[]
        for line in inp.readlines():
            ltmp=line.strip("\n"); ltmp=ltmp.split('\t');
            sds.append(ltmp)
        if np.shape(sds)[0]==1 or np.shape(sds)[1]==1:      # if the delimiter was not a tab
            sds=[]; inp = open(smpsfn,'r') # open file
            for line in inp.readlines():
                ltmp=line.strip("\n"); ltmp=ltmp.split('  ');       # trying two spaces as delimiter
                sds.append(ltmp)
            
        if np.shape(sds)[0]==1 or np.shape(sds)[1]==1 or np.shape(sds)[0]==0 or np.shape(sds)[1]==0:
            sdtmp=dict()
            sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
            sdtmp["total"]=NaN; sdtmp["time"]=NaN; sdtmp["sourcefile"]='No file found'; sdtmp["bins"]=NaN; sdtmp["data"]=NaN; sdtmp["mdp"]=NaN; sdtmp["medp"]=NaN;
        else:
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
                sdd=vstack((sdd,map(float,np.array(sds[i]))))
            nupo=np.shape(sdd)[0]; print np.shape(sdd)
            
            # size bins
            sdtmp=dict()
            for i,x in enumerate(sddttl): print ("%d- %s (%s)" % (i, x, sddunt[i]))
            Q=raw_input("Enter the first and last column of the SIZE BINS (ex: 4,24)  ")
            try: 
                [c1, c2]=Q.split(',')
                c1=int(c1); c2=int(c2)
                sdtmp["bins"]=np.array(map(float,sddunt[c1:c2+1]));
            except: sdtmp["bins"]=NaN; c1=NaN; c2=NaN;
            # concentrations of the size distribution
            if isnan(c1): sdtmp["data"]=np.ones([nupo,1])*NaN;
            elif len(np.shape(sdd))==1: 
                print("This file has only one line, we will consider it empty.")
                sdtmp["data"]=np.ones([nupo,(c2+1-c1)])*NaN;
            else: sdtmp["data"]=np.array(sdd[:,c1:c2+1]);
            # time for size distribution
            Q=raw_input("Which is(are) the time column(s)?  ")
            # MMM I didn't feel like wasting a lot of time making the time retrieval flexible, might wanna do that later.
            try: 
                Q=map(int,Q.split(','))
                if len(Q)==1: Q=Q[0]
                QQ=raw_input("Is %s an ordinal format? y/[n]  " % sdd[1,Q])
                if QQ=='y':
                    sdtmp["time"]=np.array(sdd[:,Q]);
                    if any(sdtmp["time"]>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"]<dt.datetime(2050,1,1,0,0,0).toordinal()): pass
                    elif any(sdtmp["time"][0]+xltime>dt.datetime(1950,1,1,0,0,0).toordinal()) and any(sdtmp["time"][0]+xltime<dt.datetime(2050,1,1,0,0,0).toordinal()): 
                        sdtmp["time"]=sdtmp["time"]+xltime
                    else: print("Strange ordinal...")
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
            logconc=raw_input("Are the concentrations in dN [1] format or dN/dlogDp [2]?  ")
            if logconc=="1": 
                from CloudDataSorter import dN2dNdlogDp
                sdtmp["data"]=dN2dNdlogDp(sdtmp)
            elif logconc=="2": pass
            else: print "The concentrations have been assumed to be in dN/dlogDp format"
            
            # total concentration for size distribution
            Q=raw_input("Enter the column with the total concentration:  ")
            try:
                Q=int(Q)
                sdtmp["total"]=np.array(sdd[Q])
            except: 
                try: 
                    sd4tot=dict(); sd4tot["bins"]=sdtmp["bins"]; sd4tot["data"]=sdtmp["data"].transpose(); sd4tot["time"]=sdtmp["time"];
                    from CloudDataSorter import dNdlogDp2N
                    sdtmp["total"]=dNdlogDp2N(sd4tot,nan,nan)
                    del sd4tot
                except: sdtmp["total"]=np.ones(nupo)*NaN;
                
            # mean diameter of particles
            cmdp=raw_input("Enter the column with the median volumic diameter:  ")
            try: cmdp=int(cmdp)
            except: cmdp=NaN
            if isnan(cmdp):   # assume dN/dlogDp (should be at this stage)
                try: 
                    ttmp=copy.deepcopy(sdtmp["total"])
                    ttmp[ttmp<=0]=NaN
                    sdtmp["mdp"]=np.nansum(dNdlogDp2dN(sdtmp).transpose()*sdtmp["bins"],axis=1)/ttmp
                    sdtmp["mdp"][sdtmp["mdp"]<=0]=NaN
                except: sdtmp["mdp"]=np.ones(nupo)*NaN   
            else: 
                sdtmp["mdp"]=np.array(data[cmdp])         
                sdtmp["mdp"][sdtmp["mdp"]<=0]=NaN
        
            # median diameter of particles
            cmedp=raw_input("Enter the column with the mean volumic diameter:  ")
            try: cmedp=int(cmedp)
            except: cmedp=NaN
            if isnan(cmedp):
                try: 
                    dN=dNdlogDp2dN(sdtmp)
                    ttmp=np.nansum(dN,axis=0)
                    ttmp[ttmp<=0]=NaN
                    sdtmp["medp"]=np.ones(np.shape(sdtmp["total"]))*NaN
                    for t in range(0,len(sdtmp["time"])):
                        if isnan(ttmp[t]): 
                            sdtmp["medp"][t]=NaN
                        else:
                            ntot=0; bin=-1
                            while ntot<ttmp[t]/2.0 and bin<=len(sdtmp["bins"])+1:
                                bin+=1; 
                                ntot=ntot+dN[bin,t]; 
                            sdtmp["medp"][t]=(sdtmp["bins"][bin]);
                            sdtmp["medp"][sdtmp["medp"]<=0]=NaN  
                except: sdtmp["medp"]=np.ones(nupo)*NaN   
            else: 
                sdtmp["medp"]=np.array(data[cmedp])
                sdtmp["medp"][sdtmp["medp"]<=0]=NaN      
                
            # Distribution name
            sdtmp["Distname"]=raw_input("What is the name of the size distribution? (instrument, etc.)  ")
            sdtmp["units"]=raw_input("What are the units of the size distribution? (e.g. cm-3)  ")
            sdtmp["sdtype"]=raw_input("What kind of size distribution is it? (A=aerosol, C=cloud droplet, P=precipitation)  ")
            # appending distribution to Cloud object.
            sdtmp["total"]=sdtmp["total"]
            if np.shape(sdtmp["data"])[0]==np.shape(sdtmp["time"])[0]:
                sdtmp["data"]=sdtmp["data"].transpose()
            elif np.shape(sdtmp["data"])[1]==np.shape(sdtmp["time"])[0]:
                pass
            else: print("2dc data failed. The data do not correspond to the time length.")
            sdtmp["time"]=sdtmp["time"]
            sdtmp["sourcefile"]=smpsfn
            
    return sdtmp


def todatenum(d):
    # converts to ordinal number (days since jan 1 of year 1)
    # d is a datetime.datetime or a datetime.time object
    try: d1=d.toordinal()
    except: d1=0
    d2=(d.hour+(d.minute+d.second/60.)/60.)/24.
    dn=d1+d2
    return dn
    

def dNdlogDp2N(dist,z1,z2):
    # this function return the total particle concentration when the distribution is give in dN/dlogDp
    # dN is the concentrations (in DN/dlogZ) and Z is the diameter vector
    # z1 and z2 are the boundaries within which we sum (inclusively)
    # this returns total concentration
    # SG, Feb 2008

    # modified to take in matrices dN --> lines: times columns: diameters
    # modified to make total concentrations between given diameters (z1 
    # (smaller diameter) and z2 (bigger diameter)). (z1=z2=nan if no boundaries)
    # also removes negative total concentrations and replaces it with NaN.
    if type(dist["time"])==float: return NaN
    if len(np.shape(dist["data"]))<2 and np.shape(dist["time"])[0]>1:
        N=np.ones(np.shape(dist["time"]))*NaN; 
        #print("Not enough size bins (<2).")
    else:
        if np.shape(dist["time"])[0]==1: dist["data"]=dist["data"].reshape(np.shape(dist["data"])[0],1);
        Z=copy.deepcopy(dist["bins"])
        if len(Z)==np.shape(dist["data"])[1]:
            dN=copy.deepcopy(dist["data"])
        elif len(Z)==np.shape(dist["data"])[0]:
            dN=copy.deepcopy(dist["data"].transpose())
        if isnan(z2) and isnan(z1): ix=range(len(Z));
        elif isnan(z2): ix=np.nonzero((Z>=z1))[0]
        elif isnan(z1): ix=np.nonzero((Z<=z2))[0]
        else: ix=np.nonzero((Z<=z2)*(Z>=z1))[0]
        
        t = np.shape(Z[ix])[0];
        b=np.zeros([t,1]);
        h=zeros([t,1]);
    
        # building Dmobility (dlZ) vector
        for i in range(1,t):
            b[i]=(np.log10(Z[ix[i-1]])+np.log10(Z[ix[i]]))/2.;
        for i in range(0,t-1):
            h[i]=(np.log10(Z[ix[i+1]])+np.log10(Z[ix[i]]))/2.;
        b[0]=np.log10(Z[ix[0]])+(log10(Z[ix[1]])-h[1]);
        h[t-1]=np.log10(Z[ix[t-1]])-(b[t-2]-np.log10(Z[ix[t-2]]));
    
        dlZ=h-b;
        
        N=NaN*ones([np.shape(dN)[0],1]);
        for i in range(0,np.shape(dN)[0]):
            N[i] = sum(dN[i,ix]*dlZ.transpose());       # we don't want nansum here in case the nan is not due to a negative concentration but rather to measurement problems
        N[N<0]=NaN;
        N=N.reshape(len(N),)
    return N

def dNdlogDp2dN(dist):
    # this function return the particle concentration when the distribution is given in dN/dlogDp
    # variables dN is the concentrations (in DN/dlogZ) and Z is the diameter vector
    # this returns the concentration for each size bin in dN/dDp format
    # SG, March 2012

    # also removes negative total concentrations and replaces it with NaN.
    if len(np.shape(dist["data"]))<2:
        N=N=np.ones(np.shape(dist["time"]))*NaN
    else:
        dN=copy.deepcopy(dist["data"].transpose())
        Z=copy.deepcopy(dist["bins"])
        t = np.shape(Z)[0]; # number of size bins
        b=np.zeros([t,1]);
        h=zeros([t,1]);
    
        # building Dmobility (dlZ) vector
        for i in range(1,t):
            b[i]=(np.log10(Z[i-1])+np.log10(Z[i]))/2.;
        for i in range(0,t-1):
            h[i]=(np.log10(Z[i+1])+np.log10(Z[i]))/2.;
        b[0]=np.log10(Z[0])+(log10(Z[1])-h[1]);
        h[t-1]=np.log10(Z[t-1])-(b[t-2]-np.log10(Z[t-2]));
    
        dlZ=h-b;
        
        N=NaN*ones([np.shape(dN)[0],t]);
        for i in range(0,np.shape(dN)[0]):
            N[i,:] = dN[i,:]*dlZ.transpose();
    
        N[N<0]=NaN;
        N=N.transpose()
        N=N.reshape(t,-1)
    return N
    
def dN2dNdlogDp(dist):
    # this function return the particle size dist in dN/dlogDp when the distribution (dist) is given in dN
    # variables dN is the concentrations (in dN) and Z is the diameter vector
    # this returns the concentration for each size bin in dN/dlogDp format
    # SG, March 2012

    # also removes negative total concentrations and replaces it with NaN.
    if len(np.shape(dist["data"]))<2:
        dN=np.ones(np.shape(dist["time"]))*NaN
    else:
        Z=copy.deepcopy(dist["bins"])
        if len(Z)==np.shape(dist["data"])[1]:
            dN=copy.deepcopy(dist["data"])
        elif len(Z)==np.shape(dist["data"])[0]:
            dN=copy.deepcopy(dist["data"].transpose())
                    
        t = np.shape(Z)[0]; # number of size bins
        b=np.zeros([t,1]);
        h=zeros([t,1]);
    
        # building Dmobility (dlZ) vector
        for i in range(1,t):
            b[i]=(np.log10(Z[i-1])+np.log10(Z[i]))/2.;
        for i in range(0,t-1):
            h[i]=(np.log10(Z[i+1])+np.log10(Z[i]))/2.;
        b[0]=np.log10(Z[0])+(log10(Z[1])-h[1]);
        h[t-1]=np.log10(Z[t-1])-(b[t-2]-np.log10(Z[t-2]));
    
        dlZ=h-b;
            
        NdlogDp=NaN*ones([np.shape(dN)[0],t]);
        for i in range(0,np.shape(dN)[0]):
            NdlogDp[i,:] = dN[i,:]/dlZ.transpose();
    
        NdlogDp[NdlogDp<0]=NaN;
        NdlogDp=NdlogDp.reshape(-1,t)
    return NdlogDp.transpose()
    
    
def dNdDp2dNdlogDp(dist):
    # this function converts a size distribution in dN/dDp format into a dN/dlogDp format
    # Z is a diameter
    if len(np.shape(dist["data"]))<2:
        dN=np.ones(np.shape(dist["time"]))*NaN
    else:
        Z=copy.deepcopy(dist["bins"])
        if len(Z)==np.shape(dist["data"])[1]:
            dN=copy.deepcopy(dist["data"])
        elif len(Z)==np.shape(dist["data"])[0]:
            dN=copy.deepcopy(dist["data"].transpose())
                    
        t = np.shape(Z)[0]; # number of size bins
        b=np.zeros([t,1]);
        h=zeros([t,1]);
    
        # Dividing bins linearly
        for i in range(1,t):
            b[i]=(Z[i-1]+Z[i])/2.;
        for i in range(0,t-1):
            h[i]=(Z[i+1]+Z[i])/2.;
        b[0]=Z[0]+(Z[1]-h[1]);
        h[t-1]=Z[t-1]-(b[t-2]-Z[t-2]);
        
        dZ=h-b;     # dDp
        dlZ=np.log10(h)-np.log10(b)     # dlogDp
                    
        NdlogDp=NaN*ones([np.shape(dN)[0],t]);
        for i in range(0,np.shape(dN)[0]):
            NdlogDp[i,:] = dN[i,:]*dZ.transpose()/dlZ.transpose();
    
        NdlogDp[NdlogDp<0]=NaN;
        NdlogDp=NdlogDp.reshape(-1,t)
    return NdlogDp.transpose()    


def runstats(x,n):
# Stephanie Gagne, UHel, 2010
# converted to Python, Dal, 2012
# x is an array of 1 dimension.
# n is the number of point taken in the running statistic
    """takes data, number of points for the running mean/standard deviation and returns the running mean and running standard deviation."""
    try: x.mask
    except: 
        x=ma.asanyarray(x); 
        x.mask=ones(np.shape(x))*False
    try: [ro,co]=np.shape(x)
    except: ro=np.shape(x)[0]; co=1
    if ro==1 or co==1: x=x.reshape(max(ro,co),)
    else: print("The array must be a vector (one column or row)")
    # initializing matrix
    ro=max(ro,co)
    M=ones([ro,n])*NaN;
    M=ma.asanyarray(M)
    
    # building matrix
    if n%2==1:       # if n is odd
        for j in range(int(n/2),0,-1):
            posi=int(n/2)-j       # current position
            M[0:ro-j,posi]=x[j:]
        for j in range(1,2+int(n/2),1):
            posi=int(n/2)+j-1;
            M[j-1:,posi]=x[0:(ro+1)-j]
    elif n%2==0:        # if n is even
        for j in range(n/2,0,-1):
            posi=n/2-j
            M[0:ro-j,posi]=x[j:]
        for j in range(1,n/2+1):
            posi=n/2+j-1;
            M[j-1:,posi]=x[0:(ro+1)-j]
    else: print("Well, that's pretty weird. Are you sure n is an integer?")  
    
    M.data[M.mask]=NaN
    ave=st.nanmean(M, axis=1);
    stde=st.nanstd(M, axis=1);
    return [ave, stde]      

def AdiabLWC(Temp,Pres):
    """Adiabatic lwc calculated from the temperature in Celcius and pressure in mb or hPa.
       Returns 1-Pseudoadiabatic (moist adiabatic) lapse rate [K/m], 2-height dependence of LWC (dLWCdz)."""
    # Code by Mike Earle (Env. Can.), based on Tietze et al., ACP, 2011. Adapted to Python by Stephanie Gagne
    # Constants
    epsilon = 0.622;
    R_air = 287;  # Individual gas constant for dry air in J/(kg K) 
    # Saturation vapour pressure of water [hPa]
    #p_sat = 6.1121*exp((18.678 - Temp/234.50) * Temp/(275.14 + Temp));
    p_sat = 6.112 * np.exp(17.67 * Temp / (Temp + 243.5));
    #p_sat = 2.22;
    # Saturation mixing ratio 
    w_sat = epsilon * p_sat/(Pres - p_sat);
    # Heat capacity of dry air [J/(kg K)]
    c_p = 1005;
    # Air density at a given T, P [kg/m^3]
    rho_air = (Pres*100)/(287.05 * (Temp + 273.15));
    # Heat capacity of moist air [J/(kg K)]
    c_pm = 1870;
    #c_pm = c_p * (1 + 0.9 * w_sat);
    # Heat capacity of pure water [IT cal/(g C)] 
    c_w_params = [1.000938, -2.7052e-3, -2.3235e-5, 4.3778e-6, 2.7136e-7];
    c_w = 0;
    for j in range(5):
      c_w_i = c_w_params[j] * Temp**j;
      c_w = c_w + c_w_i;
    # Convert c_w to units of J/(kg K)
    c_w = c_w * 4.1868 * 1000;
    # Calculate latent heat of vaporization for pure water
    T_0 = 0;   
    L_0 = 2.50e6;   # Latent heat of vap at T_0 in J/kg
    L = L_0 - (c_w - c_pm)*(Temp - T_0); 
    # Dry adiabatic lapse rate [K/m]
    lr_dry = 9.76e-3;
    # Pseudoadiabatic (moist adiabatic) lapse rate [K/m]
    lr_moist = lr_dry * ((1 + L * w_sat/(R_air * (Temp + 273.15)))/(1 + L**2 * epsilon * w_sat/(R_air * c_p * (Temp + 273.15)**2)));
    # Compute height dependence of LWC [g/m^4]
    dLWCdz = rho_air * c_p / L * (lr_dry - lr_moist) * 1e3;
    return [lr_moist,dLWCdz]


