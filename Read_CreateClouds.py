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


#test
# This code reads a file

from pylab import *
import numpy as np
import datetime as dt
import xlrd as xlrd
import matplotlib.pyplot as plt
import Cloud
from plotcloud import plotcloud
from plotcloud import maplwc
from CloudDataSorter import CloudDataSorter
from CloudDataSorter import dNdlogDp2N
import os
import copy as copy
from mpl_toolkits.basemap import Basemap
import Tkinter, tkFileDialog


########## File Sources ###########
dirnames=[]; pathnames=[]; fnames=[];
print("Please choose the basic aircraft directory:")
pathdef=""
while len(pathdef)==0:
    root = Tkinter.Tk()
    pathdef = tkFileDialog.askdirectory(parent=root,initialdir="/home/shared/NARE1993/Tabled_Data",title='Please select a directory')
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
    fname=fnames[Q]
except: print("Failed!")




##############################################################################################
##########################          Reading an excel file         ############################
##############################################################################################


tall=[]
sall=[]
colum=[]
colus=[]
js=0
c2=0

file = xlrd.open_workbook(fname) # open excel workbook
#for s in file.sheets(): # loop over all of the worksheets in the workbook
#   if s.name == """018F03_160805 (all + aves)""":
#s=file.sheet_by_name("""018F03_160805 (all + aves)""")  # project #1
s=file.sheet_by_index(0)       # project #2
print '\n'
i=1
while 1:
    print("line "+str(i)+" - "+str(s.row_values(i,0)[0:6]))
    L=raw_input("[t] for titles\n[u] for units\n[d] for data\n[Enter] to move on\t\t---> ")
    if L=='t':
        titles=s.row_values(i,0)
        i=i+1
    elif L=='u':
        units=s.row_values(i,0)
        i=i+1
    elif L=='d':
        dataline=i
        break
    else:
        i=i+1
for c in range(0,s.ncols):
    rtype=s.cell(dataline,c).ctype
#    print("Column %d is of type no. %d" % (c+1, rtype))
    if rtype==0 or rtype==1:  # we need to check columns with rtype==1 (string) and rtype==0 (empty) in the first row
        i=1
        while rtype==0 and dataline+i<=s.nrows-1:      # if there is at least one non-empty in the column, we are going to change the row type rtype, 
                                                # we ignore strings deliberately
            rtype=s.cell(dataline+i,c).ctype
            i+=1
#            if rtype !=0:
#                print("Column %d had empties, the first non-empty cell on row %d is of type no. %d" % (c+1, i, rtype))
    if rtype !=0:    
        tmp1=s.col_values(c,dataline) # col_values(col number,start row, optnl end row)
        for i in range(0,len(tmp1)):
            if tmp1[i]=='':
                tmp1[i]=float('NaN')       
        tmp=np.array(tmp1)       
        # although these were columns in the xls file, they are now rows.
# Notes: if I turn tmp into an array and try to remove '' (empties) --> tmp[tmp=='']=nan -->, the first line is all nan. 
        if rtype==2 or rtype==3:         # type: float (or date or time)
            if c2==0:
                mall=tmp
                c2=10
                colum.append(c)
            else:
                mall=vstack([mall,tmp])
                colum.append(c)
        else:                           # all other row types (mostly strings)
            sall.append(tmp)
            colus.append(c)

colum=np.array(colum)
colus=np.array(colus)

# Now reading the description and offering choices of columns to plot. 


for i in range(len(mall)):
    if mall[i].dtype==dtype('<U6') or mall[i].dtype==dtype('<U7') or mall[i].dtype==dtype('<U8'):
        mall[i]=mall[i].astype(float)
mall_ttl=[]; mall_unit=[];
j=0
for i in range(0,len(titles)):
    if i in colum:
        print("%d- %s (%s)" % (j, titles[i], units[i]))
        mall_ttl.append(titles[i]);   mall_unit.append(units[i]); 
        j+=1



##############################################################################################
##########################       Generating plots on demand       ############################
##############################################################################################
figi=0     # figure counter, now we have no figures opened
ion()     # interaction mode is on (no need to use show and have the code being suspended.
randomplots = raw_input("Would you like to plot data of your choice? (y/[n])  ")
if randomplots=='y':
    # prompting for choice of data to plot
    nx, ny = int(raw_input("Enter the number of the x-axis variable you wish to plot:  ")), int(raw_input("Enter the number of the y-axis variable you wish to plot:  "))
    figure(1)
    plot(mall[nx,:],mall[ny,:],'bo')
    xlabel(mall_ttl[nx])
    ylabel(mall_ttl[ny])
    # asking to add more data on the same plot
    more=raw_input("Do you wish to plot more data? (y/n)     ")
    markervec=['bo','g+','rd','cs','m*','kx','yh']
    j=1     # index of the markervec we are going to use next (j=0 was already use)
    figi=1     # index of the number of opened figures
    while more=='y':
        ng=raw_input("[Same] graph or new (n) graph?  ")
        nx, ny = int(raw_input("Enter the number of the x-axis variable you wish to plot:  ")), int(raw_input("Enter the number of the y-axis variable you wish to plot:  "))
        if ng=='n':
            figi+=1    # counting figures
            j=0      # reinitializing the marker vector markervec (since there is no confusion with the previous data
            figure(figi) # making a new figure
            xlabel(mall_ttl[nx])
            ylabel(mall_ttl[ny])
        # end of if ng loop
        plot(mall[nx,:],mall[ny,:],markervec[j])
        j+=1
        more=raw_input("Do you wish to plot more data? (y/n)     ")
        del ng
    # end of while more loop
    clfig=raw_input("Close all figures? (y/n)  ")
    if clfig=='y':
        close('all')
        figi=0     # number of open figures=0
    # end of if clfig loop
    del titles, nx, ny, more, j, markervec, clfig
# end of if randomplots=='y' loop


DCQ=raw_input("Define Cloud? [y]/n  ")
if DCQ=='n':
    pass
else:
    ###############################     Sorting the cloud data     #################################

    [Sdata, Sttl, Sunit, sd]=CloudDataSorter(mall,mall_ttl,mall_unit)


    ##############################################################################################
    #############################        Cloud Investigation       ###############################
    ##############################################################################################
    # Now we choose the columns for altitude, particle conc and CCN conc and a longitude/latitude plot
    # We then select the time during which we are "in" one cloud
    # After that we select each verticals, horizontals, below cloud, etc.
    print("\n We are now going to determine passages of the aircraft around and in the cloud (scans). \n Define only one physical cloud at the time. \n You will be asked to define the beginning and end time for a cloud and for different scans (horizontal in cloud, above cloud, below cloud and vertical). A horizontal scan with segments separated by at least 100 m must be considered as seperate scans. \n Follow the instructions closely.\n")
    # acquiring the column numbers for the needed variables.
    pos0=[i for i,x in enumerate(Sttl) if x == 'time']
    pos1=[i for i,x in enumerate(Sttl) if x == 'altitude']
    pos2=[i for i,x in enumerate(Sttl) if x == 'latitude']
    pos3=[i for i,x in enumerate(Sttl) if x == 'longitude']
    pos4=[i for i,x in enumerate(Sttl) if x == 'LWC']
    try: tim=Sdata[pos0]; alt=Sdata[pos1]; lat=Sdata[pos2]; lon=Sdata[pos3]; lwc=Sdata[pos4]
    except: print("Crucial data is missing. This needs to be corrected before we can go on.")
    ptc=[sd[i]["total"] for i,x in enumerate(sd) if 'PC'.upper() in x["Distname"].upper()]
    ccnc=[sd[i]["total"] for i,x in enumerate(sd) if 'FSSP96'.upper() in x["Distname"].upper()]
    if not ccnc: ccnc=[sd[i]["total"] for i,x in enumerate(sd) if 'FSSP100'.upper() in x["Distname"].upper()]
    try:
        c2d0=[i for i,x in enumerate(sd) if '2d'.lower() in x["Distname"].lower()][0]
        c2d=vstack((sd[c2d0]["time"],dNdlogDp2N(sd[c2d0],100,nan),dNdlogDp2N(sd[c2d0],200,nan)))
    except:
        c2d=np.ones(len(tim[0]))*NaN
    if np.shape(ptc)[0]==0:
        for i in range(len(Sttl)):
            print("%d- %s (%s)" % (i, Sttl[i], Sunit[i]))
        Q=raw_input("Is there a column corresponding to the particle number concentration? [n]/#  ")
        try: 
            Q=int(Q); ptc=Sdata[Q];
        except: ptc=np.ones(len(tim[0]))*NaN
    if np.shape(ccnc)[0]==0:
        for i in range(len(Sttl)):
            print("%d- %s (%s)" % (i, Sttl[i], Sunit[i]))
        Q=raw_input("Is there a column corresponding to the droplet number concentration? [n]/#  ")
        try: 
            Q=int(Q); ccnc=Sdata[Q];
        except: ccnc=np.ones(len(tim[0]))*NaN
    mall=vstack((tim,alt,lwc,ptc,ccnc,lat,lon))
    mcols=np.array([0, 1, 2, 3, 4, 5, 6])
    
    Adj=[1,10,100,1000,1,1]     # Default magnifying values for [0-Altitude, 1-Ptcl conc, 2-droplet conc, 3-LWC, 4-2d100 drop conc, 5-2d200 drop conc]
    Offsts=[0,0,0,0,0,0]  # default offset values for [0-Altitude, 1-Ptcl conc, 2-droplet conc, 3-LWC, 4-2d100 drop conc, 5-2d200 drop conc]
    figi+=1            # Figure number
    Zx1=np.min(tim); Zx2=np.max(tim)      # initial time span of the figure
    Zs=np.array([Zx1,Zx2,NaN,NaN,NaN,NaN,NaN,NaN])
    [Zs,Adj,Offsts,figi,fig]=plotcloud(mall,mcols,c2d,Zs,Adj,Offsts,figi,1)
    timix=np.nonzero((mall[0]>=Zx1)*(mall[0]<=Zx2))
    try:
        maplwc(mall[0][timix],mall[2][timix],mall[5][timix],mall[6][timix],1)    # appears in figure 1002
    except:
        pass
    Cnum=int(raw_input("How many clouds do you wish to enter from this file?  "))
    close(figi)
    
    #############################     Choosing the Cloud times     ###############################
    Clistmp=[]
    for Cnumi in range(Cnum):
        Zs=np.array([Zx1,Zx2,NaN,NaN,NaN,NaN,NaN,NaN])
        [Zs,Adj,Offsts,figi,fig]=plotcloud(mall,mcols,c2d,Zs,Adj,Offsts,figi,1) 
        Zsdaynum=np.floor(Zs[0])
        Times=dict()
        print("Click at the start and at the end of the cloud period.")
        L=ginput(2,timeout=0) 
        if L[0][0]<L[1][0]: I1=L[0][0]; I2=L[1][0]
        else: I1=L[1][0]; I2=L[0][0]
        Times["cloud"]=np.reshape(np.array([I1+Zsdaynum,I2+Zsdaynum]),(-1,2))
        Zs[0]=Times["cloud"][0][0]; Zs[1]=Times["cloud"][0][1];       # updating the time frame to include only the cloud we're looking at.
        malli=mall[:,nonzero((mall[0,:]>=Zs[0])*(mall[0,:]<=Zs[1]))[0]]
        
        close(figi)
        [Zs,Adj,Offsts,figi,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,0)
        close(1002)     # close flight map
        timix=np.nonzero((mall[0]>=Zs[0])*(mall[0]<=Zs[1]))
        try:
            maplwc(mall[0][timix],mall[2][timix],mall[5][timix],mall[6][timix],1)    # plot new map
        except:
            pass
        VQ=int(raw_input("How many VERTICAL scan entries do you wish to enter?  "))
        L=[]
        for j in range(VQ):
            close(figi)
            [Zzz,Azz,Offzz,iz,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,1)
            print("Vertical scan #%d. Click at the start and end of the vertical scan." % (j+1))
            Ip=ginput(2,timeout=0)
            if Ip[0][0]<Ip[1][0]: I1=Ip[0][0]; I2=Ip[1][0]
            else: I1=Ip[1][0]; I2=Ip[0][0]
            L.append(np.array([I1+Zsdaynum,I2+Zsdaynum]))
        Times["verticloud"]=np.reshape(np.array(L),(-1,2))
        
        close(figi)
        [Zs,Adj,Offsts,figi,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,0)    
        VQ=int(raw_input("How many BELOW CLOUD scan entries do you wish to enter?  "))
        L=[]
        for j in range(VQ):
            close(figi)
            [Zzz,Azz,Offzz,iz,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,1)
            print("Below cloud scan #%d. Click at the start and end of the below cloud scan." % (j+1))
            Ip=ginput(2,timeout=0)
            if Ip[0][0]<Ip[1][0]: I1=Ip[0][0]; I2=Ip[1][0]
            else: I1=Ip[1][0]; I2=Ip[0][0]
            L.append(np.array([I1+Zsdaynum,I2+Zsdaynum]))
        Times["belowcloud"]=np.reshape(np.array(L),(-1,2))
        
        close(figi)
        [Zs,Adj,Offsts,figi,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,0)    
        VQ=int(raw_input("How many ABOVE CLOUD scan entries do you wish to enter?  "))
        L=[]
        for j in range(VQ):
            close(figi)
            [Zzz,Azz,Offzz,iz,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,1)
            print("Above cloud scan #%d. Click at the start and end of the above cloud scan." % (j+1))
            Ip=ginput(2,timeout=0)
            if Ip[0][0]<Ip[1][0]: I1=Ip[0][0]; I2=Ip[1][0]
            else: I1=Ip[1][0]; I2=Ip[0][0]
            L.append(np.array([I1+Zsdaynum,I2+Zsdaynum]))
        Times["abovecloud"]=np.reshape(np.array(L),(-1,2))
        
        close(figi)
        [Zs,Adj,Offsts,figi,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,0)    
        VQ=int(raw_input("How many HORIZONTAL CLOUD scan entries do you wish to enter?  "))
        L=[]
        howmany=range(VQ)
        height_diff=100. #max difference in height before the warning
        for j in howmany:
            contin=0
            while contin==0:
                close(figi)
                [Zzz,Azz,Offzz,iz,fig]=plotcloud(malli,mcols,c2d,Zs,Adj,Offsts,figi,1)
                print("Horizontal cloud scan #%d. Click at the start and end of the horizontal cloud scan." % (j+1))
                Ip=ginput(2,timeout=0)
                if Ip[0][0]<Ip[1][0]: I1=Ip[0][0]; I2=Ip[1][0]
                else: I1=Ip[1][0]; I2=Ip[0][0]
                L.append(np.array([I1+Zsdaynum,I2+Zsdaynum]))
                minindex=(np.where(Sdata[0]>=L[-1][0])[0].data).min()
                maxindex=(np.where(Sdata[0]<=L[-1][1])[0].data).max()
                minalt=(Sdata[1][minindex:maxindex]).min()
                maxalt=(Sdata[1][minindex:maxindex]).max()
                if maxalt-minalt>height_diff:
                    erase=raw_input("WARNING: There is a difference of at least "+str(int(height_diff))+"m between the min and max height in this Horizontal cloud scan. Would you like to break it up into two horizontal scans? (y/[n])   ")
                    if erase=='y':
                        del L[-1]
                        howmany.append(len(howmany))
                        pass
                    else:
                        contin=1
                else:
                    contin=1
                
        Times["horicloud"]=np.reshape(np.array(L),(-1,2))
        
        
        ############################### Cutting the data to the cloud time ###############################
        Sdatai=Sdata[:,nonzero((Sdata[0,:]>=Times["cloud"][0][0])*(Sdata[0,:]<=Times["cloud"][0][1]))[0]]
        sdi=copy.deepcopy(sd)
        for i in range(len(sdi)):
            if type(sdi[i]["time"])==float:
                pass
            else:
                sdi[i]["total"]=sdi[i]["total"][:,nonzero((sdi[i]["time"]>=Times["cloud"][0][0])*(sdi[i]["time"]<=Times["cloud"][0][1]))[0]]
                sdi[i]["data"]=sdi[i]["data"][:,nonzero((sdi[i]["time"]>=Times["cloud"][0][0])*(sdi[i]["time"]<=Times["cloud"][0][1]))[0]]
                sdi[i]["mdp"]=sdi[i]["mdp"][:,nonzero((sdi[i]["time"]>=Times["cloud"][0][0])*(sdi[i]["time"]<=Times["cloud"][0][1]))[0]]
                sdi[i]["medp"]=sdi[i]["medp"][:,nonzero((sdi[i]["time"]>=Times["cloud"][0][0])*(sdi[i]["time"]<=Times["cloud"][0][1]))[0]]
                sdi[i]["time"]=sdi[i]["time"][:,nonzero((sdi[i]["time"]>=Times["cloud"][0][0])*(sdi[i]["time"]<=Times["cloud"][0][1]))[0]] # this line should be last in the sequence because times are used for the nonzero function of earlier lines!
                
        
    
        ###############################     Creating a Cloud instance     #################################
        
        A=Cloud.Cloud(Sdatai,Sttl,Sunit,Times,sdi)
        del Sdatai, sdi
        try:
            A["convairfile"]=fname     # name of the file where the data comes from.
            try: A["convairsheet"]=sheetname    # sheet if the file was excel or excel-like.
            except: pass
        except: print("The source is not an excel file, it doesn't have fname and/or sheetname")
        ###########################   Entering the description of the cloud   #############################
        A.describe()
        close(1002)
        close(figi)
        
        ############################## storage in temp cloud list ##############################
        Clistmp.append(A)
        
###############################     Clearing python's memory     #################################
try: del mall, Zs, tmp, colum, mcols, colus, lon, tim, lat, alt, ptc, ccnc
except: pass
try: del Sdata,Sttl, Sunit, Times, sd
except: pass
try: del  Zzz,Azz,Offzz,iz,fig
except: pass
