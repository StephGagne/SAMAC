from pylab import *
import numpy as np




#########################################################################        
############################   belowprof   ##############################
#########################################################################            
def belowprof(CL):
    """This is a returning method. This method returns a list with one dictionary per vertical scan that includes the index of what belongs to below the vertical scan, the area index (binary), and the times. This method is used to access the special profile 'belowvert'.
    Returned list format: 
    - R[0]["index"]: index of what is below cloud and within 0.0075 degrees radius (below the vertical scan)
    - R[0]["areaindex"]: index of what is in the area (within 0.0075 degrees radius)
    - R[0]["times"]: times of what is below the vertical scan."""
    BV=list()       # BV has as many list entries as the cloud has vertical scans.
    # each list element is a dictionary including times, average aerosol distribution,  ETC.?
    ctime=[i for i,x in enumerate(CL.dttl) if x == 'time']; ctime=ctime[0]
    clat=[i for i,x in enumerate(CL.dttl) if x == 'latitude']; clat=clat[0]
    clon=[i for i,x in enumerate(CL.dttl) if x == 'longitude']; clon=clon[0]
    calt=[i for i,x in enumerate(CL.dttl) if x == 'altitude']; calt=calt[0]
    res=0.0075 # resolution in degrees
    for j in range(len(CL.times["belowcloud"])):
            bctimetmp=(CL.data[ctime]>=CL.times["belowcloud"][j][0])*(CL.data[ctime]<=CL.times["belowcloud"][j][1])
            try: bctime=bctime+bctimetmp; 
            except: bctime=bctimetmp; 
    if len(CL.times["belowcloud"])==0: bctime=False*ones(np.shape(CL.data[ctime])); print("[belowprof] No below cloud scans")
    for i in range(len(CL.times["verticloud"])):      # for each vertical scan
        t=dict()
        lat=CL.data[clat][nonzero((CL.data[ctime]>=CL.times["verticloud"][i][0])*(CL.data[ctime]<=CL.times["verticloud"][i][1]))[0]]
        lon=CL.data[clon][nonzero((CL.data[ctime]>=CL.times["verticloud"][i][0])*(CL.data[ctime]<=CL.times["verticloud"][i][1]))[0]]
        scant=CL.data[ctime][nonzero((CL.data[ctime]>=CL.times["verticloud"][i][0])*(CL.data[ctime]<=CL.times["verticloud"][i][1]))[0]]
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
        arix=(np.ones(np.shape(CL.data[clat]))*0)
        for m in range(len(pt)):
            arix=arix+((CL.data[clat]-pt[m][1])**2+(CL.data[clon]-pt[m][2])**2<=res**2)
        # now finding the index of points within this area and below the cloud base in the same vertical or below cloud:
        try: Ix=nonzero( arix * ( ((CL.data[ctime]>=CL.times["verticloud"][i][0])*(CL.data[ctime]<=CL.times["verticloud"][i][1])*(CL.data[calt]<CL.props["height"][i][0])) + bctime ) )
        except: raise ValueError("[belowprof] Is the c.props['height'] defined? This may be the source of your problem.")
        
        t["index"]=Ix[0]    # index of what is below cloud and within 0.0075 degrees radius (below the vertical scan)
        t["areaindex"]=arix     # index of what is in the area (within 0.0075 degrees radius)
        t["times"]=CL.data[ctime][Ix]     # times of what is below the vertical scan
        
        BV.append(t)
    return BV
    

