
from pylab import *
import numpy as np

import samac

#########################################################################        
###############################  angles  ################################
#########################################################################
def angles(CL):
    """ This property calculates the vertical angle of the aircraft (or other real or hypothetical aerial vehicle) based on TAS and Altitude change.
    It also calculates the horizontal angle of the aircraft based on the longitude and latitude. The TAS must be in the basic data (extradata module not handled).
    Returns a list with: time, vertical angle, horizontal angle (with one point less than the time in CLOUD.data)."""
    
    # vertical angle calculations
    palt=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
    ptim=[i for i,x in enumerate(CL.dttl) if x == 'time'][0]
    ptas=[i for i,x in enumerate(CL.dttl) if 'TAS' in x.upper()]
    if len(ptas)==1: ptas=ptas[0]
    elif len(ptas)==0:
        print("[angles] True Air Speed (TAS) not found.")
        for i,title in enumerate(CL.dttl):
            print("%d - %s (%s)" % (i,title,CL.dunit[i]))
        try: ptas=int(raw_input("Enter the number of the TAS (press return if none correspond):  "))
        except: print("[angles] no TAS available: failure of method angles.")
    elif len(ptas)>=2:
        print("[angles] Multiple possible True Air Speed (TAS) were found.")
        for i in range(len(ptas)):
            print("%d - %s (%s)" % (i,CL.dttl[ptas[i]],CL.dunit[ptas[i]]))
        p=int(raw_input("Enter the number of the TAS (press return if none correspond):  "))
        try: ptas=ptas[0]
        except: print("[angles] no TAS available: failure of method.")
    else: print("[angles] Very strange indeed.")
    if type(ptas)==int: pass
    else: return "Error: could not find TAS position!"
    dalt=np.diff(samac.runstats(CL.data[palt],5)[0])
    dt=np.diff(CL.data[ptim])*24*60*60        # time in seconds
    tas=CL.data[ptas]; tas=tas[0:-1]      # removing the last point so it has the same length as the differentiated vectors
    tas=samac.runstats(tas,5)[0]      # running average of tas
    # the angle zero is when the aircraft is parallel to the ground
    VA=360*np.arcsin(dalt/(tas*dt))/(2*pi)      # plane angle, in degrees (compared to horizontal)
    
    # horizontal angle calculations
    plat=[i for i,x in enumerate(CL.dttl) if x == 'latitude'][0]
    plon=[i for i,x in enumerate(CL.dttl) if x == 'longitude'][0]
    # Let's define East as the angle zero.
    dlat=np.diff(samac.runstats(CL.data[plat],5)[0]); 
    dlon=np.diff(samac.runstats(CL.data[plon],5)[0]); 
    # assuming NA quadrant  MMM - this needs to be fixed (although right now the data has only positive latitudes/longitudes)
    # assuming parallel enough longitude lines
    HA=360*np.arctan(abs(dlat/dlon))/(2*pi)     # calculating everything as a first quadrant
    HA[((dlat>=0)*(dlon>0))]=180-HA[((dlat>=0)*(dlon>0))]   # 2nd quadrant
    HA[((dlat<0)*(dlon>0))]=180+HA[((dlat<0)*(dlon>0))]     # 3rd quadrant
    HA[((dlat<0)*(dlon<0))]=360-HA[((dlat<0)*(dlon<0))]   # 4th quadrant
    HA[((dlat<0)*(dlon==0))]=270
    HA[((dlat>0)*(dlon==0))]=90
    
    return [CL.data[ptim][0:-1], VA, HA]

