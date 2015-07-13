
from pylab import *
import numpy as np

    
########################################################################        
###############################  filler  ###############################
########################################################################
def filler(CL,ts,timestep=10.0,inst='2dc'):
    """This method returns a filled distribution (type dict), filled with zeros when no measurements were returned. 
    This was designed for 2d data: this instrument only reports data when it is detected (will not report zero measurements.
    ts: start and end time over which the data should be filled. Normally a couple of times from CloudObj.times['manoeuvre'][scan#].
    timestep: (default 10.0), in seconds (must be a float), is the time step with which the instrument returns data when measurements are non-zero.
    inst: (default '2dc') is a string format corresponding to the name of the distribution that needs to be filled."""
    newsd=dict()        # define the returning dictionary
    pos2d=[i for i,x in enumerate(CL.sd) if x["Distname"].lower()==inst.lower()]; 
    try:                # verifying that the instrument indeed exists and can be found
        pos2d=pos2d[0]
        tim2d=CL.sd[pos2d]["time"];
        prec=CL.sd[pos2d]["data"];
        ptot=CL.sd[pos2d]["total"];
        bins=CL.sd[pos2d]["bins"]; 
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

