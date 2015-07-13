
from pylab import *
import numpy as np
import scipy.stats.stats as st

#########################################################################        
################################  cinfo  ################################
#########################################################################
def cinfo(CL,param):
    """ This property returns information on the parameter in the cloud (all given in the units of the parameter). Note that the parameter is averaged over the entire cloud time at the altitude required (bottom, top or in-cloud) - not the case using vpinfo(CL,param).
        CloudObj.cinfo["bottom"]: param at the cloud base
        CloudObj.cinfo["top"]: param at the cloud top
        CloudObj.cinfo["mean"]: mean param through the cloud (in cloud)
        CloudObj.cinfo["median"]: median param through the cloud (in cloud)
        CloudObj.cinfo["stdev"]: standard deviation of the param through the cloud (in cloud)
        CloudObj.cinfo["delta"]: difference of param between the bottom and the top
        CloudObj.cinfo["slope"]: delta divided by the mean thickness
        The property can be accessed as e.g. CloudObj.cinfo["bottom"] or CloudObj.cinfo (dictionary) """
    H=dict()
    H["bottom"]=list(); H["top"]=list(); H["mean"]=list(); H["median"]=list(); H["stdev"]=list(); H["delta"]=list(); H["slope"]=list(); H["units"]=list(); 
    alt=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
    T=[i for i,x in enumerate(CL.dttl) if x == param][0]
    try:
        for i in range(len(CL.props["height"])):
            ix=nonzero((CL.data[alt]>=CL.props["height"][i][1])*(CL.data[alt]<=CL.props["height"][i][2]))
            H["bottom"].append(float(st.nanmedian(CL.data[T][nonzero((CL.data[alt]>=CL.props["height"][i][0])*(CL.data[alt]<=CL.props["height"][i][1]))])))
            H["top"].append(float(st.nanmedian(CL.data[T][nonzero((CL.data[alt]>=CL.props["height"][i][2])*(CL.data[alt]<=CL.props["height"][i][3]))])))
            H["mean"].append(float(st.nanmean(CL.data[T][ix])))
            H["median"].append(float(st.nanmedian(CL.data[T][ix])))
            H["stdev"].append(float(st.nanstd(CL.data[T][ix])))
            H["delta"].append(H["bottom"][i]-H["top"][i])
            H["slope"].append(H["delta"][i]/(np.mean([CL.props["height"][i][2], CL.props["height"][i][3]])-np.mean([CL.props["height"][i][0], CL.props["height"][i][1]])))     # units/meter
            H["units"].append(CL.dunit[T])
            del ix
    except: print("[cinfo] Height properties must be defined first using the defheight method.")
    return H
