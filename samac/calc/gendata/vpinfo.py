################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####


from pylab import *
import numpy as np
import scipy.stats.stats as st

########################################################################        
###############################  vpinfo  ###############################
########################################################################
def vpinfo(CL,param,base='bg'):
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
    altp=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
    tim=[i for i,x in enumerate(CL.dttl) if x == 'time'][0]
    T=[i for i,x in enumerate(CL.dttl) if x.lower() == param.lower()]
    if len(T)==1: 
        T=T[0]
        Td=CL.data[T]
        Tunits=CL.dunit[T]
        alt=CL.data[altp]
        ta=CL.data[tim]
    elif len(T)>1: print("[vpinfo] Parameter %s was found multiple times in the basic data." %(param)); return dict()
    elif len(T)==0:
        posx=[] 
        for i,ttl in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == param.lower()]    # check all titles matching with temperature
        if len(posx)==1: 
            Td=CL.extradata[posx[0][0]][posx[0][1]]    # loading the data
            Tunits=CL.extraunit[posx[0][0]][posx[0][1]]
            j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
            Tt=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            # adapting for too short data for interpolation
            if len(Tt)<2: Td=np.ones((2,))*NaN; Tt=np.array([CL.times["cloud"][0][0],CL.times["cloud"][0][1]]);
            # adapting the time vector to a common time vector
            ta1=np.max([CL.data[tim][0],Tt[0]]); ta2=np.min([CL.data[tim][-1],Tt[-1]]);
            ta=CL.data[tim][nonzero((CL.data[tim]>=ta1)*(CL.data[tim]<=ta2))[0]]
            alt=CL.data[altp][nonzero((CL.data[tim]>=ta1)*(CL.data[tim]<=ta2))[0]]
            fT=interpolate.interp1d(Tt,Td,kind='linear')
            Td=fT(ta)
        else: print("[vpinfo] No or multiple %s found in the basic or the extra data." %(param)); return dict()
    
    H["bottom"]=list(); H["top"]=list(); H["mean"]=list(); H["median"]=list(); H["stdev"]=list(); H["delta"]=list(); H["slope"]=list(); H["units"]=list(); H["minimum"]=list();  H["maximum"]=list();
    try:
        for i in range(len(CL.times["verticloud"])):
            if base=='4point': cb=CL.props["height"][i][1]; ct=CL.props["height"][i][2];
            else: cb=CL.props["BGheight"][i][0]; ct=CL.props["BGheight"][i][1];
            ix=nonzero((alt>=cb)*(alt<=ct)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))[0]
            if len(ix)==0:
                H["mean"].append(nan); H["median"].append(nan); H["stdev"].append(nan); H["minimum"].append(nan); H["maximum"].append(nan); H["top"].append(nan); H["bottom"].append(nan); H["delta"].append(nan); H["slope"].append(nan); H["units"].append(nan)
            else:
                H["mean"].append(float(st.nanmean(Td[ix])))
                H["median"].append(float(st.nanmedian(Td[ix])))
                H["stdev"].append(float(st.nanstd(Td[ix])))
                H["minimum"].append(float(np.nanmin(Td[ix])))
                H["maximum"].append(float(np.nanmax(Td[ix])))
                if base=='4point': 
                    if len(nonzero((alt>=ct)*(alt<=CL.props["height"][i][3])*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))[0])==0: H["top"].append(nan)
                    else: H["top"].append(float(st.nanmedian(Td[nonzero((alt>=ct)*(alt<=CL.props["height"][i][3])*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))])))
                    if len(nonzero((alt>=CL.props["height"][i][0])*(alt<=cb)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))[0])==0: H["bottom"].append(nan)
                    else: H["bottom"].append(float(st.nanmedian(Td[nonzero((alt>=CL.props["height"][i][0])*(alt<=cb)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))])))
                    H["delta"].append(H["bottom"][i]-H["top"][i])
                    H["slope"].append(H["delta"][i]/(np.mean([ct, CL.props["height"][i][3]])-np.mean([CL.props["height"][i][0], cb])))
                else: 
                    R=10     # plus/minus R meters around the cloud top
                    if len(nonzero((alt>=ct-R)*(alt<=ct+R)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))[0])==0: H["top"].append(nan)
                    else: H["top"].append(float(st.nanmedian(Td[nonzero((alt>=ct-R)*(alt<=ct+R)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))])))
                    if len(nonzero((alt>=cb-R)*(alt<=cb+R)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))[0])==0: H["bottom"].append(nan)
                    else: H["bottom"].append(float(st.nanmedian(Td[nonzero((alt>=cb-R)*(alt<=cb+R)*(ta>=CL.times["verticloud"][i][0])*(ta<=CL.times["verticloud"][i][1]))])))
                    H["delta"].append(H["bottom"][i]-H["top"][i])
                    H["slope"].append(float(H["delta"][i]/(ct-cb)))
                H["units"].append(Tunits)
            del ix
    except: 
        if base=='4point': print("[vpinfo] Height properties must be defined first using the defheight method.")
        else: print("[vpinfo] Height properties must be defined first using the defBGheight method.")
    return H     

