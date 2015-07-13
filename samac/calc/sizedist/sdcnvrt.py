################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####

from pylab import *
import numpy as np
import copy as copy    


def dNdlogDp2N(dist,z1,z2):
    # this function return the total particle concentration when the distribution is given in dN/dlogDp
    # dN is the concentrations (in dN/dlogZ) and Z is the diameter vector
    # z1 and z2 are the boundaries within which we sum (inclusively)
    # this returns the total concentration
    # SG, Feb 2008

    # modified to take in matrices dN --> lines: times columns: diameters
    # modified to make total concentrations between given diameters (z1 
    # (smaller diameter) and z2 (bigger diameter)). (z1=z2=nan if no boundaries)
    # also removes negative total concentrations and replaces it with NaN.
    if type(dist["time"])==float: return NaN
    if len(np.shape(dist["data"]))<2 and np.shape(dist["time"])[0]>1:
        N=np.ones(np.shape(dist["time"]))*NaN; 
    else:
        if np.shape(dist["time"])[0]==1: dist["data"]=dist["data"].reshape(np.shape(dist["data"])[0],1);
        Z=copy.deepcopy(dist["bins"])
        t=len(Z)    # number of bins
        if t==np.shape(dist["data"])[1]:
            dN=copy.deepcopy(dist["data"])
        elif t==np.shape(dist["data"])[0]:
            dN=copy.deepcopy(dist["data"].transpose())
        
        # predefining arrays for lower (b) and upper (h) bin limits        
        b=np.zeros([t,1]);
        h=zeros([t,1]);
        # calculating size bins lower (b) and upper (h) limits
        for i in range(1,t):
            b[i]=(np.log10(Z[i-1])+np.log10(Z[i]))/2.;
        for i in range(0,t-1):
            h[i]=(np.log10(Z[i+1])+np.log10(Z[i]))/2.;
        b[0]=np.log10(Z[0])+(log10(Z[1])-h[1]);
        h[t-1]=np.log10(Z[t-1])-(b[t-2]-np.log10(Z[t-2]));
        dlZ=h-b;    # size bin width
        
        # investigating the z1 and z2 (integration limits)
        lz1=np.log10(z1); lz2=np.log10(z2);
        if not isnan(z1):
            if lz1<b[0]: print("[dNdlogDp2N] Lower integration limit (z1=%.2f) is outside measurement range (%.2f). Integration performed from the lowest available size." %(z1,10**b[0])); z1=nan; 
        if not isnan(z2):
            if lz2>h[-1]: print("[dNdlogDp2N] Upper integration limit (z2=%.2f) is outside measurement range (%.2f). Integration performed from the highest available size." %(z2,10**h[-1])); z2=nan; 
        # finding the 'full bins': when the full size range of the bin is within z1 and z2
        if isnan(z2) and isnan(z1): ix=range(len(Z));
        elif isnan(z2): ix=np.nonzero((b>=lz1))[0]
        elif isnan(z1): ix=np.nonzero((h<=lz2))[0]
        else: ix=np.nonzero((h<=lz2)*(b>=lz1))[0]
        # now ix is the index of the full bins.
        
        # finding 'fractionned bins': when part of the size range of the bin is within z1 and z2
        pts=np.shape(dN)[0] # number of entries (points in time)
        bfix=np.nonzero((lz1>b)*(lz1<h))[0]   # smallest size fractionned bin
        hfix=np.nonzero((lz2>b)*(lz2<h))[0]   # biggest size fractionned bin
        if len(bfix)==1: 
            bfN=NaN*ones([pts,1]);
            for i in range(0,pts):
                bfN[i]=dN[i,bfix]*((h[bfix]-lz1)/(h[bfix]-b[bfix]))*dlZ[bfix].transpose()
        elif len(bfix)==0: bfN=np.zeros([pts,1]);       # if there is no such fractionned bin then the concentration is zero.
        else: raise ValueError("Lower integration limit z1 found in multiple size bins.")
        if len(hfix)==1:
            hfN=NaN*ones([pts,1]);
            for i in range(0,pts):
                hfN[i]=dN[i,hfix]*((lz2-b[hfix])/(h[hfix]-b[hfix]))*dlZ[hfix].transpose()
        elif len(hfix)==0: hfN=np.zeros([pts,1]);
        else: raise ValueError("Upper integration limit z2 found in multiple size bins.")

        N=NaN*ones([pts,1]);
        for i in range(0,pts):
            N[i] = sum(dN[i,ix]*dlZ[ix].transpose());       # we don't want nansum here in case the nan is not due to a negative concentration but rather to measurement problems
        N=N+bfN+hfN
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


