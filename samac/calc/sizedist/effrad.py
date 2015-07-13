################## Copyright 2013-2015 Stephanie Gagne ##################
#### Distributed under the terms of the GNU General Public License 3 ####

#import pylab as *
import numpy as np


########################################################################        
###############################  effrad  ###############################
########################################################################
def effrad(CL,inst,bindist='lin'):
    """ This method returns the effective radius for a given instrument for the entire cloud period. The radius is in the same units as the instrument's units (usually micrometers). Note that one might get the effective diameter if the instrument's size bins are diameters.
        example: CloudObj.effrad(inst='FSSP96',bindist='lin')
        bindist is 'lin' if the difference between bins is linearly distributed (FSSPs) and 'log' if they are logarythmically distributed (PCASP)"""
    # according to the formula in https://en.wikipedia.org/wiki/Cloud_drop_effective_radius latest access on Oct 2013.
    [pos,sd]=[[i,sd] for i,sd in enumerate(CL.sd) if sd["Distname"].lower() == inst.lower()][0]
    # building dr (dradius) vector
    R=sd["bins"]; t=len(R)
    b=np.zeros([t]);
    h=np.zeros([t]);
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

