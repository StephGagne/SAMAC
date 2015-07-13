
from pylab import *
import numpy as np
import scipy.stats.stats as st

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
