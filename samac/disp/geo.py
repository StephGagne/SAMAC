from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import copy
from math import isnan
import scipy.stats.stats as st

from mpl_toolkits.basemap import Basemap

import samac


#########################################################################        
#############################   mapcloud   ##############################
#########################################################################
def mapcloud(CL,interact=1):
    """ This method plots the trajectory of the aircraft during the whole cloud period over a map of the region. Time is colour-coded and the thickness of the line reflects the LWC.
        Use CloudObj.mapcloud() or CloudObj.mapcloud(1) to continue using ipython (interactive mode);
        Use CloudObj.mapcloud(0) if you would rather close the figure before continuing to use ipython (the plot is only plotted after the whole code is read).
        Extradata module is handled for LWC.  """
            
    # time
    pos=[i for i,x in enumerate(CL.dttl) if x == 'time']
    if len(pos)==1: tim=CL.data[pos[0]]
    elif len(pos)>1: print("[mapcloud] Multiple time entries were found. Go tidy up your data.")
    elif len(pos)==0: print("[mapcloud] No time was found in the basic data. Go tidy up your data.") 
    else: print("[mapcloud] Unknown Error loading time.")
    # LWC
    pos=[i for i,x in enumerate(CL.dttl) if x == 'LWC']
    if len(pos)>0:
        try: 
            lwc=CL.data[pos]
        except: 
            print("[mapcloud] Multiple LWC were found in the basic data. Go tidy up your data.")
            lwc=NaN
    else:
        posx=[] 
        for i,L in enumerate(CL.extrattl):     # for all extra datasets available
            posx=posx+[[i,j] for j,x in enumerate(L) if x.upper() == 'LWC']    # check all titles matching with lwc
        if len(posx)==1: 
            lwc1=CL.extradata[posx[0][0]][posx[0][1]]    # loading the lwc data
            j=[j for j,x in enumerate(CL.extrattl[i]) if x == 'time'][0]
            t1=CL.extradata[posx[0][0]][j]     # associated time stamp
            # interpolation on the general time stamp
            lwc=nan*ones([len(tim)])    # preallocating
            ix=((tim>=t1[0])*(tim<=t1[-1]))    # index of interpolating times within the limit of the original time stamp
            f=interpolate.interp1d(t1,lwc1,'linear')    
            if 'ma' not in str(type(lwc1)).lower(): lwc1=np.ma.array(lwc1,mask=False)  # creating a mask if none exists
            fma=interpolate.interp1d(t1,lwc1.mask,kind='linear')    # interpolating the mask
            lwc[ix]=f(tim[ix])
            lwc=np.ma.array(f(tim[ix]),mask=fma(tim[ix]))
            lwc=np.ma.masked_where(isnan(lwc), lwc)
        else: 
            print("[mapcloud] No LWC (or multiple) was found in the basic data, and not in the extra data either.")
            lwc=NaN
    # lat
    pos=[i for i,x in enumerate(CL.dttl) if x == 'latitude']
    try: 
        lat=CL.data[pos]
    except: 
        print("[mapcloud] Latitude not found in basic data")
        lat=NaN
    # lon
    pos=[i for i,x in enumerate(CL.dttl) if x == 'longitude']
    try: 
        lon=CL.data[pos]
    except: 
        print("[mapcloud] Longitude not found in basic data")
        lon=NaN
   
    # masked lwc become zero so the dots are still there but the size is the minimum size
    lwc[lwc.mask]=0.
    
    if interact==1: samac.maplwc(tim,lwc,lat,lon,1)     # plotting program
    else: samac.maplwc(tim,lwc,lat,lon,0)


#########################################################################        
##############################   maplwc   ###############################
#########################################################################
def maplwc(tim,lwc,lat,lon,imode,mf=0):
    """ This function maps an aircraft trajectory. """
    # this plots the cloud's lwc on a map using basemap.
    # all inputs must be the same length (except for imode, mf)
    # imode=0 => interactive mode off, imode=1 => interactive mode on
    # mf is not usually called upon. Used in if loop below. Removal?
    # Stephanie Gagne, Dal, March 2012
    lat=np.ma.masked_outside(lat, -90, 90)
    lon=np.ma.masked_outside(lon, -180, 180)
    if imode==1: plt.ion()
    elif imode==0: plt.ioff()
        
        # regional subplot
    fig2=plt.figure(1002)
    ax1=fig2.add_subplot(121)
    try:
        llon=np.floor(np.nanmin(lon)); hlon=np.ceil(np.nanmax(lon))
        llat=np.floor(np.nanmin(lat)); hlat=np.ceil(np.nanmax(lat))
    except:
        llon=np.floor(lon.min()); hlon=np.ceil(lon.max())
        llat=np.floor(lat.min()); hlat=np.ceil(lat.max())

    m = Basemap(llcrnrlon=llon-2, llcrnrlat=llat-2, urcrnrlon=hlon+2, urcrnrlat=hlat+2, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=100)
    x, y = m(lon, lat)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.2)
    m.drawmapboundary()
    m.drawcountries()
    m.drawmapboundary(fill_color='lightcyan')
    m.fillcontinents(color='tan',lake_color='lightcyan')
    ax1.plot(x,y,'ro',markeredgewidth=0)
    m.drawparallels(np.arange(llat-2,hlat+2,1),labels=[1,0,0,0])
    m.drawmeridians(np.arange(llon-2,hlon+2,1),labels=[0,0,0,1])
    ax1.set_title('Situation')

    # zoomed-on-cloud subplot
    if mf==0:
        if len(np.shape(lwc))==2:
            lwc=lwc.data[0]
        else:
            lwc=lwc.data
    ax2=fig2.add_subplot(122)
    try:
        llon=np.nanmin(lon); hlon=np.nanmax(lon)
        llat=np.nanmin(lat); hlat=np.nanmax(lat)
    except: 
        if type(lon)==np.ma.core.MaskedArray: llon=np.min(lon); hlon=np.max(lon)
        if type(lon)==np.ma.core.MaskedArray: llat=np.min(lat); hlat=np.max(lat)
    m = Basemap(llcrnrlon=llon, llcrnrlat=llat, urcrnrlon=hlon, urcrnrlat=hlat, projection='lcc', lat_1=0.5*(llat+hlat), lon_0=0.5*(llon+hlon), resolution='i', area_thresh=100)
    x, y = m(lon, lat)
    m.drawcoastlines(linewidth=0.25)
    m.drawcountries(linewidth=0.2)
    m.drawmapboundary()
    m.drawcountries()
    m.drawmapboundary(fill_color='lightcyan')
    m.fillcontinents(color='tan',lake_color='lightcyan',zorder=0)
    if hlat-llat>0.5: m.drawparallels(np.arange(np.floor(llat),np.ceil(hlat),0.5),labels=[0,1,0,0])
    elif hlat-llat<=0.5: m.drawparallels(np.arange(np.floor(llat),np.ceil(hlat),0.05),labels=[0,1,0,0])
    if hlon-llon>1.: m.drawmeridians(np.arange(np.floor(llon),np.ceil(hlon),1),labels=[0,0,0,1])
    elif hlon-llon<=1.: m.drawmeridians(np.arange(np.floor(llon),np.ceil(hlon),0.1),labels=[0,0,0,1])
    for i in range(len(lwc)):
        if isnan(lwc[i]):
            lwc[i]=0.
    ax2.scatter(x,y,c=tim,s=70*lwc+1.5, marker='o',cmap='spectral',edgecolors='none')
    ax2.set_title('Flight map: \nsize=lwc \ncolour=time (black to white)')
    
    if imode==0: plt.show()

    
########################################################################        
###############################  path3d  ###############################
########################################################################
def path3d(CL):
    """ This method will make a 3D plot of the flightpath of the plane during measurement.
    Colours correspond to the time tracking colours (colours may not work if using an old version of matplotlib)."""
    from mpl_toolkits.mplot3d import Axes3D
    
    plat=[i for i,x in enumerate(CL.dttl) if x == 'latitude'][0]
    lat=copy.deepcopy(CL.data[plat])
    plon=[i for i,x in enumerate(CL.dttl) if x == 'longitude'][0]
    lon=copy.deepcopy(CL.data[plon])
    palt=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
    # Not sure why we used to copy only the data part (and not the mask). If someone knows why, feel free to expain. 
    #alt=copy.deepcopy(CL.data[palt].data) #t=copy.defrom pylab import *epcopy(CL.data[pt].data)
    alt=copy.deepcopy(CL.data[palt])
    pt=[i for i,x in enumerate(CL.dttl) if x == 'time'][0]
    t=copy.deepcopy(CL.data[pt])

    #FIND THE QUADRANT
    if st.nanmean(lat)>0:
        quad_NS = 'N'
    else:
        quad_NS = 'S'

    if st.nanmean(lon)>0:
        quad_EW = 'E'
    else:
        quad_EW = 'W'

    # This piece used to take out anomalously deviant altitude data. We now decided to leave it to the user to go and mask the bad data.
    #M=runstats(alt,20)
    #alt=np.ma.masked_where((alt>(M[0]+M[1]*1.)+(isnan(alt))), alt)
        
    norm = matplotlib.colors.Normalize(vmin=t[0],vmax=t[-1])

    fig = figure()
    ax = Axes3D(fig)
    majorFormatter_lon = FormatStrFormatter('%.2f '+quad_EW)
    majorFormatter_lat = FormatStrFormatter('%.2f '+quad_NS)
    try:
        if int(matplotlib.__version__[0])>0:
            ax.scatter(abs(lat),abs(lon),alt,lw=0,alpha=1,cmap='spectral',norm=norm,c=t)
            ax.view_init(28,145)
            ax.yaxis.set_major_formatter(majorFormatter_lon)
            ax.xaxis.set_major_formatter(majorFormatter_lat)
            if quad_EW == 'E':
                ax.set_ylim(ax.get_ylim()[::-1])
            if quad_NS == 'S':
                ax.set_xlim(ax.get_xlim()[::-1])
        else:       # old version of matplotlib that doesn't support the color tracker
            ax.scatter(abs(lat),abs(lon),alt,lw=0,alpha=1,cmap='spectral',norm=norm)
            ax.view_init(28,145)
            ax.yaxis.set_major_formatter(majorFormatter_lon)
            ax.xaxis.set_major_formatter(majorFormatter_lat)
            if quad_EW == 'E':
                ax.set_ylim(ax.get_ylim()[::-1])
            if quad_NS == 'S':
                ax.set_xlim(ax.get_xlim()[::-1])
    except: print("[path3d] Error evaluating your version of matplotlib.")
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Longitude')
    ax.set_zlabel('Altitude')

    plt.ion()
    plt.show()
