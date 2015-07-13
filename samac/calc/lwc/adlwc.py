
from pylab import *
import numpy as np

from scipy import interpolate

import samac

#########################################################################        
################################  adlwc  ################################
#########################################################################                
def adlwc(CL,base='bg'):
    """ This method returns the adiabatic LWC based on the temperature and pressure for each vertical profile. It returns the altitude and the LWC in an array for each vertical scan. The arrays are then returned in a list. This property fully handles the extradata module.
    Example: vertical scan #1's Altitude would be R[1][0] and the Adiabatic LWC R[1][1]
    base: base=4point: to use base height selection between 4-point height average of the lower and upper estimates for the base height (defheight).
          base=bg: to use the best guess base height (defBGheight))"""
    try:
        Ti=[i for i,x in enumerate(CL.dttl) if x == 'time'][0]
        A=[i for i,x in enumerate(CL.dttl) if x == 'altitude'][0]
        t=CL.data[Ti]; # time
        alt=CL.data[A]    # altitude
        # Pressure
        P=[i for i,x in enumerate(CL.dttl) if x.lower() == 'pressure']
        if len(P)==1: P=P[0]; Pd=CL.data[P]; Pt=t; 
        elif len(P)>1: print("[adlwc] Multiple Pressure found in the basic data."); crash;
        else:       # looking for pressure in the extradata
            posx=[] 
            for i,ttl in enumerate(CL.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'pressure']    # check all titles matching with pressure
            if len(posx)==1: 
                Pd=CL.extradata[posx[0][0]][posx[0][1]]    # loading the pressure data
                j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
                Pt=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[adlwc] No Pressure found in the basic or the extra data."); crash
        # Temperature         
        Te=[i for i,x in enumerate(CL.dttl) if x.lower() == 'temperature']
        if len(Te)==1: Te=Te[0]; Ted=CL.data[Te]; Tet=t; 
        elif len(Te)>1: print("[adlwc] Multiple temperature found in the basic data."); crash;
        else:       # looking for pressure in the extradata
            posx=[] 
            for i,ttl in enumerate(CL.extrattl):     # for all extra datasets available
                posx=posx+[[i,j] for j,x in enumerate(ttl) if x.lower() == 'temperature']    # check all titles matching with temperature
            if len(posx)==1: 
                Ted=CL.extradata[posx[0][0]][posx[0][1]]    # loading the pressure data
                j=[j for j,x in enumerate(CL.extrattl[i]) if x.lower() == 'time'][0]
                Tet=CL.extradata[posx[0][0]][j]     # loading associated time stamp
            else: print("[adlwc] No temperature found in the basic or the extra data."); crash
    except: 
        print("[adlwc] One of pressure, altitude, time or temperature are missing, multiple or in the wrong format.");            
        return []
    # adapting for too short data for interpolation
    if len(Pt)<2: Pd=np.ones((2,))*NaN; Pt=np.array([CL.times["cloud"][0][0],CL.times["cloud"][0][1]]);
    if len(Tet)<2: Ted=np.ones((2,))*NaN; Tet=np.array([CL.times["cloud"][0][0],CL.times["cloud"][0][1]]);
    fP=interpolate.interp1d(Pt,Pd,kind='linear')   # f=interp(x=time of pressure, y=pressure, linear interpolation)
    fT=interpolate.interp1d(Tet,Ted,kind='linear')   # f=interp(x=time of temperature, y=temperature, linear interpolation)
    # adapting time vectors to a unique interpolable time vector.
    ta1=np.max([t[0],Pt[0],Tet[0]]); ta2=np.min([t[-1],Pt[-1],Tet[-1]]);
    ta=t[nonzero((t>=ta1)*(t<=ta2))[0]];
    alt=alt[nonzero((t>=ta1)*(t<=ta2))[0]];
    # calculations
    LWCad=[]
    dlwcdz=samac.AdiabLWC(fT(ta),fP(ta))[1]
    for i,j in enumerate(CL.times["verticloud"]):
        lwci=[]
        if base=='4point':
            cb=CL.props["height"][i][0]    # lower cloud base estimate
            ct=CL.props["height"][i][3]    # upper cloud top estimate
            if isnan(ct): ct=np.max(alt[np.nonzero((ta>=j[0])*(ta<=j[1]))])
        else:
            cb=CL.props["BGheight"][i][0]    # best guess cloud base estimate
            ct=CL.props["BGheight"][i][1]    # best guess cloud top estimate
            if isnan(ct): ct=np.max(alt[np.nonzero((ta>=j[0])*(ta<=j[1]))])
        ix=np.nonzero((ta>=j[0])*(ta<=j[1])*(alt>=cb)*(alt<=ct))
        alti=alt[ix]; dalt=np.diff(alti);
        dlwci=dlwcdz[ix];   
        dlwci=dlwci[0:-1]; alti=alti[0:-1] # cutting off the last digit
        for k in range(len(dlwci)):
            ixk=np.nonzero((alti>=cb)*(alti<=alti[k]))
            lwci.append(abs(sum(dlwci[ixk]*dalt[ixk])))
        lwci=np.array(lwci)
        LWCad.append(np.array([alti,lwci]))
    return LWCad  


#########################################################################        
#############################  AdiabLWC  ################################
#########################################################################
def AdiabLWC(Temp,Pres):
    """Adiabatic lwc calculated from the temperature in Celcius and pressure in mb or hPa.
       Returns 1-Pseudoadiabatic (moist adiabatic) lapse rate [K/m], 2-height dependence of LWC (dLWCdz)."""
    # Code by Mike Earle (Env. Can.), based on Tietze et al., ACP, 2011. Adapted to Python by Stephanie Gagne
    # Constants
    epsilon = 0.622;
    R_air = 287;  # Individual gas constant for dry air in J/(kg K) 
    # Saturation vapour pressure of water [hPa]
    #p_sat = 6.1121*exp((18.678 - Temp/234.50) * Temp/(275.14 + Temp));
    p_sat = 6.112 * np.exp(17.67 * Temp / (Temp + 243.5));
    #p_sat = 2.22;
    # Saturation mixing ratio 
    w_sat = epsilon * p_sat/(Pres - p_sat);
    # Heat capacity of dry air [J/(kg K)]
    c_p = 1005;
    # Air density at a given T, P [kg/m^3]
    rho_air = (Pres*100)/(287.05 * (Temp + 273.15));
    # Heat capacity of moist air [J/(kg K)]
    c_pm = 1870;
    #c_pm = c_p * (1 + 0.9 * w_sat);
    # Heat capacity of pure water [IT cal/(g C)] 
    c_w_params = [1.000938, -2.7052e-3, -2.3235e-5, 4.3778e-6, 2.7136e-7];
    c_w = 0;
    for j in range(5):
      c_w_i = c_w_params[j] * Temp**j;
      c_w = c_w + c_w_i;
    # Convert c_w to units of J/(kg K)
    c_w = c_w * 4.1868 * 1000;
    # Calculate latent heat of vaporization for pure water
    T_0 = 0;   
    L_0 = 2.50e6;   # Latent heat of vap at T_0 in J/kg
    L = L_0 - (c_w - c_pm)*(Temp - T_0); 
    # Dry adiabatic lapse rate [K/m]
    lr_dry = 9.76e-3;
    # Pseudoadiabatic (moist adiabatic) lapse rate [K/m]
    lr_moist = lr_dry * ((1 + L * w_sat/(R_air * (Temp + 273.15)))/(1 + L**2 * epsilon * w_sat/(R_air * c_p * (Temp + 273.15)**2)));
    # Compute height dependence of LWC [g/m^4]
    dLWCdz = rho_air * c_p / L * (lr_dry - lr_moist) * 1e3;
    return [lr_moist,dLWCdz]

