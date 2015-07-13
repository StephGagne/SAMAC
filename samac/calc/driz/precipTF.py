
from pylab import *
import numpy as np

import samac


#########################################################################        
############################  precipblwTF  ##############################
#########################################################################
def precipblwTF(CL):
    """ This property tells whether precipitation below cloud was observed during the entire cloud period: true or false for >100 & >200 microns.
        The dictionary entries can be accessed through CloudObj.precipblwTF["100"] or CloudObj.precipblwTF["200"]. """
    pos=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]       
    if len(pos)==0: pos=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()]
    if len(pos)==0: raise ValueError("[precipblwTF] No precipitation distribution found.")
    elif len(pos)==1: pos=int(pos[0])
    elif len(pos)>1:
        for p in pos:
            print("%d - %s" % (p, CL.sd[p]["Distname"]))
        pos=int(raw_input("What is the prefered precipitation distribution? Enter the number."))
    else: pos=nan   
    if type(CL.sd[pos]["time"])==float or pos==nan: 
        precb=dict(); precb["100"]=list(); precb["200"]=list();
        precb["100"].append(nan); precb["200"].append(nan);
    elif np.shape(CL.sd[pos]["time"])[0]<=1: 
        precb=dict(); precb["100"]=list(); precb["200"]=list();
        precb["100"].append(nan); precb["200"].append(nan);
    else:
        P100=samac.dNdlogDp2N(CL.sd[pos],100,nan)
        #P100=np.ma.masked_less_equal(P100,1e-2);
        P100[P100<=1e-2]=0
        P200=samac.dNdlogDp2N(CL.sd[pos],200,nan)
        P200[P200<=1e-2]=0
        #P200=np.ma.masked_less_equal(P200,1e-2);
        t100=0; t200=0; precb=dict();
        for b in range(np.shape(CL.times['belowcloud'])[0]):
            if len(P100[(CL.sd[pos]["time"]<=CL.times['belowcloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['belowcloud'][b][0])])==0: pass
            else:
                t100=t100+np.nansum(P100[(CL.sd[pos]["time"]<=CL.times['belowcloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['belowcloud'][b][0])])
                t200=t200+np.nansum(P200[(CL.sd[pos]["time"]<=CL.times['belowcloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['belowcloud'][b][0])])    
        if t100>0: precb["100"]=True
        else: precb["100"]=False
        if t200>0: precb["200"]=True
        else: precb["200"]=False
    return precb



#########################################################################        
###########################  precipblwvpTF  #############################
#########################################################################
def precipblwvpTF(CL):
    """ This property tells whether precipitation below cloud was observed for a given vertical scan and below that scan: true or false for >100 & >200 microns.
        It is important to note that the plane not detecting rain when it passed below the vertical scan doesn't mean that the shaft has not precipitated or will not precipitate.
        The dictionary entries can be accessed through CloudObj.precipblwvpTF["100"] or CloudObj.precipblwvpTF["200"]. """
    pos=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]  
    if len(pos)==0: pos=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()]  
    if len(pos)==1: pos=int(pos[0])
    else: pos=nan; print("[precipblwvpTF] Multiple precipitation distribution found, no priority is set.")
    if type(CL.sd[pos]["time"])==float or pos==nan: 
        precb=dict(); precb["100"]=list(); precb["200"]=list();
        for r in range(len(CL.times["verticloud"])):
            precb["100"].append(nan); precb["200"].append(nan);
    elif np.shape(CL.sd[pos]["time"])[0]<=1: 
        precb=dict(); precb["100"]=list(); precb["200"]=list();
        for r in range(len(CL.times["verticloud"])):
            precb["100"].append(nan); precb["200"].append(nan);
    else:
        precb=dict();
        precb["100"]=list(); precb["200"]=list();
        P100=samac.dNdlogDp2N(CL.sd[pos],100,nan)
        #P100=np.ma.masked_less_equal(P100,1e-2);   # nansum doesn't handle masked arrays.
        P100[P100<=1e-2]=0
        P200=samac.dNdlogDp2N(CL.sd[pos],200,nan)
        P200[P200<=1e-2]=0
        #P200=np.ma.masked_less_equal(P200,1e-2);
        for k in range(len(CL.times["verticloud"])):
            tsix=list()
            G=samac.belowprof(CL)[k]["index"]
            Gt=samac.belowprof(CL)[k]["times"]
            H=np.nonzero(diff(G)>5)[0]
            if np.shape(H)[0]==0: tsix.append([Gt[0],Gt[-1]])
            else:
                for c in range(len(H)+1):
                    if c==0: tsix.append([Gt[0],Gt[H[c]]])
                    elif c<=len(H)-1: tsix.append([Gt[H[c-1]+1],Gt[H[c]]])
                    elif c==len(H): tsix.append([Gt[H[c-1]+1],Gt[-1]])
            t100=0; t200=0;
            for b in range(len(tsix)):
                if len(P100[(CL.sd[pos]["time"]<=tsix[b][1])*(CL.sd[pos]["time"]>=tsix[b][0])])==0: pass
                else:
                    t100=t100+np.nansum(P100[(CL.sd[pos]["time"]<=tsix[b][1])*(CL.sd[pos]["time"]>=tsix[b][0])])
                    t200=t200+np.nansum(P200[(CL.sd[pos]["time"]<=tsix[b][1])*(CL.sd[pos]["time"]>=tsix[b][0])])
            if t100>0: precb["100"].append(True)
            else: precb["100"].append(False)
            if t200>0: precb["200"].append(True)
            else: precb["200"].append(False)
    return precb



#########################################################################        
##########################  precipincloudTF  ############################
#########################################################################
def precipincloudTF(CL,abovesize=100):
    """ This method tells whether precipitation in cloud (horicloud + verticloud) was observed during the entire cloud period: true or false. 
    Default value abovesize=100 for precipitation >100 microns. Noise level at 1e-2."""
    # Do not use to count the total amount of rain detected as some scans may be counted more than once.
    pos=[i for i,x in enumerate(CL.sd) if 'p1' in x["sdtype"].lower()]        # finding P1 size dist
    if len(pos)==0: pos=[i for i,x in enumerate(CL.sd) if 'p' in x["sdtype"].lower()]    # finding P size dist
    if len(pos)==1: pos=int(pos[0])
    elif len(pos)>1:
        for p in pos:
            print("%d - %s" % (p, CL.sd[p]["Distname"]))
        pos=int(raw_input("What is the prefered precipitation distribution? Enter the number."))
    else: pos=nan
    if type(pos)==float:
        if isnan(pos):
            precb=nan   
    elif type(CL.sd[pos]["time"])==float: 
        precb=nan; 
    elif np.shape(CL.sd[pos]["time"])[0]<=1: 
        precb=nan
    else:
        Pas=samac.dNdlogDp2N(CL.sd[pos],abovesize,nan)
        #Pas=np.ma.masked_less_equal(Pas,1e-2);     # nansum doesn't seem to handle masked arrays. 
        Pas[Pas<=1e-2]=0;
        t100=0;
        for b in range(np.shape(CL.times['horicloud'])[0]):
            if len(Pas[(CL.sd[pos]["time"]<=CL.times['horicloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['horicloud'][b][0])])==0: pass
            else: t100=t100+np.nansum(Pas[(CL.sd[pos]["time"]<=CL.times['horicloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['horicloud'][b][0])])
        for b in range(np.shape(CL.times['verticloud'])[0]):
            if len(Pas[(CL.sd[pos]["time"]<=CL.times['verticloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['verticloud'][b][0])])==0: pass
            else: t100=t100+np.nansum(Pas[(CL.sd[pos]["time"]<=CL.times['verticloud'][b][1])*(CL.sd[pos]["time"]>=CL.times['verticloud'][b][0])])
        if t100>0: precb=True
        else: precb=False
    return precb
