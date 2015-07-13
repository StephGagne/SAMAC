

import matplotlib.pyplot as plt
import samac

#########################################################################        
#############################   plotavsd   ##############################
#########################################################################            
def plotavsd(CL,prof='belowcloud',scan=0,inst='PCASP',filler=0):
    """ This method plots the average size distribution for a given instrument and scan type and number. 
    Example: CloudOjb.plotavsd(prof='abovecloud',num=1,inst='PCASP') will plot the average size distribution for the above cloud #1 for the PCASP in the current figure window.
    The defaults are prof='belowcloud',num=0,inst='PCASP'
    The Special profile 'belowvert' can also be used"""
    asd=samac.avsizedist(CL,prof=prof,scan=scan,inst=inst,filler=filler)
    if len(asd)==0: print("[plotavsd] no data available")
    else:
        if sum(asd[0])==0: plt.plot(asd[1],asd[0],'b-*')
        else: plt.loglog(asd[1],asd[0],'b-*')
        plt.xlabel("Diameter (um)")
        unitsname=[CL.sd[i]["units"] for i,x in enumerate(CL.sd) if x["Distname"].upper()==inst.upper()]
        plt.ylabel('concentration [dN/dlogDp](%s)' %(unitsname[0]))
        plt.show()
    

