
import matplotlib.pyplot as plt
import datetime as dt
import samac

#########################################################################        
############################   plotallavsd   ############################
#########################################################################
def plotallavsd(CL,prof='belowcloud',num=0,interact=None,exc='2d'):
    """ This method plots the average size distribution of all size distribution intruments (by default, except for 2d-type measurements) for a given scan. The maximum number of instruments that can be plotted is 8.
        prof: manoeuvre type (cloud, abovecloud, belowcloud, verticloud, horicloud and belowvert (special profile, see method belowprof). Default is 'belowcloud'.
        Use CloudObj.plotallavsd(0), CloudObj.plotallavsd() or CloudObj.plotallavsd(interact=0) for a plot NOT in interactive mode;   
        Use CloudObj.plotallavsd(1) or CloudObj.plotallavsd(interact=1)for a plot in interactive mode.  
        Use CloudObj.plotallavsd(exc='inst1,inst2') to exclude instruments inst1 and inst2 (includes 2d unless specified in the list).
        exc: exceptions not to be plotted. Exceptions should be stated in a single string, with instrument names (as found in c.sd[n]["Distname"]) separated by commas (e.g. exc='2dc,PCASP,FSSP100'). Default excepts all instruments that contain '2d' in their name.
        Use CloudObj.plotallavsd(exc=None) to exclude no instrument (includes 2d).
        The default is prof='belowcloud',num=0,interact=None,exc='2d'"""
    if interact==None or interact==0: interact=0
    elif interact==1: plt.ion()
    SD=list()
    specs=['b-*','r-o','c-v','m-d','g-s','k-^','b--o','r--*']
    instlist=[x["Distname"] for i,x in enumerate(CL.sd)]
    unitlist=[x["units"] for i,x in enumerate(CL.sd)]
    alllist=[[x,unitlist[i]] for i,x in enumerate(instlist)]
    if exc==None: pass
    else:
        X=exc.split(',');
        for c in range(len(X)):
            X[c]=X[c].strip()
        for c in range(len(X)):
            alllist=[x for i,x in enumerate(alllist) if X[c].lower() not in x[0].lower()]
            #alllist=[x for i,x in enumerate(alllist) if X[c].lower()!=x[0].lower()]
    for i in range(len(alllist)):
        SD.append(samac.avsizedist(CL,prof=prof,scan=num,inst=alllist[i][0]))
        plt.loglog(SD[i][1],SD[i][0],specs[i])
    instunit=list()
    for i in range(len(alllist)):
        instunit.append(alllist[i][0]+' ('+alllist[i][1]+')')
    plt.legend(instunit,loc=3)
    plt.xlabel('diameter (um)')
    plt.ylabel('concentration [dN/dlogDp]')
    post=[i for i,x in enumerate(CL.dttl) if x == 'time'] # position of the time
    plt.title("Average size dist. for %s no. %d on %s @ %s \n (%s)" % (prof, num, CL.desc["date"], (dt.datetime(1,1,1)+dt.timedelta(days=CL.data[post][0][0])).strftime('%H:%M'), CL.desc["humanplace"]))
    if interact==0: plt.show()




