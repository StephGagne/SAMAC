# top level init file

# importing most needed libraries (this way it is imported once at the beginning)
# imports in modules only import the namespace and not the packages (it will save computing time)

from pylab import *     # syntax not recommended consider changing for import pylab (and check if all names/functions still work)
import copy as copy
import os

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import math
import xlrd
import scipy.stats.stats as st
from scipy import interpolate

### SAMAC version ###
__version__="1.0.0" #
### SAMAC version ###

#############################
## importing samac package ##
#############################

# importing size distribution tools
from calc.sizedist.sdcnvrt import dNdlogDp2N         # when imported directly (like here) can be reached as samac.dNdlogDp2N
from calc.sizedist.sdcnvrt import dNdlogDp2dN
from calc.sizedist.sdcnvrt import dN2dNdlogDp
from calc.sizedist.sdcnvrt import dNdDp2dNdlogDp
from calc.sizedist.avsizedist import avsizedist
from calc.sizedist.avCDNC import avCDNC
from calc.sizedist.effrad import effrad
from calc.sizedist.filler import filler

# importing lwc-related tools
from calc.lwc.adlwc import AdiabLWC
from calc.lwc.adlwc import adlwc
from calc.lwc.lwp import lwp

# importing precipitation-related tools
from calc.driz.precipTF import precipblwTF
from calc.driz.precipTF import precipblwvpTF
from calc.driz.precipTF import precipincloudTF
from calc.driz.dalprecip import totalprecip
from calc.driz.dalprecip import totalvolprecip
from calc.driz.dalprecip import precipconc

# importing display tools
from disp.geo import maplwc
from disp.geo import mapcloud
from disp.geo import path3d
from disp.overview import overview
from disp.plotsd import plotsd
from disp.vprof import vprof
from disp.wholeprof import wholeprof
from disp.plotavsd import plotavsd
from disp.plotallavsd import plotallavsd
from disp.hrzplot import hrzplot
from disp.plotcloud import plotcloud
from disp.plotcloud import plotcloud2

# importing general tools
from writeouthdf5 import writeouthdf5
from calc.gendata.vpinfo import vpinfo
from calc.gendata.cinfo import cinfo
from calc.gendata.angles import angles
from calc.runstats import runstats
from todatenum import todatenum

# unclassified calc?
from belowprof import belowprof 
from calc.gendata.turbhrz import turbhrz
from calc.gendata.thickness import thickness


# will need to be moved in modules where needed:
# np.seterr(divide='ignore',invalid='ignore')        # ignoring warning of divisions by zero, or invalid values (such as NaNs)

