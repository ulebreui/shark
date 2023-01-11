import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
import math
#from scipy.special import erf
#from scipy.special import erfc
from math import erf
from math import erfc
plt.rcParams.update({'font.size': 13})

def round_to_n(x, n):
    " Round x to n significant figures "
    return round(x, -int(py.floor(np.sign(x) * py.log10(abs(x)))) + n)
def str_fmt(x, n=2):
    " Format x into nice Latex rounding to n"
    power = int(py.log10(round_to_n(x, 0)))
    f_SF = round_to_n(x, n) * pow(10, -power)
    if(x<0.1 or x>1000.):
        return "{}\,10^{}".format(f_SF, power)
    else:
        return str(round(x,2))

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
    

iouts=[2]
nbins=[50,100,200,300]
colors=['steelblue','coral','olive','purple','red','orange']
for ibin in range(len(nbins)):
    nbin=nbins[ibin]
    iout=iouts[0]
    #color=colors[i]
    mydata = dune.SharkData(iout,path=str(nbin)+'bins/')
    dune.plot_distribution(mydata,ir=0,show=False,normalise=True,ylim=[1e-7,1e2],label='$\mathcal{N}=$'+str(nbin),color=colors[ibin])
plt.legend()
plt.title('Convergence study for the dust coagulation')
#plt.show()
plt.savefig("distrib_dust_convergence.png")
#plt.close("all")
