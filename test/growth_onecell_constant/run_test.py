import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
import math
from scipy.special import iv
#from scipy.special import erfc
from math import erf
from math import erfc
from mpmath import *



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
    
execute=False
num_threads=8
flags="SETUP=growth_one_cell NDUST=200 ENERGY=0 MHD=0 GRAVITY=0 SPHERE=0 NHBINS=1 OPENMP=1 NGHOST=0"
#os.system("rm *png *pdf" )
if(execute):
    os.system("rm shark" )
    os.system("rm -rf output*" )
    os.chdir("../../bin")
if(execute):
    os.system("Make clean && Make "+flags+"&& cp shark ../test/growth_onecell_constant/. && Make clean " )
    os.chdir("../test/growth_onecell_constant/")    
    os.system("export OMP_NUM_THREADS="+str(num_threads)+"; ./shark growth.nml")

iouts=np.arange(2,5)

colors=['black','blue','navy','mediumblue','darkslateblue','slateblue','blueviolet','purple']
def sol_constant(m1,m2,t):
    f1= 4./(2.0+t)**2
    f2= t/(2.0+t)
    x=np.sqrt(m1*m2)
    return(m2-m1)*f1*x*np.exp(-(1.0-f2)*x)
    #return(m2-m1)*x*np.exp(-x)

for i in range(len(iouts)):
    iout=iouts[i]
    print(iout)
    color=colors[i]
    mydata = dune.SharkData(iout)
    mminus=getattr(mydata,"mminus")
    mplus=getattr(mydata,"mplus")
    m=np.sqrt(mminus*mplus)
    sol=np.zeros(len(m))
    print (getattr(mydata,"time"))
    for ii in range(len(sol)):
        sol[ii]=sol_constant(mminus[ii],mplus[ii],getattr(mydata,"time"))
    dune.plot_distribution_test(mydata,ylim=[1e-10,1],show=False,color=color,linestyle='none',marker='o',title='')
    plt.loglog(m,sol,color=color, label='t = '+str(int(getattr(mydata,"time"))) )
#plt.xscale("linear")
#plt.yscale("linear")

plt.ylim(1e-7,1.0)
plt.xlim(10,1e7)

plt.legend()
#plt.show()

plt.savefig("distrib_dust_const.png")
plt.close("all")
