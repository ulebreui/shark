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
    
execute=True
rerun=True
num_threads=8
rundir=os.getcwd()
print(rundir)
nbins=[50,100,150,200,300,400,500,1000]
for nbin in nbins:
    if(execute):
        rundir1=os.getcwd()+'/'+str(nbin)+'bins'
        if(not os.path.exists(rundir1)):
            os.mkdir(rundir1)
        if(not rerun):
            if(os.path.exists(rundir1+'/output_'+str(2).zfill(5))):
                continue
        bindir='/Users/ul264359/Documents/codes/dev/shark_master/bin'
        flags="SETUP=growth_one_cell_phy NDUST="+str(nbin)+" ENERGY=0 MHD=0 GRAVITY=0 SPHERE=0 NHBINS=2 OPENMP=0 NGHOST=0"
        os.chdir(bindir)
        os.system("Make clean && Make "+flags+"&& cp shark "+str(rundir1)+" && Make clean " )
        os.chdir(rundir1)
        os.system("rm -rf output* *pdf" )
        os.system("export OMP_NUM_THREADS="+str(num_threads)+"; ./shark ../growth.nml")
        os.chdir(rundir)
iouts=[2]

for nbin in nbins:
    for i in range(len(iouts)):
        iout=iouts[i]
        mydata = dune.SharkData(iout,path=str(nbin)+'bins/')
        dune.plot_distribution(mydata,ir=0,show=False,normalise=True,ylim=[0.1,1e9],label='$\mathcal{N}=$'+str(nbin))#,color=color)
plt.legend()
plt.show()
