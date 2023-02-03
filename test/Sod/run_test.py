import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
import glob

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
num_threads=8
NX=128
NY=128
flags="SETUP=SOD NDUST=0 SPHERE=0 NX="+str(NX)+" NY="+str(NY)+" NGHOST=2 OPENMP=0 DEBUG=1 SOLVER=2"
#os.system("rm *png *pdf" )
table=' /Users/ul264359/Documents/codes/dev/shark_master/tables/ormel_brown_no_frag_nh1d5_T10K_100bins.dat'
if(execute):
    os.system("rm shark" )
    os.system("rm -rf output*" )
    os.chdir("../../bin")
if(execute):
    os.system("Make clean && Make "+flags+"&& cp shark ../test/Sod/. && Make clean " )
    os.chdir("../test/Sod/")    
    os.system("export OMP_NUM_THREADS="+str(num_threads)+"; ./shark sod.nml"+table)

filelist = sorted(glob.glob("output*"))
number = filelist[-1].split("_")[-1]
filename= "output_"+ str(number).zfill(5)

x    = np.loadtxt(filename+"/x.dat")
y    = np.loadtxt(filename+"/y.dat")

rho  = np.loadtxt(filename+"/rho.dat")
v    = np.loadtxt(filename+"/v.dat")
P    = np.loadtxt(filename+"/P.dat")
print(P)
plt.plot(y,rho,label='Gas',color='r',marker='o',ls='none')
plt.show()
plt.plot(y,P,label='Gas',color='r',marker='o',ls='none')
plt.show()
plt.plot(y,v,label='Gas',color='r',marker='o',ls='none')
plt.show()

rho_2D= np.reshape(rho,(NY,NX),order = "C").T

xx=np.linspace(0,1,NX)
yy=np.linspace(0,1,NY)

plt.contourf(xx,yy,rho_2D,levels=100)
plt.colorbar(label='Density [code unit]')
plt.xlabel("x")
plt.ylabel("y")

plt.show()

#plt.legend()
#plt.savefig("density.png")
#plt.close("all")