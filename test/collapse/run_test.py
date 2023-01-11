import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
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
flags="SETUP=collapse_spher NDUST=10 ENERGY=0 MHD=0 GRAVITY=1 SPHERE=1 NHBINS=128 NGHOST=2 OPENMP=1"
#os.system("rm *png *pdf" )
table=' /Users/ul264359/Documents/codes/dev/shark_master/tables/ormel_brown_no_frag_nh1d5_T10K_100bins.dat'
if(execute):
    os.system("rm shark" )
    os.system("rm -rf output*" )
    os.chdir("../../bin")
if(execute):
    os.system("Make clean && Make "+flags+"&& cp shark ../test/collapse/. && Make clean " )
    os.chdir("../test/collapse/")    
    os.system("export OMP_NUM_THREADS="+str(num_threads)+"; ./shark collapse.nml"+table)
mydata = dune.SharkData(-1)
fig=plt.figure()
if(hasattr(mydata,"ndust")):
    ndust=getattr(mydata,"ndust")
    idust=np.arange(1,ndust+1)
    sdust_all=getattr(mydata,"sdust_all")
    sdust=np.array(np.log10(sdust_all[0,:]))
    colors=pl.cm.plasma(idust)
    norm = mpl.colors.Normalize(vmin=np.amin(sdust), vmax=np.amax(sdust))
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.plasma)
    cmap.set_array([])
    for i in idust:
        lab='Dust, $'+str(latex_float(sdust[i-1]))+'$ cm'
        lab=None
        dune.plot_slice(mydata,"grid","rhod"+str(i),show=False,label=lab,color=colors[i-1])
dune.plot_slice(mydata,"grid","rho",show=False,label='Gas',color='k')
if(hasattr(mydata,"ndust")):
    fig.colorbar(cmap,label='Logarithm of the dust size [cm]')
plt.legend()
plt.savefig("density_dust.png")
plt.close("all")
if(hasattr(mydata,"ndust")):
    fig=plt.figure()
    for i in idust:
        lab='Dust, $'+str(latex_float(sdust[i-1]))+'$ cm'
        lab=None
        dune.plot_slice(mydata,"grid","epsilon"+str(i),show=False,label=lab,color=colors[i-1],ylab='Dust-to-gas ratio $\\epsilon$')
    #plt.legend(ncol=2)
    fig.colorbar(cmap,label='Logarithm of the dust size [cm]')
    plt.savefig("epsilon_dust.png")
    plt.close("all")

    fig=plt.figure()

    dune.plot_distribution(mydata,radius=1,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=5,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=10,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=50,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=100,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=1000,ylim=[1e-5,0.1],show=False)
    dune.plot_distribution(mydata,radius=5000,ylim=[1e-5,0.1],show=False)

    plt.savefig("distrib_dust.png")
    plt.close("all")

    fig=plt.figure
    dune.plot_slice(mydata,"rho","eta_a",show=False,label='$\\eta_a$',color='r')
    dune.plot_slice(mydata,"rho","eta_o",show=False,label='$\\eta_o$',color='g')
    dune.plot_slice(mydata,"rho","eta_h_plus",show=False,label='$\\eta_h$',ylab='Resistivity [s$^{-1}$]',color='b')
    dune.plot_slice(mydata,"rho","eta_h_minus",show=False,label='-$\\eta_h$',ylab='Resistivity [s$^{-1}$]',color='b',linestyle='dotted')

    plt.legend(ncol=4)
    plt.savefig("resistivities.png")
    plt.close("all")

    fig=plt.figure
    dune.plot_slice(mydata,"nh","ni",show=False,label=True)
    dune.plot_slice(mydata,"nh","ne",show=False,label=True)
    plt.savefig("ionis.png")
    plt.close("all")

os.system("rm shark")
#os.system("rm -rf output*")
