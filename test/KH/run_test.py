import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
import glob
import lic
import matplotlib

def get_lic(x, y, l=30):
    lic_res = lic.lic(data_x=x, data_y=y, contrast=False, length=l)
    #amplify contrast on lic
    lim=(.2,.5)
    lic_cmap = "binary"
    lic_data_clip = np.clip(lic_res,lim[0],lim[1])
    lic_data_rgba = matplotlib.cm.ScalarMappable(norm=None, cmap=lic_cmap).to_rgba(lic_data_clip)
    lic_data_clip_rescale = (lic_data_clip-lim[0])/(lim[1]-lim[0])
    lic_data_rgba[...,3] = lic_data_clip_rescale * 1
    return lic_data_rgba

execute=True
num_threads=8

resolution=1024
NX=resolution
NY=resolution
flags="SETUP=KH NDUST=0 SPHERE=0 NX="+str(NX)+" NY="+str(NY)+" NGHOST=2 OPENMP=1 DEBUG=0"
#os.system("rm *png *pdf" )
table=' /Users/ul264359/Documents/codes/dev/shark_master/tables/ormel_brown_no_frag_nh1d5_T10K_100bins.dat'
if(execute):
    os.system("rm shark" )
    os.system("rm -rf output*" )
    os.chdir("../../bin")
if(execute):
    os.system("Make clean && Make "+flags+"&& cp shark ../test/KH/. && Make clean " )
    os.chdir("../test/KH/")    
    os.system("export OMP_NUM_THREADS="+str(num_threads)+"; ./shark kh.nml"+table)


filelist = sorted(glob.glob("output*"))
number = filelist[-1].split("_")[-1]
filename= "output_"+ str(number).zfill(5)

x    = np.loadtxt(filename+"/x.dat")
y    = np.loadtxt(filename+"/y.dat")

rho  = np.loadtxt(filename+"/rho.dat")
v    = np.loadtxt(filename+"/v.dat")
vy   = np.loadtxt(filename+"/vy.dat")
P    = np.loadtxt(filename+"/P.dat")

rho_2D = np.reshape(rho,(NY,NX),order = "C").T

xx=np.linspace(0,1,NX)
yy=np.linspace(0,1,NY)

plt.contourf(xx,yy,rho_2D,levels=100)
plt.colorbar(label='Density [code unit]')
plt.xlabel("x")
plt.ylabel("y")

plt.show()

VX = np.reshape(v,(NY,NX),order = "C").T
VY = np.reshape(vy,(NY,NX),order = "C").T

lic_res = get_lic(VX, VY, l=30)

fig,ax = plt.subplots(figsize=(8,8))
im = ax.imshow(rho_2D, origin="lower", alpha=1, cmap="inferno")
im2 = ax.imshow(lic_res, origin="lower", cmap="binary", alpha=.3)
fig.colorbar(im, label='Density [code unit]')
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.show()


#plt.legend()
#plt.savefig("density.png")
#plt.close("all")