import os
import dune
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import matplotlib as mpl
import glob
import lic
import matplotlib
from scipy.io import FortranFile
import struct
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

def get_binary_data(fmt="",ninteg=0,nlines=0,nfloat=0,nstrin=0,nquadr=0,nlongi=0,content=None,correction=0):
    offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16 + 4 + correction
    byte_size = {"i":4,"d":8,"q":8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]
    return struct.unpack(fmt, content[offset:offset+pack_size])

resolution=128
NX=resolution
NY=resolution

def read_var(filename,varname,size):
    f=open(filename+"/"+varname, "rb")
    dat = np.fromfile(f, dtype=np.float, count=size, sep='')
    #dat = np.loadtxt(filename+"/"+varname)#
    return dat
execute=True
num_threads=8


NDUST=50
flags="SETUP=KH NDUST="+str(NDUST)+" SPHERE=0 NX="+str(NX)+" NY="+str(NY)+" NGHOST=2 OPENMP=1 DEBUG=0 SOLVER=2"
print (flags)
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


x   = read_var(filename,'x',NX*NY)
y   = read_var(filename,'y',NX*NY)
v   = read_var(filename,'v',NX*NY)
vy  = read_var(filename,'vy',NX*NY)
P   = read_var(filename,'P',NX*NY)
rho = read_var(filename,'rho',NX*NY)
if(NDUST>0):
    rhod  = read_var(filename,'rhod',NX*NY*NDUST)
    vd    = read_var(filename,'vd',NX*NY*NDUST)
    vdy   = read_var(filename,'vdy',NX*NY*NDUST)

rho_2D = np.reshape(rho,(NY,NX),order = "C").T

xx=np.linspace(0,1,NX)
yy=np.linspace(0,1,NY)

plt.contourf(xx,yy,rho_2D,levels=100)
plt.colorbar(label='Density [code unit]')
plt.xlabel("x")
plt.ylabel("y")

plt.show()
if(NDUST>0):
    eps_2D = np.reshape(rhod,(NDUST,NY,NX),order = "C").T
    plt.contourf(xx,yy,eps_2D[:,:,0]/rho_2D,levels=100)
    plt.colorbar(label='Dust-to-gas ratio')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


"""
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
"""

