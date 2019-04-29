# -*- coding: utf-8 -*-
"""
Created on Sat Sep  3 11:57:22 2016

@author: florian
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:11:37 2016

@author: florian
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 10:04:54 2016

@author: florian
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
#

#Initialises an isothermal atmosphere at rest 
# with physical values c

out_dir = '/Users/florian/Desktop/test/ECLIPS3D/trunk/2D_axi/data/'


Nlatf=40
Nzf=20





r=6372000.0
Heig=80000.0
nz=100
nlat=100



dphi=np.pi/nlat
dz=Heig/nz

x_f=np.linspace(-90,90,nlat+1)*np.pi/180
z_f=np.linspace(0,Heig,nz+1)


u=np.zeros((nz+1,nlat+1))
v=np.zeros((nz+1,nlat+1))
w=np.zeros((nz+1,nlat+1))

T=250.
R=287.05
cp=1005.0
p0=1.0E5
g=9.81
omega=7.292E-5

g1=np.zeros(nz+1)
for i in range(nz+1) :
    g1[i]=g*r**2/(r+z_f[i])**2

p=np.zeros((nz+1,nlat+1))
theta=np.zeros((nz+1,nlat+1))
rho=np.zeros((nz+1,nlat+1))

for i in range(nz+1) :
    p[i]=p0*np.exp(-g1[i]*z_f[i]/(R*T))
    
#for i in range(nz+1) :
#    u[i]=2.0*omega*r*np.exp(-(x_f[:])**2/(2*0.2**2))

rho=p/(R*T)
theta=(p0/p)**(R/cp)*T

du=np.gradient(u,dz,dphi)
dv=np.gradient(v,dz,dphi)
dw=np.gradient(w,dz,dphi)
dp=np.gradient(p,dz,dphi)
dtheta=np.gradient(theta,dz,dphi)
drho=np.gradient(rho,dz,dphi)
#
#
u_phi=du[1]
v_phi=dv[1]
w_phi=dw[1]
p_phi=dp[1]
rho_phi=drho[1]
theta_phi=dtheta[1]

u_z=du[0]
v_z=dv[0]
w_z=dw[0]
p_z=dp[0]
rho_z=drho[0]
theta_z=dtheta[0]




dphi=0.5*np.pi/(Nlatf)
dz=Heig/(Nzf)

XU=(np.linspace(0,Nlatf-1,Nlatf)+.5)*dphi
XV=np.linspace(0,Nlatf,Nlatf+1)*dphi
ZU=(np.linspace(0,Nzf-1,Nzf)+0.5)*dz
ZW=(np.linspace(0,Nzf,Nzf+1))*dz




#du=np.gradient(u,dz,dphi)
#dv=np.gradient(v,dz,dphi)
#dp=np.gradient(p,dz,dphi)
#dtheta=np.gradient(theta,dz,dphi)
#drho=np.gradient(rho,dz,dphi)





fu=interpolate.interp2d(x_f,z_f,u)
fv=interpolate.interp2d(x_f,z_f,v)
fp=interpolate.interp2d(x_f,z_f,p)
fw=interpolate.interp2d(x_f,z_f,w)
ftheta=interpolate.interp2d(x_f,z_f,theta)
frho=interpolate.interp2d(x_f,z_f,rho)

fdu_dz=interpolate.interp2d(x_f,z_f,u_z)
fdv_dz=interpolate.interp2d(x_f,z_f,v_z)
fdp_dz=interpolate.interp2d(x_f,z_f,p_z)
fdw_dz=interpolate.interp2d(x_f,z_f,w)
fdtheta_dz=interpolate.interp2d(x_f,z_f,theta_z)
fdrho_dz=interpolate.interp2d(x_f,z_f,rho_z)

fdu_dphi=interpolate.interp2d(x_f,z_f,u_phi)
fdv_dphi=interpolate.interp2d(x_f,z_f,v_phi)
fdp_dphi=interpolate.interp2d(x_f,z_f,p_phi)
fdw_dphi=interpolate.interp2d(x_f,z_f,w)
fdtheta_dphi=interpolate.interp2d(x_f,z_f,theta_phi)
fdrho_dphi=interpolate.interp2d(x_f,z_f,rho_phi)



u_f=fu(XU,ZU)
u_hf=fu(XU,ZW)
u_lf=fu(XV,ZU)

v_f=fv(XU,ZU)
v_hf=fv(XU,ZW)
v_lf=fv(XV,ZU)

p_f=fp(XU,ZU)
p_hf=fp(XU,ZW)
p_lf=fp(XV,ZU)

theta_f=ftheta(XU,ZU)
theta_hf=ftheta(XU,ZW)
theta_lf=ftheta(XV,ZU)

rho_f=frho(XU,ZU)
rho_hf=frho(XU,ZW)
rho_lf=frho(XV,ZU)

w_f=fw(XU,ZU)
w_hf=fw(XU,ZW)
w_lf=fw(XV,ZU)





du_dz_f=fdu_dz(XU,ZU)
du_dz_hf=fdu_dz(XU,ZW)
du_dz_lf=fdu_dz(XV,ZU)

dv_dz_f=fdv_dz(XU,ZU)
dv_dz_hf=fdv_dz(XU,ZW)
dv_dz_lf=fdv_dz(XV,ZU)

dp_dz_f=fdp_dz(XU,ZU)
dp_dz_hf=fdp_dz(XU,ZW)
dp_dz_lf=fdp_dz(XV,ZU)

dtheta_dz_f=fdtheta_dz(XU,ZU)
dtheta_dz_hf=fdtheta_dz(XU,ZW)
dtheta_dz_lf=fdtheta_dz(XV,ZU)

drho_dz_f=fdrho_dz(XU,ZU)
drho_dz_hf=fdrho_dz(XU,ZW)
drho_dz_lf=fdrho_dz(XV,ZU)

dw_dz_f=fdw_dz(XU,ZU)
dw_dz_hf=fdw_dz(XU,ZW)
dw_dz_lf=fdw_dz(XV,ZU)




du_dphi_f=fdu_dphi(XU,ZU)
du_dphi_hf=fdu_dphi(XU,ZW)
du_dphi_lf=fdu_dphi(XV,ZU)

dv_dphi_f=fdv_dphi(XU,ZU)
dv_dphi_hf=fdv_dphi(XU,ZW)
dv_dphi_lf=fdv_dphi(XV,ZU)

dp_dphi_f=fdp_dphi(XU,ZU)
dp_dphi_hf=fdp_dphi(XU,ZW)
dp_dphi_lf=fdp_dphi(XV,ZU)

dtheta_dphi_f=fdtheta_dphi(XU,ZU)
dtheta_dphi_hf=fdtheta_dphi(XU,ZW)
dtheta_dphi_lf=fdtheta_dphi(XV,ZU)

drho_dphi_f=fdrho_dphi(XU,ZU)
drho_dphi_hf=fdrho_dphi(XU,ZW)
drho_dphi_lf=fdrho_dphi(XV,ZU)

dw_dphi_f=fdw_dphi(XU,ZU)
dw_dphi_hf=fdw_dphi(XU,ZW)
dw_dphi_lf=fdw_dphi(XV,ZU)






res=np.append(u_f.T,u_hf.T)
res=np.append(res,u_lf.T)


res=np.append(res,v_f.T)
res=np.append(res,v_hf.T)
res=np.append(res,v_lf.T)

res=np.append(res,p_f.T)
res=np.append(res,p_hf.T)
res=np.append(res,p_lf.T)

res=np.append(res,theta_f.T)
res=np.append(res,theta_hf.T)
res=np.append(res,theta_lf.T)

res=np.append(res,rho_f.T)
res=np.append(res,rho_hf.T)
res=np.append(res,rho_lf.T)

res=np.append(res,w_f.T)
res=np.append(res,w_hf.T)
res=np.append(res,w_lf.T)



dres=np.append(du_dphi_f.T,du_dz_f.T)

dres=np.append(dres,dv_dphi_f.T)
dres=np.append(dres,dv_dphi_lf.T)
dres=np.append(dres,dv_dz_lf.T)

dres=np.append(dres,dp_dphi_f.T)
dres=np.append(dres,dp_dphi_lf.T)
dres=np.append(dres,dp_dz_hf.T)

dres=np.append(dres,dtheta_dz_f.T)
dres=np.append(dres,dtheta_dphi_hf.T)
dres=np.append(dres,dtheta_dz_hf.T)

dres=np.append(dres,drho_dphi_f.T)
dres=np.append(dres,drho_dphi_hf.T)
dres=np.append(dres,drho_dphi_lf.T)


dres=np.append(dres,drho_dz_f.T)
dres=np.append(dres,drho_dz_hf.T)
dres=np.append(dres,drho_dz_lf.T)

dres=np.append(dres,dw_dphi_hf.T)
dres=np.append(dres,dw_dz_f.T)
dres=np.append(dres,dw_dz_hf.T)


file=open(out_dir+'data.dat','w')
dfile=open(out_dir+'ddata.dat','w')
#fin=np.append(fin,theta_f.T)
#file.write(w_f.tofile())
np.savetxt(out_dir+'data.dat',res, fmt='%.8e')
np.savetxt(out_dir+'ddata.dat',dres, fmt='%.8e')

file.close()
dfile.close()
