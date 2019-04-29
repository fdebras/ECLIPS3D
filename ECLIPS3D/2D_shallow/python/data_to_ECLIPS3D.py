"""
Create data for the 2D_shallow 
version of ECLIPS3D
"""

import numpy as np


# Number of points in longitude
nlong=40
#Number of points in latitude
nlat=20


#output directory
output_dir = '/Users/florian/Desktop/test/ECLIPS3D/trunk/2D_shallow/data/'


###################################################################
#Information characterizing hot Jupiters
# See Showman Polvani 2011
phi = 0.
omega = 2.06E-5
g0=9.42

rtot = 9.44E7


beta = 2.*omega*np.cos(phi)/rtot

H = 4.E6/g0

V_caract = np.sqrt(g0*H)
L_caract = np.sqrt(np.sqrt(g0*H)/beta)
T_caract = 1./(np.sqrt(np.sqrt(g0*H)*beta))
###################################################################


# Longitudinal and latitudnal grid size and points. Not used in
# this program, but necessary in ECLIPS3D


dx = (2.*np.pi*rtot/L_caract)/nlong
dy = (np.pi*rtot/(1.3*L_caract))/nlat

#Longitduinal staggered points
XU = (np.linspace(-nlong/2+1,nlong/2,nlong)-0.5)*dx
XV = (np.linspace(-nlong/2+1,nlong/2,nlong))*dx
#Latitudinal staggered points
YU = (np.linspace(1,nlat,nlat)-0.5)*dy
YV = (np.linspace(0,nlat,nlat+1))*dy




#Implementing a rest background
u=np.zeros((nlong,nlat))
v=np.zeros((nlong,nlat+1))
h=np.zeros((nlong,nlat))+1. # H_0 = 1, adimensionalized

du_dx_u = np.zeros(np.shape(u))
du_dy_u = np.zeros(np.shape(u))
v_u =np.zeros(np.shape(u))

u_v = np.zeros(np.shape(v))
dv_dx_v = np.zeros(np.shape(v))
dv_dy_v = np.zeros(np.shape(v))


u_h = np.zeros(np.shape(h))
du_dx_h = np.zeros(np.shape(h))
v_h = np.zeros(np.shape(h))
dv_dy_h = np.zeros(np.shape(h))
dh_dx_h = np.zeros(np.shape(h))
dh_dy_h = np.zeros(np.shape(h))



#creating the output files
ures=np.append(u,du_dx_u)
ures=np.append(ures,du_dy_u)
ures=np.append(ures,v_u)


vres=np.append(u_v,v)
vres=np.append(vres,dv_dx_v)
vres=np.append(vres,dv_dy_v)

hres=np.append(u_h,du_dx_h)
hres=np.append(hres,v_h)
hres=np.append(hres,dv_dy_h)

hres=np.append(hres,h)
hres=np.append(hres,dh_dx_h)
hres=np.append(hres,dh_dy_h)



#Writing up data files
ufile=open(output_dir+'u_data.dat','w')
vfile=open(output_dir+'v_data.dat','w')
hfile=open(output_dir+'h_data.dat','w')

np.savetxt(output_dir+'u_data.dat',ures, fmt='%.8e')
np.savetxt(output_dir+'v_data.dat',vres, fmt='%.8e')
np.savetxt(output_dir+'h_data.dat',hres, fmt='%.8e')


ufile.close()
vfile.close()
hfile.close()

