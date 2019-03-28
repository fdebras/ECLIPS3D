# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 18:35:42 2017

@author: florian
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:53:09 2016

@author: florian
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import gridspec

###############################################################################
###############################################################################
###############################################################################
# Careful ! With the previous set of data, 18 and 19 have to be switched 
# by 8 and 9
###############################################################################
###############################################################################
###############################################################################

#I had some troubles with NaN and stuff, this is a security in
#a way. Cleverer people might find more efficient way of using it. 
v1=18
v2=19

symmetric = True
ymax = 90

rep='/Users/florian/Desktop/MyWork/ECLIPS3D/data/3D_steady/'

infile=open(rep + 'initial_state.dat')
infile2=open(rep + 'solution_to_read.dat')



data=infile.read()
infile.close()

#Remove the blank spaces
values=data.split()

Nlong=int(values[0])
Nlat=int(values[1])
Nz=int(values[2])
nz=Nz
Heig=float(values[3])



r=9.44e7
omega=2.06E-5
Heigh_min=0
Heigh_max=1.1E7
g0=9.42
p0=2.2E7
gascons=4593.0
cp=14308.45


heigh_max=Heigh_max
heigh_min=Heigh_min

dlong=360.0/Nlong
dphi=90.0/(Nlat)
dz=(heigh_max-heigh_min)/Nz

XU=(np.linspace(0,Nlong-1,Nlong))*dlong
XV=(np.linspace(0,Nlong-1,Nlong)+0.5)*dlong
YU=(np.linspace(0,Nlat-1,Nlat)+0.5)*dphi
YV=(np.linspace(0,Nlat,Nlat+1))*dphi
ZU=heigh_min+(np.linspace(0,Nz-1,Nz)+0.5)*dz
ZW=heigh_min+(np.linspace(0,Nz,Nz+1))*dz



u_s=np.zeros((Nlong,Nlat,Nz))
v_s=np.zeros((Nlong,Nlat+1,Nz))
w_s=np.zeros((Nlong,Nlat,Nz+1))
p_s=np.zeros((Nlong,Nlat,Nz))
theta_s=np.zeros((Nlong,Nlat,Nz+1))

rho_u=np.zeros((Nlong,Nlat,Nz))
rho_v=np.zeros((Nlong,Nlat+1,Nz))
rho_w=np.zeros((Nlong,Nlat,Nz+1))
rho_p=np.zeros((Nlong,Nlat,Nz))

c_p=np.zeros((Nlong,Nlat,Nz))
n_w=np.zeros((Nlong,Nlat,Nz+1))




kill_lat=1 #number of points not to display in latitude
kill_z=1 #number of points not to display in latitude

a=4
#
for k in range(Nz) :
    for j in range(Nlat) :
        u_s[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong
            
for k in range(Nz) :
    for j in range(Nlat+1) : 
        v_s[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong

for k in range(Nz+1) :
    for j in range(Nlat) :
        w_s[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)
        a=a+Nlong
            
            
for k in range(Nz) :
    for j in range(Nlat) :
        p_s[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong

for k in range(Nz+1) :
    for j in range(Nlat) :
        theta_s[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong





for k in range(Nz) :
    for j in range(Nlat) :
        rho_u[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong
            
for k in range(Nz) :
    for j in range(Nlat+1) : 
        rho_v[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong
            
for k in range(Nz+1) :
    for j in range(Nlat) :
        rho_w[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong

for k in range(Nz) :
    for j in range(Nlat) :
        rho_p[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong

            
for k in range(Nz) :
    for j in range(Nlat) :
        c_p[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong
            
for k in range(Nz+1) :
    for j in range(Nlat) :
        n_w[:,j,k]=np.array(values[a:a+Nlong]).astype(np.float)  
        a=a+Nlong

print(a)
        
t=0
n_auto=0
stop=False

a=0


debz=10
finz=18

debl=0
finl=-1

deblong=0
finlong=-1

plot_theta=True
plot_theta=False





data2=infile2.read()
infile2.close()

#Remove the blank spaces
values2=data2.split()


u_r=np.zeros((Nlong,Nlat,Nz))
v_r=np.zeros((Nlong,Nlat+1,Nz))
p_r=np.zeros((Nlong,Nlat,Nz))
theta_r=np.zeros((Nlong,Nlat,Nz))
w_r=np.zeros((Nlong,Nlat,Nz))

#w_r and theta_r at top and bottom are zero but we add w_s and theta_s at the end


for i in range(Nlong) :
    for j in range(Nlat) :
        for k in range(Nz) :
            if (values2[a][v1]=='E' or values2[a][v2]=='E') :
                u_r[i,j,k]=np.array(values2[a]).astype(np.float)
            else :
                u_r[i,j,k]=0.0
            a=a+1
        
for i in range(Nlong) :
    for j in range(Nlat+1) :
        for k in range(Nz) :                   
            if (values2[a][v1]=='E' or values2[a][v2]=='E') :
                v_r[i,j,k]=np.array(values2[a]).astype(np.float)
            else :
                v_r[i,j,k]=0.0
            a=a+1
        
for i in range(Nlong) :
    for j in range(Nlat) :
        for k in range(Nz) :
            if (values2[a][v1]=='E' or values2[a][v2]=='E') :
                p_r[i,j,k]=np.array(values2[a]).astype(np.float)
            else :
                p_r[i,j,k]=0.0
            a=a+1
        
for i in range(Nlong) :
    for j in range(Nlat) :
        for k in range(Nz) :
            if (values2[a][v1]=='E' or values2[a][v2]=='E') :
                w_r[i,j,k]=np.array(values2[a]).astype(np.float)
            else :
                w_r[i,j,k]=0.0
            a=a+1
        
for i in range(Nlong) :
    for j in range(Nlat) :
        for k in range(Nz) :
            if (values2[a][v1]=='E' or values2[a][v2]=='E') :
                theta_r[i,j,k]=np.array(values2[a]).astype(np.float)
            else :
                theta_r[i,j,k]=0.0
            a=a+1

ftheta=interpolate.RegularGridInterpolator((XV,YU,r+ZW[1:]),theta_r)
CoordP=np.meshgrid(XV,YU,r+ZU[1:],indexing='ij')
CoordP=np.array(CoordP).T

caca = np.zeros((Nlong,Nlat,Nz))
caca[:,:,1:]=ftheta(CoordP).T


u_r/=(rho_u)
v_r/=(rho_v)        
w_r/=(rho_w[:,:,1:])       
theta_r/=rho_w[:,:,1:]/theta_s[:,:,1:]
for k in range(Nz) :
    theta_r[:,:,k]=theta_r[:,:,k]/(g0*r*r/(r+ZW[k+1])**2)
  
#
theta_r=theta_r+theta_s[:,:,1:]
p_r=p_r+p_s#
u_r=u_r+u_s
v_r=v_r+v_s
w_r=w_r+w_s[:,:,1:]


ftheta=interpolate.RegularGridInterpolator((XV,YU,r+ZW[1:]),theta_r)
fw=interpolate.RegularGridInterpolator((XV,YU,r+ZW[1:]),w_r)
CoordP=np.meshgrid(XV,YU,r+ZU[1:],indexing='ij')
CoordP=np.array(CoordP).T
theta_p=ftheta(CoordP).T

w_p = np.zeros((Nlong,Nlat,Nz))
w_p[:,:,1:] = fw(CoordP).T

T=np.sign(p_r[:,:,1:])*theta_p*(np.abs(p_r[:,:,1:])/p0)**(gascons/cp)


div=1 # the smaller, the more arrows
pole=1 # distance from the pole
time=0 #time for the frequency

uu=u_r
vv=v_r
pp=p_r

for posz in range(debz,finz+1) :
    YU2=YU[:Nlat-pole]
    plt.figure()
  

    plt.pcolormesh(XV,YU2,pp[:,:Nlat-pole,posz].T)
    plt.colorbar()
    plt.quiver(XU[::div],YU2[::div],uu[::div,:Nlat-pole:div,posz].T,vv[::div,:Nlat-pole:div,posz].T,scale=2*np.max([np.max(np.abs(uu[:,:Nlat-pole,posz])),np.max(np.abs(vv[:,:Nlat-pole,posz]))]),scale_units='inches',pivot='tail')

    plt.title('p perturbations in long-lat, z= ' + str(ZU[posz]))
    
    
if plot_theta :    
    for posz in range(debz,finz+1) :
        plt.figure()
        
        plt.pcolormesh(XV,YU2,theta_r[:,:Nlat-pole,posz].T)
        plt.colorbar()
        plt.quiver(XU[::div],YU2[::div],uu[::div,:Nlat-pole:div,posz].T,vv[::div,:Nlat-pole:div,posz].T,scale=2*np.max([np.max(np.abs(uu[:,:,posz])),np.max(np.abs(vv[:,:,posz]))]),scale_units='inches',pivot='tail')
        plt.title('Theta perturbations in long-lat, z= ' + str(ZU[posz]))
    
for posl in range(debl,finl+1)     :
    plt.figure()
    plt.pcolormesh(XV,ZU,p_r[:,posl].T)
    plt.colorbar()
    plt.quiver(XU[::div],ZU,u_r[::div,posl].T*0.005,w_r[::div,posl].T,scale=3*np.max([np.max(u_r[:,posl]*0.005),np.max(np.abs(w_r[:,posl]))]),scale_units='inches',pivot='tail')
    plt.title('p perturbations in long-height, lat= ' + str(YU[posl]))


if plot_theta :    
    for posl in range(debl,finl+1)     :   
        plt.figure()
        plt.pcolormesh(XV,ZW,theta_r[:,posl].T)
        plt.colorbar()
        plt.quiver(XU[::div],ZU,u_r[::div,posl,:].T,w_r[::div,posl].T,scale=3*np.max([np.max(u_r[::div,posl,:]),np.max(w_r[::div,posl,:-1])]),scale_units='inches',pivot='tail')
        plt.title('Theta perturbations in long-height, lat= ' + str(YU[posl]))
    
    
for poslong in range(deblong,finlong+1)     :
    plt.figure()
    plt.pcolormesh(YU2,ZU,p_r[poslong,:Nlat-pole:div].T)
    plt.colorbar()
    plt.quiver(YU2[::div],ZU,u_r[poslong,:Nlat-pole:div].T*0.005,w_r[poslong,:Nlat-pole:div].T,scale=3*np.max([np.max(u_r[poslong,:Nlat-pole,:-1]*0.005),np.max(w_r[poslong,:Nlat-pole])]),scale_units='inches',pivot='tail')
    plt.title('p perturbations in lat-height, long= ' + str(XV[poslong]))

if plot_theta :    
    for poslong in range(deblong,finlong+1)     :
        plt.figure()
        plt.pcolormesh(YU2,ZU,theta_r[poslong,:Nlat-pole:div].T)
        plt.colorbar()
        plt.quiver(YU2[::div],ZU,v_r[poslong,:Nlat-pole:div,:].T,w_r[poslong,:Nlat-pole:div].T,scale=3*np.max([np.max(v_r[poslong]),np.max(w_r[poslong])]),scale_units='inches',pivot='tail')
        plt.title('Theta perturbations in lat-height, long= ' + str(XV[poslong]))


#Mass fluxes

mass_u = rho_u*u_r*dz*r*dphi*np.pi/180
mass_v = np.zeros((Nlong,Nlat+1,Nz))
for i in range(Nlat+1) :
    mass_v[:,i] = rho_v[:,i]*v_r[:,i]*dz*r*np.cos(YV[i]*np.pi/180)*dlong*np.pi/180
    
mass_w = np.zeros((Nlong,Nlat,Nz))
for i in range(Nlat) :
    mass_w[:,i] = rho_w[:,i,1:]*w_r[:,i]*r*dlong*np.pi/180.*r*np.cos(YU[i]*np.pi/180)*dphi
    



#Eddy momentum equations
mer_transfer = np.zeros((Nlat,Nz))
ver_transfer = np.zeros((Nlat,Nz))

rho = np.zeros((Nlong,Nlat,Nz))
rho_v_prime = np.zeros((Nlong,Nlat,Nz))
rho_w_prime = np.zeros((Nlong,Nlat,Nz))

rho[:] = rho_u[:]+(p_r[:]-p_s[:])/c_p[:]-caca[:]/g0

for i in range(Nlong) :
    rho_v_prime[i] = rho[i]*v_r[i,1:] - np.mean(rho[:]*v_r[:,:-1],axis=0)
    rho_w_prime[i] = rho[i]*w_p[i] - np.mean(rho[:]*w_p[:],axis=0)



u_prime = np.zeros((Nlong,Nlat,Nz))
v_prime = np.zeros((Nlong,Nlat+1,Nz))
w_prime = np.zeros((Nlong,Nlat,Nz))
for i in range(Nlong) : 
    u_prime[i] = u_r[i]-np.mean(u_r,axis=0)
    v_prime[i] = v_r[i]-np.mean(v_r,axis=0)
    w_prime[i] = w_p[i]-np.mean(w_p,axis=0)

mer_transfer[:] = np.mean(rho_v_prime[:]*u_prime[:],axis=0)
ver_transfer[:] = np.mean(rho_w_prime[:]*u_prime[:],axis=0)

for j in range(Nlat) :
    mer_transfer[j] = mer_transfer[j]*np.cos(YU[j]*np.pi/180.)**2
    
for k in range(Nz) :
    ver_transfer[:,k] = ver_transfer[:,k]*(r+ZU[k])**3

ver_transfer[:] = np.gradient(ver_transfer[:],dphi*np.pi/180.,dz)[1]
mer_transfer[:] = np.gradient(mer_transfer[:],dphi*np.pi/180.,dz)[0]

for j in range(Nlat) : 
    for k in range(Nz) :
        mer_transfer[j,k] = mer_transfer[j,k]/((r+ZU[k])*np.cos(YU[j]*np.pi/180.)**2)
    
for k in range(Nz) : 
    ver_transfer[:,k] = ver_transfer[:,k]/((r+ZU[k])**3)


#Mean flow equations
mean_mer_transfer = np.zeros((Nlat,Nz))
mean_ver_transfer = np.zeros((Nlat,Nz))


rho_v_mean = np.zeros((Nlat,Nz))
rho_w_mean = np.zeros((Nlat,Nz))

rho_v_mean[:] = np.mean(rho[:]*v_r[:,:-1],axis=0)
rho_w_mean[:] = np.mean(rho[:]*w_p[:],axis=0)

mean_mer_transfer[:] = rho_v_mean[:]*np.mean(u_r,axis=0)[:]
mean_ver_transfer[:] = rho_w_mean[:]*np.mean(u_r,axis=0)[:]

for j in range(Nlat) :
    mean_mer_transfer[j] = mean_mer_transfer[j]*np.cos(YU[j]*np.pi/180.)**2
    
for k in range(Nz) :
    mean_ver_transfer[:,k] = mean_ver_transfer[:,k]*(r+ZU[k])**3

mean_ver_transfer[:] = np.gradient(mean_ver_transfer[:],dphi*np.pi/180.,dz)[1]
mean_mer_transfer[:] = np.gradient(mean_mer_transfer[:],dphi*np.pi/180.,dz)[0]

for j in range(Nlat) : 
    for k in range(Nz) :
        mean_mer_transfer[j,k] = mean_mer_transfer[j,k]/((r+ZU[k])*np.cos(YU[j]*np.pi/180.)**2)
    
for k in range(Nz) : 
    mean_ver_transfer[:,k] = mean_ver_transfer[:,k]/((r+ZU[k])**3) 
    
#Coriolis equations
    
coriolis_mer = np.zeros((Nlat,Nz))
coriolis_ver = np.zeros((Nlat,Nz))
for i in range(Nlat) :
    coriolis_mer[i] = -2.0*omega*rho_v_mean[i]*np.sin(YU[i]*np.pi/180.)
    
for i in range(Nlat) :
    coriolis_ver[i] = 2.0*omega*rho_w_mean[i]*np.cos(YU[i]*np.pi/180.)

tot_transfer = mer_transfer[:]+ver_transfer[:]+mean_ver_transfer[:]+ \
mean_mer_transfer[:]+coriolis_mer[:]+coriolis_ver[:]


if symmetric : 
    

    p_2 = np.zeros((Nlong,2*Nlat,Nz))
    theta_2 = np.zeros((Nlong,2*Nlat,Nz))
    u_2 = np.zeros((Nlong,2*Nlat,Nz))
    v_2 = np.zeros((Nlong,2*Nlat+1,Nz))
    w_2 = np.zeros((Nlong,2*Nlat,Nz))
    w_p_2 = np.zeros((Nlong,2*Nlat,Nz))
    caca_2 = np.zeros((Nlong,2*Nlat,Nz))
    rho_2 =  np.zeros((Nlong,2*Nlat,Nz))
    rho_3 =  np.zeros((Nlong,2*Nlat,Nz))
    
    p_2[:,Nlat:] = p_r
    theta_2[:,Nlat:] = theta_r
    u_2[:,Nlat:] = u_r
    v_2[:,Nlat:] = v_r
    w_2[:,Nlat:] = w_r
    w_p_2[:,Nlat:] = w_p
    caca_2[:,Nlat:] = caca
    rho_2[:,Nlat:] = rho
    rho_3[:,Nlat:] = rho[:]-rho_u[:]
    
    
    
    for i in range(Nlat) :
        p_2[:,i]=p_r[:,Nlat-1-i]
        theta_2[:,i]=theta_r[:,Nlat-1-i]
        u_2[:,i]=u_r[:,Nlat-1-i]
        v_2[:,i]=-v_r[:,Nlat-i] #v is antisymetric !
        w_2[:,i]=w_r[:,Nlat-1-i]
        w_p_2[:,i] = w_p[:,Nlat-1-i]
        caca_2[:,i] = caca[:,Nlat-1-i]
        rho_2[:,i] = rho[:,Nlat-1-i]
        rho_3[:,i] = rho[:,Nlat-1-i]-rho_u[:,Nlat-1-i]
        
    T_2 = np.sign(p_2[:,:,:])*theta_2[:,:,:]*(np.abs(p_2[:,:,:])/p0)**(gascons/cp)
    
    dphi_2=ymax/(Nlat)
    YUU=-ymax+(np.linspace(0,2*Nlat-1,2*Nlat)+.5)*dphi_2
    YVV=-ymax+(np.linspace(0,2*Nlat,2*Nlat+1))*dphi_2
    
    
    #Eddy momentum equations
    mer_transfer_2 = np.zeros((2*Nlat,Nz))
    ver_transfer_2 = np.zeros((2*Nlat,Nz))
    
    rho_v_prime_2 = np.zeros((Nlong,2*Nlat,Nz))
    rho_w_prime_2 = np.zeros((Nlong,2*Nlat,Nz))

    
    for i in range(Nlong) :
        rho_v_prime_2[i] = rho_2[i]*v_2[i,1:] - np.mean(rho_2[:]*v_2[:,:-1],axis=0)
        rho_w_prime_2[i] = rho_2[i]*w_2[i] - np.mean(rho_2[:]*w_2[:],axis=0)
    
    
    
    u_prime_2 = np.zeros((Nlong,2*Nlat,Nz))
    v_prime_2 = np.zeros((Nlong,2*Nlat+1,Nz))
    w_prime_2 = np.zeros((Nlong,2*Nlat,Nz))
    for i in range(Nlong) : 
        u_prime_2[i] = u_2[i]-np.mean(u_2,axis=0)
        v_prime_2[i] = v_2[i]-np.mean(v_2,axis=0)
        w_prime_2[i] = w_2[i]-np.mean(w_2,axis=0)
    
    mer_transfer_2[:] = np.mean(rho_v_prime_2[:]*u_prime_2[:],axis=0)
    ver_transfer_2[:] = np.mean(rho_w_prime_2[:]*u_prime_2[:],axis=0)
    
    for j in range(2*Nlat) :
        mer_transfer_2[j] = mer_transfer_2[j]*np.cos(YUU[j]*np.pi/180.)**2
        
    for k in range(Nz) :
        ver_transfer_2[:,k] = ver_transfer_2[:,k]*(r+ZU[k])**3
    
    ver_transfer_2[:] = np.gradient(ver_transfer_2[:],dphi_2*np.pi/180.,dz)[1]
    mer_transfer_2[:] = np.gradient(mer_transfer_2[:],dphi_2*np.pi/180.,dz)[0]
    
    for j in range(2*Nlat) : 
        for k in range(Nz) :
            mer_transfer_2[j,k] = mer_transfer_2[j,k]/((r+ZU[k])*np.cos(YUU[j]*np.pi/180.)**2)
        
    for k in range(Nz) : 
        ver_transfer_2[:,k] = ver_transfer_2[:,k]/((r+ZU[k])**3)
    
    
    #Mean flow equations
    mean_mer_transfer_2 = np.zeros((2*Nlat,Nz))
    mean_ver_transfer_2 = np.zeros((2*Nlat,Nz))
    
    
    rho_v_mean_2 = np.zeros((2*Nlat,Nz))
    rho_w_mean_2 = np.zeros((2*Nlat,Nz))
    
    rho_v_mean_3 = np.zeros((2*Nlat,Nz))
    rho_w_mean_3 = np.zeros((2*Nlat,Nz))
    
    rho_v_mean_3[:] = np.mean((rho_3[:])*v_2[:,:-1],axis=0)
    rho_w_mean_3[:] = np.mean((rho_3[:])*w_2[:],axis=0)

    rho_v_mean_2[:] = np.mean((rho_2[:])*v_2[:,:-1],axis=0)
    rho_w_mean_2[:] = np.mean((rho_2[:])*w_2[:],axis=0)
    
    mean_mer_transfer_2[:] = rho_v_mean_2[:]*np.mean(u_2,axis=0)[:]
    mean_ver_transfer_2[:] = rho_w_mean_2[:]*np.mean(u_2,axis=0)[:]
    
    for j in range(2*Nlat) :
        mean_mer_transfer_2[j] = mean_mer_transfer_2[j]*np.cos(YUU[j]*np.pi/180.)**2
        
    for k in range(Nz) :
        mean_ver_transfer_2[:,k] = mean_ver_transfer_2[:,k]*(r+ZU[k])**3
    
    mean_ver_transfer_2[:] = np.gradient(mean_ver_transfer_2[:],dphi_2*np.pi/180.,dz)[1]
    mean_mer_transfer_2[:] = np.gradient(mean_mer_transfer_2[:],dphi_2*np.pi/180.,dz)[0]
    
    for j in range(2*Nlat) : 
        for k in range(Nz) :
            mean_mer_transfer_2[j,k] = mean_mer_transfer_2[j,k]/((r+ZU[k])*np.cos(YUU[j]*np.pi/180.)**2)
        
    for k in range(Nz) : 
        mean_ver_transfer_2[:,k] = mean_ver_transfer_2[:,k]/((r+ZU[k])**3) 
        
    #Coriolis equations
        
    coriolis_mer_2 = np.zeros((2*Nlat,Nz))
    coriolis_ver_2 = np.zeros((2*Nlat,Nz))
    for i in range(2*Nlat) :
        coriolis_mer_2[i] = -2.0*omega*rho_v_mean_3[i]*np.sin(YUU[i]*np.pi/180.)
        
    for i in range(2*Nlat) :
        coriolis_ver_2[i] = 2.0*omega*rho_w_mean_3[i]*np.cos(YUU[i]*np.pi/180.)
    
    tot_transfer_2 = mer_transfer_2[:]+ver_transfer_2[:]+mean_ver_transfer_2[:]+ \
    mean_mer_transfer_2[:]+coriolis_mer_2[:]+coriolis_ver_2[:]




