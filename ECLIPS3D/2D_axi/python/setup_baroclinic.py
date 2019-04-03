# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:38:10 2017

@author: florian
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.misc import derivative
import scipy.optimize as Opt


three=False
#three = True


Nlong = 35
Nlat = 25
Nz = 12


#parameters from Ulrich et al. 2014
a = 6.371229e06
b = 2.
d0 = a/2.
g = 9.80616
k = 3.
p0 = 1.0E5
R = 287.0
T0_E = 310
T0_P = 240
V_p = 1.0
z_t = 1.5e4
Gamma = 0.005
lambda_c = np.pi/9.
phi_c = 2.*np.pi/9.
omega = 7.29212e-5

Height = 30.0E3
cp = 1005.0

A = 1./Gamma
T0 = 0.5 * (T0_E+T0_P)
B = (T0-T0_P)/(T0*T0_P)
C = (k+2.)/2. * (T0_E-T0_P)/(T0_E*T0_P)

H = R*T0/g


dlong = 2*np.pi/Nlong
dlat = 0.5*np.pi/Nlat
dz = Height/Nz

XU = np.linspace(0,Nlong-1,Nlong)*dlong
XV = XU+0.5*dlong

YU = np.linspace(0,Nlat-1,Nlat)*dlat + 0.5*dlat
YV = np.linspace(0,Nlat,Nlat+1)*dlat

ZU = np.linspace(0,Nz-1,Nz)*dz + 0.5*dz
ZW = np.linspace(0,Nz,Nz+1)*dz


#=============================================
# Intermediate quantities
#=============================================

t1 = lambda z : A*Gamma/T0*np.exp(Gamma/T0*(z-a)) + \
B*(1. - 2.*((z-a)/(b*H))**2)*np.exp(-((z-a)/(b*H))**2)

t2 = lambda z : C*(1. - 2.*((z-a)/(b*H))**2)*np.exp(-((z-a)/(b*H))**2)


t1_int = lambda z : A * (np.exp(Gamma/T0*(z-a)) - 1) + \
B*(z-a) * np.exp(-((z-a)/(b*H))**2)

t2_int = lambda z : C*(z-a) * np.exp(-((z-a)/(b*H))**2)

dt1_dz = lambda z : A*(Gamma/T0)**2*np.exp(Gamma/T0*(z-a)) + \
B*(- 6.*(z-a)/(b*H)**2 + 4.*(z-a)**3/(b*H)**4)*np.exp(-((z-a)/(b*H))**2)

dt2_dz = lambda z :C*(- 6.*(z-a)/(b*H)**2 + 4.*(z-a)**3/(b*H)**4)*np.exp(-((z-a)/(b*H))**2)


#=============================================
# Thermodynamic values
#=============================================


T = lambda y,z : (a/z)**2 * (t1(z) - t2(z)* ( (z/a * np.cos(y))**k - k/(k+2.)*(z/a * np.cos(y))**(k+2.)) ) ** (-1.)

dT_dphi = lambda y,z : (z/a)**2 * t2(z) * k*np.tan(y)* \
( -(z/a*np.cos(y))**k + (z/a*np.cos(y))**(k+2) ) * T(y,z)**2

dT_dr = lambda y,z : T(y,z) * (-2/z - (dt1_dz(z)-dt2_dz(z)* ( (z/a * np.cos(y))**k - k/(k+2.)*(z/a * np.cos(y))**(k+2.)) - \
t2(z)*k*(z**(k-1)*(np.cos(y)/a)**k - z**(k+1)*(np.cos(y)/a)**(k+2))) / \
(t1(z) - t2(z)* ( (z/a * np.cos(y))**k - k/(k+2.)*(z/a * np.cos(y))**(k+2.)) ) )


p = lambda y,z : p0 * np.exp(-g/R * t1_int(z) + g/R * t2_int(z) * \
( (z/a*np.cos(y))**k - k/(k+2.)*(z*np.cos(y)/a)**(k+2.)) )

dp_dphi = lambda y,z : k*g/R*np.tan(y)*t2_int(z)* \
( -(z/a*np.cos(y))**k + (z/a*np.cos(y))**(k+2) ) * p(y,z)

dp_dr = lambda y,z : g/R*p(y,z) * (-t1(z)+t2(z)*((z/a*np.cos(y))**k - k/(k+2.)*(z*np.cos(y)/a)**(k+2.)) \
+t2_int(z)*k*(z**(k-1)*(np.cos(y)/a)**k - z**(k+1)*(np.cos(y)/a)**(k+2)) )

theta = lambda y,z : (p0/p(y,z))**(R/cp)*T(y,z)

dtheta_dphi = lambda y,z : -(R/cp)*p0**(R/cp)/p(y,z)**(R/cp+1)*dp_dphi(y,z)*T(y,z)+ \
dT_dphi(y,z)*(p0/p(y,z))**(R/cp)

dtheta_dr = lambda y,z : -(R/cp)*p0**(R/cp)/p(y,z)**(R/cp+1)*dp_dr(y,z)*T(y,z)+ \
dT_dr(y,z)*(p0/p(y,z))**(R/cp)

rho = lambda y,z : p(y,z)/(R*T(y,z))

drho_dphi = lambda y,z : 1./(R*T(y,z))*dp_dphi(y,z) - p(y,z)/(R*T(y,z)**2)*dT_dphi(y,z)

drho_dr = lambda y,z : 1./(R*T(y,z))*dp_dr(y,z) - p(y,z)/(R*T(y,z)**2)*dT_dr(y,z)


#=============================================
# Speed
#=============================================

U = lambda y,z : g/a * k * T(y,z) *t2_int(z) * (\
(z/a*np.cos(y))**(k-1.) - (z/a*np.cos(y))**(k+1.) )

dU_dphi = lambda y,z : dT_dphi(y,z)*U(y,z)/T(y,z)+ g/a*k*T(y,z)*t2_int(z)*\
np.sin(y)*((k+1)*(z/a)**(k+1)*np.cos(y)**(k) - (k-1)*(z/a)**(k-1)*np.cos(y)**(k-2))

dU_dr = lambda y,z : dT_dr(y,z)*U(y,z)/T(y,z) + t2(z)*U(y,z)/t2_int(z) + \
g/a*k*T(y,z)*t2_int(z)*((k-1)*z**(k-2)*(np.cos(y)/a)**(k-1) - (k+1)*z**k*(np.cos(y)/a)**(k+1)) 


u = lambda y,z : -omega*z*np.cos(y) + np.sqrt((omega*z*np.cos(y))**2 + \
z*np.cos(y)*U(y,z) )

du_dphi = lambda y,z : omega*z*np.sin(y) + 1./(2*np.sqrt((omega*z*np.cos(y))**2 + \
z*np.cos(y)*U(y,z) )) * (-2*z**2*omega**2*np.sin(y)*np.cos(y) - z*np.sin(y)*U(y,z) + \
z*np.cos(y)*dU_dphi(y,z))

du_dr = lambda y,z : -omega*np.cos(y) + 1./(2*np.sqrt((omega*z*np.cos(y))**2 + \
z*np.cos(y)*U(y,z) )) * (2*z*(omega*np.cos(y))**2 + np.cos(y)*U(y,z) + \
z*np.cos(y)*dU_dr(y,z) ) 

#=============================================
# Perturbed speed
#=============================================

dzeta = lambda z : 1.-3.*(z/z_t)**2 + 2.*(z/z_t)**3

dd = lambda x,y : a*np.arccos(np.sin(phi_c)*np.sin(y) + np.cos(phi_c)*np.cos(y)*np.cos(x-lambda_c) )

u1 = lambda x,y,z : -16.*V_p/(3.*np.sqrt(3.))*dzeta(z)*np.cos(np.pi*max(0,min(d0,dd(x,y))/(2*d0)))**3* \
np.sin(np.pi*max(0,min(d0,dd(x,y))/(2*d0))) * \
(-np.sin(phi_c)*np.cos(y) + np.cos(phi_c)*np.sin(y)*np.cos(x-lambda_c))/(np.sin(dd(x,y)/a))

#u1 = lambda x,y,z : -16.*V_p/(3.*np.sqrt(3.))*dzeta(z)*np.cos(np.pi*dd(x,y)/(2*d0))**3* \
#np.sin(np.pi*dd(x,y)/(2*d0)) * \
#(-np.sin(phi_c)*np.cos(y) + np.cos(phi_c)*np.sin(y)*np.cos(x-lambda_c))/(np.sin(dd(x,y)/a))

v1 = lambda x,y,z : -16.*V_p/(3.*np.sqrt(3.))*dzeta(z)*np.cos(np.pi*max(0,min(d0,dd(x,y))/(2*d0)))**3*\
np.sin(np.pi*max(0,min(d0,dd(x,y))/(2*d0))) * \
(np.cos(phi_c)*np.sin(x-lambda_c))/(np.sin(dd(x,y)/a))



u_u = np.zeros((Nlong,Nlat,Nz))
u_v = np.zeros((Nlong,Nlat+1,Nz))
u_w = np.zeros((Nlong,Nlat,Nz+1))
u_p = np.zeros((Nlong,Nlat,Nz))

v_u = np.zeros((Nlong,Nlat,Nz))
v_v = np.zeros((Nlong,Nlat+1,Nz))
v_w = np.zeros((Nlong,Nlat,Nz+1))
v_p = np.zeros((Nlong,Nlat,Nz))

w_u = np.zeros((Nlong,Nlat,Nz))
w_v = np.zeros((Nlong,Nlat+1,Nz))
w_w = np.zeros((Nlong,Nlat,Nz+1))
w_p = np.zeros((Nlong,Nlat,Nz))


p_u = np.zeros((Nlong,Nlat,Nz))
p_v = np.zeros((Nlong,Nlat+1,Nz))
p_w = np.zeros((Nlong,Nlat,Nz+1))
p_p = np.zeros((Nlong,Nlat,Nz))

theta_u = np.zeros((Nlong,Nlat,Nz))
theta_v = np.zeros((Nlong,Nlat+1,Nz))
theta_w = np.zeros((Nlong,Nlat,Nz+1))
theta_p = np.zeros((Nlong,Nlat,Nz))

rho_u = np.zeros((Nlong,Nlat,Nz))
rho_v = np.zeros((Nlong,Nlat+1,Nz))
rho_w = np.zeros((Nlong,Nlat,Nz+1))
rho_p = np.zeros((Nlong,Nlat,Nz))



#--------------------------------------------------------
#Derivatives on u 
#--------------------------------------------------------
du_dlong_u = np.zeros((Nlong,Nlat,Nz))
du_dphi_u = np.zeros((Nlong,Nlat,Nz))
du_dz_u = np.zeros((Nlong,Nlat,Nz))
dp_dlong_u = np.zeros((Nlong,Nlat,Nz))
drho_dlong_u = np.zeros((Nlong,Nlat,Nz))
drho_dphi_u = np.zeros((Nlong,Nlat,Nz))
drho_dz_u = np.zeros((Nlong,Nlat,Nz))

#--------------------------------------------------------
#Derivatives on v
#--------------------------------------------------------
dv_dlong_v = np.zeros((Nlong,Nlat+1,Nz))
dv_dphi_v = np.zeros((Nlong,Nlat+1,Nz))
dv_dz_v = np.zeros((Nlong,Nlat+1,Nz))
dp_dphi_v = np.zeros((Nlong,Nlat+1,Nz))
drho_dlong_v = np.zeros((Nlong,Nlat+1,Nz))
drho_dphi_v = np.zeros((Nlong,Nlat+1,Nz))
drho_dz_v = np.zeros((Nlong,Nlat+1,Nz))

#--------------------------------------------------------
#Derivatives on p
#--------------------------------------------------------
du_dlong_p = np.zeros((Nlong,Nlat,Nz))
dv_dphi_p = np.zeros((Nlong,Nlat,Nz))
dp_dlong_p = np.zeros((Nlong,Nlat,Nz))
dp_dphi_p = np.zeros((Nlong,Nlat,Nz))
dw_dz_p = np.zeros((Nlong,Nlat,Nz))
dtheta_dz_p =np.zeros((Nlong,Nlat,Nz))
drho_dlong_p = np.zeros((Nlong,Nlat,Nz))
drho_dphi_p = np.zeros((Nlong,Nlat,Nz))

#--------------------------------------------------------
#Derivatives on w
#--------------------------------------------------------
dp_dz_w = np.zeros((Nlong,Nlat,Nz+1))
dw_dlong_w = np.zeros((Nlong,Nlat,Nz+1))
dw_dphi_w = np.zeros((Nlong,Nlat,Nz+1))
dw_dz_w = np.zeros((Nlong,Nlat,Nz+1))
dtheta_dlong_w = np.zeros((Nlong,Nlat,Nz+1))
dtheta_dphi_w = np.zeros((Nlong,Nlat,Nz+1))
dtheta_dz_w = np.zeros((Nlong,Nlat,Nz+1))
drho_dlong_w = np.zeros((Nlong,Nlat,Nz+1))
drho_dphi_w = np.zeros((Nlong,Nlat,Nz+1))
drho_dz_w = np.zeros((Nlong,Nlat,Nz+1))

T_u=np.zeros((Nlong,Nlat,Nz))

#=============================================
# Perturbed speed
#=============================================

for i in range(Nlat):
    for j in range(Nz) :
        u_u[:,i,j]=np.zeros(Nlong)+np.array(u(YU[i],a+ZU[j]))
        u_v[:,i,j]=np.zeros(Nlong)+np.array(u(YV[i],a+ZU[j]))
        u_w[:,i,j]=np.zeros(Nlong)+np.array(u(YU[i],a+ZW[j]))
        u_p[:,i,j]=np.zeros(Nlong)+np.array(u(YU[i],a+ZU[j]))
        
        T_u[:,i,j]=np.zeros(Nlong)+np.array(T(YU[i],a+ZU[j]))
        
        p_u[:,i,j]=np.zeros(Nlong)+np.array(p(YU[i],a+ZU[j]))
        p_v[:,i,j]=np.zeros(Nlong)+np.array(p(YV[i],a+ZU[j]))
        p_w[:,i,j]=np.zeros(Nlong)+np.array(p(YU[i],a+ZW[j]))
        p_p[:,i,j]=np.zeros(Nlong)+np.array(p(YU[i],a+ZU[j]))
        
        theta_u[:,i,j]=np.zeros(Nlong)+np.array(theta(YU[i],a+ZU[j]))
        theta_v[:,i,j]=np.zeros(Nlong)+np.array(theta(YV[i],a+ZU[j]))
        theta_w[:,i,j]=np.zeros(Nlong)+np.array(theta(YU[i],a+ZW[j]))
        theta_p[:,i,j]=np.zeros(Nlong)+np.array(theta(YU[i],a+ZU[j]))
        
        rho_u[:,i,j]=np.zeros(Nlong)+np.array(rho(YU[i],a+ZU[j]))
        rho_v[:,i,j]=np.zeros(Nlong)+np.array(rho(YV[i],a+ZU[j]))
        rho_w[:,i,j]=np.zeros(Nlong)+np.array(rho(YU[i],a+ZW[j]))
        rho_p[:,i,j]=np.zeros(Nlong)+np.array(rho(YU[i],a+ZU[j]))
        
        du_dphi_u[:,i,j]=np.zeros(Nlong)+np.array(du_dphi(YU[i],a+ZU[j]))
        du_dz_u[:,i,j] = np.zeros(Nlong)+np.array(du_dr(YU[i],a+ZU[j]))
        drho_dphi_u[:,i,j] = np.zeros(Nlong)+np.array(drho_dphi(YU[i],a+ZU[j]))
        drho_dz_u[:,i,j] = np.zeros(Nlong)+np.array(drho_dr(YU[i],a+ZU[j]))
        
        dp_dphi_v[:,i,j] = np.zeros(Nlong)+np.array(dp_dphi(YV[i],a+ZU[j]))
        drho_dphi_v[:,i,j] = np.zeros(Nlong)+np.array(drho_dphi(YV[i],a+ZU[j]))
        drho_dz_v[:,i,j] = np.zeros(Nlong)+np.array(drho_dr(YV[i],a+ZU[j]))
        
        dp_dphi_p[:,i,j] = np.zeros(Nlong)+np.array(dp_dphi(YU[i],a+ZU[j]))
        dtheta_dz_p[:,i,j] = np.zeros(Nlong)+np.array(dtheta_dr(YU[i],a+ZU[j]))
        drho_dphi_p[:,i,j] = np.zeros(Nlong)+np.array(drho_dphi(YU[i],a+ZU[j]))
        
        dp_dz_w[:,i,j] = np.zeros(Nlong)+np.array(dp_dr(YU[i],a+ZW[j]))
        dtheta_dphi_w[:,i,j] = np.zeros(Nlong)+np.array(dtheta_dphi(YU[i],a+ZW[j]))
        dtheta_dz_w[:,i,j] = np.zeros(Nlong)+np.array(dtheta_dr(YU[i],a+ZW[j]))
        drho_dphi_w[:,i,j] = np.zeros(Nlong)+np.array(drho_dphi(YU[i],a+ZW[j]))
        drho_dz_w[:,i,j] = np.zeros(Nlong)+np.array(drho_dr(YU[i],a+ZW[j]))
        
    u_w[:,i,Nz]=np.zeros(Nlong)+np.array(u(YU[i],a+ZW[Nz]))
    p_w[:,i,Nz]=np.zeros(Nlong)+np.array(p(YU[i],a+ZW[Nz]))
    theta_w[:,i,Nz]=np.zeros(Nlong)+np.array(theta(YU[i],a+ZW[Nz]))
    rho_w[:,i,Nz]=np.zeros(Nlong)+np.array(rho(YU[i],a+ZW[Nz]))
        
        
    dp_dz_w[:,i,Nz] = np.zeros(Nlong)+np.array(dp_dr(YU[i],a+ZW[Nz]))
    dtheta_dphi_w[:,i,Nz] = np.zeros(Nlong)+np.array(dtheta_dphi(YU[i],a+ZW[Nz]))
    dtheta_dz_w[:,i,Nz] = np.zeros(Nlong)+np.array(dtheta_dr(YU[i],a+ZW[Nz]))
    drho_dphi_w[:,i,Nz] = np.zeros(Nlong)+np.array(drho_dphi(YU[i],a+ZW[Nz]))
    drho_dz_w[:,i,Nz] = np.zeros(Nlong)+np.array(drho_dr(YU[i],a+ZW[Nz]))
    
for j in range(Nz) :
    u_v[:,Nlat,j]=np.zeros(Nlong)+np.array(u(YV[i],a+ZU[j]))
    p_v[:,Nlat,j]=np.zeros(Nlong)+np.array(p(YV[i],a+ZU[j]))
    theta_v[:,Nlat,j]=np.zeros(Nlong)+np.array(theta(YV[i],a+ZU[j]))
    rho_v[:,Nlat,j]=np.zeros(Nlong)+np.array(rho(YV[i],a+ZU[j]))
    
    
    
    dp_dphi_v[:,Nlat,j] = np.zeros(Nlong)+np.array(dp_dphi(YV[Nlat],a+ZU[j]))
    drho_dphi_v[:,Nlat,j] = np.zeros(Nlong)+np.array(drho_dphi(YV[Nlat],a+ZU[j]))
    drho_dz_v[:,Nlat,j] = np.zeros(Nlong)+np.array(drho_dr(YV[Nlat],a+ZU[j]))
        
        
print('perturbation ...')
#
#
for i in range(Nlong) :
    for j in range(Nlat) :
        for kk in range(Nz) :
            u_u[i,j,kk] = u_u[i,j,kk]+ u1(XU[i],YU[j],min(z_t,ZU[kk]))
            u_v[i,j,kk] += u1(XV[i],YV[j],min(z_t,ZU[kk]))
            u_w[i,j,kk] += u1(XV[i],YU[j],min(z_t,ZW[kk]))
            u_p[i,j,kk] += u1(XV[i],YU[j],min(z_t,ZU[kk]))
            
            v_u[i,j,kk] += v1(XU[i],YU[j],min(z_t,ZU[kk]))
            v_v[i,j,kk] += v1(XV[i],YV[j],min(z_t,ZU[kk]))
            v_w[i,j,kk] += v1(XV[i],YU[j],min(z_t,ZW[kk]))
            v_p[i,j,kk] += v1(XV[i],YU[j],min(z_t,ZU[kk]))
            
            du_dlong_u[i,j,kk] += derivative(lambda x : u1(x,YU[j],min(z_t,ZU[kk])),XU[i])
            du_dphi_u[i,j,kk] += derivative(lambda y : u1(XU[i],y,min(z_t,ZU[kk])),YU[j])
            du_dz_u[i,j,kk] += derivative(lambda z : u1(XU[i],YU[j],z),min(z_t,ZU[kk]))
            
            dv_dlong_v[i,j,kk] += derivative(lambda x : v1(x,YV[j],min(z_t,ZU[kk])),XV[i])
            dv_dphi_v[i,j,kk] += derivative(lambda y : v1(XV[i],y,min(z_t,ZU[kk])),YV[j])
            dv_dz_v[i,j,kk] += derivative(lambda z : v1(XV[i],YV[j],z),min(z_t,ZU[kk]))
            
            du_dlong_p[i,j,kk] += derivative(lambda x : u1(x,YU[j],min(z_t,ZU[kk])),XV[i])
            dv_dphi_p[i,j,kk] += derivative(lambda y : v1(XV[i],y,min(z_t,ZU[kk])),YU[j])
            
        u_w[i,j,Nz] += u1(XV[i],YU[j],min(z_t,ZW[Nz]))
        v_w[i,j,Nz] += v1(XV[i],YU[j],min(z_t,ZW[Nz]))
        
    for kk in range(Nz) :
        u_v[i,Nlat,kk] += u1(XV[i],YV[Nlat],min(z_t,ZU[kk]))
        v_v[i,Nlat,kk] += v1(XV[i],YV[Nlat],min(z_t,ZU[kk]))
        
        dv_dlong_v[i,Nlat,kk] += derivative(lambda x : v1(x,YV[Nlat],min(z_t,ZU[kk])),XV[i])
        dv_dphi_v[i,Nlat,kk] += derivative(lambda y : v1(XV[i],y,min(z_t,ZU[kk])),YV[Nlat])
        dv_dz_v[i,Nlat,kk] += derivative(lambda z : v1(XV[i],YV[Nlat],z),min(z_t,ZU[kk]))
#        

print ('over')
            
if three :   

    res=np.append(u_u,u_v)
    res=np.append(res,u_p)
    res=np.append(res,u_w)
    
    res=np.append(res,v_u)
    res=np.append(res,v_v)
    res=np.append(res,v_p)
    res=np.append(res,v_w)
    
    res=np.append(res,p_u)
    res=np.append(res,p_v)
    res=np.append(res,p_p)
    res=np.append(res,p_w)
    
    res=np.append(res,w_u)
    res=np.append(res,w_v)
    res=np.append(res,w_p)
    res=np.append(res,w_w)
    
    res=np.append(res,theta_u)
    res=np.append(res,theta_v)
    res=np.append(res,theta_p)
    res=np.append(res,theta_w)

    res=np.append(res,rho_u)
    res=np.append(res,rho_v)
    res=np.append(res,rho_p)
    res=np.append(res,rho_w)
        
            
            
    dures=np.append(du_dlong_u,du_dphi_u)
    dures=np.append(dures,du_dz_u)
    dures=np.append(dures,dp_dlong_u)
    dures=np.append(dures,drho_dlong_u)
    dures=np.append(dures,drho_dphi_u)
    dures=np.append(dures,drho_dz_u)
    
    
    dvres=np.append(dv_dlong_v,dv_dphi_v)
    dvres=np.append(dvres,dv_dz_v)
    dvres=np.append(dvres,dp_dphi_v)
    dvres=np.append(dvres,drho_dlong_v)
    dvres=np.append(dvres,drho_dphi_v)
    dvres=np.append(dvres,drho_dz_v)
    
    
    
    dpres=np.append(du_dlong_p,dv_dphi_p)
    dpres=np.append(dpres,dp_dlong_p)
    dpres=np.append(dpres,dp_dphi_p)
    dpres=np.append(dpres,dw_dz_p)
    dpres=np.append(dpres,dtheta_dz_p)
    dpres=np.append(dpres,drho_dlong_p)
    dpres=np.append(dpres,drho_dphi_p)
    
    
    
    dwres=np.append(dp_dz_w,dw_dlong_w)
    dwres=np.append(dwres,dw_dphi_w)
    dwres=np.append(dwres,dw_dz_w)
    dwres=np.append(dwres,dtheta_dlong_w)
    dwres=np.append(dwres,dtheta_dphi_w)
    dwres=np.append(dwres,dtheta_dz_w)
    dwres=np.append(dwres,drho_dlong_w)
    dwres=np.append(dwres,drho_dphi_w)
    dwres=np.append(dwres,drho_dz_w)
    
    
    
    file=open('/Users/florian/Desktop/MyWork/3D_sca/data/data.dat','w')
    dufile=open('/Users/florian/Desktop/MyWork/3D_sca/data/dudata.dat','w')
    dvfile=open('/Users/florian/Desktop/MyWork/3D_sca/data/dvdata.dat','w')
    dpfile=open('/Users/florian/Desktop/MyWork/3D_sca/data/dpdata.dat','w')
    dwfile=open('/Users/florian/Desktop/MyWork/3D_sca/data/dwdata.dat','w')
    
    
    np.savetxt('/Users/florian/Desktop/MyWork/3D_sca/data/data.dat',res, fmt='%.8e')
    np.savetxt('/Users/florian/Desktop/MyWork/3D_sca/data/dudata.dat',dures, fmt='%.8e')
    np.savetxt('/Users/florian/Desktop/MyWork/3D_sca/data/dvdata.dat',dvres, fmt='%.8e')
    np.savetxt('/Users/florian/Desktop/MyWork/3D_sca/data/dpdata.dat',dpres, fmt='%.8e')
    np.savetxt('/Users/florian/Desktop/MyWork/3D_sca/data/dwdata.dat',dwres, fmt='%.8e')
    
    
    file.close()
    dufile.close()
    dvfile.close()
    dpfile.close()
    dwfile.close()
    
else :
    
    u_f=u_u[0]
    u_hf=u_w[0]
    u_lf=u_v[0]
    
    v_f=v_u[0]
    v_hf=v_w[0]
    v_lf=v_v[0]
    
    w_f=w_u[0]
    w_hf=w_w[0]
    w_lf=w_v[0]
    
    p_f=p_u[0]
    p_hf=p_w[0]
    p_lf=p_v[0]
    
    theta_f=theta_u[0]
    theta_hf=theta_w[0]
    theta_lf=theta_v[0]
    
    rho_f=rho_u[0]
    rho_hf=rho_w[0]
    rho_lf=rho_v[0]

    du_dphi_f = du_dphi_u[0]
    du_dz_f = du_dz_u[0]
    
    dv_dphi_f = dv_dphi_p[0]
    dv_dphi_lf = dv_dphi_v[0]
    dv_dz_lf = dv_dz_v[0]
    
    dw_dphi_hf = dw_dphi_w[0]
    dw_dz_hf = dw_dz_w[0]   
    dw_dz_f = dw_dz_p[0]
    
    dp_dphi_f = dp_dphi_p[0]
    dp_dphi_lf = dp_dphi_v[0]
    dp_dz_hf = dp_dz_w[0]
    
    dtheta_dz_f = dtheta_dz_p[0]
    dtheta_dphi_hf = dtheta_dphi_w[0]
    dtheta_dz_hf = dtheta_dz_w[0]
    
    drho_dphi_f = drho_dphi_u[0]
    drho_dphi_hf = drho_dphi_w[0]
    drho_dphi_lf = drho_dphi_v[0]
    
    drho_dz_f = drho_dz_u[0]
    drho_dz_hf = drho_dz_w[0]
    drho_dz_lf = drho_dz_v[0]
    

    print(u_f.shape)

    res=np.append(u_f,u_hf)
    res=np.append(res,u_lf)
    
    
    res=np.append(res,v_f)
    res=np.append(res,v_hf)
    res=np.append(res,v_lf)
    
    res=np.append(res,p_f)
    res=np.append(res,p_hf)
    res=np.append(res,p_lf)
    
    res=np.append(res,theta_f)
    res=np.append(res,theta_hf)
    res=np.append(res,theta_lf)
    
    res=np.append(res,rho_f)
    res=np.append(res,rho_hf)
    res=np.append(res,rho_lf)
    
    res=np.append(res,w_f)
    res=np.append(res,w_hf)
    res=np.append(res,w_lf)
    
    
    
    dres=np.append(du_dphi_f,du_dz_f)
    
    dres=np.append(dres,dv_dphi_f)
    dres=np.append(dres,dv_dphi_lf)
    dres=np.append(dres,dv_dz_lf)
    
    dres=np.append(dres,dp_dphi_f)
    dres=np.append(dres,dp_dphi_lf)
    dres=np.append(dres,dp_dz_hf)
    
    dres=np.append(dres,dtheta_dz_f)
    dres=np.append(dres,dtheta_dphi_hf)
    dres=np.append(dres,dtheta_dz_hf)
    
    dres=np.append(dres,drho_dphi_f)
    dres=np.append(dres,drho_dphi_hf)
    dres=np.append(dres,drho_dphi_lf)
    dres=np.append(dres,drho_dz_f)
    dres=np.append(dres,drho_dz_hf)
    dres=np.append(dres,drho_dz_lf)
    
    
    dres=np.append(dres,dw_dphi_hf)
    dres=np.append(dres,dw_dz_f)
    dres=np.append(dres,dw_dz_hf)
    
    
    
        
    file=open('/Users/florian/Desktop/MyWork/2D/data_grid/steady_state/data.dat','w')
    dfile=open('/Users/florian/Desktop/MyWork/2D/data_grid/steady_state/ddata.dat','w')
    #fin=np.append(fin,theta_f.T)
    #file.write(w_f.tofile())
    np.savetxt('/Users/florian/Desktop/MyWork/2D/data_grid/steady_state/data.dat',res, fmt='%.8e')
    np.savetxt('/Users/florian/Desktop/MyWork/2D/data_grid/steady_state/ddata.dat',dres, fmt='%.8e')
    
    file.close()
    dfile.close()










