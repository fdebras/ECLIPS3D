

import numpy as np
from scipy import interpolate

#output directory
out_dir = '/Users/florian/Desktop/MyWork/ECLIPS3D/data/3D_steady/'

#Size of the output
Nlongf=20
Nlatf=15
Nzf=18


#Can be used to change boundary conditions. Be careful
no_slip = True
#no_slip = False



ymax=90.

lambda_min=0.
lambda_max=360.



#Parameter from Komacek and Showman 2016
p_top = 100
p_bottom = 1.0E6

p_rad_top = 1000
p_rad_bottom= 1.0E6
trad=10**(6.)
trad_bottom = 1.0E7


Delta_T = 100.


#Data from Mayne et al. 2014
r=9.44e7
omega=2.06E-5
Heigh_min=0
Heigh_max=1.1E7
g0=9.42
p0=2.2E7
gascons=4593.0
cp=14308.4



Nlong=200
Nlat=100
Nz=100



dlong=2*np.pi/Nlong
dphi=0.5*np.pi/(Nlat)
dz=(Heigh_max-Heigh_min)/(Nz)

long = np.linspace(0,Nlong,Nlong+1)*dlong*180/np.pi
lat = np.linspace(0,Nlat,Nlat+1)*dphi*180/np.pi
height = Heigh_min+np.linspace(0,Nz,Nz+1)*dz


u = np.zeros((Nz+1,Nlat+1,Nlong+1))
v = np.zeros((Nz+1,Nlat+1,Nlong+1))
w = np.zeros((Nz+1,Nlat+1,Nlong+1))





T = lambda x : 1696.6986 + 132.23180 *np.log10(x/1.0E5) - 174.30459 * np.log10(x/1.0E5)**2 + \
        12.579612 * np.log10(x/1.0E5)**3 + 59.513639 * np.log10(x/1.0E5)**4 \
        + 9.6706522 * np.log10(x/1.0E5)**5 - 4.1136048 * np.log10(x/1.0E5)**6 - 1.0632301 * np.log10(x/1.0E5)**7 + \
        0.064400203 * np.log10(x/1.0E5) ** 8 \
        + 0.035974396 * np.log10(x/1.0E5)**9 + 0.0025740066 * np.log10(x/1.0E5)** 10.
        
Delta_Teq = lambda x :   Delta_T * np.log(max(p_top,min(x,p_bottom))/p_bottom) / np.log(p_top/p_bottom)     



#First trad : Komacek and Showman
#Second trad : Iro et al. 2005
trad_vec = lambda x : trad_bottom * (max(p_rad_top,min(x,p_rad_bottom))/p_rad_bottom) ** (np.log(trad/trad_bottom)/ np.log(p_rad_top/p_rad_bottom) )
#trad_vec = lambda x : 10.0**(5.4659686+1.4940124*np.log10(x/1.E5)+0.66079196*np.log10(x/1.E5)**2+0.16475329*np.log10(x/1.E5)**3+0.014241552*np.log10(x/1.E5)**4)


T_eq = lambda x,y,z : T(z) + Delta_Teq(z) * np.cos((min(270,max(x,90))-180.)*np.pi/180) * np.cos(y*np.pi/180) -Delta_Teq(z)/2.


Q = lambda x,y,z : (T_eq(x,y,z)-T(z))/trad_vec(z)
#Q = lambda x,y,z : (T(z))/trad




#
#p_high=1.0E6
#p_low=100



#See Heng et al. 2011 if you need this
#for i in range(Nlong) :
#    for j in range(Nlat) :
#        for k in range(Nz) :
#            
#            if (p[i,j,k]>=p_high) :
#                pp=np.log10(p_high/1.0E5)
#                T_night=1388.2145 + 267.66586 * pp - 215.53357 * pp**2 + \
#                61.814807 * pp**3 + 135.68661 * pp**4 + 2.0149044 * pp**5 - \
#                40.907246* pp**6 - 19.015628 * pp**7 - 3.8771634 * pp**8 - \
#                0.38413901 * pp**9 - 0.015089084 * pp**10 + \
#                100.0*(1.0-np.exp(-np.log10(p[i,j,k]/p_high)))
#                
#                T_day=2149.9581+4.1395571 * pp - 186.24851 * pp**2 + \
#                135.52524 * pp **3 + 106.20433 * pp **4 - 35.851966 * pp **5 - \
#                50.022826 * pp **6 - 18.462489 * pp **7 - 3.3319965 * pp**8 - \
#                0.30295925 * pp**9 - 0.011122316*pp**10 - \
#                120.0*(1.0-np.exp(-np.log10(p[i,j,k]/p_high)))
#            elif (p[i,j,k]<p_low) :
#                pp=np.log10(p_low/1.0E5)
#                T_night=1388.2145 + 267.66586 * pp - 215.53357 * pp**2 + \
#                61.814807 * pp**3 + 135.68661 * pp**4 + 2.0149044 * pp**5 - \
#                40.907246* pp**6 - 19.015628 * pp**7 - 3.8771634 * pp**8 - \
#                0.38413901 * pp**9 - 0.015089084 * pp**10
#                
#                T_night=max(T_night*(np.exp(0.01*np.log10(p[i,j,k]/p_low))),250)      
#                
#                
#                T_day=2149.9581+4.1395571 * pp - 186.24851 * pp**2 + \
#                135.52524 * pp **3 + 106.20433 * pp **4 - 35.851966 * pp **5 - \
#                50.022826 * pp **6 - 18.462489 * pp **7 - 3.3319965 * pp**8 - \
#                0.30295925 * pp**9 - 0.011122316*pp**10
#                
#                T_day=max(T_day*(np.exp(0.015*np.log10(p[i,j,k]/p_low))),1000)
#            else : 
#                pp=np.log10(p[i,j,k]/1.0E5)
#                T_night=1388.2145 + 267.66586 * pp - 215.53357 * pp**2 + \
#                61.814807 * pp**3 + 135.68661 * pp**4 + 2.0149044 * pp**5 - \
#                40.907246* pp**6 - 19.015628 * pp**7 - 3.8771634 * pp**8 - \
#                0.38413901 * pp**9 - 0.015089084 * pp**10
#                
#                T_day=2149.9581+4.1395571 * pp - 186.24851 * pp**2 + \
#                135.52524 * pp **3 + 106.20433 * pp **4 - 35.851966 * pp **5 - \
#                50.022826 * pp **6 - 18.462489 * pp **7 - 3.3319965 * pp**8 - \
#                0.30295925 * pp**9 - 0.011122316*pp**10
#                
#            if (long[i]>=90 and long[i]<=270) :
#                T_eq[i,j,k]=(T_night**4 + (T_day**4-T_night**4)* \
#                np.cos((long[i]-180.0)*np.pi/180.0) * np.cos(lat[j]*np.pi/180.0))**0.25
#            else :
#                T_eq[i,j,k]=T_night
#if you want a perturb steady state, uncomment
#theta=T_eq/exner

p = np.zeros((Nz+1,Nlat+1,Nlong+1))
rho = np.zeros((Nz+1,Nlat+1,Nlong+1))
theta = np.zeros((Nz+1,Nlat+1,Nlong+1))


p[0] = p0
rho[0] = p[0]/(gascons*T(p[0,0,0]))
theta[0] = T(p[0,0,0])*(p0/p[0])**(gascons/cp)

exner = (p/p0)**(gascons/cp)

g=np.zeros(Nz+1)
for k in range(Nz+1) :
    g[k]=g0*r*r/(r+height[k])**2

for k in range(Nz) : 
    exner[k+1] = exner[k]*np.exp(-g[k]*dz/(cp*T(p[k,0,0])))
    p[k+1] = p0*(exner[k+1])**(cp/gascons)
    rho[k+1] = p[k+1]/(gascons*T(p[k+1,0,0]))
    theta[k+1] = T(p[k+1,0,0])*(p0/p[k+1])**(gascons/cp)



du=np.gradient(u,dz,dphi,dlong)
dv=np.gradient(v,dz,dphi,dlong)
dp=np.gradient(p,dz,dphi,dlong)
dw=np.gradient(w,dz,dphi,dlong)
dtheta=np.gradient(theta,dz,dphi,dlong)
drho=np.gradient(rho,dz,dphi,dlong)



u=u.T
v=v.T
p=p.T
w=w.T
theta=theta.T
rho=rho.T


u_z=du[0].T
v_z=dv[0].T
p_z=dp[0].T
w_z=dw[0].T
theta_z=dtheta[0].T
rho_z=drho[0].T

u_phi=du[1].T
v_phi=dv[1].T
p_phi=dp[1].T
w_phi=dw[1].T
theta_phi=dtheta[1].T
rho_phi=drho[1].T

u_long=du[2].T
v_long=dv[2].T
p_long=dp[2].T
w_long=dw[2].T
theta_long=dtheta[2].T
rho_long=drho[2].T

dlong=360./Nlongf
dphi=ymax/(Nlatf)
dz=(Heigh_max-Heigh_min)/(Nzf)


XU=np.linspace(0,Nlongf-1,Nlongf)*dlong
XV=XU+0.5*dlong
YU=(np.linspace(0,Nlatf-1,Nlatf)+.5)*dphi
YV=(np.linspace(0,Nlatf,Nlatf+1))*dphi
ZU=Heigh_min+(np.linspace(0,Nzf-1,Nzf)+0.5)*dz
ZW=Heigh_min+(np.linspace(0,Nzf,Nzf+1))*dz




#u=u.T
#v=v.T
##p=p.T
#rho=rho.T
#theta=theta.T
#
#
#u_z=u_z.T
#v_z=v_z.T
#p_z=p_z.T
#rho_z=rho_z.T
#theta_z=theta_z.T
#
#u_phi=u_phi.T
#v_phi=v_phi.T
#p_phi=p_phi.T
#rho_phi=rho_phi.T
#theta_phi=theta_phi.T



fu=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),u)
fv=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),v)
fp=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),p)
fw=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),w)
ftheta=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),theta)
frho=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),rho)



fdu_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),u_long)
fdv_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),v_long)
fdp_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),p_long)
fdw_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),w_long)
fdtheta_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),theta_long)
fdrho_dlong=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),rho_long)


fdu_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),u_phi)
fdv_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),v_phi)
fdp_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),p_phi)
fdw_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),w_phi)
fdtheta_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),theta_phi)
fdrho_dphi=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),rho_phi)


fdu_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),u_z)
fdv_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),v_z)
fdp_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),p_z)
fdw_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),w_z)
fdtheta_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),theta_z)
fdrho_dz=interpolate.RegularGridInterpolator((long[:],lat[:],height[:]),rho_z)



CoordU=np.meshgrid(XU,YU,ZU,indexing='ij')
CoordU=np.array(CoordU).T
#lol=fu(CoordU).T
CoordV=np.meshgrid(XV,YV,ZU,indexing='ij')
CoordV=np.array(CoordV).T

CoordW=np.meshgrid(XV,YU,ZW,indexing='ij')
CoordW=np.array(CoordW).T

CoordW2=np.meshgrid(XV,YU,ZW[1:-1],indexing='ij')
CoordW2=np.array(CoordW2).T

CoordP=np.meshgrid(XV,YU,ZU,indexing='ij')
CoordP=np.array(CoordP).T


u_u=fu(CoordU)
u_v=fu(CoordV)
u_w=fu(CoordW)
u_p=fu(CoordP)

v_u=fv(CoordU)
v_v=fv(CoordV)
v_w=fv(CoordW)
v_p=fv(CoordP)

p_u=fp(CoordU)
p_v=fp(CoordV)
p_w=fp(CoordW)
p_p=fp(CoordP)

w_u=fw(CoordU)
w_v=fw(CoordV)
w_w=fw(CoordW)
w_p=fw(CoordP)

theta_u=ftheta(CoordU)
theta_v=ftheta(CoordV)
theta_w=ftheta(CoordW)
theta_p=ftheta(CoordP)

rho_u=frho(CoordU)
rho_v=frho(CoordV)
rho_w=frho(CoordW)
rho_p=frho(CoordP)



du_dlong_u=fdu_dlong(CoordU)
du_dlong_v=fdu_dlong(CoordV)
du_dlong_w=fdu_dlong(CoordW)
du_dlong_p=fdu_dlong(CoordP)

dv_dlong_u=fdv_dlong(CoordU)
dv_dlong_v=fdv_dlong(CoordV)
dv_dlong_w=fdv_dlong(CoordW)
dv_dlong_p=fdv_dlong(CoordP)

dp_dlong_u=fdp_dlong(CoordU)
dp_dlong_v=fdp_dlong(CoordV)
dp_dlong_w=fdp_dlong(CoordW)
dp_dlong_p=fdp_dlong(CoordP)

dw_dlong_u=fdw_dlong(CoordU)
dw_dlong_v=fdw_dlong(CoordV)
dw_dlong_w=fdw_dlong(CoordW)
dw_dlong_p=fdw_dlong(CoordP)

dtheta_dlong_u=fdtheta_dlong(CoordU)
dtheta_dlong_v=fdtheta_dlong(CoordV)
dtheta_dlong_w=fdtheta_dlong(CoordW)
dtheta_dlong_p=fdtheta_dlong(CoordP)

drho_dlong_u=fdrho_dlong(CoordU)
drho_dlong_v=fdrho_dlong(CoordV)
drho_dlong_w=fdrho_dlong(CoordW)
drho_dlong_p=fdrho_dlong(CoordP)




du_dphi_u=fdu_dphi(CoordU)
du_dphi_v=fdu_dphi(CoordV)
du_dphi_w=fdu_dphi(CoordW)
du_dphi_p=fdu_dphi(CoordP)

dv_dphi_u=fdv_dphi(CoordU)
dv_dphi_v=fdv_dphi(CoordV)
dv_dphi_w=fdv_dphi(CoordW)
dv_dphi_p=fdv_dphi(CoordP)

dp_dphi_u=fdp_dphi(CoordU)
dp_dphi_v=fdp_dphi(CoordV)
dp_dphi_w=fdp_dphi(CoordW)
dp_dphi_p=fdp_dphi(CoordP)

dw_dphi_u=fdw_dphi(CoordU)
dw_dphi_v=fdw_dphi(CoordV)
dw_dphi_w=fdw_dphi(CoordW)
dw_dphi_p=fdw_dphi(CoordP)

dtheta_dphi_u=fdtheta_dphi(CoordU)
dtheta_dphi_v=fdtheta_dphi(CoordV)
dtheta_dphi_w=fdtheta_dphi(CoordW)
dtheta_dphi_p=fdtheta_dphi(CoordP)

drho_dphi_u=fdrho_dphi(CoordU)
drho_dphi_v=fdrho_dphi(CoordV)
drho_dphi_w=fdrho_dphi(CoordW)
drho_dphi_p=fdrho_dphi(CoordP)




du_dz_u=fdu_dz(CoordU)
du_dz_v=fdu_dz(CoordV)
du_dz_w=fdu_dz(CoordW)
du_dz_p=fdu_dz(CoordP)

dv_dz_u=fdv_dz(CoordU)
dv_dz_v=fdv_dz(CoordV)
dv_dz_w=fdv_dz(CoordW)
dv_dz_p=fdv_dz(CoordP)

dp_dz_u=fdp_dz(CoordU)
dp_dz_v=fdp_dz(CoordV)
dp_dz_w=fdp_dz(CoordW)
dp_dz_p=fdp_dz(CoordP)

dw_dz_u=fdw_dz(CoordU)
dw_dz_v=fdw_dz(CoordV)
dw_dz_w=fdw_dz(CoordW)
dw_dz_p=fdw_dz(CoordP)

dtheta_dz_u=fdtheta_dz(CoordU)
dtheta_dz_v=fdtheta_dz(CoordV)
dtheta_dz_w=fdtheta_dz(CoordW)
dtheta_dz_p=fdtheta_dz(CoordP)

drho_dz_u=fdrho_dz(CoordU)
drho_dz_v=fdrho_dz(CoordV)
drho_dz_w=fdrho_dz(CoordW)
drho_dz_p=fdrho_dz(CoordP)

if (no_slip) :
    Q_theta_w = np.zeros((Nlongf,Nlatf,Nzf))
    Q_p_p = np.zeros((Nlongf,Nlatf,Nzf))
else :
    Q_theta_w = np.zeros((Nlongf,Nlatf,Nzf-1))
    Q_p_p = np.zeros((Nlongf,Nlatf,Nzf))

for i in range(Nlongf) :
    for j in range(Nlatf) : 
        if (no_slip) :
            for k in range(Nzf) :
                Q_theta_w[i,j,k] = Q(XV[i],YU[j],p_w[k+1,0,0])*g0*r*r/(r+ZW[k+1])**2* \
                rho_w[k+1,0,0]/(theta_w[k+1,0,0]*(p_w[k+1,0,0]/p0)**(gascons/cp))*\
                np.exp(min(0,Nzf-Nzf/5.-k))
                
                Q_p_p[i,j,k] =  Q(XV[i],YU[j],p_p[k,0,0])* gascons*rho_p[k,j,i]*cp/(cp-gascons)*np.exp(min(0,Nzf-Nzf/5.-k))      
        
        else :
            for k in range(Nzf-1) :
                Q_theta_w[i,j,k] = Q(XV[i],YU[j],p_w[k+1,0,0])*g0*r*r/(r+ZW[k+1])**2* \
                rho_w[k+1,0,0]/(theta_w[k+1,0,0]*(p_w[k+1,0,0]/p0)**(gascons/cp))*\
                np.exp(min(0,Nzf-Nzf/5.-k))
                
                Q_p_p[i,j,k] =  Q(XV[i],YU[j],p_p[k,0,0])* gascons*rho_p[k,j,i]*cp/(cp-gascons)*np.exp(min(0,Nzf-Nzf/5.-k))      
            
            Q_p_p[i,j,Nzf-1] = Q(XV[i],YU[j],p_p[Nzf-1,0,0])* gascons*rho_p[Nzf-1,j,i]*cp/(cp-gascons)*np.exp(min(0,Nzf-Nzf/5.-k))





if (no_slip) :
    Qres=np.append(np.zeros(Nlongf*(Nlatf+(Nlatf+1))*Nzf),Q_p_p)
    Qres=np.append(Qres,np.zeros(Nlongf*Nlatf*(Nzf)))
    Qres=np.append(Qres,Q_theta_w)

else : 
    Qres=np.append(np.zeros(Nlongf*(Nlatf+(Nlatf+1))*Nzf),Q_p_p)
    Qres=np.append(Qres,np.zeros(Nlongf*Nlatf*(Nzf-1)))
    Qres=np.append(Qres,Q_theta_w)

res=np.append(u_u.T,u_v.T)
res=np.append(res,u_p.T)
res=np.append(res,u_w.T)

res=np.append(res,v_u.T)
res=np.append(res,v_v.T)
res=np.append(res,v_p.T)
res=np.append(res,v_w.T)

res=np.append(res,p_u.T)
res=np.append(res,p_v.T)
res=np.append(res,p_p.T)
res=np.append(res,p_w.T)

res=np.append(res,w_u.T)
res=np.append(res,w_v.T)
res=np.append(res,w_p.T)
res=np.append(res,w_w.T)

res=np.append(res,theta_u.T)
res=np.append(res,theta_v.T)
res=np.append(res,theta_p.T)
res=np.append(res,theta_w.T)

res=np.append(res,rho_u.T)
res=np.append(res,rho_v.T)
res=np.append(res,rho_p.T)
res=np.append(res,rho_w.T)
    
        
        
dures=np.append(du_dlong_u.T,du_dphi_u.T)
dures=np.append(dures,du_dz_u.T)
dures=np.append(dures,dp_dlong_u.T)
dures=np.append(dures,drho_dlong_u.T)
dures=np.append(dures,drho_dphi_u.T)
dures=np.append(dures,drho_dz_u.T)


dvres=np.append(dv_dlong_v.T,dv_dphi_v.T)
dvres=np.append(dvres,dv_dz_v.T)
dvres=np.append(dvres,dp_dphi_v.T)
dvres=np.append(dvres,drho_dlong_v.T)
dvres=np.append(dvres,drho_dphi_v.T)
dvres=np.append(dvres,drho_dz_v.T)



dpres=np.append(du_dlong_p.T,dv_dphi_p.T)
dpres=np.append(dpres,dp_dlong_p.T)
dpres=np.append(dpres,dp_dphi_p.T)
dpres=np.append(dpres,dw_dz_p.T)
dpres=np.append(dpres,dtheta_dz_p.T)
dpres=np.append(dpres,drho_dlong_p.T)
dpres=np.append(dpres,drho_dphi_p.T)



dwres=np.append(dp_dz_w.T,dw_dlong_w.T)
dwres=np.append(dwres,dw_dphi_w.T)
dwres=np.append(dwres,dw_dz_w.T)
dwres=np.append(dwres,dtheta_dlong_w.T)
dwres=np.append(dwres,dtheta_dphi_w.T)
dwres=np.append(dwres,dtheta_dz_w.T)
dwres=np.append(dwres,drho_dlong_w.T)
dwres=np.append(dwres,drho_dphi_w.T)
dwres=np.append(dwres,drho_dz_w.T)

Qdata=np.append(Q_p_p.T,Q_theta_w.T)


file=open(out_dir+'data.dat','w')
dufile=open(out_dir+'dudata.dat','w')
dvfile=open(out_dir+'dvdata.dat','w')
dpfile=open(out_dir+'dpdata.dat','w')
dwfile=open(out_dir+'dwdata.dat','w')
qfile=open(out_dir+'qdata.dat','w')
heatingfile=open(out_dir+'heating.dat','w')

np.savetxt(out_dir+'data.dat',res, fmt='%.8e')
np.savetxt(out_dir+'dudata.dat',dures, fmt='%.8e')
np.savetxt(out_dir+'dvdata.dat',dvres, fmt='%.8e')
np.savetxt(out_dir+'dpdata.dat',dpres, fmt='%.8e')
np.savetxt(out_dir+'dwdata.dat',dwres, fmt='%.8e')
np.savetxt(out_dir+'qdata.dat',Qdata, fmt='%.8e')
np.savetxt(out_dir+'heating.dat',Qres, fmt='%.8e')

file.close()
dufile.close()
dvfile.close()
dpfile.close()
dwfile.close()
qfile.close()
heatingfile.close()







