# -*- coding: utf-8 -*-
"""
Created on Mon May 16 17:00:52 2016

@author: florian
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 11 17:50:50 2016

@author: florian
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import gridspec


#Directory of files to plot
rep = '../data/'

#Name of steady initial thermodynamic profile and results.
# Defaults : rho_cs_ns.dat and selected_modes.dat

infile=open(rep + 'rho_cs_ns.dat')
infile2=open(rep + 'selected_modes.dat')

data=infile.read()
infile.close()


data2=infile2.read()
infile2.close()

#Remove the blank spaces
values2=data2.split()

#Remove the blank spaces
values=data.split()


#Read the number of points
Nlong=int(values[0])
Nlat=int(values[1])
Nz=int(values[2])
nz=Nz
Heig=float(values[3])


heigh_min=0.0
heigh_max=80.
#heigh_max=8930000.0

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




kill_lat=0 #number of points not to display in latitude
kill_z=0 #number of points not to display in latitude



a=4
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

#rho_u=rho_u[:,:,:Nz-kill_z]
#rho_v=rho_v[:,:,:Nz-kill_z]
#rho_w=rho_w[:,:,:Nz-1-kill_z]
#rho_p=rho_p[:,:,:Nz-kill_z]
#
#c_p=c_p[:,:,:Nz-kill_z]
#n_w=n_w[:,:,:Nz-1-kill_z]

        
t=0
n_auto=0
stop=False

a=0



#loop over the results, stopped by user or at the end of the file
while (a<len(values2) and stop==False) :
    t=t+1
    sigmar=float(values2[a])
    sigmai=float(values2[a+1])
    a=a+2
    
    
    long= False #If True, plot in lat-long with wavenumber zm in longitude
    tot=True# If True,plot lat-long and height-long and planet 
    tot=False
    
    mom=True #if True, display momentum
    mom=False
    
    mom_vert=True
    mom_vert=False
    
    time_ev=True #if True show the time evolution over two periods
    time_ev=False
    #longitude to display            
    debl=0
    finl=5
    
    #altitude to display
    debz = 16
    finz = 16
    
    
    posz=7
    posl=40
    num=20

    

    if (n_auto==10) :
        stop=True

    if ((t>=0) and (stop==False)) : #and (np.abs(sigmai)<1.0E-6)) :
        
        u=np.zeros((Nlong,Nlat,Nz))
        v=np.zeros((Nlong,Nlat+1,Nz))
        p=np.zeros((Nlong,Nlat,Nz))
        theta=np.zeros((Nlong,Nlat,Nz))
        w=np.zeros((Nlong,Nlat,Nz))
        
        plop=np.zeros((Nlong,Nlat,Nz))
        
        u_r=np.zeros((Nlong,Nlat,Nz))
        v_r=np.zeros((Nlong,Nlat+1,Nz))
        p_r=np.zeros((Nlong,Nlat,Nz))
        theta_r=np.zeros((Nlong,Nlat,Nz))
        w_r=np.zeros((Nlong,Nlat,Nz))
        
        u_i=np.zeros((Nlong,Nlat,Nz))
        v_i=np.zeros((Nlong,Nlat+1,Nz))
        p_i=np.zeros((Nlong,Nlat,Nz))
        theta_i=np.zeros((Nlong,Nlat,Nz))
        w_i=np.zeros((Nlong,Nlat,Nz))        
        
        print(t)
        n_auto=n_auto+1
        #characteristics of the chosen perturbation    
        print('\n')
        print("Re(Sigma) =" , sigmar, " and Im(Sigma) = ", sigmai)
        #choice=input("Study this mode ? (y for yes, n for no)")
        
        choice='y'
    
            #Some zero or -99999 frequencies also indicate the end of the file
        if ((sigmar==0 or sigmar==-99999) and (sigmai==0 or sigmai==-99999)) :
            stop=True
    
        #Go to the next frequency
        elif (choice=='n') :
            a=a+2*Nlong*(2*Nlat*Nz+(Nlat+1)*Nz+2*Nlat*(Nz))

        elif(choice=='y') :
            
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            u_r[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            u_r[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat+1) :
                    for k in range(Nz) :                   
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            v_r[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            v_r[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            p_r[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            p_r[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            w_r[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            w_r[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            theta_r[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            theta_r[i,j,k]=0.0
                        a=a+1
                    
                    
            #Imaginary part
                    
                    

            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            u_i[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            u_i[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat+1) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            v_i[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            v_i[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            p_i[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            p_i[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            w_i[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            w_i[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    for k in range(Nz) :
                        if (values2[a][8]=='E' or values2[a][9]=='E') :
                            theta_i[i,j,k]=np.array(values2[a]).astype(np.float)
                        else :
                            theta_i[i,j,k]=0.0
                        a=a+1
                    
            for i in range(Nlong) : 
                for j in range(Nlat) :
                    for k in range(Nz) :
                        u[i,j,k]=u_i[i,j,k]+u_r[i,j,k]#np.sqrt(u_r[i,j,k]**2+u_i[i,j,k]**2)
                        p[i,j,k]=p_r[i,j,k]+p_i[i,j,k]#np.sqrt(p_r[i,j,k]**2+p_i[i,j,k]**2)
                        
            for i in range(Nlong) : 
                for j in range(Nlat+1) :
                    for k in range(Nz) :
                        v[i,j,k]=v_r[i,j,k]-v_i[i,j,k]#np.sqrt(v_r[i,j,k]**2+v_i[i,j,k]**2)

                    
            for i in range(Nlong) : 
                for j in range(Nlat) :
                    for k in range(Nz) :
                        w[i,j,k]=w_r[i,j,k]-w_i[i,j,k]#np.sqrt(w_r[i,j,k]**2+w_i[i,j,k]**2)
                        theta[i,j,k]=theta_r[i,j,k]+theta_i[i,j,k]#np.sqrt(theta_r[i,j,k]**2+theta_i[i,j,k]**2)
#            plt.pcolor(u_r[5])
#            u_r=u_r[:,:,:Nz-kill_z]            
#            v_r=v_r[:,:,:Nz-kill_z]
#            p_r=p_r[:,:,:Nz-kill_z]
#            w_r=w_r[:,:,:Nz-1-kill_z]
#            theta_r=theta_r[:,:,:Nz-1-kill_z]
#            
#            u_i=u_i[:,:,:Nz-kill_z]            
#            v_i=v_i[:,:,:Nz-kill_z]
#            p_i=p_i[:,:,:Nz-kill_z]
#            w_i=w_i[:,:,:Nz-1-kill_z]
#            theta_i=theta_i[:,:,:Nz-1-kill_z]
#            
#            u=u[:,:,:Nz-kill_z]            
#            v=v[:,:,:Nz-kill_z]
#            p=p[:,:,:Nz-kill_z]
#            w=w[:,:,:Nz-1-kill_z]
#            theta=theta[:,:,:Nz-1-kill_z]            
            

            
            #Rescaling the variable as in Thuburn et al. (2002 I)
#            proutiplip
            u_r/=np.sqrt(rho_u)
            v_r/=np.sqrt(rho_v)        
            p_r/=np.sqrt(rho_p*c_p)
            w_r/=np.sqrt(rho_w[:,:,1:])       
            theta_r/=(np.sign(n_w[:,:,1:])*np.sqrt(abs(rho_w[:,:,1:]*n_w[:,:,1:])))
            
                     
            
            u_i/=np.sqrt(rho_u)
            v_i/=np.sqrt(rho_v)        
            p_i/=np.sqrt(rho_p*c_p)
            w_i/=np.sqrt(rho_w[:,:,1:])       
            theta_i/=(np.sign(n_w[:,:,1:])*np.sqrt(abs(rho_w[:,:,1:]*n_w[:,:,1:])))
            
            u/=np.sqrt(rho_u)
            v/=np.sqrt(rho_v)        
            p/=np.sqrt(rho_p*c_p)
            w/=np.sqrt(rho_w[:,:,1:])       
            theta/=(np.sign(n_w[:,:,1:])*np.sqrt(abs(rho_w[:,:,1:]*n_w[:,:,1:])))
            
            #Defining horizontal velocity maximum to be 1
            hv=2.*np.abs((u_r+1j*u_i)*(u_r-1j*u_i)+(v_r[:,1:]+1j*v_i[:,1:])*(v_r[:,1:]-1j*v_i[:,1:]))
            
            scale=np.max(hv)
            scale=np.sqrt(scale)
                     

            u_r/=scale
            v_r/=scale  
            p_r/=scale
            w_r/=scale   
            theta_r/=scale
            
            u_i/=scale
            v_i/=scale  
            p_i/=scale
            w_i/=scale   
            theta_i/=scale
            
            u/=scale
            v/=scale  
            p/=scale
            w/=scale   
            theta/=scale
            
            if mom :
                moment=momentum(posz,num)
            
            if mom_vert :
                moment2=momemtum_vert(posl,num)
            
            if time_ev :
                time_evolution(posz,num)

            #Interpolating
            
            #Interpolation factor   
            Intf=1
            
            if (Intf>1) :
            
                y=np.linspace(dphi,90,Nlat*Intf)  
                y_v=np.linspace(0,90,(Nlat+1)*Intf)
                z=np.linspace(0,Heig-dz,Nz*Intf)
                z_w=np.linspace(dz,Heig-dz,(Nz-1)*Intf)
                fu=interpolate.interp2d(YU,ZU,u[2].T)
                zu=fu(y,z)
                
                fv=interpolate.interp2d(YV,ZU,v[2].T)
                zv=fv(y_v,z)
                
                fp=interpolate.interp2d(YU,ZU,p[2].T)
                zp=fp(y,z)
                
                fw=interpolate.interp2d(YU,ZW[1:],w[2].T)
                zw=fw(y,z_w)
                
                ft=interpolate.interp2d(YU,ZW[1:],theta[2].T)
                zt=ft(y,z_w)
                
        
                
                #Plotting
                
                fig=plt.figure()
                gs=gridspec.GridSpec(3,2)      
                
                figu=fig.add_subplot(gs[0,0])
                pu=figu.pcolormesh(y,z,zu)
                plt.colorbar(pu)        
                
                plt.setp(figu,title='u')                
                
                figv=fig.add_subplot(gs[0,1])
                pv=figv.pcolormesh(y_v,z,zv)
                plt.colorbar(pv)
                
                figw=fig.add_subplot(gs[2,0])
                pw=figw.pcolormesh(y,z_w,zw)
                plt.colorbar(pw)
                
                figp=fig.add_subplot(gs[1,0])
                pp=figp.pcolormesh(y,z,zp)
                plt.colorbar(pp)
                
                figt=fig.add_subplot(gs[1,1])
                pt=figt.pcolormesh(y,z_w,zt)
                plt.colorbar(pt)
                plt.show()
                
                fig.suptitle('Freq =' + str(sigmar) ) #'   kz=' + str(zerovert) + '   klat=' + str(zerolat))
            
            
            elif (tot) :
                div=1 # the smaller, the more arrows
                pole=1 # distance from the pole
                time=0 #time for the frequency
                
                uu=u_r*np.cos(time)+u_i*np.sin(np.sign(sigmai)*time)
                vv=v_r*np.cos(time)+v_i*np.sin(np.sign(sigmai)*time)
                pp=p_r*np.cos(time)+p_i*np.sin(np.sign(sigmai)*time)
                
                for posz in range(debz,finz+1) :
                    YU2=YU[:Nlat-pole]
                    plt.figure()
      
    
                    plt.pcolormesh(XV,YU2,pp[:,:Nlat-pole,posz].T)
                    plt.colorbar()
                    plt.quiver(XU[::div],YU2[::div],uu[::div,:Nlat-pole:div,posz].T,vv[::div,:Nlat-pole:div,posz].T,scale=2*np.max([np.max(np.abs(uu[:,:,posz])),np.max(np.abs(vv[:,:,posz]))]),scale_units='inches',pivot='tail')
    
                    plt.title('p perturbations in long-lat, z= ' + str(ZU[posz]))
#                    
#                    plt.figure()
#                    plt.pcolormesh(XV,YU,theta[:,:Nlat-pole,posz].T)
#                    plt.colorbar()
#                    plt.quiver(XU[::div],YU[::div],u[::div,:Nlat-pole:div,posz].T,v[::div,:Nlat-pole:div,posz].T,scale=2*np.max([np.max(u[:,:,posz]),np.max(v[:,:,posz])]),scale_units='inches',pivot='tail')
#                    plt.title('Theta perturbations in long-lat, Freq=' + str(sigmar))
#    
#                    plt.figure()
#  
                
                
 
#                plt.pcolormesh(XV,ZU,p[:,posl].T)
#                plt.colorbar()
#                plt.quiver(XU[::div],ZU[:-1],u[::div,posl,:-1].T,w[::div,posl].T,scale=3*np.max([np.max(u[:,posl]),np.max(w[:,posl])]),scale_units='inches',pivot='tail')
#                plt.title('p perturbations in long-height, Freq=' + str(sigmar))
#                
#                
#                plt.figure()
#                plt.pcolormesh(XV,ZW,theta[:,posl].T)
#                plt.colorbar()
#                plt.quiver(XU[::div],ZU[:-1],u[::div,posl,:-1].T,w[::div,posl].T,scale=3*np.max([np.max(u[:,posl]),np.max(w[:,posl])]),scale_units='inches',pivot='tail')
#                plt.title('Theta perturbations in long-height, Freq=' + str(sigmar))
                

#                    test_planet.planet(Nlat,Nlong,XU,YU[:Nlat-2],u[:,:,posz].T,v[:,:,posz].T,p[:,:,posz].T,theta[:,:,posz].T)
                
            else : 
                for mm in [debl,finl] :
                    fig=plt.figure()
                    gs=gridspec.GridSpec(3,2)
                    gs.update(wspace=0.3, hspace=0.5)
               
                            
                    figu=fig.add_subplot(gs[0,0])
                    pu=figu.pcolormesh(YU[:Nlat-kill_lat],ZU,u[mm,:Nlat-kill_lat].T,cmap='jet')
                    plt.colorbar(pu)        
                    
                    plt.setp(figu,title='u') 
                    plt.setp(figu,xlabel='Latitude (deg)')
                    plt.setp(figu,ylabel='Height (km)')
                    
                    figv=fig.add_subplot(gs[0,1])
                    pv=figv.pcolormesh(YV[:Nlat-kill_lat],ZU,v[mm,:Nlat-kill_lat].T,cmap='jet')
                    plt.colorbar(pv)
                    
                    plt.setp(figv,title='v') 
                    plt.setp(figv,xlabel='Latitude (deg)')
                    plt.setp(figv,ylabel='Height (km)')
                    
                    
                    figw=fig.add_subplot(gs[2,0])
                    pw=figw.pcolormesh(YU[:Nlat-kill_lat],ZW[1:],w[mm,:Nlat-kill_lat].T,cmap='jet')
                    plt.colorbar(pw)
                    
                    plt.setp(figw,title='w')   
                    plt.setp(figw,xlabel='Latitude (deg)')
                    plt.setp(figw,ylabel='Height (km)')                    
                
                    figp=fig.add_subplot(gs[1,0])
                    pp=figp.pcolormesh(YU[:Nlat-kill_lat],ZU,p[mm,:Nlat-kill_lat].T,cmap='jet')
                    plt.colorbar(pp)
                    #pk=figp.quiver(ZU,XU,u,v)
                    
                    plt.setp(figp,title='p')   
                    plt.setp(figp,xlabel='Latitude (deg)')
                    plt.setp(figp,ylabel='Height (km)')                   
                    
                    figt=fig.add_subplot(gs[1,1])
                    pt=figt.pcolormesh(YU[:Nlat-kill_lat],ZW[1:],theta[mm,:Nlat-kill_lat].T,cmap='jet')
                    plt.colorbar(pt)
                    
                    plt.show()
                    
                    plt.setp(figt,title='theta')   
                    plt.setp(figt,xlabel='Latitude (deg)')
                    plt.setp(figt,ylabel='Height (km)')
                    
#                    fig.suptitle('Freq =' + str(-sigmai))# + ' and posl =' + str(mm) ) #'   kz=' + str(zerovert) + '   klat=' + str(zerolat))

                for mm in range(debz,finz+1) :
                    fig=plt.figure()
                    gs=gridspec.GridSpec(3,2)   
                    gs.update(wspace=0.3, hspace=0.5)
               
                            
                    figu=fig.add_subplot(gs[0,0])
                    pu=figu.pcolormesh(XU,YU[:Nlat-kill_lat],u[:,:Nlat-kill_lat,mm].T)
                    plt.colorbar(pu)        
                    
                    plt.setp(figu,title='u')                    
                    
                    figv=fig.add_subplot(gs[0,1])
                    pv=figv.pcolormesh(XV,YV[:Nlat-kill_lat],v[:,:Nlat-kill_lat,mm].T)
                    plt.colorbar(pv)
                    
                    plt.setp(figv,title='v')                    
                    
                    
                    figw=fig.add_subplot(gs[2,0])
                    pw=figw.pcolormesh(XV,YU[:Nlat-kill_lat],w[:,:Nlat-kill_lat,mm].T)
                    plt.colorbar(pw)
                    
                    plt.setp(figw,title='w')   
                    
                    figp=fig.add_subplot(gs[1,0])
                    pp=figp.pcolormesh(XV,YU[:Nlat-kill_lat],p[:,:Nlat-kill_lat,mm].T)
                    plt.colorbar(pp)
                    #pk=figp.quiver(ZU,XU,u,v)
                    
                    plt.setp(figp,title='p')   
                    
                    
                    figt=fig.add_subplot(gs[1,1])
                    pt=figt.pcolormesh(XV,YU[:Nlat-kill_lat],theta[:,:Nlat-kill_lat,mm].T)
                    plt.colorbar(pt)
                    plt.show()
                    
                    plt.setp(figt,title='theta')   
                    
                    fig.suptitle('Freq =' + str(-sigmai))#+ ' and posz =' + str(mm) ) #'   kz
        else :
            stop=True
    else :
        a=a+2*Nlong*(2*Nlat*Nz+(Nlat+1)*Nz+2*Nlat*(Nz))
        
        
#time evolution
# If the frequency is positive, cos + sin
# If the frequency is negative, cos - sin
        




