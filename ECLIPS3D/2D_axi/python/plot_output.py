# -*- coding: utf-8 -*-
"""
Created on Fri May 27 21:02:52 2016

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



out_dir = '/Users/florian/Desktop/ECLIPS3D/ECLIPS3D/2D_axi/data/'
#infile=open('/Users/florian/Desktop/MyWork/2D_sca/data/normal_modes/rho_cs_ns.dat')
infile=open(out_dir+'rho_cs_ns.dat')

data=infile.read()
infile.close()

#Remove the blank spaces
values=data.split()

Nlong=200
Nlat=int(values[0])
Nz=int(values[1])
nz=Nz
Heig=float(values[2])

dphi=np.pi/(2.0*(Nlat-0.5))
dz=Heig/(Nz)

XU=(np.linspace(0,Nlat-1,Nlat)+0.5)*dphi
XV=(np.linspace(0,Nlat,Nlat+1)+0.5)*dphi
ZU=(np.linspace(0,Nz-1,Nz)+0.5)*dz
ZW=(np.linspace(1,Nz-1,Nz-1))*dz


rho_u=np.zeros((Nlat,Nz))
rho_v=np.zeros((Nlat+1,Nz))
rho_w=np.zeros((Nlat,Nz+1))

c_p=np.zeros((Nlat,Nz))
n_w=np.zeros((Nlat,Nz+1))


u=np.zeros((Nlat,Nz))
v=np.zeros((Nlat+1,Nz))
p=np.zeros((Nlat,Nz))
theta=np.zeros((Nlat,Nz-1))
w=np.zeros((Nlat,Nz-1))

plop=np.zeros((Nlat,Nz))

u_r=np.zeros((Nlat,Nz))
v_r=np.zeros((Nlat+1,Nz))
p_r=np.zeros((Nlat,Nz))
theta_r=np.zeros((Nlat,Nz-1))
w_r=np.zeros((Nlat,Nz-1))

u_i=np.zeros((Nlat,Nz))
v_i=np.zeros((Nlat+1,Nz))
p_i=np.zeros((Nlat,Nz))
theta_i=np.zeros((Nlat,Nz-1))
w_i=np.zeros((Nlat,Nz-1))



a=3
for k in range(Nz) :
    for j in range(Nlat) :
        rho_u[j,k]=np.array(values[a]).astype(np.float)  
        a=a+1
            
for k in range(Nz) :
    for j in range(Nlat+1) : 
        rho_v[j,k]=np.array(values[a]).astype(np.float)  
        a=a+1
            
for k in range(Nz+1) :
    for j in range(Nlat) :
        rho_w[j,k]=np.array(values[a]).astype(np.float)  
        a=a+1
            

            
for k in range(Nz) :
    for j in range(Nlat) :
        c_p[j,k]=np.array(values[a]).astype(np.float)  
        a=a+1
            
for k in range(Nz+1) :
    for j in range(Nlat) :
        n_w[j,k]=np.array(values[a]).astype(np.float)  
        a=a+1
        
t=0
n_auto=0
stop=False

lol=360./Nlong
XXX=np.linspace(0,Nlong-1,Nlong)*lol

a=0

infile2=open(out_dir+'selected_modes.dat')
data2=infile2.read()
infile2.close()

#Remove the blank spaces
values2=data2.split()

#loop over the results, stopped by user or at the end of the file
while (a<len(values2) and stop==False) :
    t=t+1
    sigmar=float(values2[a])
    sigmai=float(values2[a+1])
    a=a+2
    
    
    long= False #If True, plot in lat-long with wavenumber zm in longitude
    tot=True# If True,plot lat-long and height-long and planet 
    tot=False
    zm=1        
    posz=0 #altitude of the point to plot in lat-long coordinate           
    posl=35 # latitude of the point to plot in height long
    
    if (n_auto==10) :
        stop=True

    if ((t>=0) and (stop==False)) :# and (np.abs(sigmai)>1.0e-5)) :
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
            a=a+2*(2*Nlat*Nz+(Nlat+1)*Nz+2*Nlat*(Nz-1))

        elif(choice=='y') :
            
            for i in range(Nlat) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        u_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        u_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat+1) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        v_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        v_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        p_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        p_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz-1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        w_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        w_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz-1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        theta_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        theta_r[i,j]=0.0
                    a=a+1
                    
                    
            #Imaginary part
                    
                    

            for i in range(Nlat) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        u_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        u_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat+1) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        v_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        v_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        p_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        p_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz-1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        w_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        w_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) :
                for j in range(Nz-1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        theta_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        theta_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlat) : 
                for j in range(Nz) :
                    u[i,j]=u_r[i,j]
                    p[i,j]=p_r[i,j]
            
            for i in range(Nlat) : 
                for j in range(Nz) :
                    v[i,j]=v_i[i,j]
                    
            for i in range(Nlat) : 
                for j in range(Nz-1) :
                    w[i,j]=w_i[i,j]
                    theta[i,j]=theta_r[i,j]

            #Rescaling the variable as in Thuburn et al. (2002 I)
            u_r/=np.sqrt(rho_u)
            v_r/=np.sqrt(rho_v)        
            p_r/=np.sqrt(rho_u*c_p)
            w_r/=np.sqrt(rho_w[:,1:-1])       
            theta_r/=np.sign(n_w[:,1:-1])*np.sqrt(abs(rho_w[:,1:-1]*n_w[:,1:-1]))
            
            u_i/=np.sqrt(rho_u)
            v_i/=np.sqrt(rho_v)        
            p_i/=np.sqrt(rho_u*c_p)
            w_i/=np.sqrt(rho_w[:,1:-1])       
            theta_i/=np.sign(n_w[:,1:-1])*np.sqrt(abs(rho_w[:,1:-1]*n_w[:,1:-1]))
            
            u/=np.sqrt(rho_u)
            v/=np.sqrt(rho_v)        
            p/=np.sqrt(rho_u*c_p)
            w/=np.sqrt(rho_w[:,1:-1])       
            theta/=np.sign(n_w[:,1:-1])*np.sqrt(abs(rho_w[:,1:-1]*n_w[:,1:-1]))
            
            #Defining horizontal velocity maximum to be 1
            hv=np.abs((u_r+1j*u_i)*(u_r-1j*u_i)+(v_r[1:]+1j*v_i[1:])*(v_r[1:]-1j*v_i[1:]))[5:-5]
            
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
            
            

            #Interpolating
            
            #Interpolation factor   
            Intf=1
            
            if (Intf>1) :
            
                x=np.linspace(dphi,90,Nlat*Intf)  
                x_v=np.linspace(0,90-dphi,Nlat*Intf)
                y=np.linspace(0,Heig-dz,Nz*Intf)
                y_w=np.linspace(dz,Heig-dz,(Nz-1)*Intf)
                fu=interpolate.interp2d(XU,ZU,u.T)
                zu=fu(x,y)
                
                fv=interpolate.interp2d(XU,ZU,v.T)
                zv=fv(x_v,y)
                
                fp=interpolate.interp2d(XU,ZU,p.T)
                zp=fp(x,y)
                
                fw=interpolate.interp2d(XU,ZW,w.T)
                zw=fw(x,y_w)
                
                ft=interpolate.interp2d(XU,ZW,theta.T)
                zt=ft(x,y_w)
                
        
                
                #Plotting
                
                fig=plt.figure()
                gs=gridspec.GridSpec(3,2)      
                
                figu=fig.add_subplot(gs[0,0])
                pu=figu.pcolormesh(x,y,zu)
                plt.colorbar(pu)        
                
                plt.setp(figu,title='u')                
                
                figv=fig.add_subplot(gs[0,1])
                pv=figv.pcolormesh(x_v,y,zv)
                plt.colorbar(pv)
                
                figw=fig.add_subplot(gs[2,0])
                pw=figw.pcolormesh(x,y_w,zw)
                plt.colorbar(pw)
                
                figp=fig.add_subplot(gs[1,0])
                pp=figp.pcolormesh(x,y,zp)
                plt.colorbar(pp)
                
                figt=fig.add_subplot(gs[1,1])
                pt=figt.pcolormesh(x,y_w,zt)
                plt.colorbar(pt)
                plt.show()
                
                fig.suptitle('Freq =' + str(sigmar) ) #'   kz=' + str(zerovert) + '   klat=' + str(zerolat))
            
            
            elif (tot) :
                XU=XU
                plt.figure()
  
                
                uu=np.zeros((Nlat,Nlong))
                for m in range(Nlong) :
                    for n in range(Nlat) :
                        uu[n,m]=u_r[n,posz]*np.cos(zm*XXX[m]*np.pi/180)-u_i[n,posz]*np.sin(zm*XXX[m]*np.pi/180)
                        
                        
                vv=np.zeros((Nlat,Nlong))
                for m in range(Nlong) :
                    for n in range(Nlat) :
                        vv[n,m]=v_r[n,posz]*np.cos(zm*XXX[m]*np.pi/180)-v_i[n,posz]*np.sin(zm*XXX[m]*np.pi/180)
                        
                ppp=np.zeros((Nlat,Nlong))
                for m in range(Nlong) :
                    for n in range(Nlat) :
                        ppp[n,m]=p_r[n,posz]*np.cos(zm*XXX[m]*np.pi/180)-p_i[n,posz]*np.sin(zm*XXX[m]*np.pi/180)
                        
                tt=np.zeros((Nlat,Nlong))
                for m in range(Nlong) :
                    for n in range(Nlat) :
                        tt[n,m]=theta_r[n,posz]*np.cos(zm*XXX[m]*np.pi/180)-theta_i[n,posz]*np.sin(zm*XXX[m]*np.pi/180)
                        
                
                plt.pcolormesh(XXX,XU*180/np.pi,ppp[:Nlat])
                plt.colorbar()
                plt.quiver(XXX[::2],XU[::2]*180/np.pi,uu[:Nlat:2,::2],vv[:Nlat:2,::2],scale=3,scale_units='inches',pivot='tail')

                plt.title('p perturbations in long-lat, Freq=' + str(sigmar))
                
                plt.figure()
                plt.pcolormesh(XXX,XU*180/np.pi,tt[:Nlat])
                plt.colorbar()
                plt.quiver(XXX[::2],XU[::2]*180/np.pi,uu[:Nlat:2,::2],vv[:Nlat:2,::2],scale=3,scale_units='inches',pivot='tail')
                plt.title('Theta perturbations in long-lat, Freq=' + str(sigmar))

                plt.figure()
  
                
                uz=np.zeros((Nz,Nlong))
                for m in range(Nlong) :
                    for n in range(Nz) :
                        uz[n,m]=u_r[posl,n]*np.cos(zm*XXX[m]*np.pi/180)-u_i[posl,n]*np.sin(zm*XXX[m]*np.pi/180)
                        
                pz=np.zeros((Nz,Nlong))
                for m in range(Nlong) :
                    for n in range(Nz) :
                        pz[n,m]=p_r[posl,n]*np.cos(zm*XXX[m]*np.pi/180)-p_i[posl,n]*np.sin(zm*XXX[m]*np.pi/180)
                        
                wz=np.zeros((Nz-1,Nlong))
                for m in range(Nlong) :
                    for n in range(Nz-1) :
                        wz[n,m]=w_r[posl,n]*np.cos(zm*XXX[m]*np.pi/180)-w_i[posl,n]*np.sin(zm*XXX[m]*np.pi/180)
                        
                tz=np.zeros((Nz-1,Nlong))
                for m in range(Nlong) :
                    for n in range(Nz-1) :
                        tz[n,m]=theta_r[posl,n]*np.cos(zm*XXX[m]*np.pi/180)-theta_i[posl,n]*np.sin(zm*XXX[m]*np.pi/180)
                        
 
                plt.pcolormesh(XXX,ZU,pz)
                plt.colorbar()
                plt.quiver(XXX[::2],ZU[:-1],uz[:-1,::2],wz[:,::2],scale=3*np.max([np.max(uz),np.max(wz)]),scale_units='inches',pivot='tail')
                plt.title('p perturbations in long-height, Freq=' + str(sigmar))
                
                
                plt.figure()
                plt.pcolormesh(XXX,ZU[:-1],tz)
                plt.colorbar()
                plt.quiver(XXX[::2],ZU[:-1],uz[:-1,::2],wz[:,::2],scale=3*np.max([np.max(uz),np.max(wz)]),scale_units='inches',pivot='tail')
                plt.title('Theta perturbations in long-height, Freq=' + str(sigmar))
                

                planet_NS.planet_NS(Nlat,Nlong,XXX,XU*180/np.pi,uu,vv,ppp,tt)
                
            else : 
                fig=plt.figure()
                gs=gridspec.GridSpec(3,2)
                gs.update(wspace=0.3, hspace=0.5)
                   
                        
                figu=fig.add_subplot(gs[0,0])
                pu=figu.pcolor(XU[5:-5]*180/np.pi,ZU/1000.,u[5:-5].T, cmap='jet')
                plt.colorbar(pu) 
                plt.setp(figu,xlabel='Latitude (deg)')
                plt.setp(figu,ylabel='Height (km)')
                
                plt.setp(figu,title='u')                    
                
                figv=fig.add_subplot(gs[0,1])
                pv=figv.pcolor(XV[5:-5]*180/np.pi,ZU/1000.,v[5:-5].T, cmap='jet')
                plt.colorbar(pv)
                
                plt.setp(figv,title='v')   
                plt.setp(figv,xlabel='Latitude (deg)')
                plt.setp(figv,ylabel='Height (km)')                 
                
                
                figw=fig.add_subplot(gs[2,0])
                pw=figw.pcolormesh(XU[5:-5]*180/np.pi,ZW/1000.,w[5:-5].T, cmap='jet')
                plt.colorbar(pw)
                
                plt.setp(figw,title='w')
                plt.setp(figw,xlabel='Latitude (deg)')
                plt.setp(figw,ylabel='Height (km)')     
                
                figp=fig.add_subplot(gs[1,0])
                pp=figp.pcolormesh(XU[5:-5]*180/np.pi,ZU[1:]/1000.,p[5:-5,1:].T, cmap='jet')
                plt.colorbar(pp)
                #pk=figp.quiver(ZU,XU,u,v)
                
                plt.setp(figp,title='p')  
                plt.setp(figp,xlabel='Latitude (deg)')
                plt.setp(figp,ylabel='Height (km)')     
                
                
                figt=fig.add_subplot(gs[1,1])
                pt=figt.pcolormesh(XU[5:-5]*180/np.pi,ZW/1000.,theta[5:-5].T, cmap='jet')
                plt.colorbar(pt)
                plt.show()
                
                plt.setp(figt,title='theta')  
                plt.setp(figt,xlabel='Latitude (deg)')
                plt.setp(figt,ylabel='Height (km)')     
                
#                fig.suptitle('Freq =' + str(sigmar) ) #'   kz=' + str(zerovert) + '   klat=' + str(zerolat))
        else :
            stop=True
    else :
        a=a+2*(2*Nlat*Nz+(Nlat+1)*Nz+2*Nlat*(Nz-1))

   