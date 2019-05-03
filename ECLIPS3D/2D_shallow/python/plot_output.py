# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 09:55:33 2017

@author: florian
"""

import numpy as np
import matplotlib.pyplot as plt


#File to plot directory
rep = '../data/'
#Name of the output file. Default : selected_modes.dat
infile2=open(rep + 'selected_modes.dat')

#First mode to plot
t0 = 20

# Number of modes to plot
number= 10

#Size of the grid, should be useless in the next version
Nlong = 60
nlong = Nlong
Nlat = 30
nlat = Nlat



#Read the file
data2=infile2.read()
infile2.close()

#Remove the blank spaces
values2=data2.split()



###################################################################
#Information characterizing hot Jupiters
# See Showman Polvani 2011
phi = 0.

omega = 2.06E-5
g0=9.42
p0=2.2E7
rtot = 9.44E7


beta = 2.*omega*np.cos(phi)/rtot

H = 4.E6/g0


V_caract = np.sqrt(g0*H)
L_caract = np.sqrt(np.sqrt(g0*H)/beta)
T_caract = 1./(np.sqrt(np.sqrt(g0*H)*beta))

###################################################################



#Initialising variables to read
t=0
n_auto=0
stop=False
a=0



u_r = np.zeros((Nlong,Nlat))
u_i = np.zeros((Nlong,Nlat))

v_r = np.zeros((Nlong,Nlat+1))
v_i = np.zeros((Nlong,Nlat+1))

h_r = np.zeros((Nlong,Nlat))
h_i = np.zeros((Nlong,Nlat))

u =  np.zeros((Nlong,Nlat))
v = np.zeros((Nlong,Nlat+1))
h = np.zeros((Nlong,Nlat))


# Staggered grid
dx = (2.*np.pi*rtot/L_caract)/nlong
dy = (np.pi*rtot/(1.3*L_caract))/nlat

XU = (np.linspace(-nlong/2+1,nlong/2,nlong)-0.5)*dx
XV = (np.linspace(-nlong/2+1,nlong/2,nlong))*dx

YU = (np.linspace(1,nlat,nlat)-0.5)*dy
YV = (np.linspace(0,nlat,nlat+1))*dy






#loop over the results, stopped by user or at the end of the file
while (a<len(values2) and stop==False) :
    t=t+1
    # Frequency and growth rate
    sigmar=float(values2[a])
    sigmai=float(values2[a+1])
    a=a+2
    if (n_auto==number) :
        stop=True
  
    if ((t>=t0) and (stop==False)): #  and (sigmai<=0)) : #and (np.abs(sigmai)>1E-6)) :# and (np.abs(sigmai)>1.0e-5)) :
        print(t)
        print(n_auto)
        n_auto=n_auto+1
        #characteristics of the chosen perturbation    
        print('\n')
        print("Re(Sigma) =" , sigmar, " and Im(Sigma) = ", sigmai)
        choice=input("Plot this mode ? (y for yes, n for no,stop for stopping)")
        
#        choice='y'
#        
        if (choice=='stop') :
            stop=True
    
            #Some zero or -99999 frequencies also indicate the end of the file
        if ((sigmar==0 or sigmar==-99999) and (sigmai==0 or sigmai==-99999)) :
            stop=True
    
        #Go to the next frequency
        elif (choice=='n') :
            a=a+2*(2*Nlat*Nlong+(Nlat+1)*Nlong)

        elif(choice=='y') :
            
            for i in range(Nlong) :
                for j in range(Nlat) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        u_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        u_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat+1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        v_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        v_r[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        h_r[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        h_r[i,j]=0.0
                    a=a+1

                    
                    
            #Imaginary part
                    
                    

            for i in range(Nlong) :
                for j in range(Nlat) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        u_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        u_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat+1) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        v_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        v_i[i,j]=0.0
                    a=a+1
                    
            for i in range(Nlong) :
                for j in range(Nlat) :
                    if (values2[a][8]=='E' or values2[a][9]=='E') :
                        h_i[i,j]=np.array(values2[a]).astype(np.float)
                    else :
                        h_i[i,j]=0.0
                    a=a+1
                    
                    
                    
            for i in range(Nlong) : 
                for j in range(Nlat) :
                    u[i,j]=u_r[i,j]#+u_i[i,j]
                    h[i,j]=h_r[i,j]#+h_i[i,j]
            
            for i in range(Nlong) : 
                for j in range(Nlat+1) :
                    v[i,j]=v_r[i,j]#+v_i[i,j]
              
              
#            if (sigmai < 0) :
#                plt.figure()
#                plt.pcolor(XV,YU,u.T*v[:,:-1].T)
#                plt.plot(4*np.mean(u.T*v[:,:-1].T,axis=1)/ \
#                np.max(np.mean(np.abs(u.T*v[:,:-1].T),axis=1)),YU)
#                plt.title(sigmai)
            
            plt.figure()
            plt.pcolor(XV,YU,h.T)
            plt.quiver(XV[::2],YU,u[::2].T,v[::2,:-1].T)
            plt.title(sigmai)
            

            
    else :
        a=a+2*(2*Nlat*Nlong+(Nlat+1)*Nlong)
        
        
        
        
        


