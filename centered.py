# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 17:49:32 2015

@author: everard
"""
import numpy as np

Layers = 5
Tstart = 0.0
Tend = 3600.0
dt = 0.5
Time = Tend/dt
Transfer = 0.0001
Theta_L = np.arange(283.15,283.10,-0.01)
G_temp = 270.15
K = .10
dn = 1.0 #1 meter
g=9.8
alpha = (np.pi)/24
gamma = 0.01

u = np.empty((Layers,Time),'float')
theta = np.empty((Layers,Time),'float')
Flux_u = np.empty((Layers+1,Time),'float')
Flux_theta = np.empty((Layers+1,Time),'float')
eddy_length = np.empty((Layers),'float')
Kr = np.empty((Layers),'float')
Layer_temp = np.empty((Layers,Time))
bouyancy = np.empty((Layers))

'''Initial conditions'''
u[:,0] = 0.0 #initial wind speed is zero
theta[:,0] = 0.0 #initial theta dviation is zero
Flux_u[:,0] = 0.0 #no initial fluxes
Flux_theta[:,0] = 0.0 # no initial fluxes

bouyancy[:] = -dt*g*np.sin(alpha)/Theta_L[:]
for t in range(0,3600):

    if t==0:
        u[:,t] = 0.0 #initial wind speed is zero
        theta[:,t] = 0.0 #initial theta dviation is zero
        Flux_u[:,t] = 0.0 #no initial fluxes
        Flux_theta[:,t] = 0.0 # no initial fluxes
        Layer_temp[:,t] = Theta_L[:]
    else:        
        Flux_u[-1,t] = 0.0 #always have zero flux top
        Flux_u[0,t] = 0.0 #no wind flux at ground
        Flux_theta[-1,t] = 0.0 #always have zero flux top
        Flux_theta[0,t] = -Transfer*(Layer_temp[0,t-1]-G_temp) #bulk aerodynamic method   
        Flux_u[1:-2,t] = K*((u[1:-1,t-1]-u[0:-2,t-1])/(dn))
        Flux_theta[1:-2,t] = K*((theta[1:-1,t-1]-theta[0:-2,t-1])/(dn))    
        u[0,t] = 0.0
        u[1:-1,t] = u[1:-1,t-1] - bouyancy[1:-1]*theta[1:-1,t-1] + (1/dn)*(Flux_u[2:-1,t]-Flux_u[1:-2,t])
        theta[0,t] = Flux_theta[0,t]
        theta[1:-1,t] = theta[1:-1,t-1] - gamma*dt*np.sin(alpha)*u[1:-1,t-1] + (1/dn)*(Flux_theta[2:-1,t]-Flux_theta[1:-2,t])
        
        
        
        
        
        
        '''defining the momentum eddy diffusivity'''
        K_h = np.zeros((6),'float')
        K_h[0] = 0.0
        K_h[1:4] = (length_MetOffice[1:4]**2)*Ri*(1/user.dn)*(y[6:9]-y[5:8])
        K_h[5] = (length_MetOffice[5]**2)*Ri*(1/user.dn)*(0.0-y[9])
#        K_h[1] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[6]-y[5])
#        K_h[2] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[7]-y[6])
#        K_h[3] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[8]-y[7])
#        K_h[4] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[9]-y[8])
#        K_h[5] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(0.0-y[9])
        
        '''defining the parameterization for turbulent stress'''
        '''There is an issue with the 9th layer, potential temperature is increasing
           thus causing the wind speeds to increase dramatically (and move upslope)
           Probably a diffusion mishap, just can't figure out exactly how to handle 
           this layer, it has some communication issues/perhaps my system has
           balance issues and this layer just happens to be the sink?
        '''
        Flux_U = np.zeros(((Varct/2)+1),'float')
        Flux_T = np.empty(((Varct/2)+1),'float')
        Flux_T[0] = user.TransferCoef*((Theta_L[0]-y[0])-G_temp)
        Flux_T[1:4] = -K_h[1:4]*(1/user.dn)*(y[1:4]-y[0:3])
        Flux_T[5] = 0.0
#        Flux_T[1] = -K_h[1]*(1/user.dn)*(y[1]-y[0])
#        Flux_T[2] = -K_h[2]*(1/user.dn)*(y[2]-y[1])
#        Flux_T[3] = -K_h[3]*(1/user.dn)*(y[3]-y[2])
#        Flux_T[4] = -K_h[4]*(1/user.dn)*(y[4]-y[3])
#        Flux_T[5] = 0.0
        Flux_U[0] = -user.TransferCoef*(y[5]**2)
        Flux_U[1:4] = -K_h[1:4]*(1/user.dn)*(y[6:9]-y[5:8])
        Flux_U[5] = 0.0
#        Flux_U[1] = -K_h[1]*(1/user.dn)*(y[6]-y[5])
#        Flux_U[2] = -K_h[2]*(1/user.dn)*(y[7]-y[6])
#        Flux_U[3] = -K_h[3]*(1/user.dn)*(y[8]-y[7])
#        Flux_U[4] = -K_h[4]*(1/user.dn)*(y[9]-y[8])

        '''Now derivatves'''
        f[0:4] = (1/user.dn)*(Flux_T[1:5]-Flux_T[0:4]) + user.gamma*np.sin(alpha)*y[5:9]
        f[5:9] = (1/user.dn)*(Flux_U[1:5]-Flux_U[0:4]) - user.g*np.sin(alpha)*y[0:4]/Theta_L[0:4]
#        f[0] = (1/user.dn)*(Flux_T[1]-Flux_T[0]) + user.gamma*np.sin(alpha)*y[5]
#        f[1] = (1/user.dn)*(Flux_T[2]-Flux_T[1]) + user.gamma*np.sin(alpha)*y[6]
#        f[2] = (1/user.dn)*(Flux_T[3]-Flux_T[2]) + user.gamma*np.sin(alpha)*y[7]
#        f[3] = (1/user.dn)*(Flux_T[4]-Flux_T[3]) + user.gamma*np.sin(alpha)*y[8]
#        f[4] = (1/user.dn)*(Flux_T[5]-Flux_T[4]) + user.gamma*np.sin(alpha)*y[9]
#        f[5] = (1/user.dn)*(Flux_U[1]-Flux_U[0]) - user.g*np.sin(alpha)*y[0]/Theta_L[0] 
#        f[6] = (1/user.dn)*(Flux_U[2]-Flux_U[1]) - user.g*np.sin(alpha)*y[1]/Theta_L[1]
#        f[7] = (1/user.dn)*(Flux_U[3]-Flux_U[2]) - user.g*np.sin(alpha)*y[2]/Theta_L[2]
#        f[8] = (1/user.dn)*(Flux_U[4]-Flux_U[3]) - user.g*np.sin(alpha)*y[3]/Theta_L[3]
#        f[9] = (1/user.dn)*(Flux_U[5]-Flux_U[4]) - user.g*np.sin(alpha)*y[4]/Theta_L[4]  
        
        
        
        
        
K_h = np.zeros((Varct/2)+1,'float')
        K_h[0] = 0.0
        K_h[1] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[6]-y[5])
        K_h[2] = (length_MetOffice[2]**2)*Ri*(1/user.dn)*(y[7]-y[6])
        K_h[3] = (length_MetOffice[3]**2)*Ri*(1/user.dn)*(y[8]-y[7])
        K_h[4] = (length_MetOffice[4]**2)*Ri*(1/user.dn)*(y[9]-y[8])
        K_h[5] = (length_MetOffice[5]**2)*Ri*(1/user.dn)*(0.0-y[9])
        K_h[6] = (length_MetOffice[6]**2)*Ri*(1/user.dn)*(y[6]-y[5])
        K_h[7] = (length_MetOffice[7]**2)*Ri*(1/user.dn)*(y[7]-y[6])
        K_h[8] = (length_MetOffice[8]**2)*Ri*(1/user.dn)*(y[8]-y[7])
        K_h[9] = (length_MetOffice[9]**2)*Ri*(1/user.dn)*(y[9]-y[8])
        K_h[10] = (length_MetOffice[10]**2)*Ri*(1/user.dn)*(0.0-y[9])
        
        '''defining the parameterization for turbulent stress'''
        '''There is an issue with the 9th layer, potential temperature is increasing
           thus causing the wind speeds to increase dramatically (and move upslope)
           Probably a diffusion mishap, just can't figure out exactly how to handle 
           this layer, it has some communication issues/perhaps my system has
           balance issues and this layer just happens to be the sink?
        '''
        Flux_U = np.zeros(((Varct/2)+1),'float')
        Flux_T = np.empty(((Varct/2)+1),'float')
        Flux_T[0] = user.TransferCoef*((Theta_L[0]-y[0])-G_temp)
        Flux_T[1] = -K_h[1]*(1/user.dn)*(y[1]-y[0])
        Flux_T[2] = -K_h[2]*(1/user.dn)*(y[2]-y[1])
        Flux_T[3] = -K_h[3]*(1/user.dn)*(y[3]-y[2])
        Flux_T[4] = -K_h[4]*(1/user.dn)*(y[4]-y[3])
        Flux_T[5] = -K_h[5]*(1/user.dn)*(y[5]-y[4])
        Flux_T[6] = -K_h[6]*(1/user.dn)*(y[6]-y[5])
        Flux_T[7] = -K_h[7]*(1/user.dn)*(y[7]-y[6])
        Flux_T[8] = -K_h[8]*(1/user.dn)*(y[8]-y[7])
        Flux_T[9] = -K_h[9]*(1/user.dn)*(y[9]-y[8])
        Flux_T[10] = 0.0
        Flux_U[0] = -user.TransferCoef*(y[10]**2)
        Flux_U[1] = -K_h[1]*(1/user.dn)*(y[11]-y[10])
        Flux_U[2] = -K_h[2]*(1/user.dn)*(y[12]-y[11])
        Flux_U[3] = -K_h[3]*(1/user.dn)*(y[13]-y[12])
        Flux_U[4] = -K_h[4]*(1/user.dn)*(y[14]-y[13])
        Flux_U[5] = -K_h[5]*(1/user.dn)*(y[15]-y[14])
        Flux_U[6] = -K_h[6]*(1/user.dn)*(y[16]-y[15])
        Flux_U[7] = -K_h[7]*(1/user.dn)*(y[17]-y[16])
        Flux_U[8] = -K_h[8]*(1/user.dn)*(y[18]-y[17])
        Flux_U[9] = -K_h[9]*(1/user.dn)*(y[19]-y[18])
        Flux_U[10] = 0.0

        '''Now derivatves'''
        f[0] = (1/user.dn)*(Flux_T[1]-Flux_T[0]) + user.gamma*np.sin(alpha)*y[10]
        f[1] = (1/user.dn)*(Flux_T[2]-Flux_T[1]) + user.gamma*np.sin(alpha)*y[11]
        f[2] = (1/user.dn)*(Flux_T[3]-Flux_T[2]) + user.gamma*np.sin(alpha)*y[12]
        f[3] = (1/user.dn)*(Flux_T[4]-Flux_T[3]) + user.gamma*np.sin(alpha)*y[13]
        f[4] = (1/user.dn)*(Flux_T[5]-Flux_T[4]) + user.gamma*np.sin(alpha)*y[14]
        f[5] = (1/user.dn)*(Flux_T[6]-Flux_T[5]) + user.gamma*np.sin(alpha)*y[15]
        f[6] = (1/user.dn)*(Flux_T[7]-Flux_T[6]) + user.gamma*np.sin(alpha)*y[16]
        f[7] = (1/user.dn)*(Flux_T[8]-Flux_T[7]) + user.gamma*np.sin(alpha)*y[17]
        f[8] = (1/user.dn)*(Flux_T[9]-Flux_T[8]) + user.gamma*np.sin(alpha)*y[18]
        f[9] = (1/user.dn)*(Flux_T[10]-Flux_T[9]) + user.gamma*np.sin(alpha)*y[19]
        
        
        
        f[10] = (1/user.dn)*(Flux_U[1]-Flux_U[0]) - user.g*np.sin(alpha)*y[0]/Theta_L[0] 
        f[11] = (1/user.dn)*(Flux_U[2]-Flux_U[1]) - user.g*np.sin(alpha)*y[1]/Theta_L[1]
        f[12] = (1/user.dn)*(Flux_U[3]-Flux_U[2]) - user.g*np.sin(alpha)*y[2]/Theta_L[2]
        f[13] = (1/user.dn)*(Flux_U[4]-Flux_U[3]) - user.g*np.sin(alpha)*y[3]/Theta_L[3]
        f[14] = (1/user.dn)*(Flux_U[5]-Flux_U[4]) - user.g*np.sin(alpha)*y[4]/Theta_L[4]
        f[15] = (1/user.dn)*(Flux_U[6]-Flux_U[5]) - user.g*np.sin(alpha)*y[5]/Theta_L[5] 
        f[16] = (1/user.dn)*(Flux_U[7]-Flux_U[6]) - user.g*np.sin(alpha)*y[6]/Theta_L[6]
        f[17] = (1/user.dn)*(Flux_U[8]-Flux_U[7]) - user.g*np.sin(alpha)*y[7]/Theta_L[7]
        f[18] = (1/user.dn)*(Flux_U[9]-Flux_U[8]) - user.g*np.sin(alpha)*y[8]/Theta_L[8]
        f[19] = (1/user.dn)*(Flux_U[10]-Flux_U[9]) - user.g*np.sin(alpha)*y[9]/Theta_L[9]          
        return f





