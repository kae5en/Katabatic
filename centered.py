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
        





