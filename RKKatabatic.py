# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 18:04:36 2015

@author: Student
"""
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
import yaml
import pandas as pd
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.colorbar as colorbar
import os,glob

def rkck_init():
    # %
    # % initialize the Cash-Karp coefficients
    # % defined in the tableau in lab 4,
    # % section "Embedded Runge Kutta"
    # %
    a = np.array([0.2, 0.3, 0.6, 1.0, 0.875])
  # c1 coefficients for the fifth order scheme
    c1 = np.array([37.0 / 378.0, 0.0, 250.0 / 621.0,
                   125.0 / 594.0, 0.0, 512.0 / 1771.0])
  # c2=c* coefficients for the fourth order schme
    c2 = np.array([2825.0 / 27648.0, 0.0, 18575.0 / 48384.0,
                   13525.0 / 55296.0, 277.0 / 14336.0, .25])
    b = np.empty([5, 5], 'float')
  # the following line is ci - ci* in lab4, \Delta_est equationl
  # this is used to calculate \Delta_est = estError for the embededd
  # Runge Kutta  \sum_^6 (c_i -c_i^*)
  #
    c2 = c1 - c2
  # this sets b values for same tableu 
    b[0, 0] = 0.2
    b[1, 0] = 3.0 / 40.0
    b[1, 1] = 9.0 / 40.0
    b[2, 0] = 0.3
    b[2, 1] = -0.9
    b[2, 2] = 1.2
    b[3, 0] = -11.0 / 54.0
    b[3, 1] = 2.5
    b[3, 2] = -70.0 / 27.0
    b[3, 3] = 35.0 / 27.0
    b[4, 0] = 1631.0 / 55296.0
    b[4, 1] = 175.0 / 512.0
    b[4, 2] = 575.0 / 13824.0
    b[4, 3] = 44275.0 / 110592.0
    b[4, 4] = 253.0 / 4096.0
    return (a, c1, c2, b)

class Integrator:

    def set_yinit(self):
        raise ValueError(
            'set_init needs to be overridden in the derived class')

    def __init__(self, coeffFileName):
        with open(coeffFileName, 'rb') as f:
            config = yaml.load(f)
        self.config = config
        # read in dt tstart tend
        timevars = namedtuple('timevars',config['timevars'].keys())
        self.timevars = timevars(**config['timevars'])
        self.rkckConsts = rkck_init()

    def __str__(self):
        out = 'integrator instance with attributes initvars, timevars,uservars, ' + \
            'adaptvars'
        return out

    def derivs(self, y, t):
        raise ValueError('derivs5 needs to be overrideen in the derived class')
        return None

    def rkckODE5(self, yold, timeStep, deltaT):

        # initialize the Cash-Karp coefficients
        # defined in the tableau in lab 4,

        a, c1, c2, b = self.rkckConsts
        i = self.initvars
        # set up array to hold k values in lab4 
        derivArray = np.empty([6, self.nvars], 'float')
        ynext = np.zeros_like(yold)
        bsum = np.zeros_like(yold)
        estError = np.zeros_like(yold)
        # vector k1 in lab4 equation 3.9
        derivArray[0, :] = self.derivs(yold, timeStep)[:]

        # calculate step
        # c1=c_i in lab 4 notation, but c2=c_i - c^*_i

        y = yold
        for i in np.arange(5):
            bsum = 0.
            for j in np.arange(i + 1):
                bsum = bsum + b[i, j] * derivArray[j, :]
            # vectors k2 through k6 in lab4 
#           pdb.set_trace()
            derivArray[i + 1, :] = self.derivs(y + deltaT * bsum, timeStep + a[i] * deltaT)[:]
            # partial sum of error in lab4 \Delta_est
            #
            #  sum the error term
            #
            estError = estError + c2[i] * derivArray[i, :]
            # print "estError: ",estError
            #
            # 5th order estimate y_{n+1}
            #
            ynext = ynext + c1[i] * derivArray[i, :]
        # final fifth order anser
        y = y + deltaT * (ynext + c1[5] * derivArray[5, :])
        # final 4th order estimate estimate
        estError = deltaT * (estError + c2[5] * derivArray[5, :])
        # print "estError final: ",estError
        timeStep = timeStep + deltaT
#       pdb.set_trace()
        return (y, estError, timeStep)
        
        def timeloop5fixed(self):
            t = self.timevars
            yold = self.yinit
            yError = np.zeros_like(yold)
            yvals = [yold]
            errorList = [yError]
            timeSteps = np.arange(t.tstart, t.tend, t.dt)
            for theTime in timeSteps[:-1]:
                yold, yError, newTime = self.rkckODE5(yold, theTime, t.dt)
                yvals.append(yold)
                errorList.append(yError)
            yvals = np.array(yvals).squeeze()
            errorVals = np.array(errorList).squeeze()
            return (timeSteps, yvals, errorVals)


class Katabatic(Integrator):
    def set_yinit(self):
        uservars = namedtuple('uservars',self.config['uservars'].keys())
        self.uservars = uservars(**self.config['uservars'])
        #read in theta, flux, and u
        initvars = namedtuple('initvars', self.config['initvars'].keys())
        self.initvars = initvars(**self.config['initvars'])
        
        '''probably a more elegant way to do this... just don't know how yet'''
        self.yinit = np.array(
            [self.initvars.Theta01,self.initvars.Theta02,self.initvars.Theta03,
             self.initvars.Theta04,self.initvars.Theta05,self.initvars.Theta06,
             self.initvars.Theta07,self.initvars.Theta08,self.initvars.Theta09,
             self.initvars.Theta10,self.initvars.U01,self.initvars.U02,
             self.initvars.U03,self.initvars.U04,self.initvars.U05,
             self.initvars.U06,self.initvars.U07,self.initvars.U08,
             self.initvars.U09,self.initvars.U10])
        self.nvars = len(self.yinit)
        timevars = namedtuple('timevars',self.config['timevars'].keys())
        self.time = timevars(**self.config['timevars'])
        return None
        
    def __init__(self, coeffFileName):
        super().__init__(coeffFileName)
        self.set_yinit()
    ###Nvars will always be an even number because there is one more flux level
    ###than layers for both theta and u... so dividing by 2 will never be an issue
   

    def derivs(self,y,t):
        user = self.uservars
        time = self.time
        alpha = (2*np.pi/360)*user.alpha
        ambient_wind = user.wind_aloft/np.sin(alpha)
        Varct = int(self.nvars)
        Ri_rad = -(user.g*user.LWO*15.0)/(user.Theta_synoptic*user.rho*user.Cp*(user.synoptic_wind**3))
        '''stability function, Ri_param'''
        #Ri = (1-5*Ri)**2 #Richardson parameter with fixed Ri
        Ri = (1-5*Ri_rad)**2 #Richardson parameter with radiation Ri
        
        f = np.empty([Varct],'float')
        G_temp = (user.LWO/(user.sigma*user.epsilon))**0.25
        zf = np.arange(0.0,user.Top,user.dn) #array with height levels

        '''defining the neutral mixing length used by the Met Office'''
        length_MetOffice = np.empty([Varct/2],'float')
        length_MetOffice[0]=0.0
        for h in np.arange(1,(Varct/2)):
            length_MetOffice[h] = 1/((user.lamb*user.k*zf[h])/(user.lamb+(user.k*zf[h])))
        
        '''defining the momentum eddy diffusivity'''
        K_h = np.empty([Varct/2],'float')
        K_h[0:8] = (length_MetOffice[0:8]**2)*Ri*(1/user.dn)*(y[((Varct/2)+1):(Varct-1)]- \
                    y[(Varct/2):(Varct-2)])
        K_h[9] = (length_MetOffice[9]**2)*Ri*(1/user.dn)*(ambient_wind-y[(Varct-1)])

        
        '''defining the parameterization for turbulent stress'''
        '''There is an issue with the 9th layer, potential temperature is increasing
           thus causing the wind speeds to increase dramatically (and move upslope)
           Probably a diffusion mishap, just can't figure out exactly how to handle 
           this layer, it has some communication issues/perhaps my system has
           balance issues and this layer just happens to be the sink?
        '''
        Flux_U = np.empty([(Varct/2)+1],'float')
        Flux_T = np.empty([(Varct/2)+1],'float')
        Flux_U[0] = user.TransferCoef*np.abs(y[Varct/2])*(1/user.dn)*(y[Varct/2]-0.0)
        Flux_U[1:((Varct/2)-1)] = K_h[1:9]*(1/user.dn)*(y[((Varct/2)+1):(Varct-1)]- \
                    y[(Varct/2):(Varct-2)])
        Flux_U[-1] = 0.0
        Flux_T[0] = (user.LWO/(user.rho*user.Cp))*(y[0]-G_temp)
        Flux_T[1] = user.rho*user.TransferCoef*(y[1]-y[0])
        Flux_T[2:((Varct/2)-1)] = K_h[2:9]*(1/user.dn)*(y[2:((Varct/2)-1)]- \
                    y[1:((Varct/2)-2)])
        Flux_T[-1] = 0.0
        
        '''Creating the ambient potential temperature profile'''
        Theta_L = np.arange(283.15,283.04,-0.01) 
        '''Creating the difference in ambient and flow potential temperature'''
        Theta_diff = np.empty([Varct/2],'float')
        Theta_diff[0:((Varct/2)-1)] = Theta_L[0:((Varct/2)-1)]- y[0:((Varct/2)-1)]
        '''The first section of the derivative is for the potential temperature profile'''
        f[0:((Varct/2)-1)] = (1/user.dn)*(Flux_T[1:-1]-Flux_T[0:-2]) + \
                    (user.gamma*np.sin(alpha)*y[(Varct/2):-1])
        '''The second section of the derivative is for the wind speed profile'''
        f[(Varct/2):-1] = (1/user.dn)*(Flux_U[1:-1]-Flux_U[0:-2]) + \
                    ((user.g)*np.sin(alpha)/Theta_L[0:((Varct/2)-1)])*(Theta_diff[0:((Varct/2)-1)]) - \
                    user.Drag*((y[(Varct/2):(Varct-1)])**2) 
        return f
    
    def timeloop5fixed(self):
        t = self.timevars
        yold = self.yinit
        yError = np.zeros_like(yold)
        yvals = [yold]
        errorList = [yError]
        timeSteps = np.arange(t.tstart, t.tend, t.dt)
        for theTime in timeSteps[:-1]:
            yold, yError, newTime = self.rkckODE5(yold, theTime, t.dt)
            yvals.append(yold)
            errorList.append(yError)
        yvals = np.array(yvals).squeeze()
        errorVals = np.array(errorList).squeeze()
        return (timeSteps, yvals, errorVals) 


theSolver=Katabatic('4LayerKflow.yaml')
TimeVals,yVals,errorVals=theSolver.timeloop5fixed()

Wind = np.empty((10,len(TimeVals)))
Temp = np.empty((10,len(TimeVals)))
for i in range(len(TimeVals)):
    Placer = yVals[i]
    Temp[:,i] = Placer[0:10]
    Wind[:,i] = Placer[10:20]




    
    






        