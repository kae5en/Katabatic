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

    def derivs(self, y, t,Fluxes,blackdar):
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
        self.yinit = np.array(
            [self.initvars.Theta1,self.initvars.Theta2,self.initvars.Theta3,
             self.initvars.Theta4,self.initvars.U1,self.initvars.U2,
             self.initvars.U3,self.initvars.U4])
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
        Ri = (1-5*user.Ri)**2
        f = np.empty([self.nvars],'float')
        G_temp = (user.LWO/(user.sigma*user.epsilon))**0.25
        
        blackdar = np.empty([(self.nvars/2)+1],'float')
        blackdar[0] = 0.0
        blackdar[-1] = 0.0
        for i in np.arange(1,self.nvars/2):
            blackdar[i] = user.lamb/(1+(user.lamb/(user.k*i*user.dn)))
        
        Fluxes = np.empty([(self.nvars)+2],'float')
        Fluxes[0] = user.LWO/(user.rho*user.Cp)
        Fluxes[1] = -user.rho*user.Cd*(y[1]-y[0])
        Fluxes[2] = (blackdar[2]**2)*Ri*(1/user.dn)*(y[6]-y[5])
        Fluxes[3] = (blackdar[3]**2)*Ri*(1/user.dn)*(y[7]-y[6])
        Fluxes[4] = 0.0
        Fluxes[5] = 0.0
        Fluxes[6] = (blackdar[1]**2)*Ri*(1/user.dn)*(y[5]-y[4])
        Fluxes[7] = (blackdar[2]**2)*Ri*(1/user.dn)*(y[6]-y[5])
        Fluxes[8] = (blackdar[3]**2)*Ri*(1/user.dn)*(y[7]-y[6])
        Fluxes[9] = 0.0
#        Fluxes[2:(self.nvars/2)] = (blackdar[2:(self.nvars/2)]**2)*Ri*(1/user.dn)*(y[(self.nvars/2)+2:-1] - \
#                y[(self.nvars/2)+1:-2])
#        Fluxes[(self.nvars/2)+1:(self.nvars)] = (blackdar[0:-1]**2)*Ri*(1/user.dn)*(y[(self.nvars/2)+1:(self.nvars-1)] - \
#                y[(self.nvars/2):(self.nvars-2)])
        Fluxes[-1] = 0.0
        f[0] = (1/(2*(user.dn**2)))*(y[0]-G_temp)*(Fluxes[1]-Fluxes[0]) 
        f[1:(self.nvars/2)-2] = (1/(2*(user.dn**2)))*(y[2:(self.nvars/2)-1] - \
                y[0:(self.nvars/2)-3])*(Fluxes[2:(self.nvars/2)-1]-Fluxes[1:(self.nvars/2)-2]) + \
                user.gamma*np.sin(alpha)*y[(self.nvars/2)+1:-2]
        f[(self.nvars/2)] = 0.0 #wind speed at ground
        f[(self.nvars/2)+1:-2] = (1/(2*(user.dn**2)))*(y[(self.nvars/2)+2:-1]-y[(self.nvars/2):-3]) + \
                ((user.g)*np.sin(alpha)/user.Theta_L)*y[1:(self.nvars/2)-2]
        f[-1] = 0 #wind speed at the top level
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









        