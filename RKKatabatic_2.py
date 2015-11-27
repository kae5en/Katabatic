# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 12:05:46 2015

@author: everard
"""

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
import time

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
        adaptvars = namedtuple('adaptvars', config['adaptvars'].keys())
        self.adaptvars = adaptvars(**config['adaptvars'])
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

    def timeloop5Err(self):
        """return errors as well as values
        """
        t = self.timevars
        a = self.adaptvars
        i = self.initvars
        nvars = self.nvars
        oldTime = t.tstart
        olddt = t.dt
        yold = self.yinit
        yerror = np.zeros_like(yold)
        num = 0
        badsteps = 0
        goodsteps = 0
        timeVals = []
        yvals = []
        errorList = []
        while(oldTime < t.tend):
            timeVals.append(oldTime)
            yvals.append(yold)
            errorList.append(yerror)
            if(num > a.maxsteps):
                raise Exception('num > maxsteps')
            # start out with goodstep false and
            # try different sizes for the next step
            # until one meets the error conditions
            # then move onto next step by setting
            # goodstep to true
            goodStep = False
            failSteps = 0
            while(not goodStep):
                # to exit this loop, need to
                # get the estimated error smaller than
                # the desired error set by the relative
                # tolerance
                if(failSteps > a.maxfail):
                    raise Exception('failSteps > a.maxfail')
                #
                # try a timestep, we may need to reverse this
                #
                ynew, yerror, timeStep = self.rkckODE5(yold, oldTime, olddt)
                # print("try a step: : ", ynew)
                #
                # lab 5 section 4.2.3
                # find the desired tolerance by multiplying the relative
                # tolerance (RTOL) times the value of y
                # compare this to the error estimate returnd from rkckODE5
                # atol takes care of the possibility that y~0 at some point
                #
                errtest = 0.
                for i in range(nvars):
                    errtest = errtest + \
                        (yerror[i] / (a.atol + a.rtol * np.abs(ynew[i])))**2.0
                errtest = np.sqrt(errtest / nvars)
                #
                # lab5 equation 4.13, S
                #
                dtchange = a.s * (1.0 / errtest)**0.2
                # print("dtchange, errtest, timeStep: ",
                #       dtchange, errtest, timeStep, ynew, yerror)
                if (errtest > 1.0):
                    # estimated error is too big so
                    # reduce the timestep and retry
                    # dtFailMax ~ 0.5, which guarantees that
                    # the new timestep is reduced by at least a
                    # factor of 2
                    # dtFailMin~0.1, which means that we don't trust
                    # the estimate to reduce the timestep by more
                    # than a factor of 10 in one loop
                    if(dtchange > a.dtfailmax):
                        olddt = a.dtfailmax * olddt
                    elif (dtchange < a.dtfailmin):
                        olddt = a.dtfailmin * olddt
                    else:
                        olddt = dtchange * olddt
                    if (timeStep + olddt == timeStep):
                        raise Exception('step smaller than machine precision')
                    failSteps = failSteps + 1
                    #
                    # undo the timestep since the error wasn't small enough
                    #
                    ynew = yold
                    timeStep = oldTime
                    # go back to top and see if this olddt produices
                    # a better yerrror
                else:
                    # errtest < 1, so we're happy
                    # try to enlarge the timestep by a factor of dtChange > 1
                    # but keep it smaller than dtpassmax
                    # try enlarging the timestep bigger for next time
                    # dtpassmin ~ 0.1 and dtpassmax ~ 5
                    if (abs((1.0 - dtchange)) > a.dtpassmin):
                        if(dtchange > a.dtpassmax):
                            dtnew = a.dtpassmax * olddt
                        else:
                            dtnew = dtchange * olddt
                    else:
                        # don't bother changing the step size if
                        # the change is less than dtpassmin
                        dtnew = olddt
                    goodStep = True
                    #
                    # overwrite the old timestep with the new one
                    #
                    oldTime = timeStep
                    yold = ynew
                    # go back up to top while(timeStep < t.tend)
                    goodsteps = goodsteps + 1
                #
                # this is number of times we decreased the step size without
                #  advancing
                #
                badsteps = badsteps + failSteps
            # special case if we're within one ortwo timesteps of the end
            # otherwise, set dt to the new timestep size
            if(timeStep + dtnew > t.tend):
                olddt = t.tend - timeStep
            elif(timeStep + 2.0 * dtnew > t.tend):
                olddt = (t.tend - timeStep) / 2.0
            else:
                olddt = dtnew
        timeVals = np.array(timeVals).squeeze()
        yvals = np.array(yvals).squeeze()
        errorVals = np.array(errorList).squeeze()
        self.timevals = timeVals
        self.yvals = yvals
        self.errorVals = errorVals
        return (timeVals, yvals, errorVals)
        
            


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
        #ambient_wind = user.wind_aloft/np.sin(alpha)
        Varct = int(self.nvars)
        #Ri_rad = (user.g*user.LWO*15.0)/(user.Theta_synoptic*user.rho*user.Cp*(user.synoptic_wind**3))
        '''Creating the ambient potential temperature profile'''
        Theta_L = np.arange(283.15,283.04,-0.01) 
        Layer_temp = np.empty((Varct/2),'float')
        Layer_temp[:] = Theta_L[0:(Varct/2)]-y[0:(Varct/2)]
        '''stability function, Ri_param'''
        Ri = (1-5*user.Ri)**2 #Richardson parameter with fixed Ri
        #Ri = (1-5*Ri_rad)**2 #Richardson parameter with radiation Ri
        
        f = np.empty((Varct),'float')
        G_temp = (user.LWO/(user.sigma*user.epsilon))**0.25
        zf = np.arange(0.0,(user.Top+1),user.dn) #array with height levels

        '''defining the neutral mixing length used by the Met Office'''
        length_MetOffice = np.empty((11),'float')
        length_MetOffice[0]=0.0
        for h in np.arange(1,(11)):
            length_MetOffice[h] = 1/((user.lamb*user.k*zf[h])/(user.lamb+(user.k*zf[h])))
        
        '''defining the momentum eddy diffusivity'''
        K_h = np.zeros((11),'float')
        K_h[0] = 0.0
        K_h[1] = (length_MetOffice[1]**2)*Ri*(1/user.dn)*(y[11]-y[10])
        K_h[2] = (length_MetOffice[2]**2)*Ri*(1/user.dn)*(y[12]-y[11])
        K_h[3] = (length_MetOffice[3]**2)*Ri*(1/user.dn)*(y[13]-y[12])
        K_h[4] = (length_MetOffice[4]**2)*Ri*(1/user.dn)*(y[14]-y[13])
        K_h[5] = (length_MetOffice[5]**2)*Ri*(1/user.dn)*(y[15]-y[14])
        K_h[6] = (length_MetOffice[6]**2)*Ri*(1/user.dn)*(y[16]-y[15])
        K_h[7] = (length_MetOffice[7]**2)*Ri*(1/user.dn)*(y[17]-y[16])
        K_h[8] = (length_MetOffice[8]**2)*Ri*(1/user.dn)*(y[18]-y[17])
        K_h[9] = (length_MetOffice[9]**2)*Ri*(1/user.dn)*(y[19]-y[18])
        K_h[10] = (length_MetOffice[10]**2)*Ri*(1/user.dn)*(0.0-y[19])
        
        '''defining the parameterization for turbulent stress'''
        '''There is an issue with the 9th layer, potential temperature is increasing
           thus causing the wind speeds to increase dramatically (and move upslope)
           Probably a diffusion mishap, just can't figure out exactly how to handle 
           this layer, it has some communication issues/perhaps my system has
           balance issues and this layer just happens to be the sink?
        '''
        Flux_U = np.zeros(((Varct/2)+1),'float')
        Flux_T = np.empty(((Varct/2)+1),'float')
        Flux_T[0] = user.TransferCoef*((Theta_L[0]+y[0])-G_temp)
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
        Flux_U[0] = user.Drag*(y[10]**3)
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

        f[10] = (1/user.dn)*(Flux_U[1]-Flux_U[0]) - user.g*np.sin(alpha)*y[0]/Theta_L[0] - 4*user.Drag*(y[10]**2)
        f[11] = (1/user.dn)*(Flux_U[2]-Flux_U[1]) - user.g*np.sin(alpha)*y[1]/Theta_L[1] - 3*user.Drag*(y[11]**2)
        f[12] = (1/user.dn)*(Flux_U[3]-Flux_U[2]) - user.g*np.sin(alpha)*y[2]/Theta_L[2] - 2*user.Drag*(y[12]**2)
        f[13] = (1/user.dn)*(Flux_U[4]-Flux_U[3]) - user.g*np.sin(alpha)*y[3]/Theta_L[3] - user.Drag*(y[13]**2)
        f[14] = (1/user.dn)*(Flux_U[5]-Flux_U[4]) - user.g*np.sin(alpha)*y[4]/Theta_L[4] - 0.5*user.Drag*(y[14]**2)
        f[15] = (1/user.dn)*(Flux_U[6]-Flux_U[5]) - user.g*np.sin(alpha)*y[5]/Theta_L[5] - 0.25*user.Drag*(y[15]**2)
        f[16] = (1/user.dn)*(Flux_U[7]-Flux_U[6]) - user.g*np.sin(alpha)*y[6]/Theta_L[6] - (1/8)*user.Drag*(y[16]**2)
        f[17] = (1/user.dn)*(Flux_U[8]-Flux_U[7]) - user.g*np.sin(alpha)*y[7]/Theta_L[7] - (1/16)*user.Drag*(y[17]**2)
        f[18] = (1/user.dn)*(Flux_U[9]-Flux_U[8]) - user.g*np.sin(alpha)*y[8]/Theta_L[8] - (1/32)*user.Drag*(y[19]**2)
        f[19] = (1/user.dn)*(Flux_U[10]-Flux_U[9]) - user.g*np.sin(alpha)*y[9]/Theta_L[9]- (1/64)*user.Drag*(y[19]**2)        
        return f       
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

'''fixed time-step solutions'''
tic_fixed = time.time()
theSolver=Katabatic('deviationFlow.yaml')
TimeVals,yVals,errorVals=theSolver.timeloop5fixed()
elapsed_fixed = time.time() - tic_fixed

'''adapted time-step solutions
tic_adapt = time.time()
theVariableSolver=Katabatic('4LayerKflow.yaml')
TimeVals_a,yVals_a,errorVals_a=theVariableSolver.timeloop5Err()
elapsed_adapt = time.time() - tic_adapt
'''
#separating the wind and temepratures for fixed time step solution    
#Wind = np.empty((10,len(TimeVals)))
#Temp = np.empty((10,len(TimeVals)))
#for i in range(len(TimeVals)):
#    Placer = yVals[i]
#    Temp[:,i] = Placer[0:10]
#    Wind[:,i] = Placer[10:20]    
#Heights = np.arange(0.0,10.0,1.0)




        