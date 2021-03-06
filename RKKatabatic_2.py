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
 

    def derivs(self,y,t):
        user = self.uservars
        time = self.time
        alpha = (2*np.pi/360)*user.alpha
#        ambient_wind = user.wind_aloft/np.sin(alpha)
        Varct = int(self.nvars)
<<<<<<< HEAD
#        Ri_rad = (user.g*user.LWO*15.0)/(user.Theta_synoptic*user.rho*user.Cp*(user.synoptic_wind**3))
        ###Creating the ambient potential temperature profile###
        Theta_L = np.arange(283.15,270.04,-0.002) 
=======
        #Ri_rad = (user.g*user.LWO*15.0)/(user.Theta_synoptic*user.rho*user.Cp*(user.synoptic_wind**3))
        '''Creating the ambient potential temperature profile'''
        Theta_L = np.arange(283.15,270.04,-0.09) 
>>>>>>> 94cfc0e6081aca704aa3f8a8690759607d5b190e
        Layer_temp = np.empty((Varct/2),'float')
        Layer_temp[:] = Theta_L[0:(Varct/2)]-y[0:(Varct/2)]
        ###stability function, Ri_param###
        Ri = (1-5*user.Ri)**2 #Richardson parameter with fixed Ri
#        Ri = (1-5*Ri_rad)**2 #Richardson parameter with radiation Ri
        
        f = np.empty((Varct),'float')
        G_temp = (user.LWO/(user.sigma*user.epsilon))**0.25
        #G_temp = 275.15
        zf = np.arange(0.0,(user.Top+1),user.dn) #array with height levels

        ###defining the neutral mixing length used by the Met Office###
        length_MetOffice = np.empty((11),'float')
        length_MetOffice[0]=0.0
        for h in np.arange(1,(11)):
            length_MetOffice[h] = 1/((user.lamb*user.k*zf[h])/(user.lamb+(user.k*zf[h])))
            
        ###defining the vineyard vegetation coverage for drag###
        A = np.empty((10),'float')
        A[0] = 3
        A[1] = 8
        A[2] = 2
        A[3:10] = 1

        ###defining the momentum eddy diffusivity###
        K_h = np.zeros((11),'float')
        K_h[0] = 0.0
        K_h[1:10] = (length_MetOffice[1:10]**2)*Ri*(1/user.dn)*(y[11:20]-y[10:19])
        K_h[10] = (length_MetOffice[10]**2)*Ri*(1/user.dn)*(user.synoptic_wind-y[19])
        
<<<<<<< HEAD
        ###defining the parameterization for turbulent stress###
        Flux_U = np.zeros(((Varct/2)+1),'float')
        Flux_T = np.empty(((Varct/2)+1),'float')
        Flux_T[0] = user.TransferCoef*((Theta_L[0]+y[0])-G_temp)
=======
        '''defining the parameterization for turbulent stress'''
     
        Flux_U = np.zeros(((Varct/2)+1),'float')
        Flux_T = np.empty(((Varct/2)+1),'float')
        Flux_T[0] = user.TransferCoef*((Theta_L[0]-y[0])-G_temp)
>>>>>>> 94cfc0e6081aca704aa3f8a8690759607d5b190e
        Flux_T[1:10] = -K_h[1:10]*(1/user.dn)*(y[1:10]-y[0:9])
        Flux_T[10] = 0.0
        
        Flux_U[0] = user.Drag*(y[10]**2)
        Flux_U[1:10] = -K_h[1:10]*(1/user.dn)*(y[11:20]-y[10:19])
        Flux_U[10] = 0.0

<<<<<<< HEAD
        ###Now derivatves###
        f[0:10] = (1/user.dn)*(Flux_T[1:11]-Flux_T[0:10]) + user.gamma*np.sin(alpha)*y[10:20]
        f[10:20] = (1/user.dn)*(Flux_U[1:11]-Flux_U[0:10]) - \
                user.g*np.sin(alpha)*y[0:10]/Theta_L[0:10] - \
                user.Drag*A[0:10]*(y[10:20]**2)
=======
        '''Now derivatves'''
        f[0:10] = (1/user.dn)*(Flux_T[1:11]-Flux_T[0:10]) + user.gamma*np.sin(alpha)*y[10:20]
        f[10:20] = (1/user.dn)*(Flux_U[1:11]-Flux_U[0:10]) - user.g*np.sin(alpha)*y[0:10]/Theta_L[0:10] - user.Drag*(y[10:20]**2)
 
>>>>>>> 94cfc0e6081aca704aa3f8a8690759607d5b190e
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

<<<<<<< HEAD
###adapted time-step solutions###
=======
'''fixed time-step solutions'''
tic_fixed = time.time()
theSolver=Katabatic('deviationFlow.yaml')
TimeVals,yVals,errorVals=theSolver.timeloop5fixed()
elapsed_fixed = time.time() - tic_fixed

'''adapted time-step solutions'''
>>>>>>> 94cfc0e6081aca704aa3f8a8690759607d5b190e
tic_adapt = time.time()
theVariableSolver=Katabatic('deviationFlow.yaml')
TimeVals_a,yVals_a,errorVals_a=theVariableSolver.timeloop5Err()
elapsed_adapt = time.time() - tic_adapt
<<<<<<< HEAD
=======

#separating the wind and temepratures for fixed time step solution    
#Wind = np.empty((10,len(TimeVals)))
#Temp = np.empty((10,len(TimeVals)))
#for i in range(len(TimeVals)):
#    Placer = yVals[i]
#    Temp[:,i] = Placer[0:10]
#    Wind[:,i] = Placer[10:20]    
#Heights = np.arange(0.0,10.0,1.0)
>>>>>>> 94cfc0e6081aca704aa3f8a8690759607d5b190e

Theta_L = np.arange(283.15,275.04,-0.002) 
Final = yVals_a[-1]
Wind = Final[10:20]
Dev = Final[0:10]
Temp = Theta_L[0:10]+Dev[0:10]
Heights = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]

figWind = plt.figure()
axWind = figWind.add_subplot(111)
axWind.plot(Wind,Heights,'k-*')
axWind.grid(True)
axWind.set_ylabel('Height [m AGL]')
axWind.set_xlabel('Wind Speed [ms-1]')
figWind.set_size_inches(6.0,8.0)
for item in ([axWind.xaxis.label,axWind.yaxis.label]):
    item.set_fontsize(12)
for item in (axWind.get_xticklabels() + axWind.get_yticklabels()):
    item.set_fontsize(10)
plt.show()

figTemp = plt.figure()
axTemp = figTemp.add_subplot(111)
axTemp.plot(Temp,Heights,'k-*')
axTemp.grid(True)
axTemp.set_ylabel('Height [m AGL]')
axTemp.set_xlabel('Potential Temperature [K]')
figTemp.set_size_inches(6.0,8.0)
for item in ([axTemp.xaxis.label,axTemp.yaxis.label]):
    item.set_fontsize(12)
for item in (axTemp.get_xticklabels() + axTemp.get_yticklabels()):
    item.set_fontsize(10)
plt.show()

###fixed time-step solutions###
#tic_fixed = time.time()
#theSolver=Katabatic('deviationFlow.yaml')
#TimeVals,yVals,errorVals=theSolver.timeloop5fixed()
#elapsed_fixed = time.time() - tic_fixed
