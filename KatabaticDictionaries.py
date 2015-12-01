# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:47:29 2015

@author: EOS314
"""

import yaml #I use a yaml file to store my dictionary
out_dict = dict() #this initialized the general dictionary itself
#below, I assign sub-dicionaries to my general one, and within those subdictionaries, I assign variables
out_dict['timevars']=dict(tstart=0.0,dt=0.25,tend=500)
out_dict['initvars']=dict(Theta01=0.0,Theta02=0.0,Theta03=0.0,
                        Theta04=0.0,Theta05=0.0, Theta06=0.0,Theta07=0.0,
                        Theta08=0.0,Theta09=0.0,Theta10=0.0,Theta11=0.0,
                        Theta12=0.0,Theta13=0.0,Theta14=0.0,Theta15=0.0,
                        Theta16=0.0,Theta17=0.0,Theta18=0.,Theta19=0.0,
                        Theta20=0.0,
                        U01=0.0,U02=0.0,U03=0.0,U04=0.0,U05=0.0,U06=0.0,
                        U07=0.0,U08=0.0,U09=0.0,U10=0.0,U11=0.0,U12=0.0,
                        U13=0.0,U14=0.0,U15=0.0,U16=0.0,U17=0.0,U18=0.0,
                        U19=0.0,U20=0.0)
out_dict['uservars']=dict(rho=1.2,g=9.81,gamma=0.002,Rd=287.,TransferCoef=0.001,
                        Cp=1004.,alpha=10.,lamb=1.,k=0.4,LWO=60.,
                        sigma=5.67*(10**-8),epsilon=0.9,cool=5.698*(10**-6),
                        dn=1.0,Ri=0.15,wind_aloft=0.0,Drag=0.01,Top=20.0,
                        synoptic_wind = 0.1, Theta_synoptic = 283.055,
                        topShear = 2.0)
out_dict['adaptvars']=dict(dtpassmin=0.1,dtfailmax=0.5,dtfailmin=0.1,s=0.9,
                            rtol=1.0e-05,atol=1.0e-05,maxsteps=2000.0,
                            maxfail=60.0,dtpassmax=5.0)
#Below, I'm writing a .yaml file called 4LayerKflow... the 'w' means write.
#I then dump my dictionary into the yaml file, and I use this yaml file in my
#katabatic flow code
with open('deviationFlow.yaml','w') as f:
    yaml.dump(out_dict,f)


    