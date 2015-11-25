# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:47:29 2015

@author: EOS314
"""

import yaml #I use a yaml file to store my dictionary
out_dict = dict() #this initialized the general dictionary itself
#below, I assign sub-dicionaries to my general one, and within those subdictionaries, I assign variables
out_dict['timevars']=dict(tstart=0.0,dt=0.25,tend=5000)
out_dict['initvars']=dict(Theta01=283.15,Theta02=283.145,Theta03=283.14,
                        Theta04=283.135,Theta05=283.13,Theta06=283.125,
                        Theta07=283.12,Theta08=283.115,Theta09=283.11,
                        Theta10=283.105, U01=0.0,U02=0.0,U03=0.0,U04=0.0,
                        U05=0.0,U06=0.0,U07=0.0,U08=0.0,U09=0.0,U10=0.0)
out_dict['uservars']=dict(rho=1.2,g=9.81,gamma=0.001,Rd=287.,TransferCoef=0.001,
                        Cp=1004.,alpha=5.,lamb=50.,k=0.4,LWO=60.,
                        sigma=5.67*(10**-8),epsilon=0.9,cool=5.698*(10**-6),
                        dn=1.0,Ri=0.1,wind_aloft=0.0,Drag=0.19,Top=10.0,
                        synoptic_wind = 2.5, Theta_synoptic = 283.055)
out_dict['adaptvars']=dict(dtpassmin=0.1,dtfailmax=0.5,dtfailmin=0.1,s=0.9,
                            rtol=1.0e-05,atol=1.0e-05,maxsteps=2000.0,
                            maxfail=60.0,dtpassmax=5.0)
#Below, I'm writing a .yaml file called 4LayerKflow... the 'w' means write.
#I then dump my dictionary into the yaml file, and I use this yaml file in my
#katabatic flow code
with open('4LayerKflow.yaml','w') as f:
    yaml.dump(out_dict,f)


deficit = dict()
deficit['timevars']=dict(tstart=0.0,dt=0.25,tend=1000)
deficit['initvars']=dict(Theta01=0.0,Theta02=0.0,Theta03=0.0,
                        Theta04=0.0,Theta05=0.0,Theta06=0.0,
                        Theta07=0.0,Theta08=0.0,Theta09=0.0,
                        Theta10=0.0, U01=0.0,U02=0.0,U03=0.0,U04=0.0,
                        U05=0.0,U06=0.0,U07=0.0,U08=0.0,U09=0.0,U10=0.0)
deficit['uservars']=dict(rho=1.2,g=9.81,gamma=0.001,Rd=287.,TransferCoef=0.001,
                        Cp=1004.,alpha=5.,lamb=50.,k=0.4,LWO=60.,
                        sigma=5.67*(10**-8),epsilon=0.9,cool=5.698*(10**-6),
                        dn=1.0,Ri=0.1,wind_aloft=0.0,Drag=0.19,Top=10.0)

with open('deficit.yaml','w') as f:
    yaml.dump(deficit,f)
    