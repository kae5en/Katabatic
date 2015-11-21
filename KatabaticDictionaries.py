# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:47:29 2015

@author: EOS314
"""

import yaml
out_dict = dict()
out_dict['timevars']=dict(tstart=0.0,dt=0.125,tend=50)
out_dict['initvars']=dict(Theta1=283.199,Theta2=283.189,Theta3=283.179,
                        Theta4=283.169,U1=0.0,U2=0.0,U3=0.0,U4=0.0)
out_dict['uservars']=dict(rho=1.2,g=9.8,gamma=-0.01,Rd=287.,Cd=0.001,
                        Cp=1004.,alpha=25.,lamb=100.,k=0.4,LWO=-60.,
                        sigma=5.67*(10**-8),epsilon=0.9,cool=5.698*(10**-6),
                        dn=1.0,Ri=0.1,Theta_L=283.169)

with open('4LayerKflow.yaml','w') as f:
    yaml.dump(out_dict,f)



    