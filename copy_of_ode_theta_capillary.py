# -*- coding: utf-8 -*-
"""
Created on Fri May 26 17:36:06 2017

@author: jlapucha
"""

#copy of ode_capillary_theta.py to modify I 
#includes arc length calculations

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import quad

# INPUT PARAMETERS
y0 = [0.5,0] # [initial position at x=0,derivative at x=0]
H = 0.2 # Distance of cantilever in meters
E = 3.95e9 # Young's modulus in Pa
R1 = 359e-6 # OD of capillary tubing
R2 = 51.8e-6 #ID of capillary tubing
P = .0001 #Force on Cantilever in N
N = 2000 # Number of points in intergration

def moment(I, R1, R2):
    I= np.pi/4(R2**4-R1**4)
    return I

def cantilever_ode_system(y,x,P,E,I):
    
    dy0_dx = y[1]
    dy1_dx = -P*x/E/I*(1+y[1]**2)**(3/2)
    
    dy_dx = [dy0_dx , dy1_dx]
    
    return dy_dx
    
def integrand(dy_dx):
    
     return np.sqrt(1+(dy_dx)**2)
     
     

     
plt.close('all')

x = np.linspace(0,H,N)
output = odeint(cantilever_ode_system,y0,x,args=(P,E,I))


y = output[:,0]
dy_dx = output[:,1]


#compute arc length via integration
output_length= quad(integrand,0,0.2)


# Compute theta in degrees
theta = np.arctan(dy_dx)/np.pi*180


# Plotting y vs. x so that it looks like a vertical beam
plt.figure()
plt.plot(y,-x)
plt.axis('equal')

plt.figure()
plt.plot(x,-theta)
#plt.axis('equal')




    


