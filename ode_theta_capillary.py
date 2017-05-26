# -*- coding: utf-8 -*-
"""
Created on Tue May 23 19:21:40 2017

@author: jlapucha
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# INPUT PARAMETERS
y0 = [0.5,0] # [initial position at x=0,derivative at x=0]
H = 0.2 # Distance of cantilever in meters
E = 3.95e9 # Young's modulus in Pa
I = 8.15499e-16 # Moment of intertia for cantilever
P = .0001 # Force on cantilever in N
N = 2000 # Number of points in intergration

def cantilever_ode_system(y,x,P,E,I):
    
    dy0_dx = y[1]
    dy1_dx = -P*x/E/I*(1+y[1]**2)**(3/2)
    
    dy_dx = [dy0_dx , dy1_dx]
    
    return dy_dx
    
    
plt.close('all')

x = np.linspace(0,H,N)
output = odeint(cantilever_ode_system,y0,x,args=(P,E,I))

y = output[:,0]
dy_dx = output[:,1]
# Compute theta in degrees
theta = np.arctan(dy_dx)/np.pi*180

# Plotting y vs. x so that it looks like a vertical beam
plt.figure()
plt.plot(y,-x)
plt.axis('equal')

plt.figure()
plt.plot(x,-theta)
#plt.axis('equal')