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
from scipy.integrate import simps

# INPUT PARAMETERS
y0 = [0.5,0] # [initial position at x=0,derivative at x=0]
H = 0.2 # Distance of cantilever in meters
E = 3.95e9 # Young's modulus in Pa
R1 = 359e-6 # OD of capillary tubing
R2 = 51.8e-6 #ID of capillary tubing
P = .0001 #Force on Cantilever in N
N = 2000 # Number of points in intergration
W= np.linspace(0,0.1,20)

def calculate_moment(r1,r2):
    I= (np.pi/4)*(r2**4-r1**4)
    return I

def cantilever_ode_system(y,x,P,E,R1,R2):
    G=calculate_moment(R1,R2)
    dy0_dx = y[1]
    dy1_dx = -P*x/E/G*(1+y[1]**2)**(3/2)
    
    dy_dx = [dy0_dx , dy1_dx]
    
    return dy_dx

#defining parameters of ode with varied load- (have to re-define constants within function?)
def cantilever_ode_system_1(y1,x,L,E,R1,R2):
    G=calculate_moment(R1,R2)
    dy1_dx1= y1[1]
    dy2_dx1= -L*x/E/G*(1+y1[1]**2)**(3/2)
    dy2_dx2= [dy1_dx1, dy2_dx1]
    return dy2_dx2

#need to define function so can call different values of i array    
def integrand(dy_dx):
     return np.sqrt(1+(dy_dx)**2)
     
plt.close('all')

#original ode function prior to looping
x = np.linspace(0,H,N)
output = odeint(cantilever_ode_system,y0,x,args=(P,E,R1,R2))
y = output[:,0]
dy_dx = output[:,1]

outputs = []
for i in range(0,len(W)-1):
    outputs.append(odeint(cantilever_ode_system,y0,x,args=(W[i],E,R1,R2)))

#compute arc length for fixed load via integration
L= integrand(dy_dx)
arc_length= quad(L,0,0.2)

#quad not working for lengths so trying simps
arc_length= simps(L,0,0.2)

#compute arc length using simps function- trying to use different elements of outputs array
#try to do loop on this for computing lengths for all outputs array result?
L2=outputs[:,0]
integrand_2= integrand(L2)
arclength_2= simps(integrand_2, 0, .2)







## Compute theta from first ODE in degrees
#theta = np.arctan(dy_dx)/np.pi*180
#
##calculating ode with varied load
#
#
## Plotting y vs. x so that it looks like a vertical beam
#plt.figure
#plt.plot(y,-x)
#plt.axis('equal')
#
#plt.figure()
#plt.plot(x,-theta)
##plt.axis('equal')




    


