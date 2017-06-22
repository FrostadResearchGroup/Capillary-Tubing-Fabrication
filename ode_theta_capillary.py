# -*- coding: utf-8 -*-
"""
Created on Tue May 23 19:21:40 2017

@author: jlapucha
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.optimize import minimize

# INPUT PARAMETERS
length = 10 # Length of cantilever in cm
thetaTarget = 10 # Desiderd angle in degrees (0 = flat tip)
hTarget = 0.9*length # Distance from clamp to surface in cm
E = 3.95e9 # Young's modulus in Pa
I = 8.15499e-16 # Moment of intertia for cantilever

def cantilever_ode_system(y,x,force,E,I):
    """
    System of odes that specifies the curvature of a cantilevered beam.
    Returns the value of the derivatives of the system as a list per the 
    requirement of the odeint function.
    y = list [values of beam position, values of beam slope]
    x = values of positions correponding to y values
    force = float: force on end of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    """
    # Define derivatives
    dy0_dx = y[1]
    dy1_dx = -force*x/E/I*(1+y[1]**2)**(3/2)
    # Put derivatives in a list as required for ode solver
    dy_dx = [dy0_dx , dy1_dx]
    
    return dy_dx
    
def solve_ode_system(h,force,E,I,N=1000):
    """
    Solve for the shape of a cantilevered beam with a force P applied at the 
    end.
    h = distance from beam clamp
    force = float: force on end of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    N = int: number of points in integration
    """
    
    # Define initial position and dependent variable (distance from clamp)
    y0 = [0,0] # [initial position at x=0,derivative at x=0]
    x = np.linspace(0,h,N)
    # Solve for shape of cantilever
    output = odeint(cantilever_ode_system,y0,x,args=(force,E,I))
    y = output[:,0]
    dy_dx = output[:,1]
    # Compute theta in degrees at the end of the beam
    theta = -np.arctan(dy_dx[-1])/np.pi*180
    # Compute arc length
    length = get_arc_length(x,dy_dx)
    
    return x, y, theta, length
    
def get_arc_length(x,dy_dx):
    """
    Compute the arc length of a curve from the value of the derivatives 
    corresponding to that curve at a set of x locations.
    x = ndarray: independent variable
    dy_dy = ndarray: 
    """
    # Define the integrand from arc length formula
    integrand = np.sqrt(1+(dy_dx)**2)
    # Reshape the arrays as required for simps()
    x = x.reshape((1,-1))
    integrand = integrand.reshape((1,-1))
    # Integrate
    arc_length = simps(integrand,x)[0]
    
    return arc_length
    
def objective_function_h(force,h,length,E,I):
    """
    Compute the squared percent difference between the length of a cantilever 
    that has a force applied to the end, which is a distance h from the clamp, 
    and the actual length of the cantilever. To be used for determining the 
    force on the cantilever tip when pressed against a surface the distance h 
    away from the clamp.
    force = float: force on end of beam
    h = distance from beam clamp
    length - float: actual length of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    """
    # Take the absolute value of the force so that negative values are not 
    # encountered.
    force = (force[0])
    # get the shape of the ode
    x,y,theta,lengthTheory = solve_ode_system(h,force,E,I)
    # compute the percentage error squared
    error = (100*(length - lengthTheory)/length)**2
    
    return error
    
def objective_function_theta(params,theta,length,E,I):
    """
    Compute the squared percent difference between the length of a cantilever 
    and the angle at the end that has a force applied to the end, which is a 
    distance h from the clamp, and the actual respective legnth and angle. 
    To be used for determining the force on the cantilever tip and the 
    distance h away from the clamp to the surface than produces a given angle
    at the end.
    params = list [force,h] = [force on end of beam, distance from beam clamp]
    theta = float: actual angle at the end
    length - float: actual length of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    """
    # Parse the parameters
    force = params[0]
    h = params[1]
    # get the shape of the ode
    x,y,thetaTheory,lengthTheory = solve_ode_system(h,force,E,I)
    # compute the percentage error squared
    error = (100*(theta - thetaTheory)/theta)**2 + \
                (100000*(length - lengthTheory)**2/length)**2
                

    return error
    
def get_beam_theta_for_h(h,length,E,I,plot=False):
    """
    Compute the angle in degrees that a cantilever will make with a surface 
    if constrained within a distance from the clamp that is less than the
    length of the beam.
    h = distance from beam clamp
    length - float: actual length of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    plot (optional) = boolean for generating figure with shape of beam
    """
    # Make an initial guess for the force in Newtons.
    forceGuess = (E*I*np.sqrt(length**2 - h**2)/length**3/2.0)
    # Run optimization routine
    output = minimize(objective_function_h,forceGuess,args=(h,length,E,I),
                      method='Nelder-Mead',tol = 1e-9)
    # parse output
    force = output.x
    # Get the shape corresponding to the actual geometry 
    x,y,theta,lengthTheory = solve_ode_system(h,force,E,I)
    # Plot if desired
    if plot:
        plt.figure()
        plt.plot(y,-x)
        x1 = np.linspace(-length,length,20)
        plt.plot(x1,-h*np.ones_like(x1))
        plt.axis('equal')
        
    if (force == forceGuess):
        print '____________________________________'
        print 'WARNING: BAD INITIAL GUESS FOR FORCE!'
        print 'SEE FUNCTION: get_beam_theta_for_h'
        print '------------------------------------'
    
    return theta, force
    
def get_beam_h_for_theta(theta,length,E,I,plot=False):
    """
    Copmute the distance from the clamp that a cantilever must be constrained 
    within that is less than the length of the beam to produce an angle in 
    degrees that it will make with a surface.
    h = distance from beam clamp
    length - float: actual length of beam
    E = float: elastic modulus of beam
    I = float: moment of inertia of beam
    plot (optional) = boolean for generating figure with shape of beam
    """
    # Make an initial guess for the force in Newtons and the distance in 
    # meters. 
    hGuess = length*np.sin((90-0.25*theta)/180.*np.pi)
    forceGuess = E*I*np.sqrt(length**2 - hGuess**2)/length**3/2.0
    guess = [forceGuess,hGuess]
    # Run optimization routine
    output = minimize(objective_function_theta,guess,args=(theta,length,E,I),
                      method='Nelder-Mead',tol = 1e-9)
    # parse output
    force = (output.x[0])
    h = (output.x[1])
    # Get the shape corresponding to the actual geometry 
    x,y,thetaTheory,lengthTheory = solve_ode_system(h,force,E,I)
    # Plot if desired
    if plot:
        plt.figure()
        plt.plot(y,-x)
        x1 = np.linspace(-length,length,20)
        plt.plot(x1,-h*np.ones_like(x1))
        plt.axis('equal')
        
    if (force == forceGuess) or (h == hGuess):
        print '______________________________________________________'
        print 'WARNING: BAD INITIAL GUESS FOR FORCE AND/OR DISTANCE!'
        print 'SEE FUNCTION: get_beam_h_for_theta'
        print '------------------------------------------------------'
    
    return h, force
    
    
if __name__ == "__main__":
    
    plt.close('all')
    
    # Convert to SI units
    length *= 0.01
    hTarget *= 0.01
    
    # Get the angle that will result for a specific distance
    theta, forceTheta = get_beam_theta_for_h(hTarget,length,E,I,plot=True)
    
    # Get the distance setting needed to produce a specific angle
    h, forceH = get_beam_h_for_theta(thetaTarget,length,E,I,plot=True)
    
    # Print the output for user
    print '**** For a distance of %g cm, the angle will be %g with a force of %g Newtons ****' \
     %(hTarget*100,theta,forceTheta)
    print '**** For an angle of %g, the distance should be %g cm with a force of %g Newtons ****' \
     %(thetaTarget,h*100,forceH)