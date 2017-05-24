# -*- coding: utf-8 -*-
"""
Created on Tue May 23 19:21:40 2017

@author: jlapucha
"""

import numpy as np
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def myfun(X,x,P):
    
    E=16666.7
    I=8.15499e-16
    
    theta=X[0]
    x2=X[1]
    
    f=[np.tan(theta),(-(P*x)/(E*I))*(1+x2**2)**(3/2)]
    return f

x0=[0.001,0.001]
x=np.linspace(0,0.194,20)
P=.0003
output=odeint(myfun,x0,x,(P,))

plt.plot(output[:,0],x)