#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 12:53:01 2018

@author: leixin
"""

import numpy as np
rho  = np.loadtxt('./RHO.csv',  delimiter='\t')
u    = np.loadtxt('./U.csv',    delimiter='\t')
v    = np.loadtxt('./V.csv',    delimiter='\t')
p    = np.loadtxt('./P.csv',    delimiter='\t')
phi_a= np.loadtxt('./PHI_a.csv',delimiter='\t')
z_a  = np.loadtxt('./Z_a.csv',  delimiter='\t')
num_x=len(rho[0])
num_y=len(rho)
if not (('d_x' in dir()) and ('d_y' in dir())):
    d_x=0.01;d_y=0.01
X = np.arange(-0.5*num_x, 0.5*num_x, 1)*d_x
Y = np.arange(-0.5*num_y, 0.5*num_y, 1)*d_y
X, Y = np.meshgrid(X, Y)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()  
ax = Axes3D(fig)
'''
plot density 
'''
ax.plot_surface(X, Y, rho, cmap=plt.cm.winter)
plt.savefig("rho.png")
#plt.show()
