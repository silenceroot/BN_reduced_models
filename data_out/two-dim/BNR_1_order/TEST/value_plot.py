#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 12:53:01 2018

@author: leixin
"""

import numpy as np
rho  = np.loadtxt('./RHO.csv',  delimiter=',')
u    = np.loadtxt('./U.csv',    delimiter=',')
v    = np.loadtxt('./V.csv',    delimiter=',')
p    = np.loadtxt('./P.csv',    delimiter=',')
phi_a= np.loadtxt('./PHI_a.csv',delimiter=',')
z_a  = np.loadtxt('./Z_a.csv',  delimiter=',')
phi_b= np.loadtxt('./PHI_b.csv',delimiter=',')
z_b  = np.loadtxt('./Z_b.csv',  delimiter=',')
num_x=len(rho[0])
num_y=len(rho)
if not ('d_x' in dir()):
    d_x=0.01
if not ('d_y' in dir()):
    d_y=0.01
X = np.arange(-0.5*num_x, 0.5*num_x, 1)*d_x
Y = np.arange(-0.5*num_y, 0.5*num_y, 1)*d_y
X, Y = np.meshgrid(X, Y)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
'''
plot density 
'''
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, rho, cmap=plt.cm.winter)
plt.savefig("rho.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, p, cmap=plt.cm.winter)
plt.savefig("p.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, u, cmap=plt.cm.winter)
plt.savefig("u.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, v, cmap=plt.cm.winter)
plt.savefig("v.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, phi_a, cmap=plt.cm.winter)
plt.savefig("phi_a.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, phi_b, cmap=plt.cm.winter)
plt.savefig("phi_b.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, z_a, cmap=plt.cm.winter)
plt.savefig("z_a.png")
fig = plt.figure()  
ax = Axes3D(fig)
ax.plot_surface(X, Y, z_b, cmap=plt.cm.winter)
plt.savefig("z_b.png")
#plt.show()
