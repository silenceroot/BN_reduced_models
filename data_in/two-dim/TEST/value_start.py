#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 10:45:31 2018

@author: leixin
"""
import numpy as np
import pandas as pd
eps=1e-9
num_x=60
num_y=200
d_y=0.1
rho_a=1
rho_b=1
rho_c=1
p_a=4.4468
p_b=1
p_c=1
u_a=0
u_b=0
u_c=0
v_a=0
v_b=0
v_c=0

size_kaikou=1
N=int(size_kaikou/d_y)

rho=np.ones((num_y,num_x))
p=np.ones((num_y,num_x))
u=np.zeros((num_y,num_x))
v=np.zeros((num_y,num_x))
phi_a=np.ones((num_y,num_x))
phi_b=np.ones((num_y,num_x))
z_a=np.ones((num_y,num_x))
z_b=np.ones((num_y,num_x))

rho[-N:,:]=rho_c
rho[:,-N:]=rho_c
rho[:-N,:-N]=rho_b
rho[-N-1,0]=rho_a

p[-N:,:]=p_c
p[:,-N:]=p_c
p[:-N,:-N]=p_b
p[-N-1,0]=p_a

u[-N-1,0]=u_a
v[-N-1,0]=v_a

phi_a=eps*phi_a
z_a=eps*z_a
phi_a[-N-1,0]=1-2*eps
z_a[-N-1,0]=1-2*eps

phi_b=eps*phi_b
phi_b[:-N,:-N]=1-2*eps
phi_b[-N-1,0]=eps
z_b=eps*z_b
z_b[:-N,:-N]=1-2*eps
z_b[-N-1,0]=eps

pd.DataFrame(rho).to_csv('RHO.csv',header=False,index=False)
pd.DataFrame(p).to_csv('P.csv',header=False,index=False)
pd.DataFrame(u).to_csv('U.csv',header=False,index=False)
pd.DataFrame(v).to_csv('V.csv',header=False,index=False)
pd.DataFrame(phi_a).to_csv('PHI_a.csv',header=False,index=False)
pd.DataFrame(phi_b).to_csv('PHI_b.csv',header=False,index=False)
pd.DataFrame(z_a).to_csv('Z_a.csv',header=False,index=False)
pd.DataFrame(z_b).to_csv('Z_b.csv',header=False,index=False)