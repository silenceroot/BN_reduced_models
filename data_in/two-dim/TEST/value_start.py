#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 10:45:31 2018

@author: leixin
"""

num_x=640
num_y=3
center=num_x//5
rho_L=2.7647
rho_R=1
p_L=4.4468
p_R=1
u_L=1.4833
u_R=0
v_L=0
v_R=0
phi_a_L=0
phi_a_R=1
z_a_L=0
z_a_R=1
import numpy as np
rho  =np.hstack((rho_L  *np.ones((num_y,center)),rho_R  *np.ones((num_y,num_x-center))))
p    =np.hstack((p_L    *np.ones((num_y,center)),p_R    *np.ones((num_y,num_x-center))))
u    =np.hstack((u_L    *np.ones((num_y,center)),u_R    *np.ones((num_y,num_x-center))))
v    =np.hstack((v_L    *np.ones((num_y,center)),v_R    *np.ones((num_y,num_x-center))))
phi_a=np.hstack((phi_a_L*np.ones((num_y,center)),phi_a_R*np.ones((num_y,num_x-center))))
z_a  =np.hstack((z_a_L  *np.ones((num_y,center)),z_a_R  *np.ones((num_y,num_x-center))))
np.savetxt("RHO.csv",    rho,delimiter='\t',fmt='%.8g')
np.savetxt("P.csv",        p,delimiter='\t',fmt='%.8g')
np.savetxt("U.csv",        u,delimiter='\t',fmt='%.8g')
np.savetxt("V.csv",        v,delimiter='\t',fmt='%.8g')
np.savetxt("PHI_a.csv",phi_a,delimiter='\t',fmt='%.8g')
np.savetxt("Z_a.csv",    z_a,delimiter='\t',fmt='%.8g')