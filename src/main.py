# -*- coding: utf-8 -*-
"""
Initial data (read)

This is a script file to read in initial data.
"""


'''
DATA IN
'''
dim = 2
while True:
    order_str=input("Enter order of the scheme (positive integer < 2):")
    try:
        order=eval(order_str)
        if (type(order)==int and order>0 and order<3):
            break
    except:
            pass
import os
src_path=os.getcwd()
while True:
    data_in_folder=input('Enter the input folder name:')
    try:
        data_in_path='../data_in/two-dim/'+data_in_folder
        if os.path.exists(data_in_path):
            break
    except:
            pass
os.chdir(data_in_path)
#Generate initial data
if os.path.exists('./value_start.py'):
    os.system("python3 value_start.py")

import sys      
import pandas as pd
conf = pd.read_table('./config.txt',header = None, comment='#',names=['idx','value'],index_col='idx')
if conf.index.is_unique==False:
    print("Configure number is not unique!")
    sys.exit(1)
global t_all,eps,d_x,d_y
t_all = conf.get_value(1,'value')
num_species = int(conf.get_value(2,'value'))
eps = conf.get_value(4,'value')
CFL = conf.get_value(7,'value')
d_x = conf.get_value(10,'value')
d_y = conf.get_value(11,'value')
num_x = int(conf.get_value(13,'value'))
num_y = int(conf.get_value(14,'value'))
gamma_a = conf.get_value(6,'value')
gamma_b = conf.get_value(106,'value')
Cv_a = conf.get_value(110,'value')
Cv_b = conf.get_value(111,'value')

rho  = (pd.read_table('./RHO.csv',  header = None).dropna(axis=1)).values
u    = (pd.read_table('./U.csv',    header = None).dropna(axis=1)).values
v    = (pd.read_table('./V.csv',    header = None).dropna(axis=1)).values
p    = (pd.read_table('./P.csv',    header = None).dropna(axis=1)).values
phi_a= (pd.read_table('./PHI_a.csv',header = None).dropna(axis=1)).values
z_a  = (pd.read_table('./Z_a.csv',  header = None).dropna(axis=1)).values
name=locals()
for flu_var in ['rho','u','v','p','phi_a','z_a']:
    if (len(name[flu_var][0])!=num_x):
        print("Columns number of ",flu_var," is incorrect!")
        sys.exit(1)
    elif (len(name[flu_var])!=num_y):
        print("Rows number of ",flu_var," is incorrect!")
        sys.exit(1) 
os.chdir(src_path)


'''
FINITE VOLUME SCHEME
'''
t=0.0
phi_b=1-phi_a
z_b=1-z_a
U_rho_a=phi_a*rho
U_rho_b=phi_b*rho
U_u=rho*u
U_e_a=phi_a*rho*u**2+z_a*p/(gamma_a-1)
U_e_b=phi_b*rho*u**2+z_b*p/(gamma_b-1)
import numpy as np
F_L=np.zeros((num_y,num_x+1))
F_R=np.zeros((num_y,num_x+1))
G_L=np.zeros((num_y+1,num_x))
G_R=np.zeros((num_y+1,num_x))
U_x_int=np.zeros((num_y,num_x+1))
U_y_int=np.zeros((num_y+1,num_x))
dx_U=np.zeros((num_y,num_x))
dy_u=np.zeros((num_y,num_x))
#Godunov-type Method
import math
from finite_volume import minmod, cons2prim, prim2cons, GRP_solver_HLLC
U_rho_a,U_rho_b,U_u,U_e_a,U_e_b=cons2prim(rho,p,u,v,z_a,phi_a)
"""
while t<t_all and not (math.isinf(t) or math.isnan(t)):
    #reconstruction (minmod limiter)
    for i in range(num_x):
        dx_U=minmod()/d_x
    for i in range(num_y):
        dy_U=minmod()/d_y
    #CFL condition
    for i in range(num_x*num_y):
        a=math.sqrt(gamma_a*p/rho)
    Smax=max(abs(u)+a)
    d_t=CFL*d_x/Smax;
    if t+d_t >= t_all:
        d_t = t_all-t+0.01*eps
    #Riemann Reoblem:compute flux
    for i in range(num_x+1):
        F_L,F_R,U_x_int=GRP_solver_HLLC(U_rho_a,U_rho_b,U_u,U_e_a,U_e_b,phi_a,d_t)
    for i in range(num_y+1):
        G_L,G_R,U_y_int=GRP_solver_HLLC(U_rho_a,U_rho_b,U_u,U_e_a,U_e_b,phi_a,d_t)       
    #compute U in next step
    for i in range(num_x*num_y):
        U_rho_a=U_rho_a+d_t/d_x*(F_R-F_L)+d_t/d_y*(G_R-G_L)
    rho,p,u,v,z_a=prim2cons(U_rho_a,U_rho_b,U_u,U_e_a,U_e_b,phi_a)
    t=t+d_t
"""

'''
DATA OUT
'''
data_out_path='../data_out/two-dim/BNR_'+order_str+'_order/'+data_in_folder
if not os.path.exists(data_out_path):
    os.makedirs(data_out_path)
os.chdir(data_out_path)
np.savetxt("RHO.csv",    rho,delimiter='\t',fmt='%.8g')
np.savetxt("P.csv",        p,delimiter='\t',fmt='%.8g')
np.savetxt("U.csv",        u,delimiter='\t',fmt='%.8g')
np.savetxt("V.csv",        v,delimiter='\t',fmt='%.8g')
np.savetxt("PHI_a.csv",phi_a,delimiter='\t',fmt='%.8g')
np.savetxt("Z_a.csv",    z_a,delimiter='\t',fmt='%.8g')
#plot
if os.path.exists('./value_plot.py'):
    with open('value_plot.py','r') as f:
        exec(f.read())
os.chdir(src_path)