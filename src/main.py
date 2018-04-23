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
import sys
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
if os.path.exists('./value_start.py'):
    with open('value_start.py','r') as f:
        exec(f.read())
        
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
    if (len(name[flu_var][0])!=num_x or len(name[flu_var])!=num_y):
        print("Rows or columns number of ",flu_var," is incorrect!")
        sys.exit(1)    
os.chdir(src_path)


'''
FINITE VOLUME SCHEME
'''
import numpy as np




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