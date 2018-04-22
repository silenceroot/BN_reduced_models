# -*- coding: utf-8 -*-
"""
Initial data (read)

This is a script file to read in initial data.
"""

from pandas import Series, DataFrame
import pandas as pd
import sys
import os

src_path=os.getcwd()
os.chdir('../data_in/two-dim/TEST')

conf = pd.read_table('./config.txt',header = None, comment='#',names=['idx','value'],index_col='idx')
if conf.index.is_unique==False:
    print("Configure number is not unique!")
    sys.exit(1)
t_all = conf.get_value(1,'value')
num_species = int(conf.get_value(2,'value'))
eps = conf.get_value(4,'value')
CFL = conf.get_value(7,'value')
d_x = conf.get_value(10,'value')
d_y = conf.get_value(11,'value')
num_x = int(conf.get_value(13,'value'))
num_y = int(conf.get_value(14,'value'))
gamma_a=conf.get_value(6,'value')
gamma_b=conf.get_value(106,'value')
Cv_a=conf.get_value(110,'value')
Cv_b=conf.get_value(111,'value')
dim = 2

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