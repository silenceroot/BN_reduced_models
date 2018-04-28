# -*- coding: utf-8 -*-
"""
Initial data (read)

This is a script file to read in initial data.
"""


'''
DATA IN
'''
import numpy as np
import sys      
import pandas as pd
import math

from Boundary import Boundary_Condition
from EOS import ComputeOfC

from finite_volume import cons2prim, prim2cons, GRP_solver_HLLC, R



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
    os.system("python value_start.py")


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
w_a = conf.get_value(601,'value')
w_b = conf.get_value(602,'value')
w_c = conf.get_value(603,'value')                  #gamma_c=w_c+1 for idea gas
Cv_a = conf.get_value(701,'value')
Cv_b = conf.get_value(702,'value')
Cv_c = conf.get_value(703,'value')
A_a = conf.get_value(201,'value')
A_b = conf.get_value(202,'value')
A_c = conf.get_value(203,'value')
B_a = conf.get_value(301,'value')
B_b = conf.get_value(302,'value')
B_c = conf.get_value(303,'value')
R1_a = conf.get_value(401,'value')
R1_b = conf.get_value(402,'value')
R1_c = conf.get_value(403,'value')
R2_a = conf.get_value(501,'value')
R2_b = conf.get_value(502,'value')
R2_c = conf.get_value(503,'value')
Q_a = conf.get_value(801,'value')
Q_b = conf.get_value(802,'value')
Q_c = conf.get_value(803,'value')
v0_a = conf.get_value(901,'value')
v0_b = conf.get_value(902,'value')
v0_c = conf.get_value(903,'value')
I = conf.get_value(101,'value')
b = conf.get_value(102,'value')
a = conf.get_value(103,'value')
x = conf.get_value(104,'value')
lambda_ig_max = conf.get_value(105,'value')
G1 = conf.get_value(111,'value')
c = conf.get_value(112,'value')
d = conf.get_value(113,'value')
y = conf.get_value(114,'value')
lambda_G1_max = conf.get_value(115,'value')
G2 = conf.get_value(121,'value')
e = conf.get_value(122,'value')
g = conf.get_value(123,'value')
z = conf.get_value(124,'value')
lambda_G2_min = conf.get_value(125,'value')



rho  = np.array(pd.read_csv('./RHO.csv',header=None,dtype=np.float64))
u    = np.array(pd.read_csv('./U.csv',header=None,dtype=np.float64))
v    = np.array(pd.read_csv('./V.csv',header=None,dtype=np.float64))
p    = np.array(pd.read_csv('./P.csv',header=None,dtype=np.float64))
phi_a= np.array(pd.read_csv('./PHI_a.csv',header=None,dtype=np.float64))
z_a  = np.array(pd.read_csv('./Z_a.csv',header=None,dtype=np.float64))
phi_b= np.array(pd.read_csv('./PHI_b.csv',header=None,dtype=np.float64))
z_b  = np.array(pd.read_csv('./Z_b.csv',header=None,dtype=np.float64))


name=locals()
for flu_var in ['rho','u','v','p','phi_a','z_a','phi_b','z_b']:
    if (name[flu_var].shape[1]!=num_x):
        print("Columns number of ",flu_var," is incorrect!")
        sys.exit(1)
    elif (name[flu_var].shape[0]!=num_y):
        print("Rows number of ",flu_var," is incorrect!")
        sys.exit(1)

#边界条件
rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,\
phi_a_with_boundary,z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary=\
Boundary_Condition(rho,p,u,v,phi_a,z_a,phi_b,z_b)


os.chdir(src_path)


'''
FINITE VOLUME SCHEME
'''

t=0.0
#Godunov-type Method
while t<t_all and not (math.isinf(t) or math.isnan(t)):

    #CFL condition
    c=np.zeros((num_y+2,num_x+2))                                                       #声速数组--加边界
    for i in range(num_y+2):
        for j in range(num_x+2):
            c[i,j]=ComputeOfC(p_with_boundary[i,j],rho_with_boundary,phi_a_with_boundary,
                                phi_b_with_boundary,z_a_with_boundary,z_b_with_boundary)
    Smax=max(np.max(abs(u_with_boundary)+c),np.max(abs(v_with_boundary)+c))

    d_t=CFL*d_x/Smax;
    if t+d_t >= t_all:
        d_t = t_all-t+0.01*eps
    #Riemann Reoblem:compute flux
    F=[]                                                                                #存储计算出来的通量
    G=[]
    U_flux=np.zeros((num_y,num_x+1))                                                    #x方向的速度通量
    V_flux=np.zeros((num_y+1,num_x))                                                    #y方向的速度通量
    z_a_x_flux=np.zeros((num_y,num_x+1))                                                #x方向的z_a通量
    z_a_y_flux=np.zeros((num_y+1,num_x))
    z_b_x_flux=np.zeros((num_y,num_x+1))
    z_b_y_flux=np.zeros((num_y+1,num_x))
    
    #守恒量
    U_rho_a=rho_with_boundary*phi_a_with_boundary
    U_rho_b=rho_with_boundary*phi_b_with_boundary
    U_rho_c=rho_with_boundary*phi_c_with_boundary
    U_u=rho_with_boundary*u_with_boundary
    U_v=rho_with_boundary*v_with_boundary
    U_E=np.zeros((num_y+2,num_x+2))
    for i in range(num_y+2):
        for j in range(num_x+2):
            U_E[i,j]=prim2cons(rho_with_boundary[i,j],p_with_boundary[i,j],u_with_boundary[i,j],
                                v_with_boundary[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],
                                phi_a_with_boundary[i,j],phi_b_with_boundary[i,j])[5]
#求x方向的各种通量
    for i in range(1,num_y+1):
        for j in range(num_x+1):
            F_computer=[]
            F_return,U_return,z_a_return,z_b_return=\
            GRP_solver_HLLC(U_rho_a[i,j],U_rho_b[i,j],U_rho_c[i,j],U_u[i,j],U_v[i,j],U_E[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],c[i,j],
                            U_rho_a[i,j+1],U_rho_b[i,j+1],U_rho_c[i,j+1],U_u[i,j+1],U_v[i,j+1],U_E[i,j+1],z_a_with_boundary[i,j+1],z_b_with_boundary[i,j+1],c[i,j+1])
            F_computer.append(F_return)
            U_flux[i-1,j]=U_return
            z_a_x_flux[i-1,j]=z_a_return
            z_b_x_flux[i-1,j]=z_b_return
        F.append(F_computer)
#求y方向上的通量
    for i in range(num_y+1):
        for j in range(1,num_x+1):
            G_computer=[]
            G_return,V_return,z_a_return,z_b_return=\
            GRP_solver_HLLC(U_rho_a[i,j],U_rho_b[i,j],U_rho_c[i,j],U_v[i,j],U_u[i,j],U_E[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],c[i,j],
                            U_rho_a[i+1,j],U_rho_b[i+1,j],U_rho_c[i+1,j],U_v[i+1,j],U_u[i+1,j],U_E[i+1,j],z_a_with_boundary[i+1,j],z_b_with_boundary[i+1,j],c[i+1,j])
            G_computer.append(G_return)
            V_flux[i,j-1]=V_return
            z_a_y_flux[i,j-1]=z_a_return
            z_b_y_flux[i,j-1]=z_b_return
#计算源项R
    R_computer=np.zeros((num_y,num_x))
    for i in range(num_y):
        for j in range(num_x):
            R_computer[i,j]=R(rho[i,j],p[i,j],phi_a[i,j],v0_b)
    R_sourse=rho*R_computer
    #compute U in next step
    #6个守恒量
    for i in range(num_y):
        for j in range(num_x):
            U_i_j=np.array([U_rho_a[i+1,j+1],U_rho_b[i+1,j+1],U_rho_c[i+1,j+1],U_u[i+1,j+1],U_v[i+1,j+1],U_E[i+1,j+1]])      #n时刻守恒量值
            R_sourse_array=np.array([R_sourse[i,j],0,0,0,0,0])
            U_i_j=U_i_j-d_t/d_x*(F[i][j+1]-F[i][j])-d_t/d_y*(G[i][j]-G[i+1][j])+d_t*R_sourse_array                          #更新下一步
            Return=cons2prim(U_i_j[0],U_i_j[1],U_i_j[2],U_i_j[3],U_i_j[4],U_i_j[5],0,0)
            rho[i,j]=Return[0]
            p[i,j]=Return[1]
            u[i,j]=Return[2]
            v[i,j]=Return[3]
            phi_a[i,j]=Return[6]
            phi_b[i,j]=Return[7]
    #2个非守恒量
    U_flux_L=U_flux[:,:-1]
    U_flux_R=U_flux[:,1:]
    V_flux_L=V_flux[1:,:]
    V_flux_R=V_flux[:-1,:]
    z_a_x_flux_L=z_a_x_flux[:,:-1]
    z_a_x_flux_R=z_a_x_flux[:,1:]
    z_b_x_flux_L=z_b_x_flux[:,:-1]
    z_b_x_flux_R=z_b_x_flux[:,1:]
    z_a=z_a-d_t/d_x*(U_flux_R*z_a_x_flux_R-U_flux_L*z_a_x_flux_L-z_a*(U_flux_R-U_flux_L))\
            -d_t/d_y*(V_flux_R*z_a_y_flux_R-V_flux_L*z_a_y_flux_L-z_a*(V_flux_R-V_flux_L))
    z_b=z_b-d_t/d_x*(U_flux_R*z_b_x_flux_R-U_flux_L*z_b_x_flux_L-z_b*(U_flux_R-U_flux_L))\
            -d_t/d_y*(V_flux_R*z_b_y_flux_R-V_flux_L*z_b_y_flux_L-z_b*(V_flux_R-V_flux_L))
            
    #边界条件
    rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,\
    phi_a_with_boundary,z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary=\
    Boundary_Condition(rho,p,u,v,phi_a,z_a,phi_b,z_b)
    
    t=t+d_t



'''
DATA OUT
'''
data_out_path='../data_out/two-dim/BNR_'+order_str+'_order/'+data_in_folder
if not os.path.exists(data_out_path):
    os.makedirs(data_out_path)
os.chdir(data_out_path)

pd.DataFrame(rho).to_csv('RHO.csv',header=False,index=False)
pd.DataFrame(p).to_csv('P.csv',header=False,index=False)
pd.DataFrame(u).to_csv('U.csv',header=False,index=False)
pd.DataFrame(v).to_csv('V.csv',header=False,index=False)
pd.DataFrame(phi_a).to_csv('PHI_a.csv',header=False,index=False)
pd.DataFrame(phi_b).to_csv('PHI_b.csv',header=False,index=False)
pd.DataFrame(z_a).to_csv('Z_a.csv',header=False,index=False)
pd.DataFrame(z_b).to_csv('Z_b.csv',header=False,index=False)

#plot
if os.path.exists('./value_plot.py'):
    with open('value_plot.py','r') as f:
        exec(f.read())
os.chdir(src_path)