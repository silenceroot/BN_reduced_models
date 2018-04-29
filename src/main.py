#!/usr/bin/env python3
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
from scipy.linalg import solve

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
ee = conf.get_value(122,'value')
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
        

os.chdir(src_path)
'''
****************************************************************************************************
'''
#Boundary
def Boundary_Condition(rho,p,u,v,phi_a,z_a,phi_b,z_b):                
    size_kaikou=1                                           #开口大小
    N=int(size_kaikou/d_y)                                  #开口占网格数
    N_start=int((num_y*d_y)//2/d_y)                           #开口开始地方的网格编号

    rho_L_R=np.hstack((rho[:,0:1],rho,rho[:,-1:]))
    rho_with_boundary=np.vstack((rho_L_R[0:1,:],rho_L_R,rho_L_R[-1:,:]))
    
    p_L_R=np.hstack((p[:,0:1],p,p[:,-1:]))
    p_with_boundary=np.vstack((p_L_R[0:1,:],p_L_R,p_L_R[-1:,:]))
    
    phi_a_L_R=np.hstack((phi_a[:,0:1],phi_a,phi_a[:,-1:]))
    phi_a_with_boundary=np.vstack((phi_a_L_R[0:1,:],phi_a_L_R,phi_a_L_R[-1:,:]))
    
    phi_b_L_R=np.hstack((phi_b[:,0:1],phi_b,phi_b[:,-1:]))
    phi_b_with_boundary=np.vstack((phi_b_L_R[0:1,:],phi_b_L_R,phi_b_L_R[-1:,:]))
    
    z_a_L_R=np.hstack((z_a[:,0:1],z_a,z_a[:,-1:]))
    z_a_with_boundary=np.vstack((z_a_L_R[0:1,:],z_a_L_R,z_a_L_R[-1:,:]))
    
    z_b_L_R=np.hstack((z_b[:,0:1],z_b,z_b[:,-1:]))
    z_b_with_boundary=np.vstack((z_b_L_R[0:1,:],z_b_L_R,z_b_L_R[-1:,:]))
        
    u_L_R=np.hstack((-u[:,0:1],u,-u[:,-1:]))
    u_L_R[N_start:N_start+N,-1:]=-u_L_R[N_start:N_start+N,-1:]
    u_with_boundary=np.vstack((u_L_R[0:1,:],u_L_R,u_L_R[-1:,:]))
    
    v_L_R=np.hstack((v[:,0:1],v,v[:,-1:]))
    v_with_boundary=np.vstack((-v_L_R[0:1,:],v_L_R,-v_L_R[-1:,:]))

    return (rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,phi_a_with_boundary,\
            z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary)

#EOS
epss=1e-8
#Mie-Gruneisen EOS 中的F函数
def FF(f_v,f_A,f_B,f_R1,f_R2,f_w):
    if f_A==0 and f_B==0:
        return 0
    result=f_A*(f_v/f_w-1/f_R1)*np.exp(-f_R1*f_v)+f_B*(f_v/f_w-1/f_R2)*np.exp(-f_R2*f_v)
    return result

#Mie-Gruneisen EOS 中的Z函数
def Z(z_v,z_A,z_B,z_R1,z_R2,z_w):
    if z_A==0 and z_B==0:
        return 0
    result=z_A*(z_v/z_w)*np.exp(-z_R1*z_v)+z_B*(z_v/z_w)*np.exp(-z_R2*z_v)
    return result

#压力to内能
def p2e(p,rho,phi_a,phi_b,z_a,z_b):
    if z_a<epss or phi_a<epss:
        e_a=0
    else:
        v=z_a/(rho*phi_a)
        e_a=p*v/w_a-FF(v,A_a,B_a,R1_a,R2_a,w_a)+FF(v0_a,A_a,B_a,R1_a,R2_a,w_a)+Q_a
    if z_b<epss or phi_b<epss:
        e_b=0
    else:
        v=z_b/(rho*phi_b)
        e_b=p*v/w_b-FF(v,A_b,B_b,R1_b,R2_b,w_b)+FF(v0_b,A_b,B_b,R1_b,R2_b,w_b)+Q_b
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_c<epss or phi_c<epss:
        e_c=0
    else:
        v=z_c/(rho*phi_c)
        e_c=p*v/w_c-FF(v,A_c,B_c,R1_c,R2_c,w_c)+FF(v0_c,A_c,B_c,R1_c,R2_c,w_c)+Q_c
    return phi_a*e_a+phi_b*e_b+phi_c*e_c

#内能to压力
def e2p(e,rho,phi_a,phi_b,z_a,z_b):
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_a<epss or phi_a<epss:
        if z_b<epss or phi_b<epss:
            v_c=z_c/(rho*phi_c)
            p=w_c/v_c*(e+Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
        elif z_c<epss or phi_c<epss:
            v_b=z_b/(rho*phi_b)
            p=w_b/v_b*(e+Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
        else:
            v_c=z_c/(rho*phi_c)
            v_b=z_b/(rho*phi_b)
            B=w_b/v_b*(Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
            C=w_c/v_c*(Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
            p=B+(e-phi_c*(B-C)*v_c/w_c)/(phi_c*v_c/w_c+phi_b*v_b/w_b)
        return p
    elif z_b<epss or phi_b<epss:
        if z_c<epss or phi_c<epss:
            v_a=z_a/(rho*phi_a)
            p=w_a/v_a*(e+Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
        else:
            v_c=z_c/(rho*phi_c)
            v_a=z_a/(rho*phi_a)
            A=w_a/v_a*(Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
            C=w_c/v_c*(Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
            p=A+(e-phi_c*(A-C)*v_c/w_c)/(phi_c*v_c/w_c+phi_a*v_a/w_a)
        return p
    elif z_c<epss or phi_c<epss:
        v_b=z_b/(rho*phi_b)
        v_a=z_a/(rho*phi_a)
        A=w_a/v_a*(Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
        B=w_b/v_b*(Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
        p=A+(e-phi_b*(A-B)*v_b/w_b)/(phi_b*v_b/w_b+phi_a*v_a/w_a)
        return p
    else:
        v_b=z_b/(rho*phi_b)
        v_a=z_a/(rho*phi_a)
        v_c=z_c/(rho*phi_c)
        A=w_a/v_a*(Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
        B=w_b/v_b*(Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
        C=w_c/v_c*(Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
        aa=np.array([[phi_a,phi_b,phi_c],[w_a/v_a,-w_b/v_b,0],[w_a/v_a,0,-w_c/v_c]])
        bb=np.array([e,B-A,C-A])
        e_a,_,_=solve(aa,bb)
        p=w_a/v_a*e_a+A
        return p
#方便计算声速定义的函数
def ZZ(zz_v,zz_A,zz_B,zz_R1,zz_R2):
    if zz_A==0 and zz_B==0:
        return 0
    result=zz_A*zz_R1*(zz_v**2)*np.exp(-zz_R1*zz_v)+zz_B*zz_R2*(zz_v**2)*np.exp(-zz_R2*zz_v)
    return result
#混合声速的计算
def ComputeOfC(p,rho,phi_a,phi_b,z_a,z_b):
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_a<epss or phi_a<epss:
        c_a_square=0                                                                   #声速的平方
    else:
        v=z_a/(rho*phi_a)
        e_a=p*v/w_a-FF(v,A_a,B_a,R1_a,R2_a,w_a)+FF(v0_a,A_a,B_a,R1_a,R2_a,w_a)+Q_a
        c_a_square=w_a*(p*v+e_a-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))+ZZ(v,A_a,B_a,R1_a,R2_a)
    if z_b<epss or phi_b<epss:
        c_b_square=0
    else:
        v=z_b/(rho*phi_b)
        e_b=p*v/w_b-FF(v,A_b,B_b,R1_b,R2_b,w_b)+FF(v0_b,A_b,B_b,R1_b,R2_b,w_b)+Q_b
        c_b_square=w_b*(p*v+e_b-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))+ZZ(v,A_b,B_b,R1_b,R2_b)
    if z_c<epss or phi_c<epss:
        c_c_square=0
    else:
        v=z_c/(rho*phi_c)
        e_c=p*v/w_c-FF(v,A_c,B_c,R1_c,R2_c,w_c)+FF(v0_b,A_c,B_c,R1_c,R2_c,w_c)+Q_c
        c_c_square=w_c*(p*v+e_c-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))+ZZ(v,A_c,B_c,R1_c,R2_c)
    c_square=phi_a*c_a_square+phi_b*c_b_square+phi_c*c_c_square                    #混合声速的平方
    return np.sqrt(c_square)                                                       #返回混合声速
#finite_volume
#原始变量to守恒量
def prim2cons(rho,p,u,v,z_a,z_b,phi_a,phi_b):
    U_rho_a=rho*phi_a
    U_rho_b=rho*phi_b
    U_rho_c=rho*(1-phi_a-phi_b)
    U_u=rho*u
    U_v=rho*v
    U_e=p2e(p,rho,phi_a,phi_b,z_a,z_b)
    U_E=rho*(U_e+0.5*(u**2+v**2))
    U_z_a=z_a
    U_z_b=z_b
    return (U_rho_a,U_rho_b,U_rho_c,U_u,U_v,U_E,U_z_a,U_z_b)

#守恒量变原始变量
def cons2prim(U_rho_a,U_rho_b,U_rho_c,U_u,U_v,U_E,U_z_a,U_z_b):
    rho=U_rho_a+U_rho_b+U_rho_c
    u=U_u/rho
    v=U_v/rho
    z_a=U_z_a
    z_b=U_z_b
    phi_a=U_rho_a/rho
    phi_b=U_rho_b/rho
    U_e=U_E/rho-0.5*(u**2+v**2)
    p=e2p(U_e,rho,phi_a,phi_b,z_a,z_b)
    return (rho,p,u,v,z_a,z_b,phi_a,phi_b)

#HLLC Riemann Solver
"""
输入6个守恒变量、2个非守恒变量和声速（左右状态）
输出6个守恒变量的通量、速度、a和b物质的体积分数（用来更新2个非守恒变量）
"""
def GRP_solver_HLLC(U_rho_a_L,U_rho_b_L,U_rho_c_L,U_u_L,U_v_L,U_E_L,U_z_a_L,U_z_b_L,c_L,\
                    U_rho_a_R,U_rho_b_R,U_rho_c_R,U_u_R,U_v_R,U_E_R,U_z_a_R,U_z_b_R,c_R):
    rho_L,p_L,u_L,v_L,z_a_L,z_b_L,phi_a_L,phi_b_L=cons2prim(U_rho_a_L,U_rho_b_L,U_rho_c_L,U_u_L,U_v_L,U_E_L,U_z_a_L,U_z_b_L)
    rho_R,p_R,u_R,v_R,z_a_R,z_b_R,phi_a_R,phi_b_R=cons2prim(U_rho_a_R,U_rho_b_R,U_rho_c_R,U_u_R,U_v_R,U_E_R,U_z_a_R,U_z_b_R)
    S_L=min(u_L-c_L,u_R-c_R)
    S_R=max(u_L+c_L,u_R+c_R)
    S_star=(p_R-p_L+rho_L*u_L*(S_L-u_L)-rho_R*u_R*(S_R-u_R))/(rho_L*(S_L-u_L)-rho_R*(S_R-u_R))
    U_L=np.array([U_rho_a_L,U_rho_b_L,U_rho_c_L,U_u_L,U_v_L,U_E_L])             #守恒量左状态
    U_R=np.array([U_rho_a_R,U_rho_b_R,U_rho_c_R,U_u_R,U_v_R,U_E_R])             #守恒量右状态
    F_L=u_L*U_L+np.array([0,0,0,p_L,0,u_L*p_L])
    F_R=u_R*U_R+np.array([0,0,0,p_R,0,u_R*p_R])
    
    if S_L>=0:
        return (F_L,u_L,z_a_L,z_b_L)
    elif S_R<=0:
        return (F_R,u_R,z_a_R,z_b_R)
    else:
        rho_compute_L=rho_L*(S_L-u_L)/(S_L-S_star)                              #为了计算方便设置的变量
        rho_compute_R=rho_R*(S_R-u_R)/(S_R-S_star)
        U_star_L=rho_compute_L*np.array([phi_a_L,phi_b_L,1-phi_a_L-phi_b_L,S_star,v_L,U_E_L/rho_L+(S_star-u_L)*(S_star+p_L/(rho_L*(S_L-u_L)))])
        U_star_R=rho_compute_R*np.array([phi_a_R,phi_b_R,1-phi_a_R-phi_b_R,S_star,v_R,U_E_R/rho_R+(S_star-u_R)*(S_star+p_R/(rho_R*(S_R-u_R)))])
        F_star_L=F_L+S_L*(U_star_L-U_L)
        F_star_R=F_R+S_R*(U_star_R-U_R)
        if S_star>=0:
            return (F_star_L,S_star,z_a_L,z_b_L)
        else:
            return (F_star_R,S_star,z_a_R,z_b_R)
#源项的函数R
def R_I(rho1,Lambda,v0_s):                                                       #Lambda是产物的质量分数\v0_s是v0_b
    if rho1*v0_s>=1+a and Lambda<lambda_ig_max:
        return I*(1-Lambda)**b*(rho1*v0_s-a-1)**x
    else:
        return 0
def R_G1(p1,Lambda):
    if Lambda>lambda_G1_max:
        return 0
    elif Lambda<=lambda_G1_max and Lambda>0:
        return G1*(1-Lambda)**c*Lambda**d*p1**y
def R_G2(p1,Lambda):
    if Lambda<lambda_G2_min:
        return 0
    else:
        return G2*(1-Lambda)**ee*Lambda**g*p1**z
def R(rho1,p1,Lambda,v0_s):
    return R_I(rho1,Lambda,v0_s)+R_G1(p1,Lambda)+R_G2(p1,Lambda)
'''
*******************************************************************************************************
'''
#边界条件
rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,\
phi_a_with_boundary,z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary=\
Boundary_Condition(rho,p,u,v,phi_a,z_a,phi_b,z_b)


'''
FINITE VOLUME SCHEME
'''

t=0.0
#Godunov-type Method
while t<t_all and not (math.isinf(t) or math.isnan(t)):
    print('time:{}'.format(t))
    #计算源项R
    R_computer=np.zeros((num_y,num_x))
#    for i in range(num_y):
#        for j in range(num_x):
#            R_computer[i,j]=R(rho[i,j],p[i,j],phi_a[i,j],v0_b)
    R_sourse=rho*R_computer
    
    #CFL condition
    c_sound=np.zeros((num_y+2,num_x+2))                                                       #声速数组--加边界
    for i in range(num_y+2):
        for j in range(num_x+2):
            c_sound[i,j]=ComputeOfC(p_with_boundary[i,j],rho_with_boundary[i,j],phi_a_with_boundary[i,j],
                                phi_b_with_boundary[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j])
    Smax=max(np.max(abs(u_with_boundary)+c_sound),np.max(abs(v_with_boundary)+c_sound))
    
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
    U_rho_c=rho_with_boundary*(1-phi_a_with_boundary-phi_b_with_boundary)
    U_u=rho_with_boundary*u_with_boundary
    U_v=rho_with_boundary*v_with_boundary
    U_E=np.zeros((num_y+2,num_x+2))
    for i in range(num_y+2):
        for j in range(num_x+2):
            U_E_computer=prim2cons(rho_with_boundary[i,j],p_with_boundary[i,j],u_with_boundary[i,j],
                                v_with_boundary[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],
                                phi_a_with_boundary[i,j],phi_b_with_boundary[i,j])
            U_E[i,j]=U_E_computer[5]
#求x方向的各种通量
    for i in range(1,num_y+1):
        F_computer=[]
        for j in range(num_x+1):
            F_return,U_return,z_a_return,z_b_return=\
            GRP_solver_HLLC(U_rho_a[i,j],U_rho_b[i,j],U_rho_c[i,j],U_u[i,j],U_v[i,j],U_E[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],c_sound[i,j],
                            U_rho_a[i,j+1],U_rho_b[i,j+1],U_rho_c[i,j+1],U_u[i,j+1],U_v[i,j+1],U_E[i,j+1],z_a_with_boundary[i,j+1],z_b_with_boundary[i,j+1],c_sound[i,j+1])
            F_computer.append(F_return)
            U_flux[i-1,j]=U_return
            z_a_x_flux[i-1,j]=z_a_return
            z_b_x_flux[i-1,j]=z_b_return
        F.append(F_computer)
#求y方向上的通量
    for i in range(num_y+1):
        G_computer=[]
        for j in range(1,num_x+1):
            G_return,V_return,z_a_return,z_b_return=\
            GRP_solver_HLLC(U_rho_a[i+1,j],U_rho_b[i+1,j],U_rho_c[i+1,j],U_v[i+1,j],-U_u[i+1,j],U_E[i+1,j],z_a_with_boundary[i+1,j],z_b_with_boundary[i+1,j],c_sound[i+1,j],
                            U_rho_a[i,j],U_rho_b[i,j],U_rho_c[i,j],U_v[i,j],-U_u[i,j],U_E[i,j],z_a_with_boundary[i,j],z_b_with_boundary[i,j],c_sound[i,j])
            G_return[3],G_return[4]=G_return[4],G_return[3]
            G_computer.append(G_return)
            V_flux[i,j-1]=V_return
            z_a_y_flux[i,j-1]=z_a_return
            z_b_y_flux[i,j-1]=z_b_return
        G.append(G_computer)
    #compute U in next step
    #6个守恒量
    for i in range(num_y):
        for j in range(num_x):
            U_i_j=np.array([U_rho_a[i+1,j+1],U_rho_b[i+1,j+1],U_rho_c[i+1,j+1],U_u[i+1,j+1],U_v[i+1,j+1],U_E[i+1,j+1]])      #n时刻守恒量值
            R_sourse_array=np.array([R_sourse[i,j],0,0,0,0,0])
            U_i_j=U_i_j-d_t/d_x*(F[i][j+1]-F[i][j])-d_t/d_y*(G[i][j]-G[i+1][j])+d_t*R_sourse_array                          #更新下一步
            Return=cons2prim(U_i_j[0],U_i_j[1],U_i_j[2],U_i_j[3],U_i_j[4],U_i_j[5],z_a[i,j],z_b[i,j])
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
    z_a_y_flux_L=z_a_y_flux[1:,:]
    z_a_y_flux_R=z_a_y_flux[:-1,:]
    z_b_x_flux_L=z_b_x_flux[:,:-1]
    z_b_x_flux_R=z_b_x_flux[:,1:]
    z_b_y_flux_L=z_b_y_flux[1:,:]
    z_b_y_flux_R=z_b_y_flux[:-1,:]
    z_a=z_a-d_t/d_x*(U_flux_R*z_a_x_flux_R-U_flux_L*z_a_x_flux_L-z_a*(U_flux_R-U_flux_L))\
            -d_t/d_y*(V_flux_R*z_a_y_flux_R-V_flux_L*z_a_y_flux_L-z_a*(V_flux_R-V_flux_L))
    z_b=z_b-d_t/d_x*(U_flux_R*z_b_x_flux_R-U_flux_L*z_b_x_flux_L-z_b*(U_flux_R-U_flux_L))\
            -d_t/d_y*(V_flux_R*z_b_y_flux_R-V_flux_L*z_b_y_flux_L-z_b*(V_flux_R-V_flux_L))
            
    #边界条件
    rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,\
    phi_a_with_boundary,z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary=\
    Boundary_Condition(rho,p,u,v,phi_a,z_a,phi_b,z_b)
    
    
    
    silence=0
    if np.isnan(rho).any():
        print('rho')
        silence=1
    if np.isnan(p).any():
        print('p')
        silence=1
    if np.isnan(u).any():
        print('u')
        silence=1
    if np.isnan(v).any():
        print('v')
        silence=1
    if np.isnan(phi_a).any():
        print('phi_a')
        silence=1
    if np.isnan(phi_b).any():
        print('phi_b')
        silence=1
    if np.isnan(z_a).any():
        print('z_a')
        silence=1
    if np.isnan(z_b).any():
        print('z_b')
        silence=1
    if silence==1:
        break
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

"""
#plot
if os.path.exists('./value_plot.py'):
    with open('value_plot.py','r') as f:
        exec(f.read())
os.chdir(src_path)
"""