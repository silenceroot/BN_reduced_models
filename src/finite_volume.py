#!/usr/bin/env python3
# -*- coding: utf-8 -*
"""
Created on Mon Apr 23 15:41:05 2018

@author: leixin
"""
import numpy as np
from EOS import p2e,e2p

def minmod(a,b,c):
    if (c>0 and a>0 and b>0):
        return min(a,b,c);
    elif (c<0 and a<0 and b<0):
        return max(a,b,c)
    else:
        return 0.0
#原始变量to守恒量
def prim2cons(rho,p,u,v,z_a,z_b,phi_a,phi_b):
    U_rho_a=rho*phi_a
    U_rho_b=rho*phi_b
    U_rho_c=rho*(1-phi_a-phi_b)
    U_u=rho*u
    U_v=rho*v
    U_e=p2e(p,rho,phi_a,phi_b,z_a,z_b)
    U_E=U_e+0.5*(u**2+v**2)
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
def R_I(rho,Lambda,v0_s):                                                       #Lambda是产物的质量分数\v0_s是v0_b
    if rho*v0_s<1+a:
        return 0
    elif rho*v0_s>=1+a and Lambda<lambda_ig_max:
        return I*np.power(1-Lambda,b)*np.power(rho*v0_s-a-1,x)
def R_G1(p,Lambda):
    if Lambda>lambda_G1_max:
        return 0
    elif Lambda<=lambda_G1_max and Lambda>0:
        return G1*np.power(1-Lambda,c)*np.power(Lambda,d)*np.power(p,y)
def R_G2(p,Lambda):
    if Lambda<lambda_G2_min:
        return 0
    else:
        return G2*np.power(1-Lambda,e)*np.power(Lambda,g)*np.power(p,z)
def R(rho,p,Lambda,v0_s):
    return R_I(rho,Lambda,v0_s)+R_G1(p,Lambda)+R_G2(p,Lambda)

if __name__=='__main__':
    print(R(2,12,0.001,1))