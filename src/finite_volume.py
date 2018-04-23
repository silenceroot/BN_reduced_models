#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 15:41:05 2018

@author: leixin
"""


def minmod(a,b,c):
    if (c>0 and a>0 and b>0):
        return min(a,b,c);
    elif (c<0 and a<0 and b<0):
        return max(a,b,c)
    else:
        return 0.0
    
def prim2cons(rho,p,u,v,z_a,phi_a):
    return (U_rho_a,U_rho_b,U_u,U_e_a,U_e_b)

def cons2prim(U_rho_a,U_rho_b,U_u,U_e_a,U_e_b,phi_a):
    return (rho,p,u,v,z_a)

def GRP_solver_HLLC(U_rho_a,U_rho_b,U_u,U_e_a,U_e_b,phi_a,d_t):
    return (F_L,F_R,U_x_int)