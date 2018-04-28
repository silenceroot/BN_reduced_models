# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 19:52:16 2018

@author: silence
"""
import numpy as np
from scipy.linalg import solve
eps=1e-8
#Mie-Gruneisen EOS 中的F函数
def F(v,A,B,R1,R2,w):
    if A==0 and B==0:
        return 0
    result=A*(v/w-1/R1)*np.exp(-R1*v)+B*(v/w-1/R2)*np.exp(-R2*v)
    return result

#Mie-Gruneisen EOS 中的Z函数
def Z(v,A,B,R1,R2,w):
    if A==0 and B==0:
        return 0
    result=A*(v/w)*np.exp(-R1*v)+B*(v/w)*np.exp(-R2*v)
    return result

#压力to内能
def p2e(p,rho,phi_a,phi_b,z_a,z_b):
    if z_a<eps or phi_a<eps:
        e_a=0
    else:
        v=z_a/(rho*phi_a)
        e_a=p*v/w_a-F(v,A_a,B_a,R1_a,R2_a,w_a)+F(v0_a,A_a,B_a,R1_a,R2_a,w_a)+Q_a
    if z_b<eps or phi_b<eps:
        e_b=0
    else:
        v=z_b/(rho*phi_b)
        e_b=p*v/w_b-F(v,A_b,B_b,R1_b,R2_b,w_b)+F(v0_b,A_b,B_b,R1_b,R2_b,w_b)+Q_b
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_c<eps or phi_c<eps:
        e_c=0
    else:
        v=z_c/(rho*phi_c)
        e_c=p*v/w_c-F(v,A_c,B_c,R1_c,R2_c,w_c)+F(v0_b,A_c,B_c,R1_c,R2_c,w_c)+Q_c
    return phi_a*e_a+phi_b*e_b+phi_c*e_c

#内能to压力
def e2p(e,rho,phi_a,phi_b,z_a,z_b):
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_a<eps or phi_a<eps:
        if z_b<eps or phi_b<eps:
            v_c=z_c/(rho*phi_c)
            p=w_c/v_c*(e+Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
        elif z_c<eps or phi_c<eps:
            v_b=z_b/(rho*phi_b)
            p=w_b/v_b*(e+Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
        else:
            v_c=z_c/(rho*phi_c)
            v_b=z_b/(rho*phi_b)
            B=w_b/v_b*(Z(v_b,A_b,B_b,R1_b,R2_b,w_b)-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))
            C=w_c/v_c*(Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
            p=B+(e-phi_c*(B-C)*v_c/w_c)/(phi_c*v_c/w_c+phi_b*v_b/w_b)
        return p
    elif z_b<eps or phi_b<eps:
        if z_c<eps or phi_c<eps:
            v_a=z_a/(rho*phi_a)
            p=w_a/v_a*(e+Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
        else:
            v_c=z_c/(rho*phi_c)
            v_a=z_a/(rho*phi_a)
            A=w_a/v_a*(Z(v_a,A_a,B_a,R1_a,R2_a,w_a)-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))
            C=w_c/v_c*(Z(v_c,A_c,B_c,R1_c,R2_c,w_c)-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))
            p=A+(e-phi_c*(A-C)*v_c/w_c)/(phi_c*v_c/w_c+phi_a*v_a/w_a)
        return p
    elif z_c<eps or phi_c<eps:
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
def ZZ(v,A,B,R1,R2):
    if A==0 and B==0:
        return 0
    result=A*R1*(v**2)*np.exp(-R1*v)+B*R2*(v**2)*np.exp(-R2*v)
    return result
#混合声速的计算
def ComputeOfC(p,rho,phi_a,phi_b,z_a,z_b):
    z_c=1-z_a-z_b
    phi_c=1-phi_a-phi_b
    if z_a<eps or phi_a<eps:
        c_a_square=0                                                                   #声速的平方
    else:
        v=z_a/(rho*phi_a)
        e_a=p*v/w_a-F(v,A_a,B_a,R1_a,R2_a,w_a)+F(v0_a,A_a,B_a,R1_a,R2_a,w_a)+Q_a
        c_a_square=w_a*(p*v+e_a-Z(v0_a,A_a,B_a,R1_a,R2_a,w_a))+ZZ(v,A_a,B_a,R1_a,R2_a)
    if z_b<eps or phi_b<eps:
        c_b_square=0
    else:
        v=z_b/(rho*phi_b)
        e_b=p*v/w_b-F(v,A_b,B_b,R1_b,R2_b,w_b)+F(v0_b,A_b,B_b,R1_b,R2_b,w_b)+Q_b
        c_b_square=w_b*(p*v+e_b-Z(v0_b,A_b,B_b,R1_b,R2_b,w_b))+ZZ(v,A_b,B_b,R1_b,R2_b)
    if z_c<eps or phi_c<eps:
        c_c_square=0
    else:
        v=z_c/(rho*phi_c)
        e_c=p*v/w_c-F(v,A_c,B_c,R1_c,R2_c,w_c)+F(v0_b,A_c,B_c,R1_c,R2_c,w_c)+Q_c
        c_c_square=w_c*(p*v+e_c-Z(v0_c,A_c,B_c,R1_c,R2_c,w_c))+ZZ(v,A_c,B_c,R1_c,R2_c)
    c_square=phi_a*c_a_square+phi_b*c_b_square+phi_c*c_c_square                    #混合声速的平方
    return np.sqrt(c_square)                                                       #返回混合声速
if __name__=='__main__':
    print(ComputeOfC(12,10,0.1,0.2,0.3,0.2))