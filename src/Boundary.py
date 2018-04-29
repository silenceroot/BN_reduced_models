# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 09:57:42 2018

@author: silence
"""
import numpy as np
global d_y
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
    u_L_R[N_start:N_start+N,:-1]=-u_L_R[N_start:N_start+N,:-1]
    u_with_boundary=np.vstack((u_L_R[0:1,:],u_L_R,u_L_R[-1:,:]))
    
    v_L_R=np.hstack((v[:,0:1],v,v[:,-1:]))
    v_with_boundary=np.vstack((-v_L_R[0:1,:],v_L_R,-v_L_R[-1:,:]))

    return (rho_with_boundary,p_with_boundary,u_with_boundary,v_with_boundary,phi_a_with_boundary,\
            z_a_with_boundary,phi_b_with_boundary,z_b_with_boundary)