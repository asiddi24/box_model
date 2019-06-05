#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 16:35:40 2019

@author: asiddi24
"""

'''Four box model'''

import numpy as np
import gsw as gsw

def fourbox(N,Kv,AI,Mek,Aredi,M_s,D0,T0s,T0n,T0l,T0d,S0s,S0n,S0l,S0d,Fws,Fwn,epsilon):
    
    Area = 3.6e14
    Area_low = 2e14
    Area_s = 1e14
    Area_n = 0.6e14
    
    Dhigh = 100
    dt = 365*86400/4
    
    inorth = 0
    isouth = 1
    ilow = 2
    ideep = 3
    
    T = np.zeros(shape=[N,4])
    S = np.zeros(shape=[N,4])
    V = np.zeros(shape=[N,4])
    M_n = np.zeros(N)
    M_upw = np.zeros(N)
    M_eddy = np.zeros(N)
    D_low = np.zeros(N)
    D_low[0] = D0 # why is this happening ?
    
    T[0,:] = np.array([T0n, T0s, T0l, T0d]) # Why are the S , T vectors not aligned
    S[0,:] = np.array([S0s, S0n, S0l, S0d])
    
    sigma0 = np.zeros(shape=[N,4])
    sigma2 = np.zeros(shape=[N,4])
   
    M_LS = np.zeros(N)
    M_LN = np.zeros(N)
    
    S_initlow = np.zeros(N)
    
    V_low = np.zeros(N)
    V_deep = np.zeros(N)
    dV_low = np.zeros(N)
    dV_deep = np.zeros(N)
    
    dS_low = np.zeros(N)
    dS_south = np.zeros(N)
    dS_deep = np.zeros(N)
    dS_north = np.zeros(N)
    dT_low = np.zeros(N)
    dT_south = np.zeros(N)
    dT_deep = np.zeros(N)
    dT_north = np.zeros(N)
    
    for j in range(len(N)):
    # Computing density using TEOS-10
        sigma0[j,:] = gsw.density.sigma0(S[j,:],T[j,:]) # returns potential density anomaly , conservative Temp is in degC and abs S in g/kg
        sigma2[j,:] = gsw.density.sigma2(S[j,:],T[j,:])            
        M_LS[j] = Aredi*2.5e7*D_low[j]/1e6 # What is happening here ?
        M_LN[j] = Aredi*5e6*D_low[j]/1e6
        
        if (sigma0[j,inorth] > sigma0[j,ilow]):
            gprime = 9.8*(sigma0[j,inorth] - sigma0[j,ilow])/sigma0[j,inorth]
            
            M_n[j] = gprime*D_low[j]**2/epsilon
            M_upw[j] = Kv*Area_low/np.min(D_low[j],3700-D_low[j]) # What is happening here ?
            M_eddy[j] = AI*D_low[j]*2.5e7/1e6 # how ? do something with 2.5e7
            
            V_deep[j] = 3700*Area-Area_n*Dhigh - Area_s*Dhigh - Area_low*D_low[j]
            V_low[j] = Area_low*D_low[j]
            
            S_initlow[j] = S[j,ilow]*Area_low*D_low[j] # What is this ?
            
            dV_low[j] = (Mek-M_eddy[j] - M_n[j] + M_upw[j] - Fws - Fwn)*dt # initialize this
            dV_deep[j] = -dV_low[j] # intialize this
            
            D_low[j+1] = D_low[j] + dV_low[j]/Area_low
            
            dS_low[j] = (Mek*S[j,isouth] - M_eddy[j]*S[j,ilow] - M_n[j]*S[j,ilow] +\
                        M_upw[j]*S[j,ideep] + M_LS[j]*(S[j,isouth] - S[j,ilow]) + \
                        M_LN[j]*(S[j,inorth] - S[j,ilow]))*dt
            dS_south[j] = ((M_eddy[j]+M_LS[j])*(S[j,ilow] - S[j,isouth]) + \
                        (Mek+M_s)*(S[j,ideep] - S[j,isouth]) - Fws*S[j,isouth])*dt
            dS_deep[j] = (M_n[j]*S[j,inorth] - (M_upw[j] + Mek + M_s)*S[j,ideep] +\
                           (M_eddy[j]+M_s)*S[j,isouth] + Fws*S[j,isouth] +\
                           Fwn*S[j,inorth])*dt
            dS_north[j] = ((M_n[j]+M_LN[j])*(S[j,ilow] - S[j,inorth]) - \
                            Fwn*S[j,inorth])*dt
            dT_low[j] = (Mek*T[j,isouth] - M_eddy[j]*T[j,ilow] - M_n[j]*T[j,ilow] +\
                        M_upw[j]*T[j,ideep] + M_LS[j]*(T[j,isouth] - T[j,ilow]) + \
                        M_LN[j]*(T[j,inorth] - T[j,ilow]) + Area_low*100*(T0l-T[j,ilow])/365/86400)*dt
            dT_south[j] = ((M_eddy[j]+M_LS)*(T[j,ilow]-T[j,isouth]) + \
                            (Mek+M_s)*(T[j,ideep]-T[j,isouth]) + \
                            Area_s*100*(T0s-T[j,isouth])/365/86400)*dt
            dT_deep[j] = ((M_n[j] + Fwn)*T[j,inorth] - (M_upw[j]+Mek+M_s)*T[j,ideep] +\
                     (M_eddy[j] + M_s + Fws)*T[j,isouth])*dt
            dT_north[j] = ((M_n[j] + M_LN[j])*(T[j,ilow] - T[j,inorth]) +\
                            Area_n*100*(T0n-T[j,inorth])/365/86400)*dt
                    
            S[j+1,inorth]=S[j,inorth]+dS_north[j]/(Dhigh*Area_n)
            S[j+1,isouth]=S[j,isouth]+dS_south[j]/(Dhigh*Area_s)
            S[j+1,ilow]=(S[j,ilow]*V_low[j]+dS_low[j])/(V_low[j]+dV_low[j])
            S[j+1,ideep]=(S[j,ideep]*V_deep[j]+dS_deep[j])/(V_deep[j]+dV_deep[j])
            T[j+1,inorth]=T[j,inorth]+dT_north[j]/(Dhigh*Area_n)
            T[j+1,isouth]=T[j,isouth]+dT_south[j]/(Dhigh*Area_s)
            T[j+1,ilow]=(T[j,ilow]*V_low[j]+dT_low[j])/(V_low[j]+dV_low[j])
            T[j+1,ideep]=(T[j,ideep]*V_deep[j]+dT_deep[j])/(V_deep[j]+dV_deep[j])      
        
    elif (sigma0[j,inorth] <= sigma0[j,ilow]):
        
            
            
            
            