#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 10:27:14 2019

@author: asiddi24
"""

'''Running the box_model'''

'''Initializing variables'''

from box_model import fourbox

N=4000
Kv=1e-5
AI=1000
Mek=25e6
Aredi=1000
M_s=15e6
D0=500
T0s=2
T0n=4
T0l=17
T0d=3
S0s=34
S0n=35
S0l=36
S0d=34.5
Fws=1e6
Fwn=1e6
epsilon=1.2e-4

(M_n, M_upw, M_eddy, D_low, T, S, sigma0) = fourbox(N,Kv,AI,Mek,Aredi,M_s,D0,T0s,T0n,T0l,T0d,S0s,S0n,S0l,S0d,Fws,Fwn,epsilon)


