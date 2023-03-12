# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:55:04 2023

@author: moureaux pierre
"""

import numpy as np
from matplotlib import pyplot as plt

"Market datas"
S0 = 100
sigma0 = 0.1
r0 = 0.05

"TRS datas"
T = 1
TRSpaymentSchedule = (25,50,75,99)
sTRS = 0.01

"Model datas"
kappa = 0.1
theta = 0.1
lamb = 0.1
nu = 0.1
etha = 0.1
alpha = 0.02
rhoSsigma = 0.1
rhoSr = 0.1
deltaT = 0.01
lT = int(T/deltaT)
nbSimul = 1000

def pathsGeneration(nbSimul,lT,S0,sigma0,r0,kappa,theta,lamb,nu,etha,alpha,rhoSsigma,rhoSr):
    paths = np.zeros(nbSimul,dtype=dict)
    for i in range(nbSimul):
        phiS=np.random.normal(0, 1, lT)
        phitild=np.random.normal(0, 1, lT)
        phihat=np.random.normal(0, 1, lT)
        S = np.zeros(lT)
        sigma = np.zeros(lT)
        r = np.zeros(lT)
        intr = np.zeros(lT)
        S[0]=S0
        sigma[0]=sigma0
        r[0]=r0
        intr[0]=r[0]
        for k in range(lT-1):
            r[k+1] = (nu-etha*r[k])*deltaT + np.sqrt(alpha*r[k])*(rhoSr*np.sqrt(deltaT)*phiS[k]+np.sqrt(1-rhoSr*rhoSr)*np.sqrt(deltaT)*phihat[k])+r[k]
            sigma[k+1] = kappa*(theta-sigma[k])*deltaT + lamb*np.sqrt(sigma[k])*(rhoSsigma*np.sqrt(deltaT)*phiS[k]+np.sqrt(1-rhoSsigma*rhoSsigma)*np.sqrt(deltaT)*phitild[k])+sigma[k]
            S[k+1] = S[k]*(r[k]*deltaT+np.sqrt(sigma[k]*deltaT)*phiS[k]) + S[k]
            intr[k+1]=intr[k]+r[k+1]
        paths[i] = {"Asset": S, "Volatility": sigma, "Interest rate": r, "Cum. interest rate": intr}
    return paths

def alphaBetaTRS(TRSpaydate,prevTRSpaydate, paths, nbSimul):
    alphak = 0
    betak = 0
    for i in range(nbSimul):
        res1 = (paths[i]['Asset'][TRSpaydate]-paths[i]['Asset'][prevTRSpaydate])/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
        res2 = (paths[i]['Asset'][prevTRSpaydate])/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
        alphak += res1
        betak += res2
    alphak /= nbSimul
    betak /= nbSimul
    return [alphak,betak]

paths = pathsGeneration(nbSimul,lT,S0,sigma0,r0,kappa,theta,lamb,nu,etha,alpha,rhoSsigma,rhoSr)

def TRS(TRSpaymentSchedule, sTRS, paths,nbSimul):
    res = 0
    alpha = 0
    beta = 0
    for s in TRSpaymentSchedule:
        alpha += alphaBetaTRS(s,TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1],paths,nbSimul)[0]
        beta += alphaBetaTRS(s,TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1],paths,nbSimul)[1]*(s-TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1])
    res = alpha - sTRS*beta
    impliedRate = alpha/beta
    return [res, impliedRate]

print("TRS market value is",TRS(TRSpaymentSchedule, sTRS, paths,nbSimul)[0])
print("TRS implied rate is",TRS(TRSpaymentSchedule, sTRS, paths,nbSimul)[1])




