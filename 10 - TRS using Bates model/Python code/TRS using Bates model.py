# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 11:13:03 2023

@author: pierre moureaux
"""

import numpy as np
from matplotlib import pyplot as plt

"Market datas"
S0 = 100
sigma0 = 0.1
r = 0.05

"TRS datas"
T = 1
TRSpaymentSchedule = (25,50,75,99)
sTRS = 0.01

"Model datas"
kappa = 0.1
theta = 0.1
gamma = 0.2
rhoSsigma = 0.1
lamb = 20
muJ = 0.1
sigmaJ = 0.1
deltaT = 0.01
lT = int(T/deltaT)
nbSimul = 1000

def pathsGeneration(nbSimul,lT,S0,sigma0,r,kappa,theta,gamma,lamb,muJ,sigmaJ,rhoSsigma):
    paths = np.zeros(nbSimul,dtype=dict)
    for i in range(nbSimul):
        phiS=np.random.normal(0, 1, lT)
        phitild=np.random.normal(0, 1, lT)
        phiJ=np.random.normal(muJ, sigmaJ, lT)
        P = np.random.poisson(deltaT*lamb, lT)
        S = np.zeros(lT)
        sigma = np.zeros(lT)
        intr = np.zeros(lT)
        S[0]=S0
        sigma[0]=sigma0
        intr[0]=r
        meanJRN = np.exp(muJ+0.5*(sigmaJ**2))-1
        if (2*kappa*theta <= gamma**2):
            kappa = (gamma**2)/(2*theta)+0.1
        for k in range(lT-1):
            sigma[k+1] = kappa*(theta-sigma[k])*deltaT + gamma*np.sqrt(sigma[k])*(rhoSsigma*np.sqrt(deltaT)*phiS[k]+np.sqrt(1-rhoSsigma*rhoSsigma)*np.sqrt(deltaT)*phitild[k])+sigma[k]
            S[k+1] = S[k]*((r-lamb*meanJRN)*deltaT+np.sqrt(sigma[k]*deltaT)*phiS[k]+(np.exp(phiJ[k])-1)*P[k]) + S[k]
            intr[k+1]=intr[k]+r
        paths[i] = {"Asset": S, "Volatility": sigma, "Cum. interest rate": intr}
    plt.plot(paths[1]['Asset'])
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

paths = pathsGeneration(nbSimul,lT,S0,sigma0,r,kappa,theta,gamma,lamb,muJ,sigmaJ,rhoSsigma)

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