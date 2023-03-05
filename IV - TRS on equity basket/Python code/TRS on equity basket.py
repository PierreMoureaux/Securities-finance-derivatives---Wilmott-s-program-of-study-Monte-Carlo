# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:55:04 2023

@author: moureaux pierre
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy
import scipy.linalg

"Market datas"
M = 3
S0 = np.arange(100,100+M,1)
sigma0 = np.full(M,0.1)
r0 = 0.01

"TRS datas"
T = 1
TRSpaymentSchedule = (25,50,75,99)
sTRS = 0.0006

"Model parameters"
kappa = np.full(M,0.1)
theta = np.full(M,0.1)
lamb = np.full(M,0.1)
nu = 0.1
etha = 0.1
alpha = 0.02
deltaT = 0.01
lT = int(T/deltaT)
nbSimul = 1000

"Model correlations"
corrSS = np.matrix([[1, 0.2,0.1], [0.2, 1,0.3],[0.1, 0.3,1]])
L = scipy.linalg.cholesky(corrSS, lower=True)
rhoSsigma = np.full(M,0.1)
rhoSr = np.full(M,0.1)

def pathsGeneration(nbSimul,lT,S0,sigma0,r0,kappa,theta,lamb,nu,etha,alpha,rhoSsigma,rhoSr):
    paths = np.zeros(nbSimul,dtype=dict)
    for i in range(nbSimul):
        x = np.random.normal(0, 1, (M,lT))
        phiS = np.dot(L,x)
        phitild=np.random.normal(0, 1, (M,lT))
        phihat=np.random.normal(0, 1, lT)
        S = np.zeros((M,lT))
        sigma = np.zeros((M,lT))
        r = np.zeros(lT)
        intr = np.zeros(lT)
        S[:,0]=S0
        sigma[:,0]=sigma0
        r[0]=r0
        intr[0]=r[0]
        for k in range(lT-1):
            for j in range(M):
                r[k+1] = (nu-etha*r[k])*deltaT + np.sqrt(alpha*r[k])*(rhoSr[j]*np.sqrt(deltaT)*phiS[j,k]+np.sqrt(1-rhoSr[j]*rhoSr[j])*np.sqrt(deltaT)*phihat[k])+r[k]
                intr[k+1]=intr[k]+r[k+1]
                sigma[j,k+1] = kappa[j]*(theta[j]-sigma[j,k])*deltaT + lamb[j]*np.sqrt(sigma[j,k])*(rhoSsigma[j]*np.sqrt(deltaT)*phiS[j,k]+np.sqrt(1-rhoSsigma[j]*rhoSsigma[j])*np.sqrt(deltaT)*phitild[j,k])+sigma[j,k]
                S[j,k+1] = S[j,k]*(r[k]*deltaT+np.sqrt(sigma[j,k]*deltaT)*phiS[j,k]) + S[j,k]
        paths[i] = {"Asset": S, "Volatility": sigma, "Interest rate": r, "Cum. interest rate": intr}
    return paths

def alphaBetaTRS(TRSpaydate,prevTRSpaydate, paths, nbSimul):
    alphak = 0
    betak = 0
    for i in range(nbSimul):
        sum1 = 0
        sum2 = 0
        for j in range(M):
            sum1 += paths[i]['Asset'][j,TRSpaydate]
            sum2 += paths[i]['Asset'][j,prevTRSpaydate]
        res1 = (sum1-sum2)/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
        res2 = (sum2)/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
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




