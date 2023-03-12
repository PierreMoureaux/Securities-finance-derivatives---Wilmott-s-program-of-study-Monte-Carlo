# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:55:04 2023

@author: moureaux pierre
"""

import numpy as np

"Market datas"
r0 = 0.05
lambda0 = 0.1
gamma0 = 0.05

"TRS datas"
T = 1
TRSpaymentSchedule = (25,50,75,99)
sTRS = 0.02

"Bond datas"
TP = 2
couponSchedule = (13,33,63)
c = 0.03
R = 0.4

"Model datas"
nu = 0.1
etha = 0.2
vola = 0.1
alpha = 0.1
beta = 0.2
sigma = 0.1
omega = 0.1
deltaT = 0.01
lT = int(T/deltaT)
lTP = int(TP/deltaT)
nbSimul = 10
nbSimulP = 10

def BondPrice(nbSimulP,tk,couponSchedule,c,lTP,r0,lambda0,gamma0,nu,etha,vola,alpha,beta,sigma,omega,R):
    Ptk = 0
    nestedPath = np.zeros(nbSimulP)
    for i in range(nbSimulP):
        phir=np.random.normal(0, 1, lTP)
        phil=np.random.normal(0, 1, lTP)
        phig=np.random.normal(0, 1, lTP)
        r = np.zeros(lTP)
        l = np.zeros(lTP)
        g = np.zeros(lTP)
        c1 = 0
        c2 = 0
        c3 = 0
        r[0] = r0
        l[0] = lambda0
        g[0] = gamma0
        for k in range(lTP-1):
            r[k+1] = etha*(nu-r[k])*deltaT + vola*np.sqrt(r[k]*deltaT)*phir[k] + r[k]
            l[k+1] = beta*(alpha-l[k])*deltaT + sigma*np.sqrt(l[k]*deltaT)*phil[k] + l[k]
            g[k+1] = omega*np.sqrt(deltaT)*phig[k] + g[k] 
        c2 = np.exp(-(np.sum(r[tk:lTP-1])+np.sum(l[tk:lTP-1])+np.sum(g[tk:lTP-1])))
        for t in range(tk+1,lTP-1):
            if t in couponSchedule:
                c1 += np.exp(-(np.sum(r[tk:t+1])+np.sum(l[tk:t+1])+np.sum(g[tk:t+1])))
            c3 += l[t]*np.exp(-(np.sum(r[tk:t+1])+np.sum(l[tk:t+1])+np.sum(g[tk:t+1])))
        nestedPath[i] = c*c1+c2+(1-R)*c3
        Ptk += nestedPath[i]
    return Ptk/nbSimulP
    
def pathsGeneration(nbSimul,nbSimulP,lT,couponSchedule,c,lTP,r0,lambda0,gamma0,nu,etha,vola,alpha,beta,sigma,omega,R):
    paths = np.zeros(nbSimul,dtype=dict)
    for i in range(nbSimul):
        phir2=np.random.normal(0, 1, lT)
        P = np.zeros(lT)
        r = np.zeros(lT)
        intr = np.zeros(lT)
        r[0]=r0
        intr[0]=r[0]
        for k in range(lT-1):
            r[k+1] = etha*(nu-r[k])*deltaT + vola*np.sqrt(r[k]*deltaT)*phir2[k] + r[k]
            intr[k+1]=intr[k]+r[k+1]
            P[k] = BondPrice(nbSimulP,k,couponSchedule,c,lTP,r0,lambda0,gamma0,nu,etha,vola,alpha,beta,sigma,omega,R);
        paths[i] = {"Bond": P, "Interest rate": r, "Cum. interest rate": intr}
    return paths

paths = pathsGeneration(nbSimul,nbSimulP,lT,couponSchedule,c,lTP,r0,lambda0,gamma0,nu,etha,vola,alpha,beta,sigma,omega,R)

def alphaBetaTRS(TRSpaydate,prevTRSpaydate, paths, nbSimul):
    alphak = 0
    betak = 0
    for i in range(nbSimul):
        res1 = (paths[i]['Bond'][TRSpaydate]-paths[i]['Bond'][prevTRSpaydate])/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
        res2 = (paths[i]['Bond'][prevTRSpaydate])/(np.exp(paths[i]['Cum. interest rate'][TRSpaydate]))
        alphak += res1
        betak += res2
    alphak /= nbSimul
    betak /= nbSimul
    return [alphak,betak]

def cTRS(TRScoupondate, paths, nbSimul):
    ck = 0
    for i in range(nbSimul):
        res3 = 1/(np.exp(paths[i]['Cum. interest rate'][TRScoupondate]))
        ck += res3
    ck /= nbSimul
    return ck


def TRS(TRSpaymentSchedule,couponSchedule, sTRS, paths,nbSimul):
    res = 0
    alpha = 0
    beta = 0
    ceta = 0
    for s in TRSpaymentSchedule:
        alpha += alphaBetaTRS(s,TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1],paths,nbSimul)[0]
        beta += alphaBetaTRS(s,TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1],paths,nbSimul)[1]*(s-TRSpaymentSchedule[TRSpaymentSchedule.index(s)-1])
    for s in couponSchedule:
        ceta += cTRS(s, paths, nbSimul)
    res = alpha +c*ceta - sTRS*beta
    impliedRate = (alpha+c*ceta)/beta
    return [res, impliedRate]

print("TRS market value is",TRS(TRSpaymentSchedule,couponSchedule, sTRS, paths,nbSimul)[0])
print("TRS implied rate is",TRS(TRSpaymentSchedule,couponSchedule, sTRS, paths,nbSimul)[1])



