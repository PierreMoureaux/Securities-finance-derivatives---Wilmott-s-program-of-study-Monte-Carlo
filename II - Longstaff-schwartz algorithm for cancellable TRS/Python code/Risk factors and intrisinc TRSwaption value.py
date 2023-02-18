import numpy as np
from numpy.random import RandomState
from matplotlib import pyplot as plt

class RepoRate:
    """Repo rate"""

    def __init__(self, gamma: float = 1.0, sigmar: float = 1.0, r0: float = 0.1):
        self.gamma = gamma
        self.sigmar = sigmar
        self.r0 = r0

    def simulate(self, t: np.array, n: int, rnd: np.random.RandomState) -> np.ndarray:
        assert t.ndim == 1, "One dimensional time vector required"
        assert t.size > 0, "At least one time point is required"
        dt = np.concatenate((t[0:1], np.diff(t)))
        assert (dt >= 0).all(), "Increasing time vector required"
        dW = (rnd.normal(size=(t.size, n)).T * np.sqrt(dt)).T
        W = np.cumsum(dW, axis=0)
        return (self.gamma+self.r0)*np.exp(self.sigmar * W.T - (self.sigmar**2 / 2) * t).T-self.gamma

class AssetPath:
    """AssetPath with dedicated r_repo linked drift"""

    def __init__(self, rrepo: RepoRate, sigma: float = 1.0, S0: float = 100):
        self.rrepo = rrepo
        self.sigma = sigma
        self.S0 = S0

    def simulate(self, t: np.array, n: int, rnd: np.random.RandomState) -> np.ndarray:
        assert t.ndim == 1
        assert t.size > 0
        dt = np.concatenate((t[0:1], np.diff(t)))
        assert (dt >= 0).all()
        dW = (rnd.normal(size=(t.size, n)).T * np.sqrt(dt)).T
        W = np.cumsum(dW, axis=0)
        return self.S0*np.exp(self.sigma * W.T + ((self.rrepo.simulate(t, n, rnd) - self.sigma**2 / 2).T) * t).T
    
class Annuity:
    """Annuity, mandatory for option pay-off, assuming constant and deterministic interest rate"""

    def __init__(self,r: float):
        self.r = r
    
    def Bond(self,t: np.array):
        return np.exp(-self.r*(t[-1] - t));

    def simulate(self, t: np.array, n: int) -> np.ndarray:
        assert t.ndim == 1
        assert t.size > 0 
        b = self.Bond(t)
        interMat = np.tile(b,(t.size,1)).T
        tenor = np.zeros(t.size)
        for i in range(1,t.size):
            tenor[i] = t[i] - t[i-1]
        return np.dot(interMat, tenor)

def fee(factor,t):
    return factor*(t[-1]-t)/t[-1]

def ExerciseValue(sTRS,a,p,S0,r,t,n,factor):
    return np.tile(sTRS*a,(n,1)).T-(((p/S0*r).T)*(t[-1]-t)).T-np.tile(fee(factor,t),(n,1)).T

"""Repo rate generation"""
gamma = 0.02
sigmar = 0.15
r0 = 0.1
rr = RepoRate(gamma=gamma, sigmar=sigmar, r0 = r0)
t = np.linspace(0, 5, 12 * 5)
rnd = RandomState(seed=1234)
X = rr.simulate(t, 50, rnd)

"""Asset generation"""
sigma = 0.15
S0 = 100
ap = AssetPath(rrepo=rr, sigma=sigma, S0 = S0)
Y = ap.simulate(t, 50, rnd)

"""Annuity generation"""
r = 0.15
a = Annuity(r=r)
Z = a.simulate(t, 50)

"""Payoff generation"""
sTRS = 0.05
po = ExerciseValue(sTRS,Z,Y,S0,X,t,50,-0.5)

"""Graph generation"""
figsize = (8, 6)
plt.figure(figsize=figsize)
plt.plot(t, X)
plt.xlabel("Time t")
plt.ylabel("Repo Rate")
plt.figure(figsize=figsize)
plt.plot(t, Y)
plt.xlabel("Time t")
plt.ylabel("Stock Price")
plt.figure(figsize=figsize)
plt.plot(t, Z)
plt.xlabel("Time t")
plt.ylabel("Annuity")
plt.figure(figsize=figsize)
plt.plot(t, po)
plt.xlabel("Time t")
plt.ylabel("Exercise Value")
