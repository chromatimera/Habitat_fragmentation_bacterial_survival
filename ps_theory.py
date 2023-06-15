import math

import variables
import numpy as np
import matplotlib.pyplot as plt
from variables import *
import scipy.special as sc


def calc_theo_survival_prob():

    b = 1
    rho_bulk = variables.initialN / variables.volume # constant in det # rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
    bigPs= np.empty([10,10])
    ps =np.empty([10,10])
    #for each AB conc;
    for A in [5,10,15,25,35,45,55,65,75,100]:
     F1 = deathrate /(b * Vmax)
     F = (A - MIC) + Km * np.log(A / MIC)  # eq4
    ##calculate rho_threshold; eq (9) from paper
     rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min (per enzyme)

     for m in [1,2,3,5,10,100, 125,250, 500,1000]:
    ## calculate N_T
       vol=10E-4 /m #ml
       lam = np.floor(rho_bulk * vol)
       N_T = np.floor(rho_T * vol)
    ## calculate the theoretical survival probability; eq. (10) from paper
       prob = np.empty([int(N_T)])
       for i in range(1,int(N_T)):
           try:
               prob[i] = (lam ** i) / math.factorial(i)
           except OverflowError:
               prob[i]=math.inf

    # ps  (prob that N(0) > nt for a given droplet)
       ps[A,m]= 1 - (np.sum(prob) * math.exp(-lam))
       bigPs[A,m]= 1- (1 - ps[A,m]) ** m
       return ps, bigPs




def upper_incomplete_gamma(a, x):  # Lentz's method
    eps = 1e-8  # Desired accuracy
    max_iter = 100  # Maximum number of iterations

    # Initialize the continued fraction
    f = 1.0
    C = f
    D = 0.0

    # Iteratively compute the continued fraction
    for i in range(max_iter):
        an = -x + (a + i) * (i + 1)
        f = an * f + 1.0
        if abs(f) < eps:
            f = eps

        D = 1.0 / (an + D)
        C = an + 1.0 / C
        delta = C * D
        f *= delta

        if abs(delta - 1.0) < eps:
            break

    # Compute the upper incomplete gamma function
    gamma = math.exp(-x) * f
    return gamma


def unsimplified_calc(a, x):  # paper eq

    nt=a
    lam=x  #pbulk*v
    probability=np.empty(nt)
    for i in range(nt):
        probability[i]= (lam**i) / math.factorial(i)

    #ps  (prob that N(0) > nt for a given droplet)
    ps=1-(np.sum(probability) *math.exp(-lam))

    # Compute the upper incomplete gamma function
    gamma = (1-ps)*math.factorial(nt)
    return ps,gamma


# RUN
a = 27
x = 10
result = upper_incomplete_gamma(a, x)
print(result)
pythonUIG=G = sc.gammaincc(a ,x)
print(pythonUIG)
result2=unsimplified_calc(a, x)
print(result2)
bigPs= 1- (1 - result2[0])**500
print(bigPs)

RES=calc_theo_survival_prob()
#save