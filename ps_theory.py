import math

import variables
import numpy as np
import matplotlib.pyplot as plt
from variables import *
import scipy.special as sc
from mpmath import *
import pandas as pd


def calc_theo_survival_prob():
    b = 1
    Ab_concs=[5, 10, 15, 25, 35,45,55,65,75,100]
    vol_fac=[1, 2, 5, 10, 100, 125, 250, 500, 1000]
    rho_bulk = variables.initialN / variables.volume # constant in det # rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
    bigPs= np.empty([len(Ab_concs),len(vol_fac)])
    ps =np.empty([len(Ab_concs),len(vol_fac)])

    #for each AB conc;
    for aa in range (0,10):
     A= Ab_concs[aa]
     F1 = deathrate /(b * Vmax)
     F = (A - MIC) + Km * np.log(A / MIC)  # eq4
    ##calculate rho_threshold; eq (9) from paper
     rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min (per enzyme)

     for mm in range (0, 9):
       m=vol_fac[mm]
       vol=1E-4 /m #ml
       lam = np.floor(rho_bulk * vol)
       N_T = np.floor(rho_T * vol)  #is this at least? OR inclusive
    ## calculate the theoretical survival probability; eq. (10) from paper
    # ps  (prob that N(0) > Nt for a given droplet)
       exp_fact = exp(-lam) # math.exp;;goes to zero when lam is too high
       #ps[aa,mm] = exp_fact * nsum(lambda j: (lam) ** j / fac(j), [N_T + 1, inf]) # from nt+1 to inf  # method='r+s+e' takes ages
       ps_int=(exp_fact * nsum(lambda j: (lam) ** j / fac(j), [0, N_T]))  # from 0 to nt
       ps[aa, mm] = 1 - ps_int
       bigPs[aa, mm] = 1 - (1 - ps[aa, mm])**m  # prob of at least 1 subvol surviving

       # save:
    np.save('prob_line.npy', bigPs)  #    pd.DataFrame(bigPs_sum).to_csv('prob_line.csv'.format(), index=None)
    plt.plot(vol_fac, bigPs[4, :])
    plt.show()
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
#a = 27
#x = 10
#result = upper_incomplete_gamma(a, x)
#print(result)
#pythonUIG=G = sc.gammaincc(a ,x)
#print(pythonUIG)
#result2=unsimplified_calc(a, x)
#print(result2)
#bigPs= 1- (1 - result2[0])**500
#print(bigPs)

RES=calc_theo_survival_prob()

