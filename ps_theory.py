import math

import variables
import numpy as np
import matplotlib.pyplot as plt
from variables import *
import scipy.special as sc
from mpmath import *
import pandas as pd

os.chdir('./output/')

def range_prod(lo,hi):
    if lo+1 < hi:
        mid = (hi+lo)//2
        return range_prod(lo,mid) * range_prod(mid+1,hi)
    if lo == hi:
        return lo
    return lo*hi


def treefactorial(n):
    if n < 2:
        return 1
    return range_prod(1,n)

def calc_theo_survival_prob(vol_fac):
    b = 1
    Ab_concs = [ 35, 55, 75]

    rho_bulk = variables.initialN / variables.volume # constant in det # rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
    bigPs= np.empty([len(Ab_concs),len(vol_fac)])
    ps =np.empty([len(Ab_concs),len(vol_fac)])

    #for each AB conc;
    for aa in range(len(Ab_concs)):
     A= Ab_concs[aa]
     F1 = deathrate /(b * Vmax)
     F = (A - MIC) + Km * np.log(A / MIC)  # eq4
    ##calculate rho_threshold; eq (9) from paper
     rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min (per enzyme)

     for mm in range(len(vol_fac)):
       m=vol_fac[mm]
       vol=1E-4 /m #ml

       lam = rho_bulk * vol
       #N_T = (rho_T * vol)
       N_T = np.ceil(rho_T * vol)


    ## calculate the theoretical survival probability; eq. (10) from paper
    # ps  (prob that N(0) > Nt for a given droplet)

       exp_fact = exp(-lam) # math.exp;;goes to zero when lam is too high


       #ps[aa,mm] = exp_fact * nsum(lambda j: (lam) ** j / fac(j), [N_T + 1 , inf],  method='r+s+e') # from nt+1 to inf  # method='r+s+e' takes ages
       ps[aa, mm]=  1 - sc.gammaincc(N_T+1, lam)
       #ps_int=(exp_fact * nsum(lambda j: (lam) ** j / fac(j), [0, N_T - 1]))  # from 0 to nt
       #ps[aa, mm] = 1 - ps_int
       bigPs[aa, mm] = 1 - (1 - ps[aa, mm])**m  # prob of at least 1 subvol surviving



       # save:
    np.save('prob_line.npy', bigPs)  #
    pd.DataFrame(bigPs).to_csv('prob_line.csv'.format(), index=None)
    #plt.plot(vol_fac, bigPs[4, :])
    #plt.show()
    return ps, bigPs


def new_theory_calc(vol_fac):
    b = 1
    Ab_concs = [75]
    rho_bulk = variables.initialN / variables.volume # constant in det # rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
    bigPs= np.empty([len(Ab_concs),len(vol_fac)])
    ps =np.empty([len(Ab_concs),len(vol_fac)])
    subvol = np.empty([len(Ab_concs), len(vol_fac)])
    #for each AB conc;
    for aa in range(len(Ab_concs)):
     A= Ab_concs[aa]
     F1 = deathrate /(b * Vmax)
     F = (A - MIC) + Km * np.log(A / MIC)  # eq4
    ##calculate rho_threshold; eq (9) from paper
     rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min (per enzyme)
     f=rho_bulk/rho_T

     for mm in range(len(vol_fac)):
       m=vol_fac[mm]
       vol=1E-4 /m #ml
       subvol[aa,mm]=vol
       lam = rho_bulk * vol
       #N_T = (rho_T * vol)
       N_T = np.ceil(rho_T * vol)


    ## calculate the theoretical survival probability; eq. (10) from paper
    # ps  (prob that N(0) > Nt for a given droplet)

       exp_fact = exp(-lam) # math.exp;;goes to zero when lam is too high


       #ps[aa,mm] = exp_fact * nsum(lambda j: (lam) ** j / fac(j), [N_T + 1 , inf],  method='r+s+e') # from nt+1 to inf  # method='r+s+e' takes ages
       ps[aa, mm]=  1 - sc.gammaincc(N_T+1, lam)
       #ps_int=(exp_fact * nsum(lambda j: (lam) ** j / fac(j), [0, N_T - 1]))  # from 0 to nt
       #ps[aa, mm] = 1 - ps_int
       bigPs[aa, mm] = 1 - (1 - ps[aa, mm])**m  # prob of at least 1 subvol surviving

      # we could now calculate C? Is that useful?


       # save:
    np.save('prob_line.npy', bigPs)
    np.save('subvol.npy', subvol) #
    pd.DataFrame(bigPs).to_csv('prob_line.csv'.format(), index=None)
    #plt.plot(vol_fac, bigPs[4, :])
    #plt.show()
    return ps, bigPs

def for_loop_sum(vol_fac):  # paper eq
    #zz = np.load('prob_line_forloop.npy')
    b = 1
    Ab_concs = [35, 55, 75]
    rho_bulk = variables.initialN / variables.volume # constant in det # rho_bulk = variables.initialN * variables.total_drop_nr/variables.volume * variables.total_drop_nr
    bigPs= np.empty([len(Ab_concs),len(vol_fac)])
    ps =np.empty([len(Ab_concs),len(vol_fac)])

    #for each AB conc;
    for aa in range(len(Ab_concs)):
     A= Ab_concs[aa]
     F1 = deathrate /(b * Vmax)
     F = (A - MIC) + Km * np.log(A / MIC)  # eq4
    ##calculate rho_threshold; eq (9) from paper
     rho_T = F1 * F ## **units: cell/vol if Vmax is ug/min (per enzyme)

     for mm in range(len(vol_fac)):
       m=vol_fac[mm]
       vol=1E-4 /m #ml

       lam = rho_bulk * vol
       #N_T = (rho_T * vol)
       N_T = np.ceil(rho_T * vol).astype('int')

    ## calculate the theoretical survival probability; eq. (10) from paper
    # ps  (prob that N(0) > Nt for a given droplet)
       exp_fact = exp(-lam) # math.exp;;goes to zero when lam is too high
       probability = np.empty(N_T-1)
       for i in range(N_T-1):
           probability[i] = ((lam ** i)*exp_fact) / math.factorial(i)
       ps[aa,mm] = 1- sum(probability)
       bigPs[aa, mm] = 1 - (1 - ps[aa, mm])**m  # prob of at least 1 subvol surviving
       np.save('prob_line_forloop.npy', bigPs)
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

vol_fac = np.arange(1,1000,2)
RES = calc_theo_survival_prob(vol_fac)

