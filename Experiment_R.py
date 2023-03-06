import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import math
import strain
from variables import Nsat, degradation, deg_type, epsilon, nr_drops_total_mass, Km, Vmax, volume
from decimal import *
from scipy.special import lambertw
import os.path
getcontext().prec = 50

global AB_Deg_rate
AB_Deg_rate = 1e-2  # rate of AB degradation per cell per min

def scientific_notation(n):
    sc_not ="{:e}".format(n)
    return print(sc_not)



class Experiment_R(object):

    def __init__(self, strain_r, AB_conc, dt, t_end):
        self.strain_r = strain_r
        self.dt = dt
        self.t_end = t_end
        self.timesteps = round(t_end/dt)
        self.AB_conc = AB_conc

    def initialise(self, init_type):

        # Zero arrays:
        self.AB_conc_array = np.empty(self.timesteps + 1)
        self.N_array = np.empty(self.timesteps + 1)#.astype(float)
        self.strain_r.N = self.strain_r.initialN
        self.prob_extinction = np.empty(self.timesteps + 1)
        self.time_extinction = np.empty(self.timesteps + 1)

        self.AB_conc_array[0] = self.AB_conc
        self.prob_extinction[0] = None
        self.time_extinction[0] = None

        if init_type == "det":
            self.N_array[0] = self.strain_r.initialN
        elif init_type == "rand":
            #randomize initial starting numbers:
            self.strain_r.N = poisson.rvs(mu=self.strain_r.initialN, size=1)
            self.N_array[0] = self.strain_r.N
        else:
            print("Error in initialisation; type 'det' or 'rand' in run()")

    def degrade_ab_1_step(self, AB_concentration, N_bact, delta_t):

        if degradation == 'MM_linear':
            ...
        elif degradation == 'exponential':
            y = 1 / Km * AB_concentration
            exp_factor_1 = math.exp(AB_concentration / Km)
            exp_factor_2 = math.exp(-Vmax * N_bact * delta_t / (volume * Km))
            AB_concentration = Km * lambertw(y * exp_factor_1 * exp_factor_2).real
        return AB_concentration

    def run(self, init_type, grow_meth):  #CHOOSE INITIALIZATION AND GROWTH METHODS

        self.initialise(init_type) #det or random?
        for i in range(1, self.timesteps + 1):

            self.AB_conc_array[i] = self.degrade_ab_1_step(self.AB_conc_array[-1], self.strain_r.N, self.dt)
            # Grow strain for dt: ###tau_grow VS grow
            if grow_meth == 'binary':
                #print('strain_r.N', self.strain_r.N)
                self.strain_r.binary_grow(self.AB_conc_array[i])
            elif grow_meth == 'tau_binary':
                #print('strain_r.N', self.strain_r.N)
                self.strain_r.tau_binary_grow(self.AB_conc_array[i])
            elif grow_meth == 'balanced':
                #print('strain_r.N', self.strain_r.N)
                self.strain_r.balanced_grow(self.AB_conc_array[i])

            # Check if growth has saturated:
            if (self.strain_r.N > Nsat):
                # Set rest of list to Nsat:
                self.N_array[i:self.timesteps] = self.strain_r.N
                self.AB_conc_array[i:self.timesteps] = self.AB_conc_array[i]
                self.prob_extinction[i:self.timesteps] = 0
                break
            # Update array:
            self.N_array[i] = self.strain_r.N

        if grow_meth == "gillespie_balanced":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.gillespie_balanced_grow(self.AB_conc)
        if grow_meth == "gillespie_binary":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.gillespie_binary_grow(self.AB_conc)
        if grow_meth == "precheck_tau_binary":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.precheck_tau_binary_grow(epsilon, self.AB_conc)
        if grow_meth == "midpoint_tau_binary":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.midpoint_tau_binary_grow(epsilon, self.AB_conc)
        if grow_meth == "midpoint_tau_balanced":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.midpoint_tau_binary_grow(epsilon, self.AB_conc)
        if grow_meth == "adaptive_tau_binary":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.adaptive_tau_binary_grow(epsilon, self.AB_conc)
        if grow_meth == "adaptive_tau_balanced":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.adaptive_tau_balanced_grow(epsilon, self.AB_conc)



experiment_script_path = __file__
experiment_script_name = os.path.basename(__file__)
