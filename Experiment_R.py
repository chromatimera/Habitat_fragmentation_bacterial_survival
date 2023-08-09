import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import math
import strain
import variables
from variables import degradation, epsilon, Km, Vmax
from decimal import *
from scipy.special import lambertw
import os.path
getcontext().prec = 50

global AB_Deg_rate
AB_Deg_rate = 1e-2 # rate of AB degradation per cell per min

def scientific_notation(n):
    sc_not ="{:e}".format(n)
    return print(sc_not)



class Experiment_R(object):

    def __init__(self, strain_r, AB_conc, nr_drops_total_mass):
        self.strain_r = strain_r
        self.dt = variables.dt
        self.t_end = variables.t_end
        self.timesteps = round(self.t_end/self.dt)
        self.AB_conc = AB_conc
        self.Nsat = variables.Nsat * nr_drops_total_mass
        self.nr_drops_total_mass = nr_drops_total_mass
        self.volume = variables.volume * nr_drops_total_mass

    def initialise(self, init_type):

        # Zero arrays:
        self.AB_conc_array = np.empty(self.timesteps)
        self.N_array = np.empty(self.timesteps).astype(float)    #astype(float)
        self.strain_r.N = self.strain_r.initialN
        self.deg_list = []

        self.AB_conc_array[0] = self.AB_conc
        self.deg_list.append(0) #when initialising there's no degradation


        if init_type == "det":
            self.N_array[0] = self.strain_r.initialN
        elif init_type == "rand":
            #randomize initial starting numbers:
            self.strain_r.N = poisson.rvs(mu=self.strain_r.initialN, size=1)
            ##print('nr drops total mass', self.nr_drops_total_mass)
            #print('initial N', self.strain_r.N)
            self.N_array[0] = self.strain_r.N
        else:
            print("Error in initialisation; type 'det' or 'rand' in run()")
        #print('initial N', self.N_array[0])


    def degrade_ab_1_step_det(self, AB_conc, N_t, delta_t, volume, nr_drops_total_mass):
        if degradation == 'MM_linear':
            new_ab = AB_conc - Vmax * AB_conc / (Km + AB_conc) * N_t / volume * delta_t
            new_ab = max(0, new_ab)
            deg_linear = Vmax * AB_conc / (Km + AB_conc) * N_t / volume * delta_t
            deg_rate = 1 - deg_linear / AB_conc

        if degradation == 'MM_exponential':
            y = 1 / Km * AB_conc
            exp_factor_1 = math.exp(AB_conc / Km)
            exp_factor_2 = math.exp(-Vmax * delta_t * N_t / (volume * Km))
            new_ab = Km * lambertw(y * exp_factor_1 * exp_factor_2).real
            deg_rate = new_ab/AB_conc
        if degradation == 'exponential_decay':
            new_ab = AB_conc * math.exp(-AB_Deg_rate * N_t / nr_drops_total_mass * delta_t)
            deg_rate = math.exp(-AB_Deg_rate * N_t / nr_drops_total_mass * delta_t)
        if degradation == 'linear_decay':
            deg_linear = AB_Deg_rate * N_t * delta_t
            deg_rate = 1 - deg_linear
            new_ab = AB_conc - deg_linear * AB_conc
            new_ab = max(0, new_ab)
        self.deg_list.append(deg_rate)

        return new_ab #deg_rate, new_ab

        # AB is degraded at rate Vmax proportional to number of resistant cells:

    def run(self, init_type, grow_meth):  #CHOOSE INITIALIZATION AND GROWTH METHODS

        self.initialise(init_type) #det or random?
        self.ts = np.arange(variables.t_start, self.t_end, self.dt)

        if grow_meth == "gillespie_binary":
            self.ts, self.N_array, self.AB_conc_array = self.strain_r.gillespie_binary_grow(self.AB_conc)
        elif  grow_meth == 'binary':
         for i in range(0, self.timesteps-1):
            # Grow strain for dt: ###tau_grow VS grow
                self.AB_conc_array[i+1] = self.degrade_ab_1_step_det(self.AB_conc_array[i], self.strain_r.N, self.dt,
                                           self.volume, self.nr_drops_total_mass) #self.deg_list[i+1],
                self.strain_r.binary_grow(self.AB_conc_array[i+1])
                self.N_array[i+1] = self.strain_r.N
            # Check if growth has saturated:
            #if (self.strain_r.N > strain.Nsat):
                # Set rest of list to Nsat:
             #   self.N_array[i+1:self.timesteps] = self.strain_r.N
              #  self.AB_conc_array[i+1:self.timesteps] = self.AB_conc_array[i+1]
            #break


experiment_script_path = __file__
experiment_script_name = os.path.basename(__file__)
