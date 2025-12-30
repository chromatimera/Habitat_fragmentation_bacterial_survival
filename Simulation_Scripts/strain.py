import math
import variables
from variables import *
import os.path
import numpy as np
import random
from scipy.special import lambertw


global AB_Deg_rate
AB_Deg_rate = 1e-2  # rate of AB degradation per cell per min

def scientific_notation(n):
    sc_not ="{:e}".format(n)
    return sc_not

class strain(object):

    def __init__(self, nr_drops_total_mass):
        self.initialN = variables.initialN * nr_drops_total_mass
        self.N = self.initialN
        self.growthrate = variables.growthrate
        self.slowrate = variables.slowrate
        self.deathrate = variables.deathrate
        self.MIC = variables.MIC
        self.clock = 1
        self.t_end = variables.t_end
        self.dt = variables.dt
        self.timesteps = round(self.t_end/self.dt)
        self.Nsat = variables.Nsat * nr_drops_total_mass
        self.nr_drops_total_mass = nr_drops_total_mass
        self.volume = variables.volume * nr_drops_total_mass
        #print(nr_drops_total_mass)
        #print('initialN strain', self.initialN)
    @staticmethod
    def set_rates_and_prop_binary(AB_conc_array, MIC, N_population_array, growthrate, deathrate):
        if AB_conc_array[-1] > MIC:
            deathrate = deathrate
            growthrate = 0
        else:
            deathrate = 0
            growthrate = growthrate

        ### state-change vector vj
        rates = np.array([growthrate, deathrate], dtype=object)

        # for determining of next reaction
        a = np.zeros(2)
        a[0] = rates[0] * N_population_array[-1]
        a[1] = rates[1] * N_population_array[-1]
        a0 = sum(a)
        return rates, a, a0

    """
       Evolve death rate as a function of antibiotic concentration.
       """
    def calculate_death_rate(self, AB_conc, E_max, EC50, k):
        """
        Calculate death rate using Hill-type model.

        Parameters:
        AB_conc (float): Antibiotic concentration (ug/mL).
        E_max (float): Maximum death rate.
        EC50 (float): Antibiotic concentration at half-maximal effect.
        n (float): Hill coefficient (steepness of the curve).

        Returns:
        float: Death rate.
        """
        death_rate = (E_max * (AB_conc ** k)) / ((AB_conc ** k) + (EC50 ** k))
        return death_rate

    def calculate_death_rate_mirrored(self, AB_conc, E_max, EC50, k):
        """
        Calculate death rate using Hill-type model.

        Parameters:
        AB_conc (float): Antibiotic concentration (ug/mL).
        E_max (float): Maximum death rate.
        EC50 (float): Antibiotic concentration at half-maximal effect.
        n (float): Hill coefficient (steepness of the curve).

        Returns:
        float: Death rate.
        """
        death_rate = (E_max * (EC50 ** k)) / ((AB_conc ** k) + (EC50 ** k))
        return death_rate
    """
   Evolve antibiotic concentration array for one tau 

   The antibiotic can be degraded in 3 ways: based on Mikaelis-Menten linear, exponential or just exponential decay - I kept this as all simulations so far with exponential decay.
    *Mikaelis-Menten linear:
    
    *Mikaelis-Menten exponential:
    
    *Simple exponential decay

   Parameters
   ----------
   
   Returns
   -------
   force: time array, ab array, N array
   """
    @staticmethod
    def degrade_ab_1_step(AB_array, N_bact_pop, delta_t, nr_drops_total_mass, volume):
        if degradation == 'MM_linear':
            new_ab = AB_array[-1] - Vmax * AB_array[-1] / (Km + AB_array[-1]) * N_bact_pop[-1] / volume * delta_t
            new_ab = max(0, new_ab)
            deg_linear = Vmax * AB_array[-1] / (Km + AB_array[-1]) * N_bact_pop[-1] / volume * delta_t
            deg_rate = 1 - deg_linear / AB_array[-1]

        if degradation == 'MM_exponential':
            exp_factor_1 = math.exp(AB_array[-1] / Km)
            exp_factor_2 = math.exp(-Vmax * delta_t * N_bact_pop[-1] / (volume * Km))
            new_ab = Km * lambertw(AB_array[-1] / Km * exp_factor_1 * exp_factor_2).real
            deg_rate = new_ab / AB_array[-1]
        if degradation == 'exponential_decay':
            new_ab = AB_array[-1] * math.exp(-AB_Deg_rate * N_bact_pop[-1] / nr_drops_total_mass * delta_t)
            deg_rate = math.exp(-AB_Deg_rate * N_bact_pop[-1] / nr_drops_total_mass * delta_t)
        if degradation == 'linear_decay':
            deg_linear = AB_Deg_rate * N_bact_pop[-1] * delta_t
            deg_rate = 1 - deg_linear
            new_ab = AB_array[-1] - deg_linear * AB_array[-1]
            new_ab = max(0, new_ab)
        AB_array = np.append(AB_array, new_ab)
        return AB_array

    """
   Evolve nr of bacteria array, antibiotic concentration array and time array for one tau. 

   The nr of bacteria N(t) can either grow or die by one bacteria.
       N(t + tau) = N(t) +/- 1, given in the code by stoichiometry[chosen reaction]

   Parameters
   ----------
   a, a0 - propensities functions (probabilities of growing/dying)
   time array , ab array, N array
   stoichiometry - [-1, 1], can add or remove one bacteria based on the chosen reaction (the index in the stoichiometry array that is chosen based on previous probabilities)

   Returns
   -------
   time array, ab array, N array
   """
    @staticmethod
    def evolve_gillespie_one_step(a, a0, time_array, AB_conc_array, N_population_array, stoichiometry, nr_drops_total_mass, volume):
        # determine time of next reaction for gillespie
        r1 = random.uniform(0, 1)
        tau = (-math.log(r1) / a0)

        # determine nature of next reaction
        r2 = random.uniform(0, 1)
        acumsum = np.divide(np.cumsum(a), a0)
        chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]

        time_array = np.append(time_array, time_array[-1] + tau)
        AB_conc_array = strain.degrade_ab_1_step(AB_conc_array, N_population_array, tau, nr_drops_total_mass, volume)
        N_population_array = np.append(N_population_array, N_population_array[-1] + stoichiometry[chosen_reaction])   ###stuck

        return time_array, AB_conc_array, N_population_array

    def binary_grow(self, AB_conc): #DET growth
        if self.MIC > AB_conc:
            self.N = self.N + self.N * self.growthrate * self.dt  ## N grows uninhibited, normalized by dt
        else:
            self.N = self.N - self.N * self.deathrate * self.dt   ## N dies with growthrate normalized by dt
        if self.N < 1:  #### #can't have less than 1 bacteria
           self.N = 0

    def resource_growth(self, AB_conc): #DET growth
        if self.MIC > AB_conc:
            self.N = self.N + self.N * self.growthrate * self.dt  ## N grows uninhibited, normalized by dt
        else:
            self.N = self.N + self.N * self.slowrate * self.dt   ## N grows slowly until resource is released to threshold density
        if self.N < 1.0:  #### #can't have less than 1 bacteria ??
           self.N = 0


    def gillespie_binary_grow(self, AB_conc):   #??
        #print('initialN', self.initialN)
        #print('Nsat', self.Nsat)
        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values
            self.t_array = np.array([0])
            self.N = np.array([round(self.N)])
            self.AB_conc_array = np.array([AB_conc])

            s = 0
            #print(MIC, AB_conc)

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):
                rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)
                self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array, self.AB_conc_array, self.N, stoichiometry, self.nr_drops_total_mass, self.volume)
                # keep counters for timesteps and AB_conc
                s += 1
           # print('final N', self.N)
            return self.t_array, self.N, self.AB_conc_array

    def balanced_grow(self, AB_conc):
        # hill function
        dynamic_deathrate = self.calculate_death_rate(AB_conc, variables.E_max, variables.EC50, variables.k)
        # mirrored hill function
        dynamic_deathrate = self.calculate_death_rate_mirrored(AB_conc, variables.E_max, variables.EC50, variables.k)
        print("dynamic_deathrate", dynamic_deathrate)
        print("self.growthrate", self.growthrate)
        self.N = self.N + self.N * self.growthrate * self.dt - self.N * dynamic_deathrate * self.dt
        if self.N < 1:  #### #can't have less than 1 bacteria
           self.N = 0




strain_script_path = __file__
strain_script_name = os.path.basename(__file__)
