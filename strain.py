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
        self.deathrate = variables.deathrate
        self.MIC = variables.MIC
        self.clock = 1
        self.t_end = variables.t_end
        self.dt = variables.dt
        self.timesteps = round(self.t_end/self.dt)
        self.Nsat = variables.Nsat * nr_drops_total_mass
        self.nr_drops_total_mass = nr_drops_total_mass
        self.volume = variables.volume * nr_drops_total_mass

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
        N_population_array = np.append(N_population_array, N_population_array[-1] + stoichiometry[chosen_reaction])

        return time_array, AB_conc_array, N_population_array

    @staticmethod
    def precheck_tau_binary_g(MIC, trial_tau, epsilon, N_t, AB_conc_t, v, growthrate, deathrate):
        # stoichiometry vector
        stoichiometry = np.array([1, -1])

        rates, a, a0 = strain.set_rates_and_prop_binary(AB_conc_t, MIC, N_t, growthrate, deathrate)

        # expectation value
        ksi = a[0] * stoichiometry[0] + a[1] * stoichiometry[1]
        exp_val = ksi * trial_tau

        if MIC > AB_conc_t:  # GROWTH
            # gradient of change in propensity function
            b_growth = v / trial_tau

            ## check dt
            lhs = abs(trial_tau * ksi * b_growth)
            rhs = epsilon * N_t * rates[0]

            if lhs <= rhs:
                return 1
            else:
                return 0

        else:  # DEATH
            b_death = -v / trial_tau

            ## check dt
            lhs = abs(trial_tau * ksi * b_death)
            rhs = epsilon * N_t * rates[1]

            if lhs <= rhs:
                return 1
            else:
                return 0

    """
       Return the force on a particle in a double well potential.

       The force is given by
           F(x) = -dV/dx = -4*a*x^3 + 2*b*x

       Parameters
       ----------
       p1: Particle1D
       a: float
       b: float

       Returns
       -------
       force: float
       """
    @staticmethod
    def set_K_values(list_boolean_crit, index_noncrit, prop_list, state_change_vector_list, tau):
        # sum of the propensity functions of the critical reactions
        a0_crit = 0
        # i takes the values of the indeces of the list having True for critial rates, False for non critical
        for i in np.where(list_boolean_crit)[0]:
            a0_crit += prop_list[i]
        if a0_crit == 0:
            t_second = float('inf')
        else:
            t_second = np.random.Generator.exponential(1 / a0_crit)

        # extracted Kj values from Poisson distrib w r.v. w mean aj(x)*tau
        K = np.zeros(len(state_change_vector_list))
        if tau < t_second:
            for i in index_noncrit:
                K[i] = np.random.poisson(prop_list[i] * tau)
        elif t_second <= tau:
            tau = t_second

            ## pick which critical reaction to fire in this leap
            # list of indeces of Critical eractions, compute associated probabilities of drawing and
            # draw the critical reaction that is filled over the leap
            aj_crit = np.where(list_boolean_crit)[0]
            prob_crit = [x / a0_crit for x in prop_list[aj_crit]]
            jc = random.choice(aj_crit, weights=prob_crit, k=1)
            for i in aj_crit:
                K[i] = 0
            K[jc] = 1

            for i in index_noncrit:
                K[i] = np.random.poisson(prop_list[i] * tau)
        return K


    def binary_grow(self, AB_conc):
        #print('initialN', self.initialN)
        #print('Nsat', self.Nsat)
        if self.MIC > AB_conc:
            self.N = self.N + self.N * self.growthrate * self.dt  ## N grows uninhibited, normalized by dt
        else:
            self.N = self.N - self.N * self.deathrate * self.dt   ## N dies with growthrate normalized by dt
        if self.N < 1.0:  #### #can't have less than 1 bacteria ??
           self.N = 0

    def midpoint_tau_binary_grow(self, epsilon, AB_conc):

        if self.initialN != 0:

            # stoichiometry vector (or state change vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values

            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])
            #print(self.dt)
            s = 0
            tau_list = []

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):

                rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array[-1], self.MIC, self.N[-1], self.growthrate, self.deathrate)

                # expectation value
                ksi = a[0] * stoichiometry[0] + a[1] * stoichiometry[1]
                exp_val = ksi * self.dt
                N_prime = self.N[-1] + exp_val / 2
                #print('N_prime', N_prime)

            ## check v
                if self.MIC > self.AB_conc_array[-1]:  # GROWTH

                    v = np.random.poisson(N_prime * self.growthrate * self.dt, 1)  # second Poisson method
                    change_in_N = v * stoichiometry[0]

                    # gradient of change in propensity function
                    check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, N_prime, self.AB_conc_array[-1], v)
                    tau = self.dt
                    while check_returned == 0:
                        tau /= 2.0
                        exp_val = ksi * tau
                        N_prime = self.N[-1] + exp_val / 2
                        #print('N_prime', N_prime)

                        v = np.random.poisson(N_prime * self.growthrate * tau, 1)
                        change_in_N = v * stoichiometry[0]
                        #print('change in N', change_in_N)

                        check_returned = self.precheck_tau_binary_g(self.MIC, rates, tau, epsilon, N_prime, self.AB_conc_array[-1], v)

                    if tau > (2/a0):

                        #update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                        #update Ab_conc for next step
                        self.AB_conc_array = self.degrade_ab_1_step(self.AB_conc_array, self.N, tau, self.volume)
                        #update N
                        self.N = np.append(self.N, self.N[-1] + change_in_N)
                        #print('actual N', self.N[-1])
                        s += 1
                        gill_trial = False

                    else:
                        gill_trial = True



                else:  # DEATH
                    # print(AB_conc)

                    v = np.random.poisson(N_prime * self.deathrate * self.dt, 1)  # second Poisson method
                    change_in_N = v * stoichiometry[1]
                    #print('change in N', change_in_N)

                    check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, N_prime, self.AB_conc_array[-1], v)
                    tau = self.dt
                    while check_returned == 0:
                        tau /= 2.0
                        exp_val = ksi * tau
                        N_prime = self.N[-1] + exp_val / 2
                        #print('N_prime', N_prime)


                        v = np.random.poisson(N_prime * self.deathrate * tau, 1)
                        change_in_N = v * stoichiometry[1]
                        #print('change in N', change_in_N)

                        check_returned = self.precheck_tau_binary_g(self.MIC, rates, tau, epsilon, N_prime, self.AB_conc_array[-1], v)

                    if self.dt > (2/a0):
                        #update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                        #update Ab_conc for next step
                        self.AB_conc_array = self.degrade_ab_1_step(self.AB_conc_array, self.N, tau, self.volume)
                        #update N
                        self.N = np.append(self.N, self.N[-1] - change_in_N)
                        #print('actual N', self.N[-1])
                        s += 1
                        gill_trial = False
                    else:
                        gill_trial = True

                if gill_trial == True:
                    self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array, self.AB_conc_array[-1], self.N, AB_Deg_rate, stoichiometry, self.volume)

    def adaptive_tau_binary_grow(self, epsilon, AB_conc):
        def compute_L_and_crit_arrays(stoichiometry, a, n_crit):
            L = np.zeros(2)
            if self.AB_conc_array[-1] <= self.MIC:
                L[0] = math.floor(self.N[-1] / abs(stoichiometry[0]))  # L_g
                L[1] = 0
            else:
                L[0] = 0
                L[1] = math.floor(self.N[-1] / abs(stoichiometry[1]))  # L_d
            #print('L0, L1,', L)

            ## compute array classifying is reaction is critical or not
            is_crit = []
            index_noncrit = []

            for i in range(len(a)):
                if a[i] > 0 and L[i] <= n_crit:
                    is_crit.append(True)
                else:
                    is_crit.append(False)
                    index_noncrit.append(i)
            return L, is_crit, index_noncrit


        if self.initialN != 0:

            # stoichiometry vector or state change vector (vij)
            stoichiometry = np.array([1, -1])

            # set up time and initial values
            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])
            #print(self.dt)
            s = 0

            # keep a record of the change in propensity fct
            #a_history = []
            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):
                rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)

                ## compute Lj the maximum nr of permitted firings
                L, is_crit, index_noncrit = compute_L_and_crit_arrays(stoichiometry, a, n_crit)

                #here are no noncritical reactions
                if not any(is_crit) == False:
                    tau = float('inf')
                else:
                    mean_array = 0
                    variance_array = 0
                    tau_less_ea0 = np.zeros(2)
                    g = 1
                    for j in range(len(rates)):
                        for j_prime in index_noncrit:
                            mean_array += stoichiometry[j_prime] * a[j_prime]
                            variance_array += stoichiometry[j_prime] ** 2 * a[j_prime]
                        tau_less_ea0[0] = max(epsilon * self.N[-1]/g, 1) / abs(mean_array)
                        tau_less_ea0[1] = max(epsilon * self.N[-1]/g, 1) ** 2 / variance_array
                        tau = min(tau_less_ea0)
                    if tau < (10/a0):
                        gill_trial = True
                    else:
                        gill_trial = False

                    if gill_trial == False:
                        K = self.set_K_values(is_crit, index_noncrit, a, stoichiometry, tau)

                        if self.N[-1] + K[0] * stoichiometry[0] >= 0 or self.N[-1] + K[1]*stoichiometry[1] >= 0:
                            reduce_t = False
                        else:
                            reduce_t = True
                        red = 0
                        while (reduce_t == True):
                            tau = tau/2  # populate with redo all steps but with time t/2
                            if tau < (10 / a0):
                                gill_trial = True
                            else:
                                gill_trial = False

                            if gill_trial == False:
                                K = self.set_K_values(is_crit, index_noncrit, a, stoichiometry, tau)

                                if self.N[-1] + K[0] * stoichiometry[0] >= 0 or self.N[-1] + K[1] * stoichiometry[1] >= 0:
                                    reduce_t = False
                                else:
                                    reduce_t = True
                            red += 1

                        if self.MIC > self.AB_conc_array[-1]:  # GROWTH
                                #update N
                                self.N = np.append(self.N, self.N[-1] + K[0] * stoichiometry[0])

                        else:  # DEATH
                                #update N
                                self.N = np.append(self.N, self.N[-1] + K[1] * stoichiometry[1])

                        # update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + tau)
                        # update Ab_conc for next step
                        self.AB_conc_array = self.degrade_ab_1_step(self.AB_conc_array, self.N, tau, self.nr_drops_total_mass, self.volume)
                        s += 1


                    elif gill_trial == True:
                            for _ in range(20):
                                while (self.t_array[s] < self.t_end) and (self.N[-1] > 0) and (self.N[-1] < self.Nsat):
                                    rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)

                                    self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array,
                                                                                                          self.AB_conc_array, self.N, stoichiometry, self.nr_drops_total_mass, self.volume)
                                    s += 1
                                break
        return self.t_array, self.N, self.AB_conc_array

    def gillespie_binary_grow(self, AB_conc):
        #print('initialN', self.initialN)
        #print('Nsat', self.Nsat)
        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values
            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])

            s = 0
            #print(MIC, AB_conc)

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):
                rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)
                self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array, self.AB_conc_array, self.N, stoichiometry, self.nr_drops_total_mass, self.volume)
                # keep counters for timesteps and AB_conc
                s += 1

            return self.t_array, self.N, self.AB_conc_array




strain_script_path = __file__
strain_script_name = os.path.basename(__file__)
