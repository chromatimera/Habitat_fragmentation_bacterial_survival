from scipy.stats import poisson
import math
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

    def __init__(self, initialN, growthrate, deathrate, MIC, dt, t_end, Nsat): #prob_ext, time_ext):
        self.initialN = initialN
        self.N = initialN
        self.growthrate = growthrate  #uninhibited growth rate in minutes
        self.deathrate = deathrate    #uninhibited deaath rate in minutes
        # MIC is given initially in ug/ml -> transformed in Droplet_R script; amp molar mass is 349.406 g/mol (ug/umol); * 1e-3 is to go from ml to ul; now MIC is in moles/ul
        self.MIC = MIC
        self.clock = 1
        self.t_end = t_end
        self.dt = dt
        self.timesteps = round(self.t_end/self.dt)
        self.Nsat = Nsat

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
        # print(rates)
        # print('last N', self.N[-1])

        # for determining of next reaction
        a = np.zeros(2)
        a[0] = rates[0] * N_population_array[-1]
        a[1] = rates[1] * N_population_array[-1]
        a0 = sum(a)
        return rates, a, a0
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
    def set_rates_and_prop_balanced(AB_conc_array, MIC, N_population_array, growthrate, deathrate):
        deathrate = AB_conc_array[-1] * 1.0e-10

        ### state-change vector vj
        rates = np.array([growthrate, deathrate], dtype=object)

        # for determining of next reaction
        a = np.zeros(2)
        a[0] = rates[0] * N_population_array[-1]
        a[1] = rates[1] * N_population_array[-1]
        a0 = sum(a)
        return rates, a, a0

    @staticmethod
    def evolve_gillespie_one_step(a, a0, time_array, AB_conc_array, N_population_array, AB_Deg_rate, stoichiometry):
        # determine time of next reaction for gillespie
        r1 = random.uniform(0, 1)
        tau = (-math.log(r1) / a0)
        time_array = np.append(time_array, time_array[-1] + tau)
        # print('tau', tau)

        # determine nature of next reaction
        r2 = random.uniform(0, 1)
        # print('r2', r2)
        # print(r2)
        acumsum = np.divide(np.cumsum(a), a0)
        # print('acumsum', acumsum)
        # print(acumsum)
        chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
        # print(chosen_reaction)

        # update AB-conc to account for degradation by the cells

        # AB is degraded at rate Vmax proportional to number of resistant cells:
        AB_conc_array = np.append(AB_conc_array,
                                  AB_conc_array[-1] * math.exp(-AB_Deg_rate * N_population_array[-1] * tau))

        # update N array
        # print('last N before append', self.N[-1])
        # print('N + stoich', self.N[-1] + stoichiometry[chosen_reaction])

        N_population_array = np.append(N_population_array, N_population_array[-1] + stoichiometry[chosen_reaction])
        # print('dN', stoichiometry[chosen_reaction])
        # print('N', self.N[-1])
        # keep counters for timesteps and AB_conc
        return time_array, AB_conc_array, N_population_array

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
    def precheck_tau_binary_g(MIC, rates, trial_tau, epsilon, N_t, AB_conc_t, v):
        # stoichiometry vector
        stoichiometry = np.array([1, -1])

        # for determining of next reaction
        a = np.zeros(2)
        a[0] = rates[0] * N_t
        a[1] = rates[1] * N_t
        a0 = sum(a)
        # print('a0', a0)

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
            # print(AB_conc)
            # gradient of change in propensity function
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
    def precheck_tau_balanced_g(MIC, rates, trial_tau, epsilon, N_t, AB_conc_t, v, stoichiometry):

        # for determining of next reaction
        a = np.zeros(2)
        a[0] = rates[0] * N_t
        a[1] = rates[1] * N_t
        a0 = sum(a)
        # print('a0', a0)

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
            # print(AB_conc)
            # gradient of change in propensity function
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

    @staticmethod
    def degrade_ab_1_step(AB_concentration_array, N_bact_pop, delta_t):
        if degradation == 'exponential':
           ...
        if degradation == 'MM_exponential':
            print('ab before deg', AB_concentration_array[-1])
            y = 1 / Km * AB_concentration_array[-1]
            exp_factor_1 = math.exp(AB_concentration_array[-1] / Km)
            exp_factor_2 = math.exp(-Vmax * N_bact_pop[-1] * delta_t / (volume * Km))
            new_ab = Km * lambertw(y * exp_factor_1 * exp_factor_2).real
            AB_concentration_array = np.append(AB_concentration_array, new_ab)
        return AB_concentration_array

    def binary_grow(self, AB_conc):
        #print(scientific_notation(self.MIC), scientific_notation(AB_conc))
        # print(self.MIC, self.AB_conc)
        #print(type(self.N))
        if self.MIC > AB_conc:
            self.N = self.N + self.N * self.growthrate * self.dt  ## N grows uninhibited, normalized by dt
        else:
            self.N = self.N - self.N * self.deathrate * self.dt   ## N dies with growthrate normalized by dt

        if self.N < 0.0:  #### #can't have less than 0 bacteria (THIS WAS 1) 01/04/2022
           self.N = 0#

    def precheck_tau_binary_grow(self, epsilon, AB_conc):

        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values

            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])

            s = 0
            tau_list = []

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):

                if self.AB_conc_array[-1] > self.MIC:
                    self.deathrate = 0.045
                    self.growthrate = 0
                else:
                    self.deathrate = 0
                    self.growthrate = 0.01

                rates = np.array([self.growthrate, self.deathrate], dtype=object)
                # print(rates)
                # print('last N', self.N[-1])

                # for determining of next reaction
                a = np.zeros(2)
                a[0] = rates[0] * self.N[-1]
                a[1] = rates[1] * self.N[-1]
                a0 = sum(a)
                # print('a0', a0)

                ## check v
                if self.MIC > self.AB_conc_array[-1]:  # GROWTH
                    v = np.random.poisson(self.N[-1] * self.growthrate * self.dt, 1)  # second Poisson method
                    # gradient of change in propensity function
                    check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, self.N[-1], self.AB_conc_array[-1], v)
                    while check_returned == 0:
                        self.dt /= 2.0
                        v = np.random.poisson(self.N[-1] * self.growthrate * self.dt, 1)
                        check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, self.N[-1], self.AB_conc_array[-1], v)

                    if self.dt > (2 / a0):
                        # update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + self.dt)

                        # update Ab_conc for next step
                        self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                            -AB_Deg_rate * self.N[-1] * self.dt))
                        # update N
                        self.N = np.append(self.N, self.N[-1] + v)

                        s += 1
                        gill_trial = False

                    else:
                        gill_trial = True

                else:  # DEATH
                    # print(AB_conc)
                    v = np.random.poisson(self.N[-1] * self.deathrate * self.dt, 1)  # second Poisson method
                    check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, self.N[-1], self.AB_conc_array[-1], v)
                    while check_returned == 0:
                        self.dt /= 2.0
                        v = np.random.poisson(self.N[-1] * self.deathrate * self.dt, 1)
                        check_returned = self.precheck_tau_binary_g(self.MIC, rates, self.dt, epsilon, self.N[-1], self.AB_conc_array[-1], v)

                    if self.dt > (2 / a0):
                        # update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + self.dt)

                        # update Ab_conc for next step
                        self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                            -AB_Deg_rate * self.N[-1] * self.dt))
                        # update N
                        self.N = np.append(self.N, self.N[-1] - v)
                        s += 1
                        gill_trial = False
                    else:
                        gill_trial = True

                if gill_trial == True:
                    # determine time of next reaction for gillespie
                    r1 = random.uniform(0, 1)
                    self.dt = (-math.log(r1) / a0)
                    tau_list.append(self.dt)
                    self.t_array = np.append(self.t_array, self.t_array[-1] + self.dt)
                    # print('tau', tau)

                    # determine nature of next reaction
                    r2 = random.uniform(0, 1)
                    # print('r2', r2)
                    # print(r2)
                    acumsum = np.divide(np.cumsum(a), a0)
                    # print('acumsum', acumsum)
                    # print(acumsum)
                    chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
                    # print(chosen_reaction)

                    # update AB-conc to account for degradation by the cells

                    # AB is degraded at rate Vmax proportional to number of resistant cells:
                    self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                        -AB_Deg_rate * self.N[-1] * self.dt))

                    # update N array
                    # print('last N before append', self.N[-1])
                    # print('N + stoich', self.N[-1] + stoichiometry[chosen_reaction])

                    self.N = np.append(self.N, self.N[-1] + stoichiometry[chosen_reaction])
                    # print('dN', stoichiometry[chosen_reaction])
                    # print('N', self.N[-1])
                    # keep counters for timesteps and AB_conc
                    s += 1
                # print(mean(tau_list))
                # print('Time array:',self.t_array)
                # print('N array:', self.N)
            return self.t_array, self.N, self.AB_conc_array

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

                if self.AB_conc_array[-1] > self.MIC:
                    self.deathrate = 0.045
                    self.growthrate = 0
                else:
                    self.deathrate = 0
                    self.growthrate = 0.01

                rates = np.array([self.growthrate, self.deathrate], dtype=object)
                #print(rates)
                # print('last N', self.N[-1])

                #for determining of next reaction
                a = np.zeros(2)
                a[0] = rates[0] * self.N[-1]
                a[1] = rates[1] * self.N[-1]
                a0 = sum(a)
                #print('a[0]', a[0])
                #print('a[1]', a[1])
                # expectation value
                ksi = a[0] * stoichiometry[0] + a[1] * stoichiometry[1]
                #print(ksi)
                #print(self.dt)
                exp_val = ksi * self.dt
                #print('exp_val', exp_val)
                #print('last N', self.N[-1])
                N_prime = self.N[-1] + exp_val / 2
                #print('N_prime', N_prime)

            ## check v
                if self.MIC > self.AB_conc_array[-1]:  # GROWTH

                    v = np.random.poisson(N_prime * self.growthrate * self.dt, 1)  # second Poisson method

                    change_in_N = v * stoichiometry[0]
                    #print('change in N', change_in_N)

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
                        self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                            -AB_Deg_rate * self.N[-1] * tau))
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
                        self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                            -AB_Deg_rate * self.N[-1] * tau))
                        #update N
                        self.N = np.append(self.N, self.N[-1] - change_in_N)
                        #print('actual N', self.N[-1])
                        s += 1
                        gill_trial = False
                    else:
                        gill_trial = True

                if gill_trial == True:

                    # determine time of next reaction for gillespie
                    r1 = random.uniform(0, 1)
                    tau = (-math.log(r1) / a0)
                    tau_list.append(self.dt)
                    self.t_array = np.append(self.t_array, self.t_array[-1] + tau)
                    #print('tau', tau)


                    # determine nature of next reaction
                    r2 = random.uniform(0, 1)
                    #print('r2', r2)
                    #print(r2)
                    acumsum = np.divide(np.cumsum(a), a0)
                    #print('acumsum', acumsum)
                    #print(acumsum)
                    chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
                    #print(chosen_reaction)

                    # update AB-conc to account for degradation by the cells

                    # AB is degraded at rate Vmax proportional to number of resistant cells:
                    self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(-AB_Deg_rate * self.N[-1] * tau))

                    # update N array
                    #print('last N before append', self.N[-1])
                    #print('N + stoich', self.N[-1] + stoichiometry[chosen_reaction])

                    self.N = np.append(self.N, self.N[-1] + stoichiometry[chosen_reaction])
                    #print('dN', stoichiometry[chosen_reaction])
                    #print('N', self.N[-1])
                    # keep counters for timesteps and AB_conc
                    s += 1
                #print(mean(tau_list))
                #print('Time array:',self.t_array)
                #print('N array:', self.N)
            return self.t_array, self.N, self.AB_conc_array

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

                #a_history.append([a[0], a[1]])

                ## compute Lj the maximum nr of permitted firings

                L, is_crit, index_noncrit = compute_L_and_crit_arrays(stoichiometry, a, n_crit)

                #here are no noncritical reactions
                if not any(is_crit) == False:
                    tau = float('inf')
                else:
                    #f_array = np.zeros((len(rates), len(index_noncrit)))
                    #mean_array = np.zeros(len(rates))
                    #variance_array = np.zeros(len(rates))
                    mean_array = 0
                    variance_array = 0
                    #tau_less_ea0 = np.zeros(len(rates),2)
                    tau_less_ea0 = np.zeros(2)
                    g = 1
                    for j in range(len(rates)):
                        for j_prime in index_noncrit:
                            #change_in_a = a_history[-1][j_prime] - a_history[-2][j_prime]
                            #f_array[j,j_prime] = change_in_a * stoichiometry[j_prime]
                            #Sec II C algorithm from Cao 2006 paper
                            #mean_array[j] = mean_array[j] + f_array[j,j_prime] * a[j_prime]
                            #variance_array[j] = variance_array[j] +  f_array[j,j_prime]**2 * a[j_prime]
                            #Sec III algorithm from Cao paper where mean depends on species that participate in
                            # reaction (i,e, aj = cj* x -> then look at all the species that are present)
                            mean_array += stoichiometry[j_prime] * a[j_prime]
                            variance_array += stoichiometry[j_prime] ** 2 * a[j_prime]
                        # Sec II C algorithm from Cao 2006 paper
                        #tau_less_ea0[j][0] = max(epsilon * a[j], rates[j]) / abs(mean_array[j])
                        #tau_less_ea0[j][1] = max(epsilon * a[j], rates[j])**2 / variance_array[j]**2
                        # Sec III algorithm from Cao paper where mean depends on species that participate in
                        # reaction (i,e, aj = cj* x -> then look at all the species that are present)
                        tau_less_ea0[0] = max(epsilon * self.N[-1]/g, 1) / abs(mean_array)
                        tau_less_ea0[1] = max(epsilon * self.N[-1]/g, 1) ** 2 / variance_array
                        tau = min(tau_less_ea0)
                    #print('tau', tau)
                    if tau < (10/a0):
                        gill_trial = True
                    else:
                        gill_trial = False
                    #print('gill_trial', gill_trial)

                    if gill_trial == False:
                        K = self.set_K_values(is_crit, index_noncrit, a, stoichiometry, tau)

                        if self.N[-1] + K[0] * stoichiometry[0] >= 0 or self.N[-1] + K[1]*stoichiometry[1] >= 0:
                            reduce_t = False
                        else:
                            reduce_t = True
                        red = 0
                        while (reduce_t == True):
                            #print('red',red)
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

                                #update time array
                                self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                                #update Ab_conc for next step
                                self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                                    -AB_Deg_rate * self.N[-1] * tau))
                                #update N
                                self.N = np.append(self.N, self.N[-1] + K[0] * stoichiometry[0])
                                s += 1

                        else:  # DEATH
                                #update time array
                                self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                                #update Ab_conc for next step
                                self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                                    -AB_Deg_rate * self.N[-1] * tau))
                                #update N
                                self.N = np.append(self.N, self.N[-1] + K[1] * stoichiometry[1])
                                #print('actual N', self.N[-1])
                                s += 1

                    elif gill_trial == True:
                            for _ in range(20):
                                while (self.t_array[s] < self.t_end) and (self.N[-1] > 0) and (self.N[-1] < self.Nsat):
                                    rates, a, a0 = self.set_rates_and_prop_binary(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)

                                    self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array,
                                                                                                          self.AB_conc_array, self.N, AB_Deg_rate, stoichiometry)
                                    s += 1
                                    #print('s', s)
                                    #print('time', self.t_array)
                                    #print('N', self.N)
                                    #print('AB', self.AB_conc_array)
                                break
        #print('time', self.t_array)
        #print(len(self.t_array))
        return self.t_array, self.N, self.AB_conc_array

    def gillespie_binary_grow(self, AB_conc):
        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values
            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])

            s = 0
            tau_list = []

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):

                # setting the division, death rates
                if self.AB_conc_array[-1] > self.MIC:
                    self.deathrate = 0.045
                    self.growthrate = 0
                else:
                    self.deathrate = 0
                    self.growthrate = 0.01
                rates = np.array([self.growthrate, self.deathrate], dtype=object)

                #propensities function
                a = np.zeros(2)
                a[0] = rates[0] * self.N[-1]
                a[1] = rates[1] * self.N[-1]
                a0 = sum(a)

                # determine time of next reaction
                r1 = random.uniform(0, 1)
                tau = (-math.log(r1) / a0)
                tau_list.append(tau)

                # determine nature of next reaction
                r2 = random.uniform(0, 1)
                acumsum = np.divide(np.cumsum(a), a0)
                chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
                print(chosen_reaction)

                #update time array
                self.t_array = np.append(self.t_array, self.t_array[-1] + tau)
                # AB is degraded at rate Vmax proportional to number of resistant cells:
                self.AB_conc_array = self.degrade_ab_1_step(self.AB_conc_array, self.N, tau)
                # update N array
                self.N = np.append(self.N, self.N[-1] + stoichiometry[chosen_reaction])
                # keep counters for timesteps and AB_conc
                s += 1

            return self.t_array, self.N, self.AB_conc_array


    def balanced_grow(self, AB_conc):
        if self.N != 0:  # if there is bacteria to grow...
            #print(self.N)
            #### Returns: out ndarray or scalar after random.poisson
            #print(type(self.N))

            self.deathrate = AB_conc * 1.1e-10 #/self.N

            self.N = self.N + self.N * self.growthrate * self.dt - self.N * self.deathrate * self.dt  ## N grows as a balance
            #print(self.growthrate, scientific_notation(self.deathrate))
            # between growth and death
            x = self.deathrate/self.growthrate
            # print("normalized death rate x".format(x))

            # ## compute probability of extinction
            #
            # if x == 1:
            #     self.prob_ext = (Nsat - self.N) / Nsat
            #     # print("extinction probability: {}".format(self.prob_ext))
            # else:
            #     self.prob_ext = (pow(x, Nsat) - pow(x, self.N)) / (pow(x, Nsat) - 1)
            #     # print("extinction probability: {}".format(self.prob_ext))
            # if math.isnan(self.prob_ext):
            #     self.prob_ext = 1 - pow(x, self.N - Nsat)
            #     # print("extinction probability: {}".format(self.prob_ext))

            # ### compute extinction time
            # ### compute extinction time if x = 1
            # if x == 1:
            #     ssum1 = 0
            #     ssum2 = 0
            #     for k in range(1, Nsat - 1):
            #         m = pow(Nsat - k, 2) / k
            #         ssum1 += m
            #         # print(k, m, ssum1)
            #     for j in range(1, self.N - 1):
            #         o = (Nsat - j) * (self.N - j) / j
            #         ssum2 += o
            #         # print(j,o, ssum2)
            #     self.time_ext = 1 / (self.growthrate * (Nsat - self.N)) * ((self.N / Nsat * ssum1) - ssum2)
            #     print("extinction time if x is 1: {}".format(self.time_ext))
            # # elif x == 0:
            # #     self.time_ext = 1 / self.growthrate * (1 / self.N + x / (self.N + 1))
            # #     print(self.time_ext)
            # else:
            #     ssum1 = 0
            #     ssum2 = 0
            #     for k in range(1, Nsat - 1):
            #         m = (x ** Nsat - x ** k) * (x ** (Nsat - k) - 1) / k
            #         ssum1 += m
            #         # print(k, m, ssum1)
            #     for j in range(1, int(self.N) - 1):
            #         o = (x ** Nsat - x ** j) * (x ** (self.N - j) - 1) / j
            #         ssum2 += o
            #         # print(j,o, ssum2)
            #     self.time_ext = 1 / (self.growthrate * (x ** Nsat - x ** self.N) * (x - 1)) * \
            #                     (((1 - x ** self.N) / (1 - x ** Nsat)) * ssum1 - ssum2)
            #     print(self.time_ext)
            #
            # # if self.time_ext == float("inf"):
            # #     ssum1 = 0
            # #     ssum2 = 0
            # #     for k in range(1, self.N - 1):
            # #         for j in range(0, self.N - k - 1):
            # #             o = x ** j
            # #             ssum2 = ssum2 + o
            # #             # print(j,o, ssum2)
            # #         m = (1 / k) * ssum2
            # #         ssum1 = ssum1 + m
            # #         # print(k, m, ssum1)
            # #     self.time_ext = 1 / self.growthrate * ((x ** self.N - 1) / (x - 1) * math.log(x / (x - 1)) - ssum1)
            # #     print(self.time_ext)

        if self.N < 0.0:  #### #can't have less than 0 bacteria
            self.N = 0

    def midpoint_tau_balanced_grow(self, epsilon, AB_conc): ### not finished, need to see how to do the tau leap evolution of the system
        def precheck_tau(trial_tau, epsilon, stoichiometry, old_prop):
            rates, a, a0 = self.set_rates_and_prop_balanced(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)

            # expectation value
            ksi = a[0] * stoichiometry[0] + a[1] * stoichiometry[1]

            # gradient of change in propensity function
            b_growth = (a[0] - old_prop[-1]) / trial_tau
            b_death = (a[1] - old_prop[-1]) / trial_tau

            lhs_growth = trial_tau * abs(ksi * b_growth)
            lhs_death = trial_tau *  abs(ksi * b_death)

            ## check dt
            rhs = epsilon * a0

            if lhs_growth <= rhs and lhs_death <= rhs:
                return 1
            else:
                return 0

        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values

            self.t_array = np.array([0])
            self.N = np.array([self.N])
            self.AB_conc_array = np.array([AB_conc])
            #print(self.dt)
            s = 0
            tau_list = []
            a_history = []

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):

                rates, a, a0 = self.set_rates_and_prop_balanced(self.AB_conc_array, self.MIC, self.N, self.growthrate, self.deathrate)

                a_history.append([a[0], a[1]])
                #print('a[0]', a[0])
                #print('a[1]', a[1])
                # expectation value
                ksi = a[0] * stoichiometry[0] + a[1] * stoichiometry[1]
                #print(ksi)
                #print(self.dt)
                exp_val = ksi * self.dt
                #print('exp_val', exp_val)
                #print('last N', self.N[-1])
                N_prime = self.N[-1] + exp_val / 2
                #print('N_prime', N_prime)

            ## check v with precheck
                K_growth = np.random.poisson(N_prime * self.dt * self.growthrate, 1)
                K_death = np.random.poisson(N_prime * self.dt * self.deathrate, 1)

                change_in_N = K_growth * stoichiometry[0] + K_death * stoichiometry[1]
                #print('change in N', change_in_N)

                # gradient of change in propensity function
                check_returned = precheck_tau(self.dt, epsilon, stoichiometry, a_history)
                tau = self.dt
                while check_returned == 0:
                    tau /= 2.0
                    exp_val = ksi * tau
                    N_prime = self.N[-1] + exp_val / 2
                    #print('N_prime', N_prime)

                    K_growth = np.random.poisson(N_prime * self.dt * self.growthrate, 1)
                    K_death = np.random.poisson(N_prime * self.dt * self.deathrate, 1)

                    change_in_N = K_growth * stoichiometry[0] + K_death * stoichiometry[1]
                    #print('change in N', change_in_N)

                    check_returned = precheck_tau(self.dt, epsilon, stoichiometry, a_history)

                if self.dt > (2/a0):
                    #update time array
                    self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                    #update Ab_conc for next step
                    self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                        -AB_Deg_rate * self.N[-1] * tau))
                    #update N
                    self.N = np.append(self.N, self.N[-1] + change_in_N)
                    #print('actual N', self.N[-1])
                    s += 1
                    gill_trial = False
                else:
                    gill_trial = True

            ## if tau is less than 2/a0  -> evolve system with Gillespie
                if gill_trial == True:

                    # determine time of next reaction for gillespie
                    r1 = random.uniform(0, 1)
                    tau = (-math.log(r1) / a0)
                    tau_list.append(self.dt)
                    self.t_array = np.append(self.t_array, self.t_array[-1] + tau)
                    #print('tau', tau)

                    # determine nature of next reaction
                    r2 = random.uniform(0, 1)
                    #print('r2', r2)
                    #print(r2)
                    acumsum = np.divide(np.cumsum(a), a0)
                    #print('acumsum', acumsum)
                    #print(acumsum)
                    chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
                    #print(chosen_reaction)

                    # update AB-conc to account for degradation by the cells

                    # AB is degraded at rate Vmax proportional to number of resistant cells:
                    self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(-AB_Deg_rate * self.N[-1] * tau))

                    # update N array
                    #print('last N before append', self.N[-1])
                    #print('N + stoich', self.N[-1] + stoichiometry[chosen_reaction])

                    self.N = np.append(self.N, self.N[-1] + stoichiometry[chosen_reaction])
                    #print('dN', stoichiometry[chosen_reaction])
                    #print('N', self.N[-1])
                    # keep counters for timesteps and AB_conc
                    s += 1
                #print(mean(tau_list))
                #print('Time array:',self.t_array)
                #print('N array:', self.N)
            return self.t_array, self.N, self.AB_conc_array

    def adaptive_tau_balanced_grow(self, epsilon, AB_conc):
        def set_rates_and_prop_balanced(N_population_array, growthrate):
            self.deathrate = self.AB_conc_array[-1] * 1.0e-10

            ### state-change vector vj

            rates = np.array([growthrate, deathrate], dtype=object)

            # for determining of next reaction
            a = np.zeros(2)
            a[0] = rates[0] * N_population_array[-1]
            a[1] = rates[1] * N_population_array[-1]
            a0 = sum(a)
            return rates, a, a0

        def compute_L_and_crit_arrays(stoichiometry, a, n_crit):
            L = np.zeros(2)
            L[0] = math.floor(self.N[-1]/abs(stoichiometry[0]))
            L[1] = math.floor(self.N[-1]/abs(stoichiometry[1]))

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
                rates, a, a0 = set_rates_and_prop_balanced(self.N, self.growthrate)

                #a_history.append([a[0], a[1]])
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


                        #update time array
                        self.t_array = np.append(self.t_array, self.t_array[-1] + tau)

                        #update Ab_conc for next step
                        self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(
                            -AB_Deg_rate * self.N[-1] * tau))
                        #update N
                        self.N = np.append(self.N, self.N[-1] + K[0] * stoichiometry[0] + K[1] * stoichiometry[1])
                        s += 1


                    elif gill_trial == True:
                            for _ in range(20):
                                while (self.t_array[s] < self.t_end) and (self.N[-1] > 0) and (self.N[-1] < self.Nsat):
                                    rates, a, a0 = set_rates_and_prop_balanced( self.N, self.growthrate)

                                    self.t_array, self.AB_conc_array, self.N = self.evolve_gillespie_one_step(a, a0, self.t_array,
                                                                                                          self.AB_conc_array, self.N, AB_Deg_rate, stoichiometry)
                                    s += 1
                                    #print('s', s)
                                    #print('time', self.t_array)
                                    #print('N', self.N)
                                    #print('AB', self.AB_conc_array)
                                break
        return self.t_array, self.N, self.AB_conc_array

    def gillespie_balanced_grow(self, AB_conc):
        if self.initialN != 0:

            # stoichiometry vector
            stoichiometry = np.array([1, -1])

            # set up time and initial values

            self.t_array = np.array([0])
            self.N = np.array([self.N])
            # print("self.N", self.N)
            self.AB_conc_array = np.array([AB_conc])

            s = 0
            tau_list = []

            while (self.t_array[s] < self.t_end) and (self.N[-1] != 0) and (self.N[-1] < self.Nsat):

                # propensities (growth rates)
                #print('AB',self.AB_conc_array[-1])
                self.deathrate = self.AB_conc_array[-1] * 1.0e-10

                ### For changing growthrate, deathrate
                # growthrate_bin = abs(np.random.normal(self.growthrate, 0.1))
                # deathrate_bin = abs(np.random.normal(self.deathrate, 0.1))
                # rates = np.array([growthrate_bin, deathrate_bin], dtype=object)

                rates = np.array([self.growthrate, self.deathrate], dtype=object)
                #print(rates)
                # print('last N', self.N[-1])
                #for determining of next reaction
                a = np.zeros(2)
                a[0] = rates[0] * self.N[-1]
                a[1] = rates[1] * self.N[-1]
                a0 = sum(a)
                #print('a0', a0)



                # determine time of next reaction
                r1 = random.uniform(0, 1)
                tau = (-math.log(r1) / a0)
                tau_list.append(tau)
                self.t_array = np.append(self.t_array, self.t_array[-1] + tau)
                #print('tau', tau)


                # determine nature of next reaction
                r2 = random.uniform(0, 1)
                #print('r2', r2)
                #print(r2)
                acumsum = np.divide(np.cumsum(a), a0)
                #print('acumsum', acumsum)
                #print(acumsum)
                chosen_reaction = np.argwhere(np.array(acumsum) >= r2)[0]
                #print(chosen_reaction)

                # update AB-conc to account for degradation by the cells
                # AB is degraded at rate Vmax proportional to number of resistant cells:
                self.AB_conc_array = np.append(self.AB_conc_array, self.AB_conc_array[-1] * math.exp(-AB_Deg_rate * self.N[-1] * tau))


                # update N array
                #print('last N before append', self.N[-1])
                #print('N + stoich', self.N[-1] + stoichiometry[chosen_reaction])

                self.N = np.append(self.N, self.N[-1] + stoichiometry[chosen_reaction])
                #print('dN', stoichiometry[chosen_reaction])
                #print('N', self.N[-1])
                # keep counters for timesteps and AB_concs
                s += 1
            #print(mean(tau_list))
            #print('Time array:', self.t_array)
            #print('N array:', self.N)
            return self.t_array, self.N, self.AB_conc_array





strain_script_path = __file__
strain_script_name = os.path.basename(__file__)
