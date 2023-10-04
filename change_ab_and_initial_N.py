from DropTest_fct import DropTest
import variables
import os
import numpy as np
import importlib
variables.initialN = 3000

START_AB=0
END_AB= 61
STEP_AB=1

START_N=1
END_N=10
STEP_N=1

path = os.getcwd()
results_path = path + '/output/'
print(results_path)

ab_list=np.arange(START_AB, END_AB+1, STEP_AB)
print(ab_list)
init_N_list=np.arange(START_N, END_N+1, START_N)
print(init_N_list)


for Ni in init_N_list:
    print(Ni)
    variables.initialN = Ni
    for ab in ab_list:
        print(ab)
        variables.AB_conc = ab

        simulate = DropTest()
        simulate.calc_survival_prob_total_nr_bact_diff_part(variables.step, variables.spec_time, variables.total_sim)
        os.chdir('..')
        os.chdir('..')
        print(os.getcwd())

print('All processes complete')
variables.AB_conc = 15
variables.initialN = 5