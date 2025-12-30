from DropTest_fct import DropTest
import variables
import os
import numpy as np
import importlib
variables.initialN = 3000

START_AB=10
END_AB= 30
STEP_AB=5

path = os.getcwd()
results_path = path + '/output/'
print(results_path)

ab_list=np.arange(START_AB, END_AB+1, STEP_AB)
#print(ab_list)

for ab in ab_list:
   # print(ab)
    print(os.getcwd())
    variables.AB_conc = ab

    simulate = DropTest()
    simulate.calc_survival_prob_total_nr_bact_diff_part(variables.step, variables.spec_time, variables.total_sim)
    os.chdir('..')
    os.chdir('..')

print('All processes complete')
variables.AB_conc = 15
variables.initialN = 5