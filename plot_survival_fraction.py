import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt


rootdir = 'output/'
ab = [10, 25, 30, 35, 40, 55, 70]

os.chdir(rootdir)
print(os.getcwd())


plt.figure(1)

for i in ab:
    os.chdir('dropnr_1000_loading_rand_growth_gillespie_binary_initialN_10_abconc_{}'.format(i))
    print(os.getcwd())
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    print(onlyfiles)

    surv_fraction = pd.read_csv(onlyfiles[3])
    print(surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])
    print(part_fact)

    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'
    print(surv_fraction_transpose)


    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)


# plt.grid(True)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_transpose["Surv frac"].plot.line(yerr = surv_fraction_errors)#, color = 'orange')

plt.title('Fraction of droplets surviving a round of treatment')
plt.ylabel('Fraction of droplets surviving')
plt.xlabel('Partitioning factor')
plt.legend(ab, title='Antibiotic conc',  loc='upper right')
plt.savefig('./output/Survival fraction Gillespie + errors.png')
plt.show()