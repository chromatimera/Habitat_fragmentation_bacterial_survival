import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt

growth1 = 'binary'
growth2 = 'gillespie_binary'
rootdir = 'output/'
ab = [35, 45]

os.chdir(rootdir)
print('current dir', os.getcwd())


plt.figure(1)
label_list = []
for antib in ab:
    ### plot the deterministic value + error
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth1, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)

    surv_fraction1 = pd.read_csv(onlyfiles[3])
    part_fact1 = np.loadtxt(onlyfiles[2])

    os.chdir('..')

    ### plot the stochastic value + error
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth2, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)

    surv_fraction2 = pd.read_csv(onlyfiles[3])
    part_fact2 = np.loadtxt(onlyfiles[2])

    os.chdir('..')

    ### transpose of dataframe
    surv_fraction_transpose1 = surv_fraction1.T
    surv_fraction_transpose1.index.name = 'Part_fact'

    surv_fraction_transpose2 = surv_fraction2.T
    surv_fraction_transpose2.index.name = 'Part_fact'

    surv_fraction_transpose1.columns = ['Surv frac']
    surv_fraction_transpose1['Error95'] = surv_fraction_transpose1.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose1['Error99'] = surv_fraction_transpose1.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)

    surv_fraction_transpose2.columns = ['Surv frac']
    surv_fraction_transpose2['Error95'] = surv_fraction_transpose2.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose2['Error99'] = surv_fraction_transpose2.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)

    # plt.grid(True)
    surv_fraction_errors1 = surv_fraction_transpose1.Error95.to_frame('Surv frac')
    # reorder m - nr of subvolumes index in both dataframes - surv frac df + error df
    surv_fraction_errors1.index = surv_fraction_errors1.index.map(int)
    surv_fraction_errors1 = surv_fraction_errors1.sort_index(ascending=True)

    surv_fraction_transpose1.index = surv_fraction_transpose1.index.map(int)
    surv_fraction_transpose1 = surv_fraction_transpose1.sort_index(ascending=True)
    surv_fraction_transpose1["Surv frac"].plot.line(yerr=surv_fraction_errors1)  # , color = 'orange')
    label_list.append('{} deterministic'.format(antib))

    # plt.grid(True)
    surv_fraction_errors2 = surv_fraction_transpose2.Error95.to_frame('Surv frac')
    surv_fraction_errors2.index = surv_fraction_errors2.index.map(int)
    surv_fraction_errors2 = surv_fraction_errors2.sort_index(ascending=True)

    surv_fraction_transpose2.index = surv_fraction_transpose2.index.map(int)
    surv_fraction_transpose2 = surv_fraction_transpose2.sort_index(ascending=True)
    surv_fraction_transpose2["Surv frac"].plot.line(yerr=surv_fraction_errors2, linestyle='dashed')  # , color = 'orange')
    label_list.append('{} stochastic'.format(antib))


plt.title('Fraction of repeats with any bacteria surviving', fontsize=text_size)
plt.ylabel('Fraction of droplets surviving', fontsize=text_size)
plt.ylim(-0.1,1.1)
plt.xlabel('m (number of subvolumes)', fontsize=text_size)
plt.legend(label_list, title='Antibiotic conc and type of growth', loc='upper right')
plt.savefig('Survival fraction {} + {} + errors diff ab.png'.format(growth1, growth2))
plt.show()