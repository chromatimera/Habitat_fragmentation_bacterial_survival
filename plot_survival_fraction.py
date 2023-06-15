import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

growth = 'binary'
rootdir = 'output/'
ab = [15, 25, 45]

zz=np.load('prob_line.npy')
zzz=np.load('prob_line_finitesum.npy')
os.chdir(rootdir)
#print('current dir', os.getcwd())


plt.figure(1)
color = iter(cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
for antib, c in zip(ab, color):
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)

    surv_fraction = pd.read_csv(onlyfiles[3])
    #print('surf fraction df', surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])
    theory = pd.read_csv(onlyfiles[4], index_col='0')
    print(theory)
    #theory.set_index('0')
    theory.index.name = 'Part_fact'
    theory.columns = ['small_ps', 'big_Ps']

    theory = theory.sort_values(by = "Part_fact",ascending = False)
    print(theory)

    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'
    #print('transpose', surv_fraction_transpose)

    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    print(surv_fraction_transpose)

    plt.figure(1)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_errors.index = surv_fraction_errors.index.map(int)
    surv_fraction_errors = surv_fraction_errors.sort_index(ascending=False)
    #print(surv_fraction_errors)
    surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)
    surv_fraction_transpose = surv_fraction_transpose.sort_index(ascending=True)
    surv_fraction_transpose["Surv frac"].plot.line(yerr = surv_fraction_errors, c=c)#, color = 'orange')
    theory["big_Ps"].plot.line(c=c, linestyle='dashed')#, color = 'orange')


    os.chdir('..')
    #print(os.getcwd())

#plt.title('Fraction of repeats with any bacteria surviving', fontsize=text_size)
plt.ylabel('Survival fraction', fontsize=text_size)
plt.xlabel('m (number of subvolumes)', fontsize=text_size)
plt.legend(ab, title='Antibiotic conc', loc='upper right')
plt.savefig('Survival fraction {} + errors diff ab.png'.format(growth))
plt.show()