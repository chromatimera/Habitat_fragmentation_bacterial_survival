import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt

from ps_theory import vol_fac

BIGGER_SIZE = 16

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

#rootdir = './output/'
ab = [15, 35, 55, 75]

zz=np.load('prob_line.npy')
#os.chdir(rootdir)
zzz= zz.T
#os.chdir(rootdir)

print('current dir', os.getcwd())


plt.figure(figsize=(7, 7))
color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
label_list = []
print(color)

for antib, c, ind in zip(ab, color, range(len(ab))):
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    print(onlyfiles)

    if '.DS_Store' in onlyfiles:
        onlyfiles.remove('.DS_Store')
    else:
        pass

    surv_fraction = pd.read_csv(onlyfiles[3])
    print('surf fraction df', surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])

    theory_line_df = pd.DataFrame(zzz[:, ind], columns=['big_Ps'], index=vol_fac)
    theory_line_df.index.name = 'Vol_fac'
    theory_line_df = theory_line_df.sort_values(by="Vol_fac", ascending=True)


    print('THEORY DF', theory_line_df)


    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'
    print('transpose', surv_fraction_transpose)

    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    print(surv_fraction_transpose)

    plt.figure(1)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_errors.index = surv_fraction_errors.index.map(int)
    #surv_fraction_errors = surv_fraction_errors.sort_index(ascending=True)
    print('errors',surv_fraction_errors)

    surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)
    #surv_fraction_transpose = surv_fraction_transpose.sort_index(ascending=True)
    print('trp', surv_fraction_transpose)
    theory_line_df["big_Ps"].plot.line(c=c, linestyle='dashed', label='_nolegend_')#, color = 'orange')
    surv_fraction_transpose["Surv frac"].plot.line(yerr=surv_fraction_errors, c=c)#, color = 'orange')
    label_list.append('{}'.format(antib))


    os.chdir('..')
    #print(os.getcwd())

plt.ylabel(r'\bf{Probability of survival}')
plt.xlabel(r'\bf{m (number of subvolumes)}')

plt.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)
plt.savefig('Survival fraction {} + errors diff ab+ legend.png'.format(growth))
plt.show()