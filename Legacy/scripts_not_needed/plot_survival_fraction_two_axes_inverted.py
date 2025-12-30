import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from collections import Counter

from ps_theory import vol_fac

BIGGER_SIZE = 32

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

#rootdir = './output/'
ab = [15, 35, 55, 75]
rho = 5e7 ## this is the initial density of cells/mL; for sim starting with lamda = 5; change accordingly

zz=np.load('../prob_line.npy')
#os.chdir(rootdir)
zzz= zz.T
#os.chdir(rootdir)

print('current dir', os.getcwd())


fig = plt.figure(figsize=(11,11))
ax1 = fig.add_subplot(111)
color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
label_list = []
print(color)

for antib, c, ind in zip(ab, color, range(len(ab))):
    print('ab conc', antib)
    print(c)

    if ind == 2:
        c = next(color)
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    #print(onlyfiles)

    if '.DS_Store' in onlyfiles:
        onlyfiles.remove('.DS_Store')
    else:
        pass

    surv_fraction = pd.read_csv(onlyfiles[3])
    #print('surf fraction df', surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])

    theory_line_df = pd.DataFrame(zzz[:, ind], columns=['big_Ps'], index=vol_fac)
    theory_line_df.index.name = 'Vol_fac'
    theory_line_df = theory_line_df.sort_values(by="Vol_fac", ascending=True)
    theory_line_df['M'] = theory_line_df.index.astype(int)
    theory_line_df['RhoV'] = theory_line_df.apply(lambda x: rho * 1e-4 / x['M'], axis=1)

    plt.figure(1)

    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'

    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['M'] = surv_fraction_transpose.index.astype(int)
    surv_fraction_transpose['RhoV'] = surv_fraction_transpose.apply(lambda x: rho * 1e-4 / x['M'], axis=1)


    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_errors.index = surv_fraction_errors.index.map(int)
    surv_fraction_errors.index = surv_fraction_transpose['RhoV']
    #print('errors',surv_fraction_errors)

    surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)

    surv_fraction_transpose = surv_fraction_transpose.set_index('RhoV', drop=True)
    theory_line_df = theory_line_df.set_index('RhoV', drop=True)

    theory_line_df["big_Ps"].plot.line(ax=ax1, c=c, linestyle='dashed', label='_nolegend_', logx=True)#, color = 'orange'
    surv_fraction_transpose["Surv frac"].plot.line(ax=ax1, yerr=surv_fraction_errors, c=c, logx=True)#, color = 'orange')
    label_list.append('{}'.format(antib))

    os.chdir('../..')


plt.xlim(10**(0.6), 10**(2.3))

xticks = [5, 10, 25, 50, 100, 200]
xticks_label = [5, 10, 25, 50, 100, 5000]
second_ticks = [1000, 500, 200, 100, 50, 1]

ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks_label)

ax2 = ax1.twiny()
ax2.set_xticks(np.linspace(ax1.get_xticks()[0], ax1.get_xticks()[-1], len(ax1.get_xticks())))
ax2.set_xticklabels(second_ticks[::-1])

plt.xlabel(r'\bf{m (number of subvolumes)}')
ax1.set_xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
ax1.set_ylabel(r'\bf{Probability of survival}')

ax1.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.45), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)
plt.tight_layout()
plt.savefig('Survival fraction RhoV_ double axes'.format(growth), dpi=600)

plt.show()