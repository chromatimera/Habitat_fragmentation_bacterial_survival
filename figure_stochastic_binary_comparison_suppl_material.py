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

BIGGER_SIZE = 22

#### FIGURE SUPPL MATERIAL COMPARISON BETWEEN RAND/DET LOADING + STOCH/DET GROWTH


plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

#rootdir = './output/'
ab = [35]
rho = 5e7 ## this is the initial density of cells/mL; for sim starting with lamda = 5; change accordingly
total_sim = 1000 # number of simulation repeats

zz=np.load('prob_line.npy')
#os.chdir(rootdir)
zzz= zz.T
#os.chdir(rootdir)

# markers = ['-', '--', ':', '-.']
linestyles = [ '--', ':','-', '-.'] #['dotted', 'dashed', 'dashdot', 'long dash with offset']
print('current dir', os.getcwd())


fig = plt.figure(figsize=(8,9))
ax1 = fig.add_subplot(111)
color = plt.cm.rainbow(np.linspace(0, 1,5))

label_list = []
ind = 0
for loading in (['det', 'rand']):
    for growth in (['binary', 'gillespie_binary']):

        os.chdir('dropnr_1000_loading_{}_growth_{}_initialN_5_abconc_35'.format(loading, growth))
        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)
        #print(onlyfiles)

        if '.DS_Store' in onlyfiles:
            onlyfiles.remove('.DS_Store')
        else:
            pass

        surv_fraction = pd.read_csv(onlyfiles[3])
        print(onlyfiles[3])
        #print('surf fraction df', surv_fraction)
        part_fact = np.loadtxt(onlyfiles[2])


        plt.figure(1)

        ### transpose of dataframe
        surv_fraction_transpose = surv_fraction.T
        surv_fraction_transpose.index.name = 'Part_fact'

        surv_fraction_transpose.columns = ['Surv frac']
        surv_fraction_transpose['M'] = surv_fraction_transpose.index.astype(int)
        surv_fraction_transpose['RhoV'] = surv_fraction_transpose.apply(lambda x: rho * 1e-4 / x['M'], axis=1)


        surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(total_sim), axis=1)
        surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(total_sim), axis=1)

        surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
        surv_fraction_errors.index = surv_fraction_errors.index.map(int)
        surv_fraction_errors.index = surv_fraction_transpose['RhoV']
        #print('errors',surv_fraction_errors)

        surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)

        surv_fraction_transpose = surv_fraction_transpose.set_index('RhoV', drop=True)

        surv_fraction_transpose["Surv frac"].plot.line(ax=ax1, yerr=surv_fraction_errors, c = color[0], linestyle = linestyles[ind], logx=True)#, linestyle=next(markers))#, color = 'orange')
        label_list.append('{}, {}'.format(loading, growth))
        print('{}, {}'.format(loading, growth))
        os.chdir('..')
        print(surv_fraction_transpose)
        ind +=1




ax1.set_xlim(10**(2.3), 10**(0.6))
#plt.grid(True)
xticks = [200, 100, 50, 25, 10, 5]
xticks_label = [200, 100, 50, 25, 10, 5]
second_ticks = [25, 50, 100, 200, 500, 1000]
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks_label)

ax2 = ax1.secondary_xaxis("top")
ax2.xaxis.set_ticks(xticks[::-1], labels=second_ticks[::-1])


plt.xlabel(r'\bf{m (number of subvolumes)}')
ax1.set_xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
ax1.set_ylabel(r'\bf{Probability of survival}')
ax2.set_xlabel(r'\bf{m (number of subvolumes)}')
plt.xlim([200,5])
ax1.legend(label_list, title=r'\bf{Type of simulation}', loc='upper center', bbox_to_anchor=(0.5, 1.45), ncol=2, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)
plt.tight_layout()
plt.savefig('Prob of survival comparison det to stoc'.format(growth), dpi=600)
plt.show()