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
from matplotlib.lines import Line2D
from collections import Counter
from matplotlib import rc

#### SCRIPT TO COMPARE VARIANCE BETWEEN GILLESPIE SIMULATIONS - NOT WORKING YET


BIGGER_SIZE = 22

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsiz

growth1 = 'binary'
growth2 = 'gillespie_binary'
loading1 = 'det'
loading2 = 'rand'
rootdir = 'output/'
ab = [35]
rho = 5e7

os.chdir(rootdir)
print('current dir', os.getcwd())

colors = plt.cm.rainbow(np.linspace(0, 1, 5))
color_list = []
label_list = []
ind = 0

#print(colors)
fig = plt.figure(figsize=(9, 10))
plt.subplots_adjust(hspace=0.46)
# loop through two types of plots
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

for loading in (['det','rand']):
    for growth in (['binary', 'gillespie_binary']):
        if loading == 'rand': ### if using the DET loading - only 10 repeats that are identical, so the total_sim is 10
            total_sim = 1000
        else:
            total_sim = 10

        print(loading, growth)
        os.chdir('dropnr_1000_loading_{}_growth_{}_initialN_5_abconc_35'.format(loading, growth))
        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)
        # print(onlyfiles)

        if '.DS_Store' in onlyfiles:
            onlyfiles.remove('.DS_Store')
        else:
            pass

        surv_fraction = pd.read_csv(onlyfiles[3])
        # print('surf fraction df', surv_fraction)

        part_fact1 = np.loadtxt(onlyfiles[2])
        df_total_mass1 = pd.read_csv(onlyfiles[1])
        m_list1 = [round(1 / x) for x in part_fact1]


        ### DO THE SURVIVAL FRACTION DFS
        surv_fraction_transpose = surv_fraction.T
        surv_fraction_transpose.index.name = 'Part_fact'

        surv_fraction_transpose.columns = ['Surv frac']
        surv_fraction_transpose['M'] = surv_fraction_transpose.index.astype(int)
        surv_fraction_transpose['RhoV'] = surv_fraction_transpose.apply(lambda x: rho * 1e-4 / x['M'], axis=1)

        surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(
            lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
        surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(
            lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)

        rhov = list(surv_fraction_transpose['RhoV'])

        ## DO THE N(T) DFS

        ##start building the average N(t) and SD(t) over simulations
        total_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))
        error_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))
        var_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))

        ##start building the average N(t) and SD(t) over simulations
        #total_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))
        #error_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))
        #var_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))

        ## for each partition factor, calculate the sum over the simulation of N(t)


        for i in range(0, len(part_fact1), 1):
            for j in range(0, total_sim, 1):
                k = i * total_sim + j
                #print(i,j,k)

                total_nr_bact1[:, i] += df_total_mass1.iloc[:, k]
                #total_nr_bact2[:, i] += df_total_mass2.iloc[:, k]

        ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
        nr_simu = np.array(total_sim)
        avg_nr_bact1 = np.divide(total_nr_bact1, nr_simu)
        #avg_nr_bact2 = np.divide(total_nr_bact2, nr_simu)

        for i in range(0, len(part_fact1), 1):
            for j in range(0, total_sim, 1):
                k = j + i * total_sim
                var_nr_bact1[:, i] += (df_total_mass1.iloc[:, k] - avg_nr_bact1[:, i]) ** 2
                #var_nr_bact2[:, i] += (df_total_mass2.iloc[:, k] - avg_nr_bact2[:, i]) ** 2

        var_nr_bact1 = np.sqrt(np.divide(var_nr_bact1, nr_simu - 1))
        #var_nr_bact2 = np.sqrt(np.divide(var_nr_bact2, nr_simu - 1))
        print('VAR', var_nr_bact1)
        ### standard deviation of the mean which is
        error_nr_bact_N_tot1 = np.divide(var_nr_bact1, np.sqrt(nr_simu))
        error_nr_bact_N_tot1 = pd.DataFrame(error_nr_bact_N_tot1, columns=rhov)
        avg_nr_bact1 = pd.DataFrame(avg_nr_bact1, columns=rhov)
        var_nr_bact1 = pd.DataFrame(var_nr_bact1)

        #error_nr_bact_N_tot2 = np.divide(var_nr_bact2, np.sqrt(nr_simu))
        #error_nr_bact_N_tot2 = pd.DataFrame(error_nr_bact_N_tot2, columns=rhov)
        #avg_nr_bact2 = pd.DataFrame(avg_nr_bact2, columns=rhov)
        #var_nr_bact2 = pd.DataFrame(var_nr_bact2)

        #print(avg_nr_bact1.to_string())
        #pd.DataFrame(error_nr_bact_N_tot1).to_csv('sd_error_mean_growth_{}.csv'.format(growth1), index=None)
        #pd.DataFrame(error_nr_bact_N_tot2).to_csv('sd_error_mean_growth_{}.csv'.format(growth2), index=None)


        avg_nr_bact1.iloc[-1, 0:(len(part_fact1))].plot(yerr=error_nr_bact_N_tot1.iloc[-1, 0:len(part_fact1)].tolist(), c=colors[0], ax=ax1, logx=True)  ### plot initial nr of bacteria
        #avg_nr_bact2.iloc[-1, 0:(len(part_fact2))].plot(yerr=error_nr_bact_N_tot2.iloc[-1, 0:len(part_fact2)].tolist(),linestyle='dashed', c=c, ax=ax1, logx=True, label='_nolegend_')

        label_list.append('{}, {}'.format(loading, growth))
        print('{}, {}'.format(loading, growth))
        os.chdir('..')
        ind += 1

# chart formatting
xticks = [200, 100, 50, 25, 10, 5]
xticks_label = [5000, 100, 50, 25, 10, 5]
second_ticks = [1, 50, 100, 200, 500, 1000]

ax1.set_xlim(10**(2.3), 10**(0.6))
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks_label)
ax1.set_xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
ax1.set_ylabel(r'\bf{N(300 mins)}')

ax3 = ax1.secondary_xaxis("top")
ax3.xaxis.set_ticks(xticks[::-1], labels=second_ticks[::-1])
ax3.set_xlabel(r'\bf{m (number of subvolumes)}')

ax2.set_xlim(10**(2.3), 10**(0.6))
ax2.set_xticks(xticks)
ax2.set_xticklabels(xticks_label)
ax2.set_xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
ax2.set_ylabel(r'\bf{$P_{s}$}')

ax4 = ax2.secondary_xaxis("top")
ax4.xaxis.set_ticks(xticks[::-1], labels=second_ticks[::-1])
ax4.set_xlabel(r'\bf{m (number of subvolumes)}')

#ax1.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.35), ncol=4, fancybox=True, shadow=True, fontsize= BIGGER_SIZE)
plt.tight_layout()
#plt.savefig('Survival fraction and Nf_vs_part_fact side by side.png', dpi=600)
plt.show()