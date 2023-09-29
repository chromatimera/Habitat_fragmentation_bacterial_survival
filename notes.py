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

BIGGER_SIZE = 16

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsiz

growth1 = 'binary'
growth2 = 'gillespie_binary'
rootdir = 'output/'
ab = [35]

os.chdir(rootdir)
print('current dir', os.getcwd())

colors = plt.cm.rainbow(np.linspace(0, 1, 5))
color_list = []
label_list = []

rho = 5e7

fig = plt.figure(figsize=(11,13))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
plt.subplots_adjust(hspace=0.5)

# loop through two types of plots
for n in range(0, 2, 1):
    # add a new subplot iteratively
    ax = plt.subplot(2, 1, n + 1)

    for antib, c in zip(ab, colors):
        #print(ind)
        if Counter(c) == Counter(colors[2]):
           c = colors[3]
        ### plot the deterministic value + error
        os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth1, antib))
        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)

        surv_fraction1 = pd.read_csv(onlyfiles[3])
        part_fact1 = np.loadtxt(onlyfiles[2])
        df_total_mass1 = pd.read_csv(onlyfiles[1])
        m_list1 = [round(1 / x) for x in part_fact1]

        os.chdir('..')
        ### plot the stochastic value + error
        os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth2, antib))
        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)

        surv_fraction2 = pd.read_csv(onlyfiles[3])
        part_fact2 = np.loadtxt(onlyfiles[2])
        df_total_mass2 = pd.read_csv(onlyfiles[1])
        m_list2 = [round(1 / x) for x in part_fact2]

        os.chdir('..')

        ### DO THE SURVIVAL FRACTION DFS
        surv_fraction_transpose1 = surv_fraction1.T
        surv_fraction_transpose1.index.name = 'Part_fact'

        surv_fraction_transpose2 = surv_fraction2.T
        surv_fraction_transpose2.index.name = 'Part_fact'

        surv_fraction_transpose1.columns = ['Surv frac']
        surv_fraction_transpose1['M'] = surv_fraction_transpose1.index.astype(int)
        surv_fraction_transpose1['RhoV'] = surv_fraction_transpose1.apply(lambda x: rho * 1e-4 / x['M'], axis=1)
        surv_fraction_transpose1['Error95'] = surv_fraction_transpose1.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
        surv_fraction_transpose1['Error99'] = surv_fraction_transpose1.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)

        surv_fraction_transpose2.columns = ['Surv frac']
        surv_fraction_transpose2['M'] = surv_fraction_transpose2.index.astype(int)
        surv_fraction_transpose2['Error95'] = surv_fraction_transpose2.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
        surv_fraction_transpose2['Error99'] = surv_fraction_transpose2.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac'])) / math.sqrt(variables.total_sim), axis=1)
        surv_fraction_transpose2['RhoV'] = surv_fraction_transpose2.apply(lambda x: rho * 1e-4 / x['M'], axis=1)

        # plt.grid(True)
        surv_fraction_errors1 = surv_fraction_transpose1.Error95.to_frame('Surv frac')
        surv_fraction_errors1.index = surv_fraction_errors1.index.map(int)
        surv_fraction_errors1.index = surv_fraction_transpose1['RhoV']
        #label_list.append('{} deterministic'.format(antib))

        # plt.grid(True)
        surv_fraction_errors2 = surv_fraction_transpose2.Error95.to_frame('Surv frac')
        surv_fraction_errors2.index = surv_fraction_errors2.index.map(int)
        surv_fraction_errors2.index = surv_fraction_transpose2['RhoV']
        print('antib', antib)
        #print('surv frac errors 1', surv_fraction_errors1)
        #print('surv frac errors 2', surv_fraction_errors2)

        surv_fraction_transpose2.index = surv_fraction_transpose2.index.map(int)
        surv_fraction_transpose2 = surv_fraction_transpose2.sort_index(ascending=True)

        ## DO THE N(T) DFS

        ##start building the average N(t) and SD(t) over simulations
        total_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))
        error_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))
        var_nr_bact1 = np.zeros((df_total_mass1.shape[0], len(part_fact1)))

        ##start building the average N(t) and SD(t) over simulations
        total_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))
        error_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))
        var_nr_bact2 = np.zeros((df_total_mass2.shape[0], len(part_fact2)))

        ## for each partition factor, calculate the sum over the simulation of N(t)
        for i in range(0, len(part_fact1), 1):
            for j in range(0, total_sim, 1):
                k = i * total_sim + j
                total_nr_bact1[:, i] += df_total_mass1.iloc[:, k]
                total_nr_bact2[:, i] += df_total_mass2.iloc[:, k]

        ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
        nr_simu = np.array(total_sim)
        avg_nr_bact1 = np.divide(total_nr_bact1, nr_simu)
        avg_nr_bact2 = np.divide(total_nr_bact2, nr_simu)

        for i in range(0, len(part_fact1), 1):
            for j in range(0, total_sim, 1):
                k = j + i * total_sim
                var_nr_bact1[:, i] += (df_total_mass1.iloc[:, k] - avg_nr_bact1[:, i]) ** 2
                var_nr_bact2[:, i] += (df_total_mass2.iloc[:, k] - avg_nr_bact2[:, i]) ** 2

        var_nr_bact1 = np.sqrt(np.divide(var_nr_bact1, nr_simu - 1))
        var_nr_bact2 = np.sqrt(np.divide(var_nr_bact2, nr_simu - 1))

        ### standard deviation of the mean which is
        error_nr_bact_N_tot1 = np.divide(var_nr_bact1, np.sqrt(nr_simu))
        error_nr_bact_N_tot1 = pd.DataFrame(error_nr_bact_N_tot1, columns=part_fact1)
        avg_nr_bact1 = pd.DataFrame(avg_nr_bact1, columns=part_fact1)
        var_nr_bact1 = pd.DataFrame(var_nr_bact1)

        error_nr_bact_N_tot2 = np.divide(var_nr_bact2, np.sqrt(nr_simu))
        error_nr_bact_N_tot2 = pd.DataFrame(error_nr_bact_N_tot2, columns=part_fact2)
        avg_nr_bact2 = pd.DataFrame(avg_nr_bact2, columns=part_fact2)
        var_nr_bact2 = pd.DataFrame(var_nr_bact2)

        pd.DataFrame(error_nr_bact_N_tot1).to_csv('sd_error_mean_growth_{}.csv'.format(growth1), index=None)
        pd.DataFrame(error_nr_bact_N_tot2).to_csv('sd_error_mean_growth_{}.csv'.format(growth2), index=None)

        #os.chdir('..')

        if n==0:

            avg_nr_bact1.iloc[-1, 0:(len(part_fact1))].plot(
                yerr=error_nr_bact_N_tot1.iloc[-1, 0:len(part_fact1)].tolist(),
                logy=False, c=c, ax=ax2)  ### plot initial nr of bacteria
            avg_nr_bact2.iloc[-1, 0:(len(part_fact2))].plot(
                yerr=error_nr_bact_N_tot2.iloc[-1, 0:len(part_fact2)].tolist(),
                linestyle='dashed', logy=False, c=c, ax=ax2)
            color_list.append(c)
            #label = '_nolegend_'

        else:
            ind = ab.index(antib)
            ## plot survival fraction
            surv_fraction_transpose1["Surv frac"].plot.line(yerr=surv_fraction_errors1, c=color_list[ind], ax=ax1, logx = True)  # , color = 'orange')
            surv_fraction_transpose2["Surv frac"].plot.line(yerr=surv_fraction_errors2, linestyle='dashed',  c=color_list[ind], ax=ax1, logx = True)  # color = 'orange')
        label_list.append('{}'.format(antib))
    # chart formatting
    ax2.set_xlabel(r'\bf{m (number of subvolumes)}')
    if n==0:
        plt.ylabel(r'\bf{N(300 mins)}')

    else:

        plt.xlim(10 ** (2.3), 10 ** (0.6))

        xticks = [200, 100, 50, 25, 10, 5]
        xticks_label = [5000, 100, 50, 25, 10, 5]
        second_ticks = [1, 50, 100, 200, 500, 1000]

        ax1.set_xticks(xticks)
        ax1.set_xticklabels(xticks_label)

        ax3 = ax1.twiny()
        ax3.set_xticks(np.linspace(ax1.get_xticks()[0], ax1.get_xticks()[-1], len(ax1.get_xticks())))
        ax3.set_xticklabels(second_ticks[::-1])

        plt.xlabel(r'\bf{m (number of subvolumes)}')
        ax1.set_xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
        ax1.set_ylabel(r'\bf{Probability of survival}')

        ax1.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center',
                   bbox_to_anchor=(0.5, 1.45), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)

    #ax.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=4, fancybox=True, shadow=True, fontsize= 'large')
plt.savefig('Survival fraction and Nf_vs_part_fact side by side.png')
plt.show()