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

growth1 = 'binary'
growth2 = 'gillespie_binary'
rootdir = 'output/'
ab = [15, 25]

os.chdir(rootdir)
print('current dir', os.getcwd())

label_list = []
colors = plt.cm.rainbow(np.linspace(0, 1, 5))
color_list = []

plt.figure(figsize=(14, 7))
plt.subplots_adjust(hspace=0.5)

# loop through two types of plots
for n in range(0, 2, 1):
    # add a new subplot iteratively
    ax = plt.subplot(1, 2, n + 1)

    for antib, c in zip(ab, colors):
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
        #label_list.append('{} deterministic'.format(antib))

        # plt.grid(True)
        surv_fraction_errors2 = surv_fraction_transpose2.Error95.to_frame('Surv frac')
        surv_fraction_errors2.index = surv_fraction_errors2.index.map(int)
        surv_fraction_errors2 = surv_fraction_errors2.sort_index(ascending=True)

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
            ## plot survival fraction
            surv_fraction_transpose1["Surv frac"].plot.line(yerr=surv_fraction_errors1, color=c, ax=ax)  # , color = 'orange')
            surv_fraction_transpose2["Surv frac"].plot.line(yerr=surv_fraction_errors2, linestyle='dashed', color=c,
                                                        label='_nolegend_', ax=ax)  # color = 'orange')
            color_list.append(c)

        else:
            ind = ab.index(antib)
            avg_nr_bact1.iloc[-1, 0:(len(part_fact1))].plot(yerr=error_nr_bact_N_tot1.iloc[-1, 0:len(part_fact1)].tolist(),
                                                        logy=False, c=color_list[ind], ax=ax)  ### plot initial nr of bacteria
            avg_nr_bact2.iloc[-1, 0:(len(part_fact2))].plot(yerr=error_nr_bact_N_tot2.iloc[-1, 0:len(part_fact2)].tolist(),
                                                        linestyle='dashed', logy=False, c=color_list[ind], label='_nolegend_', ax=ax)
        label_list.append('{}'.format(antib))
    # chart formatting
    ax.set_xlabel('m (number of subvolumes)', fontsize=text_size)
    if n==0:
       plt.ylim(-0.1, 1.1)
       plt.ylabel('Ps', fontsize=text_size)
    else:
        plt.ylabel('N(300)', fontsize=text_size)
    ax.legend(label_list, title='Antibiotic concentration in μg/mL', loc='upper center', bbox_to_anchor=(0.5, 1.05),
              ncol=3, fancybox=True, shadow=True)

# plt.legend(label_list, title='Antibiotic concentration in μg/mL', loc="lower right", ncol=3)
# #plt.legend(label_list, title='Antibiotic conc', fontsize='medium', loc='upper left')
print(label_list)
# Make room on top now
#plt.legend(label_list,bbox_to_anchor=(-0.9, 1.2), title='Antibiotic concentration in μg/mL', loc="upper center")
plt.savefig('Survival fraction and Nf_vs_part_fact side by side.png')
plt.show()