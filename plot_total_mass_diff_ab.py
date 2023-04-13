from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import total_sim
import matplotlib.pyplot as plt

growth = 'binary'
rootdir = 'output/'
ab = [40, 55, 70]

os.chdir(rootdir)
print(os.getcwd())
for i in ab:
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_10_abconc_{}'.format(growth, i))
    print(os.getcwd())
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    print(onlyfiles)

    df_total_mass = pd.read_csv(onlyfiles[1])
    print(df_total_mass)
    part_fact = np.loadtxt(onlyfiles[2])
    print(part_fact)

    ##start building the average N(t) and SD(t) over simulations
    total_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))
    error_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))
    var_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))

    ## for each partition factor, calculate the sum over the simulation of N(t)
    for i in range(0, len(part_fact), 1):
        for j in range(0, total_sim, 1):
            k = j + i * total_sim
            total_nr_bact[:, i] += df_total_mass.iloc[:, k]

    ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
    nr_simu = np.array(total_sim)
    avg_nr_bact = np.divide(total_nr_bact, nr_simu)

    for i in range(0, len(part_fact), 1):
        for j in range(0, total_sim, 1):
            k = j + i * total_sim
            var_nr_bact[:, i] += (df_total_mass.iloc[:, k] - avg_nr_bact[:, i]) ** 2

    var_nr_bact = np.sqrt(np.divide(var_nr_bact, nr_simu - 1))

    ### standard deviation of the mean which is
    error_nr_bact = np.divide(var_nr_bact, np.sqrt(nr_simu))
    error_nr_bact = pd.DataFrame(error_nr_bact, columns=part_fact)
    print(error_nr_bact)

    avg_nr_bact = pd.DataFrame(avg_nr_bact, columns=part_fact)
    var_nr_bact = pd.DataFrame(var_nr_bact)

    #pd.DataFrame(error_nr_bact).to_csv('output/sd_error_mean.csv', index=None)

    ### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
    ## to plot Nf, Ni  as a function of partition factors

    print(error_nr_bact.iloc[-1, 1:len(part_fact) + 1].tolist())

    #avg_nr_bact.iloc[0, 1:(len(part_fact) + 1)].plot(
    #    yerr=error_nr_bact.iloc[0, 1:len(part_fact) + 1].tolist())  ### plot initial nr of bacteria
    avg_nr_bact.iloc[-1, 1:(len(part_fact) + 1)].plot(yerr=error_nr_bact.iloc[-1, 1:len(part_fact) + 1].tolist())  ### plot final nr of bacteria

    os.chdir('..')

plt.grid(True)
plt.title('Total mass at time 300 versus partitioning factor')
plt.ylabel('Total mass (nr of bacteria)')
plt.xlabel('Partition factor')
plt.legend(ab, title='Antibiotic conc',  loc='upper right')
plt.savefig('Nf_vs_part_fact {} + error.png'.format(growth))
plt.show()







