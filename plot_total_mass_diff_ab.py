from os import listdir
import os
from os.path import isfile, join
from variables import text_size
import pandas as pd
import numpy as np
from variables import total_sim
import matplotlib.pyplot as plt


### FIGURE 4 in paper
growth = 'binary'
rootdir = 'output/'
ab = [35, 45, 55]

os.chdir(rootdir)
print(os.getcwd())

plt.figure(figsize=(7.5,5))

for i in ab:
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, i))
    #print(os.getcwd())
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    #print(onlyfiles)

    df_total_mass = pd.read_csv(onlyfiles[1])
    #print(df_total_mass)
    part_fact = np.loadtxt(onlyfiles[2])
    #print('part fact', part_fact)
    m_list = [round(1/x) for x in part_fact]
    #print('m list', m_list)

    ##start building the average N(t) and SD(t) over simulations
    total_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))
    error_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))
    var_nr_bact = np.zeros((df_total_mass.shape[0], len(part_fact)))

    ## for each partition factor, calculate the sum over the simulation of N(t)
    for i in range(0, len(part_fact), 1):
        for j in range(0, total_sim, 1):
            k = i * 5 + j
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
    #print(error_nr_bact)

    avg_nr_bact = pd.DataFrame(avg_nr_bact, columns=part_fact)
    var_nr_bact = pd.DataFrame(var_nr_bact)
    print('avg_nr_bact', avg_nr_bact)
    print('var_nr_bact', var_nr_bact)

    print(avg_nr_bact.iloc[-1, 0:(len(part_fact))])
    #pd.DataFrame(error_nr_bact).to_csv('output/sd_error_mean.csv', index=None)

    ### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
    ## to plot Nf, Ni  as a function of partition factors

    print(error_nr_bact.iloc[-1, 0:len(part_fact)].tolist())

    #avg_nr_bact.iloc[0, 0:(len(part_fact))].plot(
    #    yerr=error_nr_bact.iloc[0, 1:len(part_fact) + 1].tolist())  ### plot initial nr of bacteria

    ## IF PLOTTING NORMALIZED FRACTIONAL INCREASE
    Nt_over_N0 = (avg_nr_bact.iloc[-1, 0:(len(part_fact))] / avg_nr_bact.iloc[0, 0:len(part_fact)])
    norm_Nt_over_N0 = Nt_over_N0.apply(lambda x: (x - min(Nt_over_N0))/(max(Nt_over_N0)-min(Nt_over_N0)))
    print(norm_Nt_over_N0)
    (norm_Nt_over_N0).plot()

    ##IF PLOTTING FRACTIONAL INCREASE
    (avg_nr_bact.iloc[-1, 0:(len(part_fact))]/avg_nr_bact.iloc[0, 0:len(part_fact)]).plot(yerr=error_nr_bact.iloc[-1, 0:len(part_fact)].tolist())  ### plot final nr of bacteria


    os.chdir('..')

plt.grid(True)
plt.title('Total mass at time 300 versus partitioning factor', fontsize=text_size)
plt.ylabel('Surviving nr of bacteria)', fontsize=text_size)
plt.xlabel('m (number of subvolumes)', fontsize=text_size)
plt.legend(ab, title='Antibiotic conc', fontsize='large',  loc='upper right')
plt.savefig('Nf_vs_part_fact {} + error.png'.format(growth))
plt.show()
