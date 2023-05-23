import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

import variables
from variables import *

df_total_mass = pd.read_csv('./output/dropnr_1000_loading_rand_growth_gillespie_binary_initialN_10_abconc_57/df_growth_{}_starting_nr_drops_{}_ABconc{}.csv'.format(growth, variables.total_drop_nr, AB_conc))

part_fact = np.loadtxt("output/dropnr_1000_loading_rand_growth_gillespie_binary_initialN_10_abconc_57/part_fact.txt", delimiter=",", unpack=False)

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

pd.DataFrame(error_nr_bact).to_csv('output/sd_error_mean_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc), index=None)

### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
## to plot Nf, Ni  as a function of partition factors


print(error_nr_bact.iloc[-1, 1:len(part_fact)+1].tolist())

avg_nr_bact.iloc[0, 1:(len(part_fact) + 1)].plot(yerr = error_nr_bact.iloc[0, 1:len(part_fact)+1].tolist())  ### plot initial nr of bacteria
avg_nr_bact.iloc[-1, 1:(len(part_fact) + 1)].plot(yerr = error_nr_bact.iloc[-1, 1:len(part_fact)+1].tolist()) ### plot final nr of bacteria
plt.grid(False)
plt.title('Total mass versus partitioning factor, ab conc {}'.format(AB_conc), fontsize=text_size)
plt.ylabel('Total mass (nr of bacteria) Nf, Ni',fontsize=text_size)
plt.xlabel('m (number of subvolumes)', fontsize=text_size)
plt.legend(['Ni', 'Nf'], loc='upper left')
plt.savefig('./output/Nf_Ni_vs_part_fact_startin_nr_of_drops_{} +error.png'.format(variables.total_drop_nr))
plt.show()

