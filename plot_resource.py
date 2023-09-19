from os import listdir
import os
from os.path import isfile, join
from variables import text_size
import pandas as pd
import numpy as np
from variables import total_sim
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

ab = [75]
growth = 'resource'
rootdir = 'output/'


os.chdir(rootdir)
color = iter(cm.rainbow(np.linspace(0, 1, 5)))

plt.figure(figsize=(7.5,5))

for antib, c in zip(ab, color):
    os.chdir('RESOURCE_dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)

    df_bact_count = pd.read_csv(onlyfiles[1])
    part_fact = np.loadtxt(onlyfiles[2]) #m
    m_list = [round(1/x) for x in part_fact]
    ##start building the average N(t) and SD(t) over simulations
    total_nr_bact = np.zeros((df_bact_count.shape[0], len(part_fact)))
    error_nr_bact = np.zeros((df_bact_count.shape[0], len(part_fact)))
    var_nr_bact = np.zeros((df_bact_count.shape[0], len(part_fact)))


    #TIME== x  axis

    plt.plot(df_bact_count['0 1'], label='m=1')
    plt.plot(df_bact_count['0 1000'], label='m=1000')




    ## for each partition factor, calculate the sum over the simulation of N(t)
  #  for i in range(0, len(part_fact), 1):
  #      for j in range(0, total_sim, 1):
  #          k = i * 5 + j  ##??
   #         total_nr_bact[:, i] += df_total_mass.iloc[:, k]

    ## at this stage we have a np array with a sum of N(t) across all iterations -> we need to divide it by the nr of sim
 #   nr_simu = np.array(total_sim)
  #  avg_nr_bact = np.divide(total_nr_bact, nr_simu)
#
 #   for i in range(0, len(part_fact), 1):
  #      for j in range(0, total_sim, 1):
  #          k = j + i * total_sim
 #           var_nr_bact[:, i] += (df_total_mass.iloc[:, k] - avg_nr_bact[:, i]) ** 2

   # var_nr_bact = np.sqrt(np.divide(var_nr_bact, nr_simu - 1))

    ### standard deviation of the mean which is
  #  error_nr_bact = np.divide(var_nr_bact, np.sqrt(nr_simu))
   # error_nr_bact = pd.DataFrame(error_nr_bact, columns=part_fact)
    # print('error nr bact', error_nr_bact)
   # avg_nr_bact = pd.DataFrame(avg_nr_bact, columns=part_fact)
   # var_nr_bact = pd.DataFrame(var_nr_bact)

   # pd.DataFrame(error_nr_bact).to_csv('sd_error_mean.csv', index=None)

    ### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
    ## to plot Nf, Ni  as a function of partition factors

    # print(error_nr_bact.iloc[0, 1:len(part_fact)+1].tolist())

    ### old way of plotting
   # avg_nr_bact.iloc[-1, 0:(len(part_fact))].plot(yerr=error_nr_bact.iloc[-1, 0:len(part_fact)].tolist(), logy=True,
                                           #       c=c)  ### plot initial nr of bacteria     is this intial or final?

    # ## IF PLOTTING NORMALIZED FRACTIONAL INCREASE
    # Nt_over_N0 = (avg_nr_bact.iloc[-1, 0:(len(part_fact))] / avg_nr_bact.iloc[0, 0:len(part_fact)])
    # norm_Nt_over_N0 = Nt_over_N0.apply(lambda x: (x - min(Nt_over_N0))/(max(Nt_over_N0)-min(Nt_over_N0)))
    # print(norm_Nt_over_N0)
    # (norm_Nt_over_N0).plot()
    #
    # ##IF PLOTTING FRACTIONAL INCREASE
    # (avg_nr_bact.iloc[-1, 0:(len(part_fact))]/avg_nr_bact.iloc[0, 0:len(part_fact)]).plot(yerr=error_nr_bact.iloc[-1, 0:len(part_fact)].tolist())  ### plot final nr of bacteria
    #

    os.chdir('..')

plt.grid(False)
plt.ylabel('Number of bacteria', fontsize=text_size)
plt.xlabel('Time (min)', fontsize=text_size)
plt.legend( title='Number of subvolumes', fontsize='large', loc='upper right')
plt.show()
