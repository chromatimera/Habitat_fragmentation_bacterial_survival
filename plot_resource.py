from os import listdir
import os
from os.path import isfile, join
from variables import text_size
import pandas as pd
import numpy as np
from variables import total_sim
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

BIGGER_SIZE = 22
LEGEND_SIZE = 20

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
plt.rcParams['figure.figsize'] = [8, 6.5]



ab = [100]
sim_reps= 100
dt=1
growth = 'resource'
rootdir = 'output/'


os.chdir(rootdir)
color = iter(cm.rainbow(np.linspace(0, 1, 5)))


for antib, c in zip(ab, color):
    os.chdir('RESOURCE_dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)

    df_bact_count = pd.read_csv(onlyfiles[1])
    part_fact = np.loadtxt(onlyfiles[2]) #m
    m_list = [round(1/x) for x in part_fact]
    df_average= pd.read_csv(onlyfiles[0])
    error_sem = np.zeros((df_bact_count.shape[0], len(part_fact)))
    ### standard deviation of the mean;
 #   error_sem = np.divide(np.std(df_bact_count,1), np.sqrt(sim_reps))
   #
    j=0
    k=0 #better way to do this; cant use col name in data frame..?
    for i in range(0, len(part_fact), 1):
            k = k+ sim_reps #end col
            error_sem[:,i] = np.divide(np.std(df_bact_count.iloc[:,j:k],1), np.sqrt(sim_reps))
            j = j + sim_reps #start col

    error_sem = pd.DataFrame(error_sem, columns=part_fact)
    #TIME== x  axis
    plt.figure(1)
    timee = np.arange(0, 100, dt)
    plt.plot(timee, df_bact_count['0 1'], label='m=1')
    plt.plot(timee,df_bact_count['0 1000'], label='m=1000')

    for k in range (0,sim_reps):  #for each sim repeat
     kk1=str(k)+' 1'
     kk500 = str(k) + ' 500'
     kk1000=str(k)+' 1000'

     plt.figure(2)
     timee_cut=np.arange(0, 100, dt)
     indx=int(100/dt)
     first_100= df_bact_count.iloc[0:indx]/df_bact_count.iloc[0]
     plt.plot( timee_cut,first_100[ kk1], label='m=1', color='g')
     plt.plot( timee_cut,first_100[ kk500], label='m=500', color='b')
     plt.plot(timee_cut,first_100[kk1000], label='m=1000', color='m')

     plt.figure(3)
     plt.plot( timee,df_bact_count[ kk1], label='m=1', color='g')
     plt.plot(timee,df_bact_count[kk1000], label='m=1000', color='m')

    os.chdir('..')


plt.grid(False)
plt.figure(2)
plt.ylabel('', )
plt.ylabel(r'\bf{N/N0}')
plt.xlabel(r'\bf{Time (min)}')
plt.legend([r'm=1', r'm=500', r'm=1000'], title=r'\bf{Number of subvolumes:}',fontsize = LEGEND_SIZE, title_fontsize=LEGEND_SIZE, loc='upper left')
plt.tight_layout()
plt.savefig('100_ugml-resource.png', dpi=600)
plt.show()

###Average;
plt.figure(4)
#plt.plot(timee, df_average['0 1'], label='m=1')
#plt.plot(timee, df_average['1'], label='m=1000')
plt.fill_between(timee, df_average['1'] - error_sem.iloc[:,0], df_average['1'] + error_sem.iloc[:,0],
                 color='gray', alpha=0.2)
plt.fill_between(timee, df_average['1000'] - error_sem.iloc[:,2], df_average['1000'] + error_sem.iloc[:,2],
                 color='gray', alpha=0.2)
plt.show()



## range??


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
