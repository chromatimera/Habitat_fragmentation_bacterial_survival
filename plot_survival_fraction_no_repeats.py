import math
from variables import *
import pandas as pd
import os
import matplotlib.pyplot as plt

BIGGER_SIZE = 28


### FIGURE 5 survival fraction droplets MAIN TEXT

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

print(os.getcwd())

#g1 = 'binary'
#g2 = 'gillespie_binary'
nr_drops = 1000

#droplet_list = np.arange(0, 20001, 400)
droplet_list = np.arange(0, 1001, 20)  ##nv
droplet_list[0] = 1
initialN = 5
antib = [10,15,20,30]
#color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))
color = plt.cm.rainbow(np.linspace(0, 1,5))
color_2 = plt.cm.viridis(np.linspace(0, 100,5))

m_list = []
survival_outcomes = []
color_list = []

colors = ['lightseagreen', color[0], 'deeppink', 'darkblue']
ind= 0

plt.figure(figsize=(8,9))

for ab, index in zip(antib, range(len(antib))):
    #if ind == 2:
    #    c = next(color)
    for i in range(len(droplet_list)):

        ### read dataframe and ignore first column as it's the index column
        filename = 'output/dropnr_{}_loading_rand_growth_binary_initialN_{}_abconc_{}_gr_0.01/initialN{}_growthrate0.01_MIC1_totaldropnr{}_ABconc{}_dt1_loadingrand_growthbinary.csv'.format(
            droplet_list[-1], initialN, ab, initialN, droplet_list[i], ab)
        with open(filename) as x:
            ncols = len(x.readline().split(','))
        binary_df = pd.read_csv(filename, usecols=range(1,ncols))

        binary_df_trp = binary_df.T

        last_row = np.array(binary_df_trp.iloc[-1:])
        first_row = np.array(binary_df_trp.iloc[0])

        surviving_last_row = np.count_nonzero(last_row)
        first_row_nonzero = np.count_nonzero(first_row)
        if first_row_nonzero == 0:
            percentage = 0
        else:

            percentage = surviving_last_row/first_row_nonzero
        m_subvol= len(last_row[0,:])

        m_list.append(m_subvol)
        survival_outcomes.append(percentage)


    list_of_tuples = list(zip(m_list, survival_outcomes))
    df = pd.DataFrame(list_of_tuples, columns = ['M', 'Surv frac'])

    df['Error95'] = df.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(x['M']), axis=1)
    df['Error99'] = df.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(x['M']), axis=1)
    df = df.set_index('M')




    df["Surv frac"].plot.line(yerr = df['Error95'], c = colors[ind])
    ind += 1
    print(df)
    m_list.clear()
    survival_outcomes.clear()
    plt.ylabel(r'\bf{Subpopulation survival probability}' + "\n" + r'\bf{p_{s} for occupied subvolumes}',
               multialignment='center')
    #second_ticks = [1, 5000, 10000, 15000, 20000]
    #plt.xticks(second_ticks)

plt.ylabel(r'Subpopulation survival probability $p_s$')
plt.xlabel(r'm (number of subvolumes)')
plt.tight_layout()
#plt.legend(antib, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.18), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE-5)
plt.savefig('./output/Survival fraction {} no legend .png'.format(droplet_list[-1]), dpi=600)
#plt.show()
