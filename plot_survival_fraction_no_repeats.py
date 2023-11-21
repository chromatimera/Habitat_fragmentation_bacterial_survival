import math
from variables import *
import pandas as pd
import os
import matplotlib.pyplot as plt

BIGGER_SIZE = 22


### FIGURE 5 survival fraction droplets

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

droplet_list = np.arange(0, 1001, 10)
droplet_list[0] = 1
antib = [10, 15, 20, 30]
color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))

m_list = []
survival_outcomes = []
color_list = []

plt.figure(figsize=(8,7))

for ab, c, ind in zip(antib, color, range(len(antib))):
    if ind == 2:
        c = next(color)

    for i in range(len(droplet_list)):

        ### read dataframe and ignore first column as it's the index column
        filename = 'output/dropnr_1000_loading_rand_growth_binary_initialN_5_abconc_{}_gr_0.01/initialN5_growthrate0.01_MIC1_totaldropnr{}_ABconc{}_dt1_loadingrand_growthbinary.csv'.format(ab,droplet_list[i],ab)
        with open(filename) as x:
            ncols = len(x.readline().split(','))
        binary_df = pd.read_csv(filename, usecols=range(1,ncols))

        binary_df_trp = binary_df.T

        last_row = np.array(binary_df_trp.iloc[-1:])
        first_row = np.array(binary_df_trp.iloc[0])

        surviving_last_row = np.count_nonzero(last_row)
        first_row_nonzero = np.count_nonzero(first_row)
        percentage = surviving_last_row/first_row_nonzero
        m_subvol= len(last_row[0,:])

        m_list.append(m_subvol)
        survival_outcomes.append(percentage)


    list_of_tuples = list(zip(m_list, survival_outcomes))
    df = pd.DataFrame(list_of_tuples, columns = ['M', 'Surv frac'])

    df['Error95'] = df.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(x['M']), axis=1)
    df['Error99'] = df.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(x['M']), axis=1)
    df = df.set_index('M')




    df["Surv frac"].plot.line(yerr = df['Error95'], c=c)
    m_list.clear()
    survival_outcomes.clear()
plt.ylabel(r'\bf{Fraction of droplets surviving}')
plt.xlabel(r'\bf{m (number of subvolumes)}')
plt.legend(antib, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.18), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE-5)
plt.savefig('./output/Survival fraction 20 ab more points zoomed in.png', dpi=600)
plt.show()
