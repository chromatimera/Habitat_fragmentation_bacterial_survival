import pandas as pd
import os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import variables
from variables import *


def read_csv():
    df = pd.read_csv('./output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr))
    return df

### PLS IGNORE THIS BIT - FOR ADDING ERROR BARS
#def add_error95_bar(df):
#    df['Error95_{}'.format(x)] = df.apply(lambda x: 2 * math.sqrt(df.iloc[:, x] * (1 - df.iloc[:, x])) / math.sqrt(variables.total_drop_nr * ), axis=1)
#    return df
print(os.getcwd())
df = read_csv()
#df = add_error95_bar()
print(df)

new_df = pd.DataFrame()
fig = plt.figure()
ax = plt.axes(projection='3d')

# for m in range(abmin, abmax, step):
#     for i in range(part_min, part_max, step):
#         for j in range(part_min, part_max, step):
#             new_i = 5 ** i * 2 ** j
#             new_nr_drops_total_mass = new_i
#             total_drop_nr = round(variables.total_drop_nr / new_nr_drops_total_mass)
#             part = 1/total_drop_nr
#             new_df.append(df['p_{}'.format(part)])
ax.plot_surface(x, y, df.iloc[-1, 1:], cmap='viridis', edgecolor='none')
ax.set_title('Surface plot')
plt.show()



### PLS IGNORE THIS BIT - FOR ADDING ERROR BARS
#adaptive_tau_binary_transposed['Error95'] = adaptive_tau_binary_transposed.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
#adaptive_tau_binary_transposed['Error99'] = adaptive_tau_binary_transposed.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)

## to plot Nf, Ni  as a function of partition factors
df.iloc[0, 1:].plot()  ### plot initial nr of bacteria
df.iloc[-1, 1:].plot() ### plot final nr of bacteria
plt.grid(True)
plt.title('Total mass versus partitioning factor, ab conc {}'.format(AB_conc))
plt.ylabel('Total mass (nr of bacteria) Nf, Ni')
plt.xlabel('Partition factor')
plt.legend(['Ni', 'Nf'], loc='upper left')
plt.savefig('./output/Nf_Ni_vs_part_fact_startin_nr_of_drops_{}.png'.format(variables.total_drop_nr))
plt.show()

