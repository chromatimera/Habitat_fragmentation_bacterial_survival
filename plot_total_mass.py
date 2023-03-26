import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

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


part_fact = np.loadtxt("output/part_fact.txt", delimiter=",", unpack=False)
print(part_fact)
print(len(part_fact))
for i in range(0, len(part_fact)):
    print(part_fact[i])
    df['Error95_{}'.format(part_fact[i])] = df.apply(lambda x: 200 / math.sqrt(variables.total_drop_nr/(variables.total_drop_nr * part_fact[i])), axis=1)

    #df['Error95_{}'.format(part_fact[i])] = df.apply(lambda x: 2 * math.sqrt(x['{}'.format(part_fact[i])] * (1 - x['{}'.format(part_fact[i])])) / math.sqrt(total_drop_nr * part_fact[i]), axis=1)
    #df['Error99_{}'.format(part_fact[i])] = df.apply(lambda x: 2.6 * math.sqrt(x['{}'.format(part_fact[i])] * (1 - x['{}'.format(part_fact[i])])) / math.sqrt(total_drop_nr * part_fact[i]), axis=1)

print(df)
### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
## to plot Nf, Ni  as a function of partition factors
errors = list(df.iloc[0, (len(part_fact)+1):])
print(errors)

df.iloc[0, 1:(len(part_fact)+1)].plot(yerr = errors)  ### plot initial nr of bacteria
df.iloc[-1, 1:(len(part_fact)+1)].plot(yerr = errors) ### plot final nr of bacteria
plt.grid(True)
plt.title('Total mass versus partitioning factor, ab conc {}'.format(AB_conc))
plt.ylabel('Total mass (nr of bacteria) Nf, Ni')
plt.xlabel('Partition factor')
plt.legend(['Ni', 'Nf'], loc='upper left')
plt.savefig('./output/Nf_Ni_vs_part_fact_startin_nr_of_drops_{}.png'.format(variables.total_drop_nr))
plt.show()

