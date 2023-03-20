import pandas as pd
import os
import matplotlib.pyplot as plt

import variables
from variables import *

print(os.getcwd())
df = pd.read_csv('./output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr))


### This plot plots Nr of bacteria as a function of time; each line is a different multiplication factor ### this is only for troubleshooting purposes
df.iloc[1:,1:].plot()
plt.grid(True)
plt.title('Total mass versus partitioning factor, growth {}'.format(growth))
plt.ylabel('Total mass (nr of bacteria)')
plt.xlabel('Partition factor')
plt.legend(loc = 'upper left')
#plt.savefig('./output/Total mass growth {}, nr of droplets {}, difference.png'.format(g, 10000))
plt.show()

## to plot Nf  as a function of partition factors
