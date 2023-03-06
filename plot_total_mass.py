import pandas as pd
import os
import matplotlib.pyplot as plt
from variables import *

print(os.getcwd())
g = 'gillespie_binary'
drops_10 = pd.read_csv('./output/df_growth_gillespie_binary_nr_drops_1000_big_droplet_nr_1.csv')
big_drop = pd.read_csv('./output/df_growth_gillespie_binary_nr_drops_1_big_droplet_nr_1000.csv')

difference = drops_10 - big_drop
print(difference)
# plt.grid(True)
difference.iloc[0:20000,10:15].plot()
plt.title('Total mass versus time, growth {}'.format(g))
plt.ylabel('Total mass (nr of bacteria)')
plt.xlabel('Time (min)')
#plt.legend(loc = 'upper left')
plt.savefig('./output/Total mass growth {}, nr of droplets {}, difference.png'.format(g, 10000))
plt.show()

# ax = drops_10.iloc[:,1:2].plot()
# big_drop.iloc[:,1:2].plot()
# plt.title('Total mass versus time, growth {}'.format(g))
# plt.ylabel('Total mass (nr of bacteria)')
# plt.xlabel('Time (min)')
# #plt.legend(loc = 'upper left')
# plt.savefig('./output/Total mass growth {}, nr of droplets {}.png'.format(g, 10000))
# plt.show()
