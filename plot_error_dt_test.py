import pandas as pd
import os
import matplotlib.pyplot as plt
from decimal import *
from variables import nr_timesteps

getcontext().prec = 15
decay1 = 'MM_exponential'
decay2 = 'MM_linear'

print(os.getcwd())
os.chdir("./output/df_test_dt_nr_timesteps_{}.csv".format(nr_timesteps))
print(os.getcwd())

df1 = pd.read_csv('df_degradation_{}_nr_timesteps_{}.csv'.format(decay1, nr_timesteps))#, usecols=range(0, 11))
df2 = pd.read_csv('df_degradation_{}}_nr_timesteps_{}.csv'.format(decay2, nr_timesteps))#, usecols=range(1, 11))

difference_df = pd.DataFrame(df1 - df2)
print('Difference df', difference_df)

column_names = list(difference_df.columns.values)
pd.DataFrame(difference_df).to_csv('difference_df_{}-{}.csv'.format(decay1, decay2), index=None)

plt.figure(1)
#plt.grid(True)
for i in range(0, 5, 1): #len(expo_df.axes[1]), 1):
    ## ignoring first row as it's populated with 0 - no degradation present when inputing the values for initial N and ab conc
    plt.plot(difference_df.iloc[1:,i], label=column_names[i])
#plt.show()
# plt.grid(True)
# #plt.title('Concentration of antibiotic over time in each droplet.')
# plt.ylabel('Error in rate (log scale)')
# plt.xlabel('Nr of dts (timesteps)')
# plt.legend()
plt.yscale('log')
# #plt.ylim(bottom=0)
plt.show()
# plt.savefig('Difference between linear and exponential degradation {} steps.svg'.format(nr_timesteps))
