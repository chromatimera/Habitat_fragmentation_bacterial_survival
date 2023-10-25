import pandas as pd
import os
import matplotlib.pyplot as plt
from decimal import *
from variables import nr_timesteps


### IGNORE THIS SCRIPT FOR NOW

getcontext().prec = 15
decay1 = 'exponential_decay'
decay2 = 'linear_decay'
#nr_timesteps = int(3e3)
print(nr_timesteps)
def read_csv(decay1, decay2, nr_timesteps):
    os.chdir("./output/df_test_dt_nr_timesteps_{}.csv".format(nr_timesteps))

    df1 = pd.read_csv('df_degradation_{}_nr_timesteps_{}.csv'.format(decay1, nr_timesteps))#, usecols=range(0, 11))
    #print(df1)
    df2 = pd.read_csv('df_degradation_{}_nr_timesteps_{}.csv'.format(decay2, nr_timesteps))#, usecols=range(1, 11))
    #print(df2)
    difference_df = pd.DataFrame(df1 - df2)
    #print(difference_df)

    pd.DataFrame(difference_df).to_csv('difference_df_{}-{}.csv'.format(decay1, decay2), index=None)
    return difference_df

def plot_difference(difference_df):
    column_names = list(difference_df.columns.values)
    plt.figure(1)
    for i in range(0, len(difference_df.axes[1]), 1):
        ## ignoring first row as it's populated with 0 - no degradation present when inputing the values for initial N and ab conc
        plt.plot(difference_df.iloc[1:, i], label=column_names[i])
    plt.grid(True)
    plt.title('Comparison between {} and {} over time'.format(decay1, decay2))
    plt.ylabel('Error in rate (log scale)')
    plt.xlabel('Nr of dts (timesteps)')
    plt.legend(loc = 1)
    plt.yscale('log')
    #plt.show()
    plt.savefig('Difference between {} and {} degradation {} steps.png'.format(decay1, decay2, nr_timesteps))

def main():
    difference = read_csv(decay1, decay2, nr_timesteps)
    #difference = pd.read_csv('./output/df_test_dt_nr_timesteps_{}.csv/difference_df_{}-{}.csv'.format(nr_timesteps, decay1, decay2))
    plot_difference(difference)

main()