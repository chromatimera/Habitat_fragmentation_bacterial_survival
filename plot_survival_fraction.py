import math
from variables import *
import pandas as pd
import os
import matplotlib.pyplot as plt

print(os.getcwd())

g1 = 'adaptive_tau_binary'
g2 = 'gillespie_balanced'
nr_drops = 1000


### read dataframes
#tau_df = pd.read_csv('output/df_growth_{}_count_survival_nr_drop_10000.csv'.format(g1))
gillespie_df = pd.read_csv('output/df_growth_gillespie_balanced_count_survival_nr_drop_1000_ab_range_0_50.csv')
midpoint_tau_binary = pd.read_csv('output/df_growth_midpoint_tau_balanced_count_survival_nr_drop_1000_ab_range_0_50.csv')
adaptive_tau_binary = pd.read_csv('output/df_growth_adaptive_tau_balanced_count_survival_nr_drop_1000_ab_range_0_50.csv')


### transpose of dataframe
adaptive_tau_binary_transposed = adaptive_tau_binary.T
gillespie_transposed = gillespie_df.T
midpoint_tau_binary_transposed = midpoint_tau_binary.T

print(adaptive_tau_binary_transposed)

adaptive_tau_binary_transposed.columns = ['Surv frac']
adaptive_tau_binary_transposed['Error95'] = adaptive_tau_binary_transposed.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
adaptive_tau_binary_transposed['Error99'] = adaptive_tau_binary_transposed.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)


gillespie_transposed.columns = ['Surv frac']
gillespie_transposed['Error95'] = gillespie_transposed.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
gillespie_transposed['Error99'] = gillespie_transposed.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)

midpoint_tau_binary_transposed.columns = ['Surv frac']
midpoint_tau_binary_transposed['Error95'] = midpoint_tau_binary_transposed.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
midpoint_tau_binary_transposed['Error99'] = midpoint_tau_binary_transposed.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)


if len(adaptive_tau_binary.columns) == len(gillespie_df.columns):



    plt.figure(1)
    # plt.grid(True)
    adaptive_tau_binary_errors = adaptive_tau_binary_transposed.Error95.to_frame('Surv frac')
    adaptive_tau_binary_transposed["Surv frac"].plot.line(yerr = adaptive_tau_binary_errors)#, color = 'orange')

    gillespie_errors = gillespie_transposed.Error95.to_frame('Surv frac')
    gillespie_transposed["Surv frac"].plot.line(yerr = gillespie_errors)#, color = 'blue')

    midpoint_tau_binary_errors = midpoint_tau_binary_transposed.Error95.to_frame('Surv Frac')
    midpoint_tau_binary_transposed['Surv frac'].plot.line(yerr = midpoint_tau_binary_errors, color = 'green')
    plt.title('Fraction of droplets surviving a round of treatment with same loading')
    plt.ylabel('Fraction of droplets surviving')
    plt.xlabel('Antibiotic concentration (um/ml)')
    plt.legend([ 'Adaptive_tau binary epsilon = 0.03', 'Gillespie binary', 'Midpoint tau binary epsilon = 0.03'],
               loc='upper right', title='Growth')
    plt.savefig('./output/Survival fraction adaptive vs gillespie balanced g.png')
    plt.show()

else:
    print('Dataframes do not have the same nr of columns')