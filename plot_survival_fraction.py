import math

import variables
from variables import *
import pandas as pd
import os
import matplotlib.pyplot as plt

print(os.getcwd())

### read dataframes
surv_fraction = pd.read_csv('output/survival_fraction_growth_{}_loading_{}_ABconc{}.csv'.format(growth, loading, AB_conc))
print(print(surv_fraction.columns.tolist()))

### transpose of dataframe
surv_fraction_transpose = surv_fraction.T
surv_fraction_transpose.index.name = 'Part_fact'
print(surv_fraction_transpose)


surv_fraction_transpose.columns = ['Surv frac']
surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)


plt.figure(1)
# plt.grid(True)
surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
surv_fraction_transpose["Surv frac"].plot.line(yerr = surv_fraction_errors)#, color = 'orange')
plt.title('Fraction of droplets surviving a round of treatment')
plt.ylabel('Fraction of droplets surviving')
plt.xlabel('Partitioning factor')
plt.savefig('./output/Survival fraction Gillespie + errors.png')
plt.show()