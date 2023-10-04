import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import os
from os.path import isfile, join
import matplotlib.ticker as tkr
import math

BIGGER_SIZE = 32

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)  ## https://matplotlib.org/stable/tutorials/text/usetex.html
plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize

growth = 'binary'
total_sim = 1
m = 1

Ni = np.arange(1, 9, 1).tolist()
antib = np.arange(0, 26, 1).tolist()
print(Ni)
print(antib)

os.chdir('./output/')
print(os.getcwd())
# for ab in antib:
os.chdir('./survival_fraction_heatmap_{}/'.format(m))
print(os.getcwd())

colnames = ['lambda']
df_heatmap_survival = pd.DataFrame(Ni)
df_heatmap_survival.columns = colnames
print(df_heatmap_survival)

df_heatmap_survival.set_index('lambda', inplace=True)

colnames = [str(x) for x in antib]
df_heatmap_survival[colnames] = np.random.randint(10, size=(len(Ni), len(antib)))

print(df_heatmap_survival)

df_heatmap_survival.insert = antib
s = 0
for l in Ni:
    survival_outcomes = []

    for ab in antib:
        path = os.getcwd()
        os.chdir('./dropnr_{}_loading_rand_growth_binary_initialN_{}_abconc_{}_gr_0.01/'.format(m, l, ab))
        ### calculate probability of survival

        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)

        if '.DS_Store' in onlyfiles:
            onlyfiles.remove('.DS_Store')
        else:
            pass

        surv_fraction = pd.read_csv(onlyfiles[1])
        os.chdir('..')

        #
        binary_df_trp = surv_fraction.T
        last_row = np.array(binary_df_trp.iloc[-1:])

        surviving = np.count_nonzero(last_row)

        percentage = surviving / len(last_row[0, :])
        m_subvol = len(last_row[0, :])

        survival_outcomes.append(percentage)

        #     ## append last row from the surv_fraction dataframe
    df_heatmap_survival.iloc[s] = survival_outcomes
    s = s + 1

print(df_heatmap_survival)

df_heatmap_survival['Rho'] = df_heatmap_survival.index * 1e7
df_heatmap_survival['Rho'] = df_heatmap_survival['Rho'].astype(int)
df_heatmap_survival.set_index('Rho', inplace=True, drop=True)
print(df_heatmap_survival)
rho_list = list(df_heatmap_survival.index)

# Define the plot
fig, ax = plt.subplots(figsize=(15, 8))
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
os.chdir('..')

# Use the heatmap function from the seaborn package
sns.heatmap(df_heatmap_survival, annot=True)
ax.invert_yaxis()
plt.xlabel(r'\bf{Initial antibiotic concentration ($\mu$g/ml)}')
plt.ylabel(r'\bf{$\rho$ (initial cells/ml)}')
plt.tight_layout()
plt.savefig('heatmap {}.png'.format(m), dpi=600)
plt.show()
