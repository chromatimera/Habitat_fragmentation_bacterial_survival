import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import os
from os.path import isfile, join
import math

BIGGER_SIZE = 22

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) ## https://matplotlib.org/stable/tutorials/text/usetex.html
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize


growth = 'binary'
total_sim = 10
m=1

Ni = np.arange(1, 10, 1).tolist()
antib = np.arange(0, 56, 5).tolist()
print(antib)



os.chdir('./output/')
print(os.getcwd())
#for ab in antib:
os.chdir('./changing lamda for binary/')
print(os.getcwd())

colnames=['lambda']

df_heatmap_survival=pd.DataFrame(Ni)
df_heatmap_survival.columns=colnames
print(df_heatmap_survival)

df_heatmap_survival.set_index('lambda', inplace=True)

colnames = [str(x) for x in antib]
df_heatmap_survival[colnames] = np.random.randint(10, size=(len(Ni), len(antib)))


print(df_heatmap_survival)

df_heatmap_survival.insert = antib
s=0
for l in Ni:
    surv_diff_ab_same_Ni = []

    for ab in antib:
        print(os.getcwd())

        os.chdir('./dropnr_1000_loading_rand_growth_binary_initialN_{}_abconc_{}/'.format(l, ab))
        ### calculate probability of survival

        path = os.getcwd()

        onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
        onlyfiles = sorted(onlyfiles)

        surv_fraction = pd.read_csv(onlyfiles[3])
        part_fact = np.loadtxt(onlyfiles[2])
        os.chdir('..')

        surv_diff_ab_same_Ni.append(surv_fraction.iloc[0,0])
        print(surv_diff_ab_same_Ni)

    ## append last row from the surv_fraction dataframe
    print(surv_diff_ab_same_Ni,'surv_diff_ab_same_Ni')
    df_heatmap_survival.iloc[s] = surv_diff_ab_same_Ni
    s = s + 1


df_heatmap_survival['Rho'] = df_heatmap_survival.index * 1e7
df_heatmap_survival['Rho'] = df_heatmap_survival['Rho'].astype(int)
df_heatmap_survival.set_index('Rho', inplace=True, drop=True)
print(df_heatmap_survival)
rho_list = list(df_heatmap_survival.index)

# Define the plot
fig, ax = plt.subplots(figsize=(8,5))

# Add title to the Heat map

os.chdir('..')

# Use the heatmap function from the seaborn package
sns.heatmap(df_heatmap_survival, annot=True)
ax.invert_yaxis()
plt.xlabel(r'\bf{Initial antibiotic concentration ($\mu$g/ml)}')
plt.ylabel(r'\bf{$\rho$ (initial cells/ml)}')
plt.tight_layout()
plt.savefig('heatmap {}.png'.format(m), dpi=600)
plt.show()

#for each simulation out of the total nr of simulations and for each
# loading factor calculate the probability of survival



