import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import os
from os.path import isfile, join
import math

growth = 'binary'
total_sim = 5

Ni = np.arange(5, 106, 5).tolist()
antib = np.arange(5, 101, 5).tolist()
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

print(df_heatmap_survival)

# Define the plot
fig, ax = plt.subplots(figsize=(13,7))

# Add title to the Heat map
title = "Probability of survival with random loading, gillespie_binary"

# Set the font size and the distance of the title from the plot
plt.title(title,fontsize=18)
ttl = ax.title
ttl.set_position([0.5,1.05])

# Use the heatmap function from the seaborn package
sns.heatmap(df_heatmap_survival, annot=True)
# Display the Pharma Sector Heatmap
plt.xlabel('antib')
plt.ylabel('initial N')
plt.show()

#for each simulation out of the total nr of simulations and for each
# loading factor calculate the probability of survival



