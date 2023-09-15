import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

growth = 'binary'
total_sim = 50

Ni = [5,10,15]
antib = np.arange(5, 101, 5).tolist()
print(antib)


os.chdir('./output/')
print(os.getcwd())
#for ab in antib:
#os.chdir('./changing lamda for gillespie_binary/')
##print(os.getcwd())

colnames=['lambda']
df_heatmap_survival=pd.DataFrame(Ni)
df_heatmap_survival.columns=colnames
print(df_heatmap_survival)

df_heatmap_survival.set_index('lambda', inplace=True)
for l in Ni:
    for ab in antib:
        print(os.getcwd())

        os.chdir('./dropnr_1000_loading_rand_growth_binary_initialN_{}_abconc_{}/'.format(l, ab))
        ### calculate probability of survival
        #df = pd.read_csv('Prob_survival.csv')
        #print(df)

        df_heatmap_survival['ab {}'.format(ab)]=df.at[df.index[3],'Surv frac']
        os.chdir('..')

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
plt.show()

#for each simulation out of the total nr of simulations and for each
# loading factor calculate the probability of survival



