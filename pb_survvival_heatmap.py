import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt


growth = 'gillespie_binary'
antib = 35
total_sim = 100
antib = [35,55,75]
os.chdir('./output/')
print(os.getcwd())
#for ab in antib:
os.chdir('./changing lamda for gillespie_binary/')
print(os.getcwd())
lambda_list = list(pd.read_csv('Lambda.txt', header=None))
print(lambda_list)
colnames=['lambda']
df_heatmap_survival=pd.read_csv('Lambda.txt', header=None, names=colnames)
print(df_heatmap_survival)
df_heatmap_survival.set_index('lambda', inplace=True)
for l in lambda_list:
    for ab in antib:
        print(os.getcwd())

        os.chdir('./dropnr_1000_loading_rand_growth_gillespie_binary_initialN_5_abconc_35/')
        df = pd.read_csv('Prob_survival.csv')
        print(df)

        df_heatmap_survival['ab {}'.format(ab)]=df.at[df.index[3],'Duration']
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



