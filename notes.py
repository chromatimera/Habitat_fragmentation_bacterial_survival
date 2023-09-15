import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np



growth = 'gillespie_binary'
antib = 35
total_sim = 100

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
df_heatmap_survival.set_index('lambda', inplace=True)
print(df_heatmap_survival)