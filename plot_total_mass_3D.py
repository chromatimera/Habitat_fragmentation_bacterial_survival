import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

import variables
from variables import *
from numpy import loadtxt


def read_csv():
    df = pd.read_csv('./output/df_growth_{}_starting_nr_drops_{}.csv'.format(growth, variables.total_drop_nr))
    return df

print(os.getcwd())
df = read_csv()
print(df)

part_fact = loadtxt("output/part_fact.txt", delimiter=",", unpack=False)
print(part_fact)
print(len(part_fact))

ab_conc = loadtxt("output/ab_conc.txt", delimiter=",", unpack=False)
print(ab_conc)
print(len(ab_conc))

nf_diff_ab = np.zeros((len(part_fact),len(ab_conc)))
for j in range(len(ab_conc)):
    for i in range(len(part_fact)):
        nf_diff_ab[i][j] = df.iloc[-1, j*4 + i + 1]

Z = nf_diff_ab

X, Y = np.meshgrid(ab_conc, part_fact)


fig,ax=plt.subplots(1,1)
cp = ax.contourf(X, Y, Z)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
#ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
plt.show()


fig = plt.figure(figsize=(6.4*2.5,4.8*2.5))
ax = plt.axes(projection='3d')
ax.view_init(25,25)
ax.plot_surface(X, Y, nf_diff_ab, rstride=1, cstride=1, cmap='viridis', edgecolor='none', linewidth=0, antialiased=False)
#ax.set_title("Fraction of atoms in molecules vs P & T", fontsize='large', fontweight='bold')
ax.set_xlabel('Antibiotic concentration', fontsize='large', fontweight='bold')
ax.set_ylabel('Partitioning factors', fontsize='large', fontweight='bold')
ax.set_zlabel("Final nr of bacteria", fontsize='large', fontweight='bold')
ax.grid()
#fig.colorbar(surf, shrink=0.3, aspect=5)
#plt.savefig('3D-xTP-low_mmm_delta_1_2_2_0_logP.png')
plt.show()

### Plot p-T projection of phase diagram - x colour coded
plt.figure(figsize=(6.4*2.5,4.8*2.5))
plt.xlabel('Pressure', fontsize=text_size, fontweight='bold')
plt.ylabel('Temperature', fontsize=text_size, fontweight='bold')
plt.contourf(X, Y, nf_diff_ab, 20, rstride=1,  cstride=1, cmap='viridis', linewidth=0, antialiased=False)
#im = plt.imshow(x_array_bdu, cmap='viridis', origin='lower')
cbar = plt.colorbar(orientation='vertical',ticks=[0,0.2,0.4,0.6,0.8,1])
cbar.set_label('x (fraction of atoms in molecules)',size=18)
plt.show()


### PLS IGNORE THIS BIT - FOR ADDING ERROR BARS
#adaptive_tau_binary_transposed['Error95'] = adaptive_tau_binary_transposed.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
#adaptive_tau_binary_transposed['Error99'] = adaptive_tau_binary_transposed.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(nr_drops), axis=1)
#


### PLOT FOR 2D DATAFRAME; if dataframe is associated with one antibiotic concentration
## to plot Nf, Ni  as a function of partition factors
df.iloc[0, 1:].plot()  ### plot initial nr of bacteria
df.iloc[-1, 1:].plot() ### plot final nr of bacteria
plt.grid(True)
plt.title('Total mass versus partitioning factor, ab conc {}'.format(AB_conc))
plt.ylabel('Total mass (nr of bacteria) Nf, Ni')
plt.xlabel('Partition factor')
plt.legend(['Ni', 'Nf'], loc='upper left')
plt.savefig('./output/Nf_Ni_vs_part_fact_startin_nr_of_drops_{}.png'.format(variables.total_drop_nr))
plt.show()