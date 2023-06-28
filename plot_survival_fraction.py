import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm, rc
from ps_theory import vol_fac

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

growth = 'binary'
rootdir = 'output/'
ab = [15, 35, 55, 75]

#os.chdir(rootdir)
zz=np.load('prob_line.npy')
zzz= zz.T
#print('current dir', os.getcwd())


plt.figure(figsize=(7, 6))
color = iter(cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
for antib, c, ind in zip(ab, color, range(len(ab))):
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()

    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    print(onlyfiles)

    if '.DS_Store' in onlyfiles:
        onlyfiles.remove('.DS_Store')
    else:
        pass

    surv_fraction = pd.read_csv(onlyfiles[3])
    print('surf fraction df', surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])

    theory_line_df = pd.DataFrame(zzz[:, ind], columns=['big_Ps'], index=vol_fac)
    theory_line_df.index.name = 'Vol_fac'
    theory_line_df = theory_line_df.sort_values(by="Vol_fac", ascending=True)


    print('THEORY DF', theory_line_df)


    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'
    print('transpose', surv_fraction_transpose)

    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    print(surv_fraction_transpose)

    plt.figure(1)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_errors.index = surv_fraction_errors.index.map(int)
    #surv_fraction_errors = surv_fraction_errors.sort_index(ascending=True)
    print('errors',surv_fraction_errors)

    surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)
    #surv_fraction_transpose = surv_fraction_transpose.sort_index(ascending=True)
    print('trp', surv_fraction_transpose)
    surv_fraction_transpose["Surv frac"].plot.line(yerr=surv_fraction_errors, c=c)#, color = 'orange')
    theory_line_df["big_Ps"].plot.line(c=c, linestyle='dashed')#, color = 'orange')


    os.chdir('..')
    #print(os.getcwd())

plt.ylabel(r'\bf{Probability of survival}', fontsize=text_size)
plt.xlabel(r'\bf{m (number of subvolumes)}', fontsize=text_size)

plt.legend(ab, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, fancybox=True, shadow=True, fontsize= 'large')
plt.savefig('Survival fraction {} + errors diff ab+ legend.png'.format(growth))
plt.show()