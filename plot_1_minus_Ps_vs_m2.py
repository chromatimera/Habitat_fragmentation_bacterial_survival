import math
import variables
from os import listdir
from os.path import isfile, join
import pandas as pd
from variables import *
import matplotlib.pyplot as plt
from ps_theory import vol_fac

BIGGER_SIZE = 22

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

ab = [35, 55, 75 ]

zz=np.load('prob_line.npy')
zzz= zz.T

print('current dir', os.getcwd())

plt.figure(figsize=(11, 9))
color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
label_list = []
print(color)

for antib, c, ind in zip(ab, color, range(len(ab))):

    if ind == 2:
        c = next(color)

    theory_line_df = pd.DataFrame(zzz[:, ind], columns=['big_Ps'], index=vol_fac)
    theory_line_df.index.name = 'Vol_fac'
    theory_line_df = theory_line_df.sort_values(by="Vol_fac", ascending=True)
    print(theory_line_df)

    ## for plot of log (1-Ps) vs 1/m2

    one_minus_Ps = 1 - theory_line_df
    print('one minus ps first', one_minus_Ps)
    one_minus_Ps['m2']=one_minus_Ps.index **2

    one_minus_Ps['log']=np.log(one_minus_Ps['big_Ps'])
    log_list = list(one_minus_Ps['log'])

    for i in range(len(log_list)):
        if str(log_list[i]) == '-inf':
           log_list[i] = -50

    one_minus_Ps['log']=log_list
    print('1-Ps: ', one_minus_Ps)

    plt.figure(1)

    #theory_line_df["big_Ps"].plot.line(c=c, linestyle='dashed', label='_nolegend_')#, color = 'orange')
    label_list.append('{}'.format(antib))

    ## plot 1-Ps versus 1/m2
    plt.plot(one_minus_Ps['m2'], one_minus_Ps['log'])
    os.chdir('..')
    #print(os.getcwd())

plt.ylabel(r'\bf{Probability of survival}')
plt.xlabel(r'\bf{$m^{2}$ (number of subvolumes)}')

plt.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}', loc='upper center', bbox_to_anchor=(0.5, 1.17), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)
plt.savefig('1-ps'.format(growth), dpi=600)
plt.show()