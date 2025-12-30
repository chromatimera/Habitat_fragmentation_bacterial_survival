import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from collections import Counter

from ps_theory import vol_fac

BIGGER_SIZE = 25

#### FIGURE 4 SURVIVAL PROBABILITY OF A POPULATION OF BACTERIA

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

#rootdir = './output/'
ab = [15, 35, 55, 75]
rho = 5e7 ## this is the initial density of cells/mL; for sim starting with lamda = 5; change accordingly

zz=np.load('prob_line.npy')
#os.chdir(rootdir)
zzz= zz.T
#os.chdir(rootdir)
subvol_list = 1e-4 / vol_fac
rhoV = 5E7 * subvol_list
a = [x - MIC for x in ab]
a_m=[x / MIC for x in ab]
F=(a )+6.7* np.log (a_m)
rhoT=(deathrate/ Vmax)*F

print('current dir', os.getcwd())


fig = plt.figure(figsize=(10,10))
ax1 = fig.add_subplot(111)
color = iter(plt.cm.rainbow(np.linspace(0, 1, 5)))
color_list = []
label_list = []
print(color)
g=0
for antib, c, ind in zip(ab, color, range(len(ab))):
    print('ab conc', antib)

    if ind == 2:
        c = next(color)
    os.chdir('dropnr_1000_loading_rand_growth_{}_initialN_5_abconc_{}'.format(growth, antib))
    path = os.getcwd()
    #print(path)
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    #print(onlyfiles)

    if '.DS_Store' in onlyfiles:
        onlyfiles.remove('.DS_Store')
    else:
        pass

    surv_fraction = pd.read_csv(onlyfiles[3])
    #print('surf fraction df', surv_fraction)
    part_fact = np.loadtxt(onlyfiles[2])

    theory_line_df = pd.DataFrame(zzz[:, ind], columns=['big_Ps'], index=vol_fac)
    theory_line_df.index.name = 'Vol_fac'
    theory_line_df = theory_line_df.sort_values(by="Vol_fac", ascending=True)
    theory_line_df['M'] = theory_line_df.index.astype(int)
    theory_line_df['RhoV'] = theory_line_df.apply(lambda x: rho * 1e-4 / x['M'], axis=1)
    theory_line_df['logPs'] = theory_line_df.apply(lambda x:np.log(x["big_Ps"]), axis=1)
    theory_line_df['subvol'] = theory_line_df.apply(lambda x: 1e-4/x['M'], axis=1)

    #theory_line_df['f(rho,rho*)RV'] = theory_line_df.apply(lambda x: (rhoT[g] / 5E7 - (rhoT[g] / 5E7) * np.log(1 + ((rhoT[g] - 5E7) / 5E7)) - 1) * x['RhoV'], axis=1)
    theory_line_df['f(rho,rho*)RV'] = theory_line_df.apply(lambda x: rhoT[g]*x['subvol'] * (1+np.log(5E7/rhoT[g]) -5E7/rhoT[g]) +np.log(1E-4/x['subvol'])-0.5*np.log(2*math.pi *  rhoT[g] *x['subvol']), axis=1)

    ### transpose of dataframe
    surv_fraction_transpose = surv_fraction.T
    surv_fraction_transpose.index.name = 'Part_fact'

    surv_fraction_transpose.columns = ['Surv frac']
    surv_fraction_transpose['M'] = surv_fraction_transpose.index.astype(int)
    surv_fraction_transpose['subvol'] = surv_fraction_transpose.apply(lambda x: 1e-4/x['M'], axis=1)

    surv_fraction_transpose['RhoV'] = surv_fraction_transpose.apply(lambda x: rho * 1e-4 / x['M'], axis=1)
    surv_fraction_transpose['logPs'] = surv_fraction_transpose.apply(lambda x:np.log(x["Surv frac"]), axis=1)
    #Original (email calc);
   # surv_fraction_transpose['f(rho,rho*)RV'] = surv_fraction_transpose.apply(lambda x: (rhoT[g] / 5E7 - (rhoT[g] / 5E7) * np.log(1 + ((rhoT[g] - 5E7) / 5E7)) - 1) * x['RhoV'], axis=1)
    #email calc minus 3/2log factor;
    #surv_fraction_transpose['f(rho,rho*)RV'] = surv_fraction_transpose.apply(lambda x: (rhoT[g] / 5E7 - (rhoT[g] / 5E7) * np.log(1 + ((rhoT[g] - 5E7) / 5E7)) - 1) * x['RhoV'] -np.log(x['subvol'])*(3/2), axis=1)
    #Current paper calc;
    #surv_fraction_transpose['f(rho,rho*)RV'] = surv_fraction_transpose.apply(lambda x: rhoT[g]*x['subvol'] * (1+np.log(5E7/rhoT[g]) -5E7/rhoT[g]) -np.log(x['subvol'])*(3/2), axis=1)
#new eq 06/03/24
    surv_fraction_transpose['f(rho,rho*)RV'] = surv_fraction_transpose.apply(lambda x: rhoT[g]*x['subvol'] * (1+np.log(5e7/rhoT[g]) -5e7/rhoT[g]) +np.log(1e-4/x['subvol'])-0.5*np.log(2*math.pi *  rhoT[g] *x['subvol']), axis=1)

    g=g+1

    surv_fraction_transpose['Error95'] = surv_fraction_transpose.apply(lambda x: 2 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error95_log'] = surv_fraction_transpose.apply(lambda x: abs(x['Error95']/x['Surv frac']), axis=1)
    surv_fraction_transpose['Error95_SEM'] = surv_fraction_transpose.apply(lambda x: 2 *np.std(x['Surv frac'] )/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_transpose['Error99'] = surv_fraction_transpose.apply(lambda x: 2.6 * math.sqrt(x['Surv frac'] * (1 - x['Surv frac']))/ math.sqrt(variables.total_sim), axis=1)
    surv_fraction_errors = surv_fraction_transpose.Error95.to_frame('Surv frac')
    surv_fraction_errors.index = surv_fraction_errors.index.map(int)
    surv_fraction_errors.index = surv_fraction_transpose['RhoV']
    #print('errors',surv_fraction_errors)
    surv_fraction_transpose.index = surv_fraction_transpose.index.map(int)
    surv_fraction_transpose = surv_fraction_transpose.set_index('RhoV', drop=True)
    theory_line_df = theory_line_df.set_index('RhoV', drop=True)

    plt.figure(1)
    theory_line_df["big_Ps"].plot.line(ax=ax1, c=c, linestyle='dashed', label='_nolegend_', logx=True)#, color = 'orange'
    surv_fraction_transpose["Surv frac"].plot.line(ax=ax1, marker='o',    linestyle='None',yerr=surv_fraction_errors, c=c, logx=True)#, color = 'orange')
    label_list.append('{}'.format(antib))

    #fig 2
    logPs=np.log(theory_line_df)
    logPs_sim=np.log(surv_fraction_transpose)

    #plt.figure(2)
    #logPs.plot.line(y='big_Ps', c=c)
   # plt.plot( rhoV,logPs, c=c)
    #logPs['big_Ps'].plot.line( c=c,linestyle='dashed') #theory
    #surv_fraction_transpose['logPs'].plot.line(marker='o',linestyle='None' ,c=c)
#=========================
    theory_line_df = theory_line_df.set_index('f(rho,rho*)RV', drop=True)
    surv_fraction_transpose = surv_fraction_transpose.set_index('f(rho,rho*)RV', drop=True)

    #plt.figure(3)
    #theory_line_df['logPs'].plot.line( style='--', c=c,label='_nolegend_')
    #surv_fraction_transpose['logPs'].plot.line(marker='o', c=c, linestyle='None',yerr=surv_fraction_transpose['Error95'])

    plt.figure(4)
    surv_fraction_transpose['logPs'].plot.line(marker='o', c=c, linestyle='None',yerr=surv_fraction_transpose['Error95_log'])
 #   theory_line_df['logPs'].plot( style='--',c=c)  #x=theory_line_df['f(rho,rho*)RV']


    surv_fraction_transpose.to_csv('../Survival_fraction_transpose_{}'.format(antib))
    os.chdir('..')

ax1.set_xlim(10**(2.3), 10**(0.6))
#plt.grid(True)
xticks = [200, 100, 50, 25, 10, 5]
xticks_label = [5000, 100, 50, 25, 10, 5]
second_ticks = [1, 50, 100, 200, 500, 1000]
ax1.set_xticks(xticks)
ax1.set_xticklabels(xticks_label)

ax2 = ax1.secondary_xaxis("top")
ax2.xaxis.set_ticks(xticks[::-1], labels=second_ticks[::-1])

plt.figure(1)
plt.xlabel(r'm (number of subvolumes)')
ax1.set_xlabel(r'$\rho$v (average number of cells per subvolume)')
ax1.yaxis.tick_right()
ax1.set_ylabel(r'$P_{s}$')
ax1.yaxis.set_label_position("right")
ax2.set_xlabel(r'm (number of subvolumes)')
plt.xlim([100,5])
ax1.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/ml}', loc='upper center', bbox_to_anchor=(0.5, 1.55), ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE)
plt.tight_layout()
plt.savefig('Survival fraction {} y vs logx square'.format(growth), dpi=600)
plt.savefig('Survival fraction {} y vs logx square.svg'.format(growth),format='svg')

#plt.figure(2)
#plt.xlim([1,100])
#plt.ylim([-30,1])
#plt.tight_layout()
#plt.ylabel('Log (Ps) /// log(sim_surv_fract)')
#plt.xlabel(r'\bf{$\rho$v (number of cells in droplet)}')
#plt.gca().invert_xaxis()
#plt.savefig('Log_Ps_plus sim'.format(growth), dpi=600)


#plt.figure(3)
#plt.ylabel('Log ($P_s$)')  #/// log(sim_surv_fract)
#plt.xlabel(r'\bf{$\rho$v f($\rho$ , $\rho$*)}')
#plt.xlim([-15,0.5])
#plt.ylim([-15,0.5])
#plt.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}',  loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE-5)
#plt.tight_layout()
#plt.savefig('Log_Ps_V_f(rho,rhoT) plus sim_zoom', dpi=600)

#
# plt.figure(4)
# plt.axline((0, 0), slope=1,linestyle='--' , label='_nolegend_', c='k')
# plt.ylabel('Log ($P_s$)')  #/// log(sim_surv_fract)
# plt.xlabel(r'\bf{f($v$, $\rho$*)}')
# plt.xlim([-12,0.5])
# plt.ylim([-12,0.5])
# plt.legend(label_list,title=r'\bf{Antibiotic concentration in $\mu$g/mL}',ncol=4,prop={'size':10},loc='lower right')
# #plt.legend(label_list, title=r'\bf{Antibiotic concentration in $\mu$g/mL}',  loc='upper center', bbox_to_anchor=(0.5, 1.4),ncol=4, fancybox=True, shadow=True, title_fontsize=BIGGER_SIZE-5, borderpad=0.5, labelspacing=0.55)
# plt.tight_layout()
# plt.savefig('Log_Ps_V_f(v,rho)_onlysim', dpi=600)
# plt.show()
# fig 3;
# plt.plot(rhoV, logPs_sim,marker='o',linestyle='None' ,c=c)
# logPs['RhoV'] = logPs.apply(lambda x: rho * 1e-4 / x['M'], axis=1)
# logPs['f(rho,rho*)RV'] = logPs.apply(lambda x: (rhoT[g] / 5E7 - (rhoT[g] / 5E7) * np.log(1 + ((rhoT[g] - 5E7) / 5E7)) - 1) * x['RhoV'], axis=1)
# logPs = logPs.set_index('f(rho,rho*)RV', drop=True)  #new xaxis
# logPs_sim = logPs_sim.set_index('f(rho,rho*)RV', drop=True)
# logPs_sim['RhoV'] = logPs_sim.apply(lambda x: rho * 1e-4 / x['M'], axis=1)
# logPs_sim['f(rho,rho*)RV'] = logPs_sim.apply(lambda x: (rhoT[g] / 5E7 - (rhoT[g] / 5E7) * np.log(1 + ((rhoT[g] - 5E7) / 5E7)) - 1) * x['RhoV'], axis=1)
# logPs_sim = logPs_sim.set_index('f(rho,rho*)RV', drop=True)  #new xaxis