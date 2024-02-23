import numpy as np
import matplotlib.pyplot as plt
from variables import *
from matplotlib import rc
import matplotlib.pyplot as plt
#pylatex
BIGGER_SIZE = 28

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True) ## https://matplotlib.org/stable/tutorials/text/usetex.html
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize


color = plt.cm.rainbow(np.linspace(0, 1,5))
color_2 = plt.cm.viridis(np.linspace(0, 100,5))

color_list = []

colors = ['deeppink', color[0], 'blue', 'teal']
ind= 0

#Equation 9: rhoT
medkm = 6.7  #UNITS: ug/mL
lowkm = 1
highkm = 15
AB_x = np.arange(0, 95, 0.5)
#calculate rho values;
F=(AB_x-MIC )+medkm* np.log (AB_x/MIC)
rhoT=(deathrate/ Vmax)*F
rhoT_low=((AB_x-MIC )+ lowkm* np.log (AB_x/MIC) )* (deathrate/ Vmax)
rhoT_high= ((AB_x-MIC )+ highkm* np.log (AB_x/MIC)) * (deathrate/ Vmax)
#example calc;
RHOt_eg =((30-MIC )+ medkm* np.log (30/MIC) )* (deathrate/ Vmax)
print('{:.5E}'.format(RHOt_eg))
y =rhoT
y[0]=0
x=AB_x
(ax1, ax2) = plt.subplots(sharex=True, sharey=True)

#plt.plot(x,rhoT_high,'--', linewidth=4, label="Large KM", color='black')
#plt.plot(x,rhoT_low,'.', linewidth=4,label="Small KM", color='black')
plt.plot(x,rhoT, linewidth=4,label="Small KM", color='black')

plt.rcParams['hatch.color'] = 'g'
ax2.fill_between(x, y, facecolor='red', alpha=0.1)
ax2.fill_between(x, y, plt.ylim()[1], facecolor='green', hatch = '.', alpha=0.1 )

ab_conc = [10,15,20,30]

for antib, c, ind in zip(ab_conc, color, range(len(ab_conc))):

    plt.plot(antib, 5e7, marker='*', ls='none', color = colors[ind], ms=20)
    ind += 1



xcoord = x[int((x.size/2.5)*2)]
ycoord = y[int((x.size/6)*2)] / 2
plt.xlabel(r'$a_{init}$($\mu$g/ml)')
plt.ylabel(r'$\rho$ (initial cells/ml)')
plt.xlim(0,35)
plt.ylim(0,0.9e8)
plt.tight_layout()
plt.xticks([0,5,10,15,20,30,35])
plt.tight_layout()
plt.savefig('km_inset_figure_5.png', dpi=800)
plt.savefig('km_inset_figure_5.svg', format="svg",dpi=600)
plt.show()
