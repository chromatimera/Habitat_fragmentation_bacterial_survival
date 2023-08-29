import numpy as np
import matplotlib.pyplot as plt
from variables import *
from matplotlib import rc
import matplotlib.pyplot as plt
#pylatex
BIGGER_SIZE = 20

#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('text', usetex=True) ## https://matplotlib.org/stable/tutorials/text/usetex.html
#sudo apt-get install latexmk  ??
plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize

#Equation 9: rhoT
#plt.rc('font', size=18)  # controls default text size
#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('text', usetex=True)
#font = {'family':'serif', 'serif': ['computer modern roman']}
#plt.rc('font',**font)
#plt.rcParams["font.family"] = "serif"



#plt.rc('font',**{'family':'serif','serif':['Times']})
#plt.rc('text', usetex=True)

medkm = 6.7  #UNITS: ug/mL
lowkm = 1
highkm = 15
AB_x = np.arange(0, 25,0.5)
#calculate rho values;
F=(AB_x-MIC )+medkm* np.log (AB_x/MIC)
rhoT=(deathrate/ Vmax)*F
rhoT_low=((AB_x-MIC )+ lowkm* np.log (AB_x/MIC) )* (deathrate/ Vmax)
rhoT_high= ((AB_x-MIC )+ highkm* np.log (AB_x/MIC)) * (deathrate/ Vmax)

y =rhoT
y[0]=0
x=AB_x
(ax1, ax2) = plt.subplots(sharex=True, sharey=True)

plt.plot(x,rhoT_high,'--', linewidth=4, label="Large KM", color='black')
plt.plot(x,rhoT_low,'.', linewidth=4,label="Small KM", color='black')
plt.plot(x,rhoT, linewidth=4,label="Small KM", color='black')

plt.rcParams['hatch.color'] = 'r'
ax2.fill_between(x, y, facecolor='red', alpha=0.2, hatch="//")
ax2.fill_between(x, y, plt.ylim()[1], facecolor='green', alpha=0.2 ) #hatch="."  //,hatch="o" +


xcoord = x[int((x.size/2.5)*2)]
ycoord = y[int((x.size/6)*2)] / 2
#plt.text(xcoord,ycoord,"DEATH", weight='bold')
#plt.text(2,8e7,"SURVIVAL", weight='bold')
plt.xlabel(r'\bf{Initial antibiotic concentration ($\mu$g/ml)}')
plt.ylabel(r'\bf{$\rho$ (initial cells/ml)}')
plt.xlim(0,24)
plt.ylim(0,0.9e8)
plt.tight_layout()
plt.show()
