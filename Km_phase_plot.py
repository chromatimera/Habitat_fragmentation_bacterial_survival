import numpy as np
import matplotlib.pyplot as plt
import variables
from variables import *

#Equation 9: rhoT
plt.rc('font', size=14)  # controls default text size
#plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rc('text', usetex=True)


medkm = 6.7  #UNITS: ug/mL
lowkm = 1
highkm = 15
AB_x = np.arange(0, 25,0.5)
#calculate rho values;
F= (AB_x-MIC )+ medkm* np.log (AB_x/MIC)
rhoT=    (deathrate/ Vmax) *F
rhoT_low=((AB_x-MIC )+ lowkm* np.log (AB_x/MIC) )* (deathrate/ Vmax)
rhoT_high= ((AB_x-MIC )+ highkm* np.log (AB_x/MIC)) * (deathrate/ Vmax)

#



y =rhoT
y[0]=0
x=AB_x
(ax1, ax2) = plt.subplots(sharex=True, sharey=True)

plt.plot(x,rhoT_high,'--', linewidth=4, label="Large KM", color='black')
plt.plot(x,rhoT_low,'.', linewidth=4,label="Small KM", color='black')
plt.plot(x,rhoT, linewidth=4,label="Small KM", color='black')
ax2.fill_between(x, y, facecolor='red')
ax2.fill_between(x, y, plt.ylim()[1], facecolor='green')


xcoord = x[int((x.size/2.5)*2)]
ycoord = y[int((x.size/6)*2)] / 2
plt.text(xcoord,ycoord,"DEATH", weight='bold')
plt.text(2,8e7,"SURVIVAL", weight='bold')
plt.xlabel('Initial antibiotic concentration ($\u03BC$g/mL)')
plt.ylabel(r'$\rho$ (initial cells per mL)' )
plt.xlim(0,24)
plt.ylim(0,0.9e8)
#plt.legend()
plt.show()