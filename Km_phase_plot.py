import numpy as np
import matplotlib.pyplot as plt
import variables
from variables import *

#Equation 9: rhoT

medkm = 6.7  #UNITS: ug/mL
lowkm = 1
highkm = 15
AB_x = np.arange(0, 25, 1)
#calculate rho values;
F= (AB_x-MIC )+ medkm* np.log (AB_x/MIC)
rhoT=    (deathrate/ Vmax) *F
rhoT_low=((AB_x-MIC )+ lowkm* np.log (AB_x/MIC) )* (deathrate/ Vmax)
rhoT_high= ((AB_x-MIC )+ highkm* np.log (AB_x/MIC)) * (deathrate/ Vmax)

#



y =rhoT
x=AB_x
(ax1, ax2) = plt.subplots(sharex=True, sharey=True)

plt.plot(rhoT_high, linewidth=4, label="Large KM")
plt.plot(rhoT_low, linewidth=4,label="Small KM")
ax2.fill_between(x, y, facecolor='red')
ax2.fill_between(x, y, plt.ylim()[1], facecolor='green')


xcoord = x[int((x.size/3)*2)]
ycoord = y[int((x.size/4)*2)] / 2
plt.text(xcoord,ycoord,"Death")
plt.text(2,5e7,"Survival")
plt.xlabel('Initial AB (ug/mL)')
plt.ylabel(r'$\rho$ (initial cells per mL)' )

plt.legend()
plt.show()