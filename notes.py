import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
X = np.arange(0, self.t_end, self.dt).tolist()
XX= np.tile(X, (self.total_drop_number, 1))

fig, ax= plt.subplots(figsize=(8,6))
plt.plot (XX.T, self.N_r_array.T)
plt.grid(False)
#plt.title('Growth of resistant strain')
#plt.title ('{}_growth'.format(growth))
plt.ylabel(r'\bf{Number of bacteria}')
plt.xlabel(r'\bf{Time (min)}')
plt.xlim(0,self.t_end)
plt.xlim(0,self.t_end)
#ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))
plt.ylim(bottom=0)
plt.fill_between(x, y1, y2, where=(y1 > y2), color='C0', alpha=0.3,
                 interpolate=True)
plt.savefig('plot_Nbact_loading_{}_growth_{}ab_conc_{}'.format(loading, growth, AB_conc))
plt.show()


plt.figure(2, figsize=(8,6))
plt.plot(XX.T, self.AB_conc_array.T)
plt.grid(False)
#plt.title('Concentration of antibiotic over time in each droplet.')
plt.ylabel(r'\bf{Antibiotic ($\mu$g/mL)}')
plt.xlabel(r'\bf{Time (min)}')
plt.xlim(0, self.t_end)
plt.ylim(bottom=0)
#plt.title ('{}_growth'.format(growth))
plt.savefig('plot_ABconc_{}_loading_{}_growth_{}'.format(self.AB_conc, loading, growth))
plt.show()
