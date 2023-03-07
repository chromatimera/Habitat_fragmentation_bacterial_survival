import math
import numpy as np
from scipy.special import lambertw
from variables import *

## set all variables here and settings
import os.path
from decimal import *

getcontext().prec = 50

total_drop_nr = 5

## not sure about the ime units

t_start = 0
t_end = 300
dt = 1

n_crit = 2
## nr of drops to simulate for total mass
nr_drops_total_mass = 1
volume = 1e-4 * nr_drops_total_mass   # volume of droplets ~100pL; 1pL is 1e-6 ul; UNITS: ul

#check_total_mass = False


Nsat = 100 * nr_drops_total_mass

initialN = 2 * nr_drops_total_mass
growthrate = 0.01 # per minute
deathrate =  0.045  # per minute

AB_molar_mass = 349.406 #g/mol (ug/umol)
# Starting from ug/mL (per droplet) -> #converts MIC into UNITS: umoles/ul
MIC = 1 /AB_molar_mass * 1e-3
AB_conc = 6 /AB_molar_mass * 1e-3 #ug/ml (per droplet)#\frac{ug}/{ml}/\frac{ug}{umole}-> UNITS: umoles/ul

### Mikaelis Menten parameters
#Km given as 6.7e-3 μg/uL -> transformed to umoles/uL
Km = 6.7e-3 /AB_molar_mass  #UNITS: umoles/uL
#Vmax given as 6e-4 μg/cell/min
Vmax = 6e-4 /AB_molar_mass  #UNITS: ug/cell/min

##ab conc values range for testing ab
abmin = 0
abmax = 50
step = 1

### drop number for testing convergence or surv frac
dropmin = 1
dropmax = 5000
nr_points = 25

## epsilon parameter for tau precheck
epsilon = 0.03

##type of loading and growth
loading = "det"
#growth = "midpoint_tau_binary"
#growth = "adaptive_tau_binary"
growth = 'gillespie_binary'
degradation = 'MM_exponential'
deg_type = 'rate'

### r timesteps for dt test
nr_timesteps = 400

print('Vmax', Vmax)
print('Km', Km)
print('volume', volume)

variables_script_path = __file__
variables_script_name = os.path.basename(__file__)



