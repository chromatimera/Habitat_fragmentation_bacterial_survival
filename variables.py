## set all variables here and settings
import os.path
from decimal import *
import numpy as np

getcontext().prec = 50

### nr of droplets and the power of i and j in the partitioning loop are related i.e. for 100 droplets, the loop goes from 0 to 3.
#droplet_list = [1]
droplet_list = np.arange(0, 1001, 5)
droplet_list[0] = 1

# setting up time values
t_start = 0
t_end = 300
dt = 1
spec_time = 299  # this is a specific time at which we calculate the number of bacteria left alive; it is set to 299 as it is the last timepoint in the for loop
total_sim = 1 # number of simulation repeats
step = 1


### nr timesteps for dt test; ignore dt for this test
nr_timesteps = int(3e6)

# n_crit used in one of the tau_leaping algorithms
n_crit = 2
volume = 1e-7   # volume of droplets from experiments ~100pL; 1pL is 1e-6 ul; 100pL - 1e-4 ul UNITS: mL

Nsat = 1e8
initialN = 5
growthrate = 0.01 # per minute from experimental data Nia thesis
deathrate  = 0.045  # per minute from Gore 2013
slowrate =0.0 #per min; for resource model

#AB_molar_mass = 349.406 #g/mol (ug/umol)
MIC = 1 # ug/mL
AB_conc = 30 #ug per mL

### Mikaelis Menten parameters
Km = 6.7  #UNITS: ug/mL
Vmax = 3.5e-8 #ug/cell/min


## epsilon parameter for tau precheck
epsilon = 0.03

##type of loading and growth
loading = "rand"  # rand #det
#growth = "binary"
growth = 'binary'
degradation = 'MM_exponential'


# Parameters for plotting
text_size = 'large'

variables_script_path = __file__
variables_script_name = os.path.basename(__file__)
