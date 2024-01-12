## set all variables here and settings
import os.path
from decimal import *
import numpy as np

getcontext().prec = 50

### nr of droplets and the power of i and j in the partitioning loop are related i.e. for 100 droplets, the loop goes from 0 to 3.
#droplet_list = [2000]
droplet_list = np.arange(0, 20001, 400)
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
#volume_exp = 1e-7   # volume of droplets from experiments ~100pL; 1pL is 1e-6 ul; 100pL - 1e-4 ul UNITS: mL

### The volume of the simulation is set at 1e-4 mL for the big droplet and it is kept constant no matter
# the nr of droplets simulated (i.e. if we simulate 10, 100 or 10 000 droplets); see README file for more info.
volume_big_drop = 1e-4  ## should be 1e-4
param_small_droplets = 1/droplet_list[-1]
volume = round(volume_big_drop * param_small_droplets, 12)

Nsat = 1e12

## Cchange N to be equivalent to 5 bact for 1000 droplets - to link with the other simulations
initialN = 0.25 # lambda value defined for the highest/max m value simulated
growthrate = 0.01 # per minute from experimental data Nia thesis
deathrate  = 0.045  # per minute from Gore 2013
slowrate =0.0 #per min; for resource model

#AB_molar_mass = 349.406 #g/mol (ug/umol)
MIC = 1 # ug/mL
AB_conc = 10 #ug/mL

### Mikaelis Menten parameters
Km = 6.7  #UNITS: ug/mL
Vmax = 3.5e-8 #ug/cell/min


## epsilon parameter for tau precheck
epsilon = 0.03

##type of loading and growth
loading = "rand"  # rand #det
growth = 'binary' #gillespie_binary
degradation = 'MM_exponential'


# Parameters for plotting
text_size = 'large'

variables_script_path = __file__
variables_script_name = os.path.basename(__file__)
