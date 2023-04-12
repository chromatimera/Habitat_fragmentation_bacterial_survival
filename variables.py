### To do list:

## 1. do countTotalMass for all AB conc and add it to survival fracion - get survival fraction and total mass at the same time,
##do the for loop ptr
## 2. N(t), survival fraction for deterministic case look weird
## why changing Nsat changes results??? -

## set all variables here and settings
import os.path
from decimal import *
import math
getcontext().prec = 50

### nr of droplets and the power of i and j in the partitioning loop are related i.e. for 100 droplets, the loop goes from 0 to 3.
total_drop_nr = 1000
## don't change the 2 lines below
part_min = 0
part_max = int(math.log10(total_drop_nr) + 1)

t_start = 0
t_end = 300
dt = 1
spec_time = 240
total_sim = 3

### nr timesteps for dt test; ignore dt for this test
nr_timesteps = int(3e6)

# n_crit used in one of the tau_leaping algorithms
n_crit = 2
volume = 1e-7   # volume of droplets ~100pL; 1pL is 1e-6 ul; 100pL - 1e-4 ul UNITS: mL

##simulating multiple fake droplets (this a factor ie. 2 means 2 * initial N, but total_droplet_nr/2)
nr_drop_min = 1
nr_drop_max = 10
step_drop = 5

Nsat = 1e8
initialN = 10
growthrate = 0.01 # per minute from experimental data Nia thesis
deathrate  = 0.045  # per minute from Gore 2013

#AB_molar_mass = 349.406 #g/mol (ug/umol)
MIC = 1 # ug/mL
AB_conc = 10 #ug/mL

### Mikaelis Menten parameters
Km = 6.7  #UNITS: ug/mL
Vmax = 3.5e-8 #ug/cell/min


##ab conc values range for testing ab
abmin = 0
abmax = 50
step = 1

## epsilon parameter for tau precheck
epsilon = 0.03

##type of loading and growth
loading = "det"  # rand #det
#growth = "midpoint_tau_binary" # for troubleshooting -- needs updating
growth = "binary"
#growth = 'gillespie_binary'
degradation = 'MM_exponential'



variables_script_path = __file__
variables_script_name = os.path.basename(__file__)