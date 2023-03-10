## set all variables here and settings
import os.path
from decimal import *

getcontext().prec = 50

total_drop_nr = 1

t_start = 0
t_end = 300
dt = 0.001

### nr timesteps for dt test; ignore dt for this test
nr_timesteps = int(3e6)

# n_crit used in one of the tau_leaping algorithms
n_crit = 2

volume = 1e-7   # volume of droplets ~100pL; 1pL is 1e-6 ul; 100pL - 1e-4 ul UNITS: mL

##simulating multiple fake droplets (this a factor ie. 2 means 2 * initial N, but total_droplet_nr/2)
nr_drop_min = 1
nr_drop_max = 10
step_drop = 5

Nsat = 100
initialN = 10
growthrate = 0.01 # per minute from experimental data Nia thesis
deathrate  = 0.045  # per minute

#AB_molar_mass = 349.406 #g/mol (ug/umol)
MIC = 1 # ug/mL
AB_conc = 4 #ug/mL

### Mikaelis Menten parameters
Km = 6.7  #UNITS: ug/mL
Vmax = 6.0e-4 #ug/cell/min


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
growth = "binary"
#growth = 'gillespie_binary'
degradation = 'MM_linear'



variables_script_path = __file__
variables_script_name = os.path.basename(__file__)