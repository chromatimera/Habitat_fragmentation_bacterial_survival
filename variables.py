## set all variables here and settings
import os.path
from decimal import *

getcontext().prec = 50

total_drop_nr = 100

t_start = 0
t_end = 300
dt = 1

# n_crit used in one of the tau_leaping algorithms
n_crit = 2

volume = 1e-7   # volume of droplets ~100pL; 1pL is 1e-6 ul; UNITS: ul

##simulating multiple fake droplets (this a factor ie. 2 means 2 * initial N, but total_droplet_nr/2)
nr_drop_min = 1
nr_drop_max = 10
step_drop = 5

Nsat = 1e2
initialN = 10
growthrate = 0.0066 # per minute from experimental data Nia thesis
deathrate =  0.045  # per minute

#AB_molar_mass = 349.406 #g/mol (ug/umol)
# Starting from ug/mL (per droplet) -> #converts MIC into UNITS: umoles/ul
MIC = 1 # ug/ml
AB_conc = 6 #ug/ml (per droplet)-> UNITS: umoles/ul

### Mikaelis Menten parameters
#Km given as 6.7e-3 μg/uL -> transformed to umoles/uL
Km = 6.7e-1  #UNITS: umoles/uL
#Vmax given as 6e-4 μg/cell/min
Vmax = 6e-4 #6e-4 ug/cell/min /AB_molar_mass  #UNITS: umoles/cell/min


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
degradation = 'MM_linear'
deg_type = 'rate'

### r timesteps for dt test
nr_timesteps = 400

variables_script_path = __file__
variables_script_name = os.path.basename(__file__)