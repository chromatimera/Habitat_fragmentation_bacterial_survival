### To do list:
### When running for loop different partitions, the bacteria seem to grow - I believe there is something wrong in either degradation or volume or something else
##  deg is faster for bigger droplets-- why -- partially because nsat too low (make a concentration rather than fixed value?)
#check TOTAL initial bact (this should be the same for all loops in det loading)-- plot Ni vs partition ( i think this is ok?)
#change units for AB in plots ;; ug /ml???
#  changed from 0 to 1  ;; if self.N < 1.0: for det growth  #### #can't have less than 1 bacteria ; is this the same in Gillespie?
#QUESTIONS
#can you explain i and j in the partition loop?
#How does poission loading work for the big droplet?
#==> ie what kind of repeats does your loop do
#Final N plots (with error bars from repeats)



## set all variables here and settings
import os.path
from decimal import *

getcontext().prec = 50

total_drop_nr = 500  ###

t_start = 0
t_end = 300
dt = 0.1

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
initialN = 5
growthrate = 0.01 # per minute from experimental data Nia thesis
deathrate  = 0.045  # per minute from Gore 2013

#AB_molar_mass = 349.406 #g/mol (ug/umol)
MIC = 1 # ug/mL
AB_conc = 25 #ug/mL

### Mikaelis Menten parameters
Km = 6.7  #UNITS: ug/mL
Vmax = 3.5e-8 #ug/cell/min


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
loading = "det"  # rand #det
#growth = "midpoint_tau_binary" # for troubleshooting -- needs updating
growth = "binary"
#growth = 'gillespie_binary'
#degradation = 'MM_linear'
degradation = 'MM_exponential'



variables_script_path = __file__
variables_script_name = os.path.basename(__file__)