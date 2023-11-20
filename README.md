# droplets_vs_bulk_survival_paper
 This code simulates the growth of beta-lactamase producing bacteria, inside droplets or in bulk and the associated survival rates.

There are 4 types of growth: binary - deterministic, Gillespie - stochastic simulation algorithm, becomes prohibitive as N is larger. Midpoint and adaptive Tau are two Tau-leaping algorithms, one proposed by Gillespie, one more advanced propose by Cao. They both map very similar to Gillespie, thus can be used to simulate large populations.

Notes on choice of volume:

In the variables.py file the volume of the big droplet is kept constant at 1e-4 mL, which translates to 0.1 uL. The  This big droplet is loaded with 5e3 bacteria which is consistent with experimental values of overnight E. coli cultures going up to a maximum of 1e6 - 1e7 CFU/uL. For simulating smaller droplets (i.e. 1000), the volume decreases to 100 pL, which is equal to 1e-4 uL and 1e-7 mL.








