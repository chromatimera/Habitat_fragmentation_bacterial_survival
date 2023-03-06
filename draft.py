import math
import numpy as np
from scipy.special import lambertw
from variables import *

N = np.linspace(0, 11, 12)
print(N)

t = 1
y = -Vmax * t
w = Km * lambertw(math.exp(y/Km)/Km)
print(w.real)

tau = [0,1,2,3,4,5,6,7,8,9]
print(sum(tau))