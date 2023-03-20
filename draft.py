import numpy as np
from variables import *
import pandas as pd

x = pd.read_csv('./output/df_growth_gillespie_binary_starting_nr_drops_500.csv')
part_fact = [0, 0.0008, 0.0016, 0.004, 0.008, 0.02, 0.04, 0.1, 0.2, 0.5, 1.0]
pf = pd.DataFrame(part_fact)
print(pf)
x.loc[-1] = [0, 0.0008, 0.0016, 0.004, 0.008, 0.02, 0.04, 0.1, 0.2, 0.5, 1.0]  # adding a row
x.index = x.index + 1  # shifting index
x.sort_index(inplace=True)
print(x)