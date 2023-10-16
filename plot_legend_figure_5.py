import pylab
import math
import variables
from os import listdir
import os
from os.path import isfile, join
import pandas as pd
import numpy as np
from variables import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
from collections import Counter
from matplotlib import rc


##### IGNORE THIS SCRIPT FOR NOW
BIGGER_SIZE = 22

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsiz


fig = pylab.figure()
figlegend = pylab.figure(figsize=(3,2))
ax = fig.add_subplot(111)
lines = ax.plot(range(10), pylab.randn(10), range(10), pylab.randn(10))
figlegend.legend(lines, (35, 55, 75), 'center')
fig.show()
figlegend.show()
figlegend.savefig('legend.png')