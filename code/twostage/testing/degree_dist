#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:48:17 2023

@author: jsheen
"""

# Overdispersed
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
mean_degree = 15
k_overdispersion = 1
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)

N_cluster = 1000
z = []
for i in range(N_cluster):
    deg = 0
    while (deg == 0):
        deg = np.random.negative_binomial(k_overdispersion, p)
    z.append(deg)
if (sum(z) % 2 == 0):
    continue_loop = False
    
plt.hist(z, bins=20, weights=np.ones(len(z)) / len(z))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Degree', fontsize=18)
plt.ylabel('Percentage', fontsize=18)

#Poisson
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
mean_degree = 15

N_cluster = 10000
z = []
for i in range(N_cluster):
    deg = 0
    while (deg == 0):
        deg = np.random.poisson(mean_degree)
    z.append(deg)
if (sum(z) % 2 == 0):
    continue_loop = False
    
plt.hist(z, bins=15, weights=np.ones(len(z)) / len(z))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.xlabel('Degree', fontsize=18)
plt.ylabel('Percentage', fontsize=18)