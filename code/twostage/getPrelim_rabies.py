#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 09:33:36 2022

@author: jsheen

@description: script used to get preliminary data: 
              1. day of intervention
              2. beta (transmission rate) of R0
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
#import matplotlib.pyplot as plt
import statistics
import multiprocessing as mp
import networkx as nx
from collections import defaultdict, Counter
import EoN
from pathlib import Path
home = str(Path.home())

# Set parameter sets ----------------------------------------------------------
Ns = [1000, 10000]
overdispersions = [0.1, 0.7, 1]
R0s = [0.25, 1.1]
morts = [0.85]
param_sets = []
for i in Ns:
    for j in overdispersions:
        for k in R0s:
            for l in morts:
                param_sets.append([i, j, k, l])
                   
# Function to find day of intervention and transmission rates for R0s ---------
def getPrelim(param_set):
    N_cluster = param_set[0]
    k_overdispersion = param_set[1]
    R0 = param_set[2]
    mort = param_set[3]
    eit = 0.005
    mean_degree = 15
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    if N_cluster == 1000:
        initial_infections_per_cluster = 4
    else:
        initial_infections_per_cluster = 40
    ave_inc_period = 5
    ave_inf_period_recover = 10
    ave_inf_period_dead = 4
    ave_inf_period_rate = (mort * (1 / ave_inf_period_dead)) + ((1-mort) * (1 / ave_inf_period_recover))
    if eit <= initial_infections_per_cluster / N_cluster:
        raise NameError("Script assumes expected It / N strictly < initial infections per cluster / N.")
    # Joel C Miller's methods to estimate_R0 ----------------------------------
    def get_Pk(G):
        Nk = Counter(dict(G.degree()).values())
        Pk = {x:Nk[x]/float(G.order()) for x in Nk.keys()}
        return Pk
    def get_PGFPrime(Pk):
        maxk = max(Pk.keys())
        ks = np.linspace(0,maxk, maxk+1)
        Pkarray = np.array([Pk.get(k,0) for k in ks])
        return lambda x: Pkarray.dot(ks*x**(ks-1))
    def get_PGFDPrime(Pk):
        maxk = max(Pk.keys())
        ks = np.linspace(0,maxk, maxk+1)
        Pkarray = np.array([Pk.get(k,0) for k in ks])
        return lambda x: Pkarray.dot(ks*(ks-1)*x**(ks-2))
    def estimate_R0(G, tau = None, gamma = None):
        transmissibility = tau/(tau+gamma)
        Pk = get_Pk(G)
        psiDPrime = get_PGFDPrime(Pk)
        psiPrime = get_PGFPrime(Pk)
        return transmissibility * psiDPrime(1.)/psiPrime(1.)
    # Find median beta that leads to desired R0 (2000 sims) ----------------
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    beta_lst = []
    for i in range(2000):
        #if (i % 100 == 0):
            #print(i)
        continue_loop = True
        while (continue_loop):
            z = []
            for i in range(N_cluster):
                deg = 0
                deg = np.random.negative_binomial(k_overdispersion, p)
                z.append(deg)
            if (sum(z) % 2 == 0):
                continue_loop = False
        G=nx.configuration_model(z)
        G=nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))
        est_R0=3.3
        beta=0.04
        while est_R0 > R0:
            beta = beta - 0.0001
            est_R0 = estimate_R0(G, tau=beta, gamma=ave_inf_period_rate) 
        beta_lst.append(beta)
        #print(beta)
    #plt.hist(beta_lst)
    #print("Median beta value for beta_R0: " + str(statistics.median(beta_lst)))
    beta_R0 = statistics.median(beta_lst)
    
    # Specify transitions and transmissions -----------------------------------
    H = nx.DiGraph()
    H.add_node('S')
    H.add_node('V')
    H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
    H.add_edge('I', 'R', rate = (1 - mort) * (1 / ave_inf_period_recover))
    H.add_edge('I', 'D', rate = (mort) * (1 / ave_inf_period_dead))
    return_statuses = ('S', 'E', 'I', 'R', 'D', 'V')
    J = nx.DiGraph()
    J.add_edge(('I', 'S'), ('I', 'E'), rate = beta_R0, weight_label='transmission_weight')
    # There are no vaccinated nodes in this simulation so the following two transmissions are irrelevant
    J.add_edge(('I', 'V'), ('I', 'E'), rate = 0, weight_label='transmission_weight')
    J.add_edge(('V', 'S'), ('V', 'V'), rate = 0, weight_label='transmission_weight')
    
    # Find day on average when expected_It_N of active infections (2000 sims) -
    nsim = 2000
    I_series = []
    while (len(I_series) < nsim):
        #if (len(I_series) % 200 == 0):
            #print(len(I_series))
        continue_loop = True
        while (continue_loop):
            z = []
            for i in range(N_cluster):
                deg = np.random.negative_binomial(k_overdispersion, p)
                z.append(deg)
            for i in range(len(z)):
                if (z[i] == 0):
                    z[i] == 1
            if (sum(z) % 2 == 0):
                continue_loop = False
        G=nx.configuration_model(z)
        G=nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))
        node_attribute_dict = {node: 1 for node in G.nodes()}
        edge_attribute_dict = {edge: 1 for edge in G.edges()}
        nx.set_node_attributes(G, values=node_attribute_dict, name='expose2infect_weight')
        nx.set_edge_attributes(G, values=edge_attribute_dict, name='transmission_weight')
        IC = defaultdict(lambda: 'S')
        for node in range(initial_infections_per_cluster):
            IC[node] = 'I'
        t, S, E, I, R, D, V = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = 200)
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t)):
            if t[t_dex] >= next_t:
                to_add_row.append(I[t_dex])
                next_t += 1
        I_series.append(to_add_row)
    interrupt_t = None
    
    # Find first day of sim where the ave. num. of infects >= expected_It_N ---
    for day_dex in range(nsim):
        focal_dist = []
        for I_series_dex in range(len(I_series)):
            if len(I_series[I_series_dex]) > day_dex:
                focal_dist.append(I_series[I_series_dex][day_dex] / N_cluster)
        if len(focal_dist) <= 200:
            print("Not enough simulations (<10%) to get average number of infections on this day.")
            interrupt_t = 'na'
            break
        else:
            #print(len(focal_dist))
            #print(statistics.mean(focal_dist))
            if statistics.mean(focal_dist) >= eit:
                interrupt_t = day_dex
                break
        
    # Write output ------------------------------------------------------------    
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0) + "_mort" + str(mort) + "_eit" + str(eit) + '.csv'
    with open(filename, 'w') as out_f:
        out_f.write(str(beta_R0))
        out_f.write(",")
        out_f.write(str(interrupt_t))
    
if __name__ == '__main__':
    pool = mp.Pool(mp.cpu_count() - 1) # Don't use all CPUs
    pool.map(getPrelim, param_sets)
    pool.close()
    pool.join()
    
    
    