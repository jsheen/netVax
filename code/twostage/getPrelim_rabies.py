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
import matplotlib.pyplot as plt
import statistics
import multiprocessing as mp
import networkx as nx
from collections import defaultdict
import EoN
from pathlib import Path
home = str(Path.home())
plot_overdispersion = False
run_sims = False

# Fixed variables -------------------------------------------------------------
alpha = 0.1
tau1 = 1 / 21
tau2 = 1 / 5.78
epsilon = 1 / 135
delta = 1 / 6

# Set parameter sets ----------------------------------------------------------
Ns = [1000]
overdispersions = [1]
R0s = [0.25, 1.1, 3]
param_sets = []
for i in Ns:
    for j in overdispersions:
        for k in R0s:
            param_sets.append([i, j, k])
            
# Show overdispersion ---------------------------------------------------------
if plot_overdispersion:
    k_overdispersion = 1
    mean_degree = 15
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    draws = []
    for i in range(1000):
        draws.append(np.random.negative_binomial(k_overdispersion, p))
    plt.hist(draws, bins=100)
            
# Test that epidemics are functioning properly --------------------------------
if run_sims:
    for param_set in param_sets:
        N_cluster = param_set[0]
        k_overdispersion = param_set[1]
        R0 = param_set[2]
        eit = 0.005
        mean_degree = 15
        p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
        if N_cluster == 1000:
            initial_infections_per_cluster = 4
        elif N_cluster == 10000:
            initial_infections_per_cluster = 40
        if eit <= initial_infections_per_cluster / N_cluster:
            raise NameError("Script assumes expected It / N strictly < initial infections per cluster / N.")
        # Method to estimate_beta from Blackwood estimate of R0 -------------------
        def estimate_beta(R0):
            beta = R0 / (alpha * ((1 / tau2) + (1 / delta)))
            return beta
        beta_R0 = estimate_beta(R0)
        # Specify transitions and transmissions -------------------------------
        H = nx.DiGraph()
        H.add_node('S')
        H.add_node('V')
        H.add_edge('E', 'I_N', rate = alpha * tau1)
        H.add_edge('E', 'T', rate = (1 - alpha) * tau1)
        H.add_edge('I_N', 'I_R', rate = tau2)
        H.add_edge('I_R', 'D', rate = delta)
        H.add_edge('T', 'S', rate = epsilon)
        return_statuses = ('S', 'E', 'T', 'I_N', 'I_R', 'D', 'V')
        J = nx.DiGraph()
        J.add_edge(('I_N', 'S'), ('I_N', 'E'), rate = beta_R0)
        J.add_edge(('I_R', 'S'), ('I_R', 'E'), rate = beta_R0)
        # There are no vaccinated nodes in this simulation so the following transmissions are irrelevant
        J.add_edge(('I_N', 'V'), ('I_N', 'E'), rate = 0)
        J.add_edge(('I_R', 'V'), ('I_R', 'E'), rate = 0)
        J.add_edge(('V', 'S'), ('V', 'V'), rate = 0)
        # Run simulations -----------------------------------------------------
        nsim = 100
        cc_greater_1 = []
        for sim_num in range(nsim):
            continue_loop = True
            while (continue_loop):
                z = []
                for i in range(N_cluster):
                    deg = 0
                    while (deg == 0):
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
            # Remove singletons
            list_of_deg0 = [node for node in G.nodes if G.degree(node) == 0]
            for node_deg0 in list_of_deg0:
                G.add_edge(node_deg0, np.random.choice(G.nodes()))
            degree_sequence = [d for n, d in G.degree()]
            if len(np.where(degree_sequence == 0)[0]) > 0:
                raise NameError('There are singletons in this graph.')
            if (nx.number_connected_components(G) != 1): # Check that number of connected components is not large
                cc_greater_1.append(nx.number_connected_components(G))
            IC = defaultdict(lambda: 'S')
            for node in range(initial_infections_per_cluster):
                IC[node] = 'I_N'
            t, S, E, T, I_N, I_R, D, V = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = float('Inf'))
            plt.plot(t, I_N + I_R)
            plt.plot(t, D)        
        # Should show graphically that the number of connected components is small
        print(str(len(cc_greater_1) / nsim))
                   
# Function to find day of intervention and transmission rates for R0s ---------
def getPrelim(param_set):
    N_cluster = param_set[0]
    k_overdispersion = param_set[1]
    R0 = param_set[2]
    eit = 0.005
    mean_degree = 15
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    if N_cluster == 1000:
        initial_infections_per_cluster = 4
    elif N_cluster == 10000:
        initial_infections_per_cluster = 40
    if eit <= initial_infections_per_cluster / N_cluster:
        raise NameError("Script assumes expected It / N strictly < initial infections per cluster / N.")
    # Method to estimate_beta from Blackwood estimate of R0 -------------------
    def estimate_beta(R0):
        beta = R0 / (alpha * ((1 / tau2) + (1 / delta)))
        return beta
    beta_R0 = estimate_beta(R0)
    
    # Specify transitions and transmissions -----------------------------------
    H = nx.DiGraph()
    H.add_node('S')
    H.add_node('V')
    H.add_edge('E', 'I_N', rate = alpha * tau1)
    H.add_edge('E', 'T', rate = (1 - alpha) * tau1)
    H.add_edge('I_N', 'I_R', rate = tau2)
    H.add_edge('I_R', 'D', rate = delta)
    H.add_edge('T', 'S', rate = epsilon)
    return_statuses = ('S', 'E', 'T', 'I_N', 'I_R', 'D', 'V')
    J = nx.DiGraph()
    J.add_edge(('I_N', 'S'), ('I_N', 'E'), rate = beta_R0)
    J.add_edge(('I_R', 'S'), ('I_R', 'E'), rate = beta_R0)
    # There are no vaccinated nodes in this simulation so the following two transmissions are irrelevant
    J.add_edge(('I_N', 'V'), ('I_N', 'E'), rate = 0)
    J.add_edge(('I_R', 'V'), ('I_R', 'E'), rate = 0)
    J.add_edge(('V', 'S'), ('V', 'V'), rate = 0)
    
    # Find day on average when expected_It_N_N of active infections (2000 sims)
    nsim = 2000
    I_series = []
    while (len(I_series) < nsim):
        continue_loop = True
        while (continue_loop):
            z = []
            for i in range(N_cluster):
                deg = 0
                while (deg == 0):
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
        # Remove singletons
        list_of_deg0 = [node for node in G.nodes if G.degree(node) == 0]
        for node_deg0 in list_of_deg0:
            G.add_edge(node_deg0, np.random.choice(G.nodes()))
        degree_sequence = [d for n, d in G.degree()]
        if len(np.where(degree_sequence == 0)[0]) > 0:
            raise NameError('There are singletons in this graph.')
        IC = defaultdict(lambda: 'S')
        for node in range(initial_infections_per_cluster):
            IC[node] = 'I_N'
        t, S, E, T, I_N, I_R, D, V = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = 200)
        next_t = 0
        to_add_row = []
        for t_dex in range(len(t)):
            if t[t_dex] >= next_t:
                to_add_row.append(I_N[t_dex] + I_R[t_dex])
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
            if statistics.mean(focal_dist) >= eit:
                interrupt_t = day_dex
                break
        
    # Write output ------------------------------------------------------------    
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0) + "_eit" + str(eit) + '_rabies.csv'
    with open(filename, 'w') as out_f:
        out_f.write(str(beta_R0))
        out_f.write(",")
        out_f.write(str(interrupt_t))

if __name__ == '__main__':
    pool = mp.Pool(mp.cpu_count() - 1) # Don't use all CPUs
    pool.map(getPrelim, param_sets)
    pool.close()
    pool.join()
    
    
    