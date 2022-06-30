#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:17:07 2022

@author: jsheen

@description: example plots
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
np.random.seed(0)
import networkx as nx
from collections import defaultdict
import EoN
import math
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())

# Set parameter sets ----------------------------------------------------------
Ns = [1000]
overdispersions = [1]
R0_wts = [2.5]
vaxs = ['R0=0_treat=0.5']
morts = [0.85]
vax_effs = [0.6]
sim_num = 100
param_sets = []
for i in Ns:
    for j in overdispersions:
        for k in R0_wts:
            for l in vaxs:
                for m in morts:
                    for n in vax_effs:
                        for o in range(sim_num):
                            param_sets.append([i, j, k, l , m, n, o])


# For each parameter set, create 3,000 simulations ----------------------------
store_t = []
store_I = []
def runSim(param_set):
    N_cluster = param_set[0]
    k_overdispersion = param_set[1]
    R0_wt = param_set[2]
    vax = param_set[3]
    mort = param_set[4]
    vax_eff = param_set[5]
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
    
    # Get pre-set beta values and time of intervention value ------------------
    # Wild-type
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_wt) + "_mort" + str(mort) + "_eit" + str(eit) + '.csv'
    prelim_data_wt = pd.read_csv(filename, header=None)
    beta_R0_wt = prelim_data_wt[0][0]
    interrupt_t = prelim_data_wt[1][0]
    # Vax
    R0_vax = float(vax.split('_')[0].replace('R0=', ''))
    R0_vax = ('%f' % R0_vax).rstrip('0').rstrip('.')
    if float(R0_vax) > 0:
        filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_vax) + "_mort" + str(mort) + "_eit" + str(eit) + '.csv'
        prelim_data_vax = pd.read_csv(filename, header=None)
        beta_R0_vax = prelim_data_vax[0][0]
    else:
        beta_R0_vax = 0
    vax_treat = float(vax.split('_')[1].replace('treat=', ''))
    
    # Specify transitions and transmissions -----------------------------------
    H = nx.DiGraph()
    H.add_node('S')
    H.add_node('V')
    H.add_edge('E', 'I', rate = 1 / ave_inc_period, weight_label='expose2infect_weight')
    H.add_edge('I', 'R', rate = (1 - mort) * (1 / ave_inf_period_recover))
    H.add_edge('I', 'D', rate = (mort) * (1 / ave_inf_period_dead))
    return_statuses = ('S', 'E', 'I', 'R', 'D', 'V')
    J = nx.DiGraph()
    J.add_edge(('I', 'S'), ('I', 'E'), rate = beta_R0_wt, weight_label='transmission_weight')
    J.add_edge(('I', 'V'), ('I', 'E'), rate = (1 - vax_eff) * beta_R0_wt, weight_label='transmission_weight')
    J.add_edge(('V', 'S'), ('V', 'V'), rate = beta_R0_vax, weight_label='transmission_weight')
    
    # Set threshold value of number of infections at time t -------------------
    threshold = 1
    
    # Simulate epidemics with vaccination -------------------------------------
    continue_loop = True
    while (continue_loop):
        z = []
        for i in range(N_cluster):
            deg = np.random.negative_binomial(k_overdispersion, p)
            z.append(deg)
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
    full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = math.ceil(interrupt_t), return_full_data=True) 
    t_first_half = full_first_half.t()
    I_first_half = full_first_half.I()
    if I_first_half[-1] >= threshold:
        curr_IC = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        suscep_nodes = [k for k,v in curr_IC.items() if v == 'S']
        suscep_nodes_half = np.random.choice(suscep_nodes, size=int(np.ceil(0.5 * len(suscep_nodes))), replace=False)
        assign_treat = list(np.random.choice(suscep_nodes_half, size=int(np.ceil(vax_treat * len(suscep_nodes_half))), replace=False))
        curr_IC.update(curr_IC.fromkeys(assign_treat, 'V'))
        assign_control = list(set(suscep_nodes_half).difference(set(assign_treat)))
        if len(assign_treat) + len(assign_control) != len(suscep_nodes_half):
            raise NameError("Unequal treatment and control groups.")
        if set(assign_treat).union(set(assign_control)) != set(suscep_nodes_half):
            raise NameError("Unequal sets of treatment and control groups.")
        test_vax_cnt = 0
        for k,v in curr_IC.items():
            if v == 'V':
                test_vax_cnt += 1
        if test_vax_cnt != len(assign_treat):
            raise NameError("Error in treatment assignment.")
        full_second_half = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float('Inf'), return_full_data=True)    
        store_t.append(full_second_half.t())
        store_I.append(full_second_half.I())
        """
        tot_con_inf = 0
        tot_trt_inf = 0
        for node in G.nodes():
            if 'I' in full_second_half.node_history(node)[1]:
                if node in assign_treat:
                    tot_trt_inf += 1
                elif node in assign_control:
                    tot_con_inf += 1
        if tot_con_inf <= tot_trt_inf:
            print('More or equal infections in treatment than control.')
            print(str(tot_con_inf - tot_trt_inf))
            plt.plot(full_second_half.t(), full_second_half.I())
        """
        
for param_set in param_sets:
    runSim(param_set)

plt.plot(store_t[0], store_I[0])
for sim_dex in range(1, len(store_t)):
    plt.plot(store_t[sim_dex], store_I[sim_dex])
            
