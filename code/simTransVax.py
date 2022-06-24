#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 19:30:09 2022

@author: jsheen

@description: script used to create village simulations, where roughly half of
              individual nodes are assigned to transmissible vaccine treatment 
              and the other half are assigned to control
              
@assumptions: we assume that vaccine is protective for entire observation period
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
import networkx as nx
from collections import defaultdict
import EoN
import math
import pandas as pd
import multiprocessing as mp
from pathlib import Path
home = str(Path.home())

# Set parameter sets ----------------------------------------------------------
Ns = [1000, 10000]
overdispersions = [1]
R0_wts = [2]
R0_vaxs = [1.5]
morts = [0.85]
sim_num = 3000
param_sets = []
for i in Ns:
    for j in overdispersions:
        for k in R0_wts:
            for l in R0_vaxs:
                for m in morts:
                    for n in range(sim_num):
                        param_sets.append([i, j, k, l , m, n])


# For each parameter set, create 3,000 simulations ----------------------------
def runSim(param_set):
    N_cluster = param_set[0]
    k_overdispersion = param_set[1]
    R0_wt = param_set[2]
    R0_vax = param_set[3]
    mort = param_set[4]
    sim_num = param_set[5]
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
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0wt" + str(R0_wt) + "_R0vax" + str(R0_vax) + "_mort" + str(mort) + "_eit" + str(eit) + '.csv'
    prelim_data = pd.read_csv(filename, header=None)
    beta_R0_wt = prelim_data[0][0]
    beta_R0_vax = prelim_data[1][0]
    interrupt_t = prelim_data[2][0]
    
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
        nodes_first_half_final = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        curr_IC = defaultdict(lambda: 'S')
        to_treat = True
        assign_treat = []
        assign_control = []
        for node in G.nodes():
             status = nodes_first_half_final[node]
             # Assign half of susceptible individuals to treatment, the other half to control
             if (status == 'S'):
                 if (to_treat):
                     curr_IC[node] = 'V'
                     to_treat = False
                     assign_treat.append(node)
                 else:
                     curr_IC[node] = status
                     to_treat = True
                     assign_control.append(node)
             else:
                 curr_IC[node] = status
        if abs(len(assign_treat) - len(assign_control)) > 1:
            raise NameError("Unequal treatment and control groups.")
        full_second_half = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float('Inf'), return_full_data=True)    
        t_no_inf = full_second_half.t()[np.where((full_second_half.I() == 0) & (full_second_half.summary()[1]['E'] == 0))[0][0]]
        surv_inf = dict.fromkeys(G.nodes(), t_no_inf)
        surv_dead = dict.fromkeys(G.nodes(), t_no_inf)
        for node in G.nodes():
            node_hist = full_second_half.node_history(node)
            if 'I' in node_hist[1]:
                surv_inf[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'I')[0][0]]
            if 'D' in node_hist[1]:
                surv_dead[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'D')[0][0]]
                if node_hist[1][0] == 'D':
                    surv_inf[node] = -1
        with open(home + '/netVax/code_output/sim_results/N' + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0wt" + str(R0_wt) + "_R0vax" + str(R0_vax) + "_mort" + str(mort) + "_eit" + str(eit) + '_sim' + str(sim_num) + '.csv', 'w') as out_f:
            out_f.write('node,assignment,time2inf,time2death\n')
            for node in G.nodes():
                out_f.write(str(node))
                out_f.write(',')
                if node in assign_treat and node not in assign_control:
                    out_f.write('t')
                elif node not in assign_treat and node in assign_control:
                    out_f.write('c')
                elif node not in assign_treat and node not in assign_control:
                    out_f.write('na')
                out_f.write(',')
                out_f.write(str(surv_inf.get(node)))
                out_f.write(',')
                out_f.write(str(surv_dead.get(node)))
                out_f.write('\n')
            
if __name__ == '__main__':
    pool = mp.Pool(mp.cpu_count() - 1) # Don't use all CPUs
    pool.map(runSim, param_sets)
    pool.close()
            
            
