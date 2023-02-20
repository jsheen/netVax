#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 14:32:32 2022

@author: jsheen

@description: script used to run a two stage experiment to identify the spillover
              effect of a transmissible vaccine
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
np.random.seed(0)
import networkx as nx
from collections import defaultdict
import EoN
import math
import pandas as pd
import multiprocessing as mp
from pathlib import Path
home = str(Path.home())

# Fixed variables -------------------------------------------------------------
alpha = 0.1
tau1 = 1 / 21
tau2 = 1 / 5.78
epsilon = 1 / 135
delta = 1 / 6

# Set parameter sets ----------------------------------------------------------
Ns = [1000]
overdispersions = [1]
R0_wts = [3]
vaxs = [0, 0.25, 1.1]
vax_effs = [0.8]
assigns = [0, 0.1, 0.2]
sim_num = 1000
param_sets = []
for i in Ns:
    for j in overdispersions:
        for k in R0_wts:
            for l in vaxs:
                for m in vax_effs:
                    for n in assigns:
                        for o in range(sim_num):
                            param_sets.append([i, j, k, l , m, n, o])


# For each parameter set, create simulations ----------------------------
def runSim(param_set):
    N_cluster = param_set[0]
    k_overdispersion = param_set[1]
    R0_wt = param_set[2]
    vax = param_set[3]
    vax_eff = param_set[4]
    assign = param_set[5]
    sim_num = param_set[6]
    eit = 0.005
    mean_degree = 15
    p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
    if N_cluster == 1000:
        initial_infections_per_cluster = 4
    elif N_cluster == 10000:
        initial_infections_per_cluster = 40
    
    # Get pre-set beta values and time of intervention value ------------------
    # Wild-type
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_wt) + "_eit" + str(eit) + '_rabies.csv'
    prelim_data_wt = pd.read_csv(filename, header=None)
    beta_R0_wt = prelim_data_wt[0][0]
    interrupt_t = prelim_data_wt[1][0]
    # Vax
    R0_vax = vax
    if float(R0_vax) > 0:
        filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_vax) + "_eit" + str(eit) + '_rabies.csv'
        prelim_data_vax = pd.read_csv(filename, header=None)
        beta_R0_vax = prelim_data_vax[0][0]
    else:
        beta_R0_vax = 0
    
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
    J.add_edge(('I_N', 'S'), ('I_N', 'E'), rate = beta_R0_wt)
    J.add_edge(('I_R', 'S'), ('I_R', 'E'), rate = beta_R0_wt)
    J.add_edge(('I_N', 'V'), ('I_N', 'E'), rate = (1 - vax_eff) * beta_R0_wt)
    J.add_edge(('I_R', 'V'), ('I_R', 'E'), rate = (1 - vax_eff) * beta_R0_wt)
    J.add_edge(('V', 'S'), ('V', 'V'), rate = beta_R0_vax)
    
    # Set threshold value of number of infections at time t -------------------
    threshold = 1
    
    # Simulate epidemics with vaccination -------------------------------------
    continue_loop = True
    while (continue_loop):
        z = []
        for i in range(N_cluster):
            deg = 0
            while (deg == 0):
                deg = np.random.negative_binomial(k_overdispersion, p)
            z.append(deg)
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
    full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = math.ceil(interrupt_t), return_full_data=True) 
    t_first_half = full_first_half.t()
    I_N_first_half = full_first_half.summary()[1]['I_N']
    I_R_first_half = full_first_half.summary()[1]['I_R']
    if I_N_first_half[-1] + I_R_first_half[-1] >= threshold:
        curr_IC_con = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        curr_IC = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        suscep_nodes = [k for k,v in curr_IC.items() if v == 'S']
        enrolled_nodes = np.random.choice(suscep_nodes, size=int(np.ceil(assign * len(suscep_nodes))), replace=False)
        curr_IC.update(curr_IC.fromkeys(enrolled_nodes, 'V'))
        test_vax_cnt = 0
        for k,v in curr_IC.items():
            if v == 'V':
                test_vax_cnt += 1
        if test_vax_cnt != len(enrolled_nodes):
            raise NameError("Error in assignment.")
        full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC_con, return_statuses, tmax = float(500), return_full_data=True)    
        full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float(500), return_full_data=True)    
        # Control
        I_N_second_half_con = full_second_half_con.summary()[1]['I_N']
        I_R_second_half_con = full_second_half_con.summary()[1]['I_R']
        if len(np.where((I_N_second_half_con == 0) & (I_R_second_half_con == 0) & (full_second_half_con.summary()[1]['E'] == 0))[0]) == 0:
            t_no_inf_con = 500
        else:
            t_no_inf_con = full_second_half_con.t()[np.where((I_N_second_half_con == 0) & (I_R_second_half_con == 0) & (full_second_half_con.summary()[1]['E'] == 0))[0][0]]
        surv_inf_con = dict.fromkeys(G.nodes(), t_no_inf_con)
        surv_dead_con = dict.fromkeys(G.nodes(), t_no_inf_con)
        for node in G.nodes():
            node_hist = full_second_half_con.node_history(node)
            if 'E' in node_hist[1]:
                surv_inf_con[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'E')[0][0]]
            if 'D' in node_hist[1]:
                surv_dead_con[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'D')[0][0]]
                if node_hist[1][0] == 'D': # This corrects for deaths in the first time step, which should not be counted
                    surv_inf_con[node] = -1
        # Treatment
        I_N_second_half_trt = full_second_half_trt.summary()[1]['I_N']
        I_R_second_half_trt = full_second_half_trt.summary()[1]['I_R']
        if len(np.where((I_N_second_half_trt == 0) & (I_R_second_half_trt == 0) & (full_second_half_trt.summary()[1]['E'] == 0))[0]) == 0:
            t_no_inf_trt = 500
        else:
            t_no_inf_trt = full_second_half_trt.t()[np.where((I_N_second_half_trt == 0) & (I_R_second_half_trt == 0) & (full_second_half_trt.summary()[1]['E'] == 0))[0][0]]
        surv_inf_trt = dict.fromkeys(G.nodes(), t_no_inf_trt)
        surv_dead_trt = dict.fromkeys(G.nodes(), t_no_inf_trt)
        for node in G.nodes():
            node_hist = full_second_half_trt.node_history(node)
            if 'E' in node_hist[1]:
                surv_inf_trt[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'E')[0][0]]
            if 'D' in node_hist[1]:
                surv_dead_trt[node] = node_hist[0][np.where(np.array(node_hist[1]) == 'D')[0][0]]
                if node_hist[1][0] == 'D':
                    surv_inf_trt[node] = -1
        # Write results
        with open(home + '/netVax/code_output/twostage/sims/2stg_N' + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0wt" + str(R0_wt) + 
                  "_R0vax" + str(R0_vax) + "_eit" + str(eit) +
                  '_vaxEff' + str(vax_eff) + '_assign' + str(assign) + '_sim' + str(sim_num) + '_rabies.csv', 'w') as out_f:
            out_f.write('node, assignment, time2inf_trt, time2death_trt, time2inf_con, time2death_con\n')
            for node in G.nodes():
                out_f.write(str(node))
                out_f.write(',')
                if node in enrolled_nodes:
                    out_f.write('e')
                elif node not in enrolled_nodes:
                    out_f.write('na')
                out_f.write(',')
                out_f.write(str(surv_inf_trt.get(node)))
                out_f.write(',')
                out_f.write(str(surv_dead_trt.get(node)))
                out_f.write(',')
                out_f.write(str(surv_inf_con.get(node)))
                out_f.write(',')
                out_f.write(str(surv_dead_con.get(node)))
                out_f.write('\n')
    else:
        with open(home + '/netVax/code_output/twostage/sims/2stg_N' + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0wt" + str(R0_wt) + 
          "_R0vax" + str(R0_vax) + "_eit" + str(eit) +
          '_vaxEff' + str(vax_eff) + '_assign' + str(assign) + '_sim' + str(sim_num) + '_rabies.csv', 'w') as out_f:
            out_f.write('node, assignment, time2inf_trt, time2death_trt, time2inf_con, time2death_con\n')
            out_f.write('na\n')

if __name__ == '__main__':
    pool = mp.Pool(mp.cpu_count() - 1) # Don't use all CPUs
    pool.map(runSim, param_sets)
    pool.close()
    pool.join()
            
            
