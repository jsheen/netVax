#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 11:57:51 2023

@author: jsheen
"""

# Import libraries and set seeds ----------------------------------------------
import numpy as np
np.random.seed(0)
import networkx as nx
from collections import defaultdict
import EoN
import math
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
home = str(Path.home())

# Fixed variables -------------------------------------------------------------
sigma = 1 / 5
gamma = 1 / 10
psi = 1 / 120
mean_degree = 15
nsim = 100
cutoff = 150

# Set parameters --------------------------------------------------------------
param_set = [1000, 2, 1.1, 0.8, 0.05]
N_cluster = param_set[0]
R0_wt = param_set[1]
vax = param_set[2]
vax_eff = param_set[3]
assign = param_set[4]
eit = 0.005
if N_cluster == 1000:
    initial_infections_per_cluster = 4
elif N_cluster == 10000:
    initial_infections_per_cluster = 40

# Get pre-set beta values and time of intervention value ------------------
# Wild-type
filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_R0" + str(R0_wt) + "_eit" + str(eit) + '_SEIR_Pois.csv'
prelim_data_wt = pd.read_csv(filename, header=None)
beta_R0_wt = prelim_data_wt[0][0]
interrupt_t = prelim_data_wt[1][0]
# Vax
R0_vax = vax
if float(R0_vax) > 0:
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_R0" + str(R0_vax) + "_eit" + str(eit) + '_SEIR_Pois.csv'
    prelim_data_vax = pd.read_csv(filename, header=None)
    beta_R0_vax = prelim_data_vax[0][0]
else:
    beta_R0_vax = 0

# Specify transitions and transmissions -------------------------------
H = nx.DiGraph()
H.add_node('S')
H.add_edge('V_E', 'V_I', rate=sigma)
H.add_edge('V_I', 'V_R', rate=gamma)
H.add_edge('E', 'I', rate=sigma)
H.add_edge('I', 'R', rate=gamma)
return_statuses = ('S', 'E', 'I', 'R', 'V_E', 'V_I', 'V_R')
J = nx.DiGraph()
J.add_edge(('I', 'S'), ('I', 'E'), rate = beta_R0_wt)
J.add_edge(('I', 'V_R'), ('I', 'E'), rate = (1 - vax_eff) * beta_R0_wt)
J.add_edge(('V_I', 'S'), ('V_I', 'V_E'), rate = beta_R0_vax)

# Set threshold value of number of infections at time t -------------------
threshold = 1

count_burnout = 0

# Simulate epidemics with vaccination -------------------------------------
for sim in range(nsim):
    continue_loop = True
    while (continue_loop):
        z = []
        for i in range(N_cluster):
            deg = 0
            while (deg == 0):
                deg = np.random.poisson(mean_degree)
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
        IC[node] = 'I'
    full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = math.ceil(interrupt_t), return_full_data=True) 
    t_first_half = full_first_half.t()
    I_first_half = full_first_half.I()
    R_first_half = full_first_half.R()
    if I_first_half[-1] >= threshold:
        curr_IC_con = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        curr_IC = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
        suscep_nodes = [k for k,v in curr_IC.items() if v == 'S']
        enrolled_nodes = np.random.choice(suscep_nodes, size=int(np.ceil(assign * len(suscep_nodes))), replace=False)
        curr_IC.update(curr_IC.fromkeys(enrolled_nodes, 'V_I'))
        test_vax_cnt = 0
        for k,v in curr_IC.items():
            if v == 'V_I':
                test_vax_cnt += 1
        if test_vax_cnt != len(enrolled_nodes):
            raise NameError("Error in assignment.")
        full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC_con, return_statuses, tmax = float(500), return_full_data=True)    
        full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float(500), return_full_data=True)    
        
        # Vaccination increase
        new_t_con = np.concatenate((t_first_half, np.array(full_second_half_con.t()) + t_first_half[-1]))
        new_R_con = np.concatenate((R_first_half, np.array(full_second_half_con.R())))
        new_t_trt = np.concatenate((t_first_half, np.array(full_second_half_trt.t()) + t_first_half[-1]))
        new_R_trt = np.concatenate((R_first_half, np.array(full_second_half_trt.R())))

        if new_t_con[-1] < 150:
            new_t_con = np.append(new_t_con, 150)
            new_R_con = np.append(new_R_con, new_R_con[-1])
        if new_t_trt[-1] < 150:
            new_t_trt = np.append(new_t_trt, 150)
            new_R_trt = np.append(new_R_trt, new_R_trt[-1])
        plt.plot(new_t_con, new_R_con / 1000, 'grey')
        plt.plot(new_t_trt, new_R_trt / 1000, 'orange')
        plt.xlim(0, cutoff)
        plt.ylim(-0.05, 1)
        plt.xlabel('', fontdict={'size':20})
        plt.ylabel('', fontdict={'size':20})
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(4.5, 0.95, 'Reactionary Trial\nNo Overdispersion', fontsize=17,
        verticalalignment='top', bbox=props)
plt.savefig(home + "/netVax/code_output/figs/pan2.png", format="png", dpi=400)
plt.show()
    
        