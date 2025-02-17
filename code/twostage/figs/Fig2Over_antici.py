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
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
home = str(Path.home())

# Fixed variables -------------------------------------------------------------
sigma = 1 / 5
gamma = 1 / 10
psi = 1 / 120
nsim = 100
cutoff = 150
anticipatory = 220

# Set parameters --------------------------------------------------------------
param_set = [1000, 1, 2, 1.1, 0.8, 0.05]
N_cluster = param_set[0]
k_overdispersion = param_set[1]
R0_wt = param_set[2]
vax = param_set[3]
vax_eff = param_set[4]
assign = param_set[5]
eit = 0.005
mean_degree = 15
p = 1.0 - mean_degree / (mean_degree + k_overdispersion)
if N_cluster == 1000:
    initial_infections_per_cluster = 4
elif N_cluster == 10000:
    initial_infections_per_cluster = 40

# Get pre-set beta values and time of intervention value ------------------
# Wild-type
filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_wt) + "_eit" + str(eit) + '_SEIR.csv'
prelim_data_wt = pd.read_csv(filename, header=None)
beta_R0_wt = prelim_data_wt[0][0]
interrupt_t = prelim_data_wt[1][0]
# Vax
R0_vax = vax
if float(R0_vax) > 0:
    filename = home + "/netVax/code_output/prelim/N" + str(N_cluster) + "_k" + str(k_overdispersion) + "_R0" + str(R0_vax) + "_eit" + str(eit) + '_SEIR.csv'
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
J.add_edge(('V_I', 'S'), ('V_I', 'V_I'), rate = beta_R0_vax)

# Set threshold value of number of infections at time t -------------------
threshold = 1

# Simulate epidemics with vaccination -------------------------------------
for sim in range(nsim):
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
    # First half
    IC = defaultdict(lambda: 'S')
    enrolled_nodes = np.random.choice(G.nodes, size=int(np.ceil(assign * len(G.nodes))), replace=False)
    for node in enrolled_nodes:
        IC[node] = 'V_I'
    full_first_half = EoN.Gillespie_simple_contagion(G, H, J, IC, return_statuses, tmax = anticipatory, return_full_data=True) 
    t_first_half = full_first_half.t()
    # Second half
    curr_IC_con = defaultdict(lambda: 'S')
    for node in enrolled_nodes:
        curr_IC_con[node] = 'V_R'
    for inf_node in np.random.choice(np.array(list(set(G.nodes()).difference(enrolled_nodes))), initial_infections_per_cluster, replace=False):
        curr_IC_con[inf_node] = 'I'
    curr_IC = full_first_half.get_statuses(list(G.nodes()), t_first_half[-1])
    suscep_nodes = [k for k,v in curr_IC.items() if v == 'S']
    for inf_node in np.random.choice(suscep_nodes, initial_infections_per_cluster, replace=False):
        curr_IC[inf_node] = 'I'
    full_second_half_con = EoN.Gillespie_simple_contagion(G, H, J, curr_IC_con, return_statuses, tmax = float(cutoff), return_full_data=True)    
    full_second_half_trt = EoN.Gillespie_simple_contagion(G, H, J, curr_IC, return_statuses, tmax = float(cutoff), return_full_data=True)    
    new_t_con = full_second_half_con.t()
    new_R_con = np.array(full_second_half_con.R())
    new_t_trt = full_second_half_trt.t()
    new_R_trt = np.array(full_second_half_trt.R())
    if new_t_con[-1] < 150:
        new_t_con = np.append(new_t_con, 150)
        new_R_con = np.append(new_R_con, new_R_con[-1])
    if new_t_trt[-1] < 150:
        new_t_trt = np.append(new_t_trt, 150)
        new_R_trt = np.append(new_R_trt, new_R_trt[-1])
    # Vaccination increase
    plt.plot(new_t_con, new_R_con / 1000, 'grey')
    plt.plot(new_t_trt, new_R_trt / 1000, 'orange')
    plt.xlim(0, cutoff)
    plt.ylim(-0.05, 1)
    plt.xlabel('', fontdict={'size':20})
    plt.ylabel('', fontdict={'size':20})
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    plt.text(4.5, 0.95, 'Anticipatory Trial\nOverdispersion', fontsize=17,
        verticalalignment='top', bbox=props)
plt.savefig(home + "/netVax/code_output/figs/pan3.png", format="png", dpi=400)
plt.show()
