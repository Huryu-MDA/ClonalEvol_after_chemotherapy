#!/usr/bin/env python
# coding: utf-8

import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt


#phylo_fig_output = sys.argv[1]

PID_prefix = sys.argv[1]
phylo_fig_output = PID_prefix + "_phylo_fig_REFSEQpos_NoMod.pdf"

# root_len means the length between the germline and the toppest bifurcation. This time, the length are set as zero. But usually, there are some length between germline and the toppest node.
root_len = 0.0015

edge_df = pd.read_csv("tr_edge_mtx.tsv", sep = "\t")
len_edge = [float(line.split("\n")[0]) for line in open("edge_length.txt")]
edge_df["len"] = len_edge

tip_dict = {ind + 1 : line.split("\n")[0] for ind, line in enumerate(open("tip_label.txt"))}

tip_list = np.array(sorted(tip_dict.keys()))
node_max = max(edge_df[["prox", "dist"]].max())
node_list = np.arange(len(tip_list) + 1, node_max + 1)
#print(node_max, node_list)

tip_node_list = np.append(tip_list, node_list)
#tip_node_list

#edge_df["dist"]

# # This will return the parent of the node
# 
# parent_ser = edge_df["prox"]
# child_ser = edge_df["dist"]
# 
# parent_dict = {num: int(parent_ser.loc[child_ser == num]) for num in tip_node_list if num in list(child_ser)}
# child_dict = {num: list(child_ser[parent_ser == num]) for num in tip_node_list if num in list(parent_ser)}

parent_list = list(edge_df["prox"])
child_list = list(edge_df["dist"])
edge_len_list = list(edge_df["len"])

C_to_P_dict = {c_node: p_node for c_node, p_node in zip(child_list, parent_list)}
P_to_paired_C_dict = {num: list(edge_df.loc[edge_df["prox"] == num, "dist"]) for num in parent_list}
C_P_len_dict = {node: edge_len for node, edge_len in zip(child_list, edge_len_list)}

# node_num would be int in my first trial

def node_coordinate(node_num):
    x = node_x(node_num)
    y = node_y(node_num)
    return (x, y)

# child_list = list(child_ser)
def node_x(node_num, 
           root_len = root_len, 
           child_list = child_list, 
           C_to_P_dict = C_to_P_dict,
           C_P_len_dict = C_P_len_dict):
    
    if node_num in child_list:
        p_node = C_to_P_dict[node_num]
        p_node_x = node_x(p_node)
        c_to_p_len = C_P_len_dict[node_num]
        x = p_node_x + c_to_p_len
        return x
    else:
        # in this case, node_num is the toppest node in the tree.
        x = root_len
        return x

def node_y(node_num, tip_list = tip_list,
           P_to_C_dict = P_to_paired_C_dict):
    if node_num in tip_list:
        y = node_num
        return y
    else:
        child_pair = P_to_paired_C_dict[node_num]
        y_child1 = node_y(child_pair[0])
        y_child2 = node_y(child_pair[1])
        y = (y_child1 + y_child2) / 2
        return y



def write_tip_label(x_max, tip_list = tip_list, tip_dict = tip_dict):
    # print(list(tip_list))
    # print(tip_dict)
    for node in tip_list:
        tiplabel=tip_dict[node]
        #plt.text(x = x_max + 0.01, y = node - 0.2, s = tiplabel)
        plt.text(x = x_max * 1.05, y = node - 0.2, s = tiplabel)


x_max_node = max([node_x(node) for node in tip_node_list])

# plot the tip label on the phylo
tip_labelling = "on"
if tip_labelling == "on":
    write_tip_label(x_max_node)

for node in child_list:
    p_node = C_to_P_dict[node]
    n_x = node_x(node)
    n_y = node_y(node)
    np_x = node_x(p_node)
    np_y = node_y(p_node)
    #print(node, p_node)
    #print(node, p_node)
    #plt.plot([n_x, np_x], [n_y, np_y])
    #plt.show()
    # C_P_hori
    #plt.plot([n_x, np_x], [n_y, n_y], c = "blue", linewidth = 2)
    plt.plot([n_x, np_x], [n_y, n_y], c = "black", linewidth = 2)
    # ---> This part would be replaced with the stacked_bar_graph
    
    # C_P_vert
    #plt.plot([np_x, np_x], [n_y, np_y], c = "red", linewidth = 2)
    plt.plot([np_x, np_x], [n_y, np_y], c = "black", linewidth = 2)

plt.title(PID_prefix)
plt.xticks([])
plt.yticks([])


plt.savefig(phylo_fig_output, transparent=True, bbox_inches='tight')
#plt.show()
    #break

