#!/usr/bin/env python
# coding: utf-8

import sys, os, glob, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser("args_handling in phylo_tree plot with sig_counts.")
parser.add_argument("arg1", help="PID_prefix. e.g. PID0005_1stACST ")
parser.add_argument("arg2", help="the path to the file of signature relative contribution. e.g.: PHYML/ML_mmSig/mmsig_relative_20230327140128.tsv")
#parser.add_argument("arg3", help="the path to the matrix (col1: node_label, col2: variant counts in each node.). If this code is used in the pipeline, the value is 'node_x.tsv'")
parser.add_argument("--lw_h", type = int, default= 15, help = "line width of horizontal branches." )
parser.add_argument("--lw_v", type = int, default= 3,  help = "line width of vertical branches." )
parser.add_argument("--fg_h", type = int, default= 35, help = "horizontal figsize of plot" )
parser.add_argument("--fg_v", type = int, default= 30, help = "vertical figsize of plot" )
parser.add_argument("--fg_f", type = str, default= "pdf", help = "output figure format. (default = pdf)")
args = parser.parse_args()

PID_prefix = args.arg1
sig_rel_count = args.arg2
#VarCounts = args.arg3

# plot_setting (phylo_tree horiz-vert size, linewidth of branch (horiz or vert))
linewidth_horiz = args.lw_h
linewidth_vert  = args.lw_v
figsize_horiz   = args.fg_h
figsize_vert    = args.fg_v
figfmt          = args.fg_f
# root_len means the length between the germline and the toppest bifurcation. This time, the length are set as zero. But usually, there are some length between germline and the toppest node.
# root_len = 0.01

# option of decoration, sig_MM2 can be turned off!
tip_labelling = "on"
tip_labelling = input("tip_labelling: on or off ? ")
node_labelling = "off"
node_labelling = input("node_labelling: on or off ? ")
sig_mm2_signal = "off"
sig_mm2_signal = input("sig_mm2_signal: on or off ? ")

# the figure format (e.g. .pdf) will be attached at the timing of the savefig().
phylo_fig_out_prefix = PID_prefix + "_REFSEQtipRemoved_SigCountsOnBranch_" + "SBSMM2_BG" + sig_mm2_signal
#phylo_fig_output = PID_prefix + "_REFSEQtipRemoved_SigCountsOnBranch_" + "SBSMM2_BG" + sig_mm2_signal + ".pdf"


# Data loading
edge_df = pd.read_csv("tr_edge_mtx.tsv", sep = "\t")
len_edge = [float(line.split("\n")[0]) for line in open("edge_length.txt")]
edge_df["len"] = len_edge

tip_dict = {ind + 1 : line.split("\n")[0] for ind, line in enumerate(open("tip_label.txt"))}

tip_list = np.array(sorted(tip_dict.keys()))
node_max = max(edge_df[["prox", "dist"]].max())
node_list = np.arange(len(tip_list) + 1, node_max + 1)
#print(node_max, node_list)

tip_node_list = np.append(tip_list, node_list)

parent_list = list(edge_df["prox"])
child_list = list(edge_df["dist"])
edge_len_list = list(edge_df["len"])

C_to_P_dict = {c_node: p_node for c_node, p_node in zip(child_list, parent_list)}
P_to_paired_C_dict = {num: list(edge_df.loc[edge_df["prox"] == num, "dist"]) for num in parent_list}
C_P_len_dict = {node: edge_len for node, edge_len in zip(child_list, edge_len_list)}


#### node_num to nodelabel/tiplabel
tip_node_label_dict = {}
for node in tip_node_list:
    if node in tip_list:
        tip_node_label_dict[node] = tip_dict[node]
    else:
        tip_node_label_dict[node] = "n-" + str(node)

sig_df = pd.read_csv(sig_rel_count, sep = "\t")
if set(tip_node_label_dict.values()) - set(sig_df.Samples) == {"REFSEQ"}:
    sig_df.loc["refseq_row"] = 0
    sig_df.loc["refseq_row", "Samples"] = "REFSEQ"
else:
    print("Discrepancy between signature_rel_matrix and the tr_edge_mtx.tsv.")
    print("Exitting")
    sys.exit()

var_sig_df = sig_df.rename(columns={"mutations": "var_counts"}).set_index("Samples").copy()

#del sig_df["mutations"]
#var_count_df = pd.read_csv(VarCounts, sep = "\t", names=["Samples", "var_counts"])
#var_sig_df = pd.merge(left = var_count_df, right = sig_df, how = "outer", on = "Samples").fillna(0)
#var_sig_df = var_sig_df.set_index("Samples").copy()

# sig_mm2 bg_noise: On or Off
if sig_mm2_signal == "off":
    var_sig_df["SBS-MM2"] = 0

var_sig_df = var_sig_df.T

#sig_list = list(sig_df.columns[1:])
sig_list = [idx for idx in var_sig_df.index if "sbs" in idx.lower()]

# The colors are picked up with the function of PowerPoint (EyeDropper function), then convert the HEX style with "https://www.rgbtohex.net/"
sig_colors_dict = {'SBS1':"#ACF2D0", 'SBS2':"#1D76B3", 'SBS5':"#63D69E",
                   'SBS8':"#DFC4F5", 'SBS9':"#EBF5BC", 'SBS13':"#F39332", 
                   'SBS18':"#CF5901", 'SBS-MM2':"#E377C2"}


nodelabel_count_dict = {col:int(var_sig_df[col]["var_counts"]) for col in var_sig_df.columns}
#nodelabel_relsig_dict = {col:np.array(var_sig_df[col][var_sig_df.index[1:]]) for col in var_sig_df.columns}
nodelabel_relsig_dict = {col:np.array(var_sig_df[col][sig_list]) for col in var_sig_df.columns}

for node in parent_list:
    C_nodes = P_to_paired_C_dict[node]
    for c_node in C_nodes:
        x_p = nodelabel_count_dict[tip_node_label_dict[node]]
        x_c = nodelabel_count_dict[tip_node_label_dict[c_node]]
        if (x_p >= x_c) & (tip_node_label_dict[c_node] != "REFSEQ"):
            #print(node, c_node, x_p, x_c)
            nodelabel_count_dict[tip_node_label_dict[c_node]] = nodelabel_count_dict[tip_node_label_dict[node]] + 1
    #print(node, C_nodes)
    #break

def node_x(node_num,
           tip_node_label_dict = tip_node_label_dict,
           nodelabel_count_dict = nodelabel_count_dict):
    nodelabel = tip_node_label_dict[node_num]
    var_counts = nodelabel_count_dict[nodelabel]
    return var_counts

def node_coordinate(node_num):
    x = node_x(node_num)
    y = node_y(node_num)
    return (x, y)


# This part is the node_x function in phy_tree_plot.py
# def node_x_deprecated(node_num, 
#            root_len = root_len, 
#            child_list = child_list, 
#            C_to_P_dict = C_to_P_dict,
#            C_P_len_dict = C_P_len_dict):
#    
#     if node_num in child_list:
#        p_node = C_to_P_dict[node_num]
#        p_node_x = node_x(p_node)
#        c_to_p_len = C_P_len_dict[node_num]
#        x = p_node_x + c_to_p_len
#        return x
#    else:
#        # in this case, node_num is the toppest node in the tree.
#        x = root_len
#        return x

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

# TIP LABEL
def write_tip_label(x_max, tip_list = tip_list, tip_dict = tip_dict, x_pos_fine_tune = 0.005):
    # print(list(tip_list))
    # print(tip_dict)
    for node in tip_list:
        tiplabel=tip_dict[node]
        plt.text(x = x_max + x_pos_fine_tune, y = node - 0.2, s = tiplabel)

# NODE LABEL
def write_node_label(node_list = node_list, tip_node_label_dict = tip_node_label_dict):
    # print(list(tip_list))
    # print(tip_dict)
    for node in node_list:
        nodelabel=tip_node_label_dict[node]
        #plt.text(x = node_x(node), y = node_y(node) - 0.5, s = str(node))
        plt.text(x = node_x(node), y = node_y(node) - 0.5, s = nodelabel)

# This was the function of phy_tree_plot.py
# def plot_C_to_P(node):
#     p_node = C_to_P_dict[node]
#     n_x = node_x(node)
#     n_y = node_y(node)
#     np_x = node_x(p_node)
#     np_y = node_y(p_node)
     #print(node, p_node)
     #print(node, p_node)
     #plt.plot([n_x, np_x], [n_y, np_y])
     #plt.show()

    # C_P_horiz
#     plt.plot([n_x, np_x], [n_y, n_y], c = "black", linewidth = 1.5)
    # C_P_vert
#     plt.plot([np_x, np_x], [n_y, np_y], c = "black", linewidth = 1.5)

# This was the function of phy_tree_plot.py (REFSEQ(+) tree, no signature decoration. etc)
# def plot_top_to_root(node, root_len = root_len):
#     n_x = node_x(node)
#     n_y = node_y(node)
#     root_x = n_x - root_len
#     root_y = n_y
#     plt.plot([n_x, root_x], [n_y, root_y], c = "black", linewidth = 1)


# if node in child_list
def C_to_P_abs_sig(node, 
                   child_list = child_list, tip_dict = tip_dict, C_to_P_dict = C_to_P_dict,
                   tip_node_label_dict = tip_node_label_dict, 
                   nodelabel_count_dict = nodelabel_count_dict, nodelabel_relsig_dict = nodelabel_relsig_dict):
    if node not in child_list:
        #print(node, " is the top node.")
        #print("p_node is the REFSEQ")
        p_node = [node for node in tip_list if tip_dict[node] == "REFSEQ"][0]
        # use abs 
    else:
        p_node = C_to_P_dict[node]
    
    node_lab, p_node_lab = tip_node_label_dict[node], tip_node_label_dict[p_node]
    
    node_abs_sig = nodelabel_count_dict[node_lab] * nodelabel_relsig_dict[node_lab]
    p_node_abs_sig = nodelabel_count_dict[p_node_lab] * nodelabel_relsig_dict[p_node_lab]
    C_P_sig_diff = node_abs_sig - p_node_abs_sig
    
    # dicrepancy handling (the case where the parent node has more count of a certain signature than the child node.)
    min_discrep = C_P_sig_diff.min()
    if (min_discrep <= -16) and (abs(min_discrep) >= nodelabel_count_dict[node_lab] * 0.05): 
        print("!!! Warning !!!")
        print(p_node_abs_sig)
        print(node_abs_sig)
        print(C_P_sig_diff)
        print("C_P disparity exists")
        print("#####################\n")
    
    C_P_sig_diff_no_negative = np.where(C_P_sig_diff < 0, 0, C_P_sig_diff)
    return C_P_sig_diff_no_negative


def plot_sig_len(sig_len, x = 0, y = 0, linewidth = 5):
    start_x = x
    for num,sig in zip(sig_len, sig_list):
        color = sig_colors_dict[sig]
        #print(start_x, num, sig, color)
        plt.plot([start_x, start_x + num] ,[y, y], color = color, linewidth = linewidth, solid_capstyle="butt")
        start_x += num

def plot_C_to_P_sig_horiz(node, linewidth = 10):
    p_node = C_to_P_dict[node]
    X, X_p = node_x(node), node_x(p_node)
    y = node_y(node)
    sig_len = C_to_P_sig_len(node)
    #print("X, X_p, y", X, X_p, y)
    plot_sig_len(sig_len, x = X_p, y = y, linewidth = linewidth)


def C_to_P_sig_len(node):
#for node in [73]:
    #print("\n")
    #print(node)
    #continue
    p_node = C_to_P_dict[node]
    X, X_p = node_x(node), node_x(p_node)
    C_P_len = X - X_p
    nodelabel, p_nodelabel = tip_node_label_dict[node], tip_node_label_dict[p_node]
    #Y = node_y(node)
 
    

    #node_lab = tip_node_label_dict[node]
    #p_node_lab = tip_node_label_dict[p_node]

    #node_abs_sig = nodelabel_count_dict[node_lab] * nodelabel_relsig_dict[node_lab]
    #p_node_abs_sig = nodelabel_count_dict[p_node_lab] * nodelabel_relsig_dict[p_node_lab]


    #print(Y, X, X_p, node_lab, p_node_lab, node_abs_sig, p_node_abs_sig)
    #print("node ", node, tip_node_label_dict[node])
    #print("p_node ", p_node, tip_node_label_dict[p_node])
    #print(Y, X, X_p, X - X_p)
    #print(sig_list)
    #print([sig_colors_dict[sig] for sig in sig_list])
    C_to_P_abs_sig_arr = C_to_P_abs_sig(node)
    #print(X, X_p, C_P_len, C_to_P_abs_sig_arr)
    
    
    if (C_P_len == 0):
        return np.zeros(len(sig_list))
    elif nodelabel == "REFSEQ":
        return np.zeros(len(sig_list))
    elif C_to_P_abs_sig_arr.sum() == 0:
        return np.zeros(len(sig_list))
    #    print("Error raised.")
    #    print("node: ", node, ", parent_node: ", p_node )
    #    print("nodelabel: ", nodelabel, ", parent_node: ", p_nodelabel )
    #    raise ValueError("The sig count of the branch is zero. Check the code or the orignal tree shape.")
 
    C_to_P_rel_sig_arr = C_to_P_abs_sig_arr / C_to_P_abs_sig_arr.sum()
    C_P_sig_len = C_to_P_rel_sig_arr * C_P_len
    return C_P_sig_len

# This is for the children node (not top node!)
def plot_C_to_P_sig(node, linewidth_horiz = 10, linewidth_vert = 1.5):
    p_node = C_to_P_dict[node]
    n_x = node_x(node)
    n_y = node_y(node)
    np_x = node_x(p_node)
    np_y = node_y(p_node)
    #print(node, n_x, n_y)
    #print(p_node, np_x, np_y)
    
    # C_P_horiz
    plot_C_to_P_sig_horiz(node, linewidth = linewidth_horiz)
    # C_P_vert
    plt.plot([np_x, np_x], [n_y, np_y], c = "gray", linewidth = linewidth_vert)

# TIP_LABEL without REFSEQ
def write_tip_label_no_refseq(x_max, tip_list = tip_list, tip_dict = tip_dict, x_pos_fine_tune = 0.005):
    # print(list(tip_list))
    # print(tip_dict)
    for node in tip_list:
        tiplabel=tip_dict[node]
        if tiplabel != "REFSEQ":
            plt.text(x = x_max + x_pos_fine_tune, y = node - 0.2, s = tiplabel)

# The parent of the top_node is the conception timing. x_start -> 0
def plot_top_node_sig_horiz(top_node, x_start = 0, linewidth = 10, 
                            tip_node_label_dict = tip_node_label_dict,
                            nodelabel_count_dict = nodelabel_count_dict,
                            nodelabel_relsig_dict = nodelabel_relsig_dict):
    node_lab = tip_node_label_dict[top_node]
    y = node_y(top_node)
    sig_len = nodelabel_count_dict[node_lab] * nodelabel_relsig_dict[node_lab]
    plot_sig_len(sig_len, x = x_start, y = y, linewidth = linewidth)



#for node in child_list:
    #p_node = C_to_P_dict[node]
    #X, X_p = node_x(node), node_x(p_node)
    #sig_len = C_to_P_sig_len(node)
    #print(node, p_node, X, X_p, sig_len)
    #break
    #print(node, C_to_P_dict[node], C_to_P_sig_len(node), file = open("stderr", "a"))

# This part is the confirmation point of REFSEQ removal.
print("REFSEQ_problem")
print("Going to remove REFSEQ and find the good way for the tree plot.")

print("REFSEQ node index", [node for node in tip_list if tip_dict[node] == "REFSEQ"])
refseq_node = [node for node in tip_list if tip_dict[node] == "REFSEQ"][0]
print("Parent of the REFSEQ (index): ", C_to_P_dict[refseq_node])

P_refseq_node = C_to_P_dict[refseq_node]
print("Children of the REFSEQ_Parent (index): ", P_to_paired_C_dict[P_refseq_node])

# The refseq_parent_node was set as top with figtree. This part is the check for that status.
print("P_refseq_node in child_list: ", P_refseq_node in child_list)

if ((P_refseq_node in child_list) is False) and len(P_to_paired_C_dict[P_refseq_node]) >= 3:
    print(P_refseq_node, " is top node.")
    top_node = P_refseq_node
elif ((P_refseq_node in child_list) is False) and len(P_to_paired_C_dict[P_refseq_node]) == 2:
    print("Children of the REFSEQ_Parent (index): ", P_to_paired_C_dict[P_refseq_node])
    print("Children of the REFSEQ_Parent (index) other than REFSEQ is only one. The top node is ", [elm for elm in P_to_paired_C_dict[P_refseq_node] if elm != refseq_node])
    top_node = [elm for elm in P_to_paired_C_dict[P_refseq_node] if elm != refseq_node][0]


# When we try to remove refseq (this was assumed as one of the tips in this way.)
# --> REFSEQ will have a parent (would be the top node. In this tree [53]) 
# --> Top node have paired children (or three???) [REFSEQ, another, one more] (This can be found with P_to_paired_C_dict)

# If the top node has two children other than REFSEQ, top node can be considered as top as it is after removal of REFSEQ.

# If the top node has only one children other than REFSEQ, the only left children would be the next top after removal or REFSEQ.
# In this situation, the X count from the ROOT (variant counts) is the same as it is. the ROOT_TOP line should be plotted for the next top_root line. The older top node should be removed. 


x_max_node = max([node_x(node) for node in tip_node_list])


plt.figure(figsize=(figsize_horiz, figsize_vert))

if tip_labelling == "on":
    #write_tip_label(x_max_node, x_pos_fine_tune = 0.005)
    write_tip_label_no_refseq(x_max_node, x_pos_fine_tune = 5)
    
if node_labelling == "on":
    write_node_label()

#for node in child_list:
for node in tip_node_list:
    if node in [refseq_node, top_node]:
        print("This is refseq or top node: ", node, tip_node_label_dict[node])
        pass
    else:
        plot_C_to_P_sig(node, linewidth_horiz = linewidth_horiz, linewidth_vert = linewidth_vert)

plot_top_node_sig_horiz(top_node, linewidth = linewidth_horiz)


plt.title(PID_prefix, fontsize = 80)
plt.xticks([])
plt.yticks([])


ax = plt.gca()
#Hide all four spines
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)


legend_elements = []

# Fill the list with patches 
for item, color in sig_colors_dict.items():
    legend_elements.append(mpatches.Patch(color=color, label=item))

# Create the legend
plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1.05), fontsize = 30)

plt.savefig(phylo_fig_out_prefix + "." + figfmt, format = figfmt, transparent=True, bbox_inches='tight')
#plt.show()
    #break
