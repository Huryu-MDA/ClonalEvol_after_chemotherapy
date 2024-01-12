import sys, os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

root_len = 0.00
print("root_len was set as " + str(root_len) + ". This should be changed in the future setting.")


def retree(tree):
    edge_df = tree.edge_df
    tip_list = tree.tip_list
    tip_node_label_dict = tree.tip_node_label_dict
    tip_label = [tip_node_label_dict[node] for node in tip_list]
    return PhyloTree(edge_df, tip_label)

def REFSEQ_set_as_root(tree):
    edge_df = tree.edge_df
    tip_list = tree.tip_list
    tip_node_list = tree.tip_node_list
    tip_node_label_dict = tree.tip_node_label_dict
    tip_label = [tip_node_label_dict[node] for node in tip_list]
    C_to_P_dict         = tree.C_to_P_dict
    C_P_len_dict        = tree.C_P_len_dict
    
    refseq_node = [node for node in tip_node_list if tip_node_label_dict[node].upper() == "REFSEQ"][0]
    P_refseq_node = C_to_P_dict[refseq_node]
    refseq_to_topnode_len = C_P_len_dict[refseq_node]
    top_refseq_raw = edge_df[(edge_df["prox"] == P_refseq_node) & (edge_df["dist"] == refseq_node)]
    to_be_dropped_index = top_refseq_raw.index
    top_to_refseq_drop = edge_df.drop(to_be_dropped_index).copy()
    root_to_top_raw = pd.DataFrame([[refseq_node, P_refseq_node, refseq_to_topnode_len]],
                                   columns=["prox", "dist", "length"])

    refseq_root_edge_df = pd.concat([top_to_refseq_drop, root_to_top_raw], ignore_index=True).copy()
    if "REFSEQ" in tip_label:
        renewed_tip_label = [elm for elm in tip_label if elm != "REFSEQ"]
    return PhyloTree(refseq_root_edge_df, renewed_tip_label)

class PhyloTree:
    def __init__(self, edge_df, tip_label = False):
        edge_df.columns= ["prox", "dist", "length"]
        # The order in the edge_df is considered.
        self.edge_df       = edge_df
        self.parent_list   = list(edge_df["prox"])
        self.child_list    = list(edge_df["dist"])
        self.edge_len_list = list(edge_df["length"])
        self.tip_list      = [elm for elm in self.child_list if elm not in self.parent_list]
        
        # The order in the edge_df is not considered.
        self.tip_node_list = sorted(set(set(edge_df["prox"]) | set(edge_df["dist"])))

        # TIP_NODE_LABEL_DICT
        if tip_label == False:
            self.tip_node_label_dict = {tip: "t-" + str(tip) for tip in self.tip_list}
        else:
            self.tip_node_label_dict = {key: value for key, value in zip(self.tip_list, tip_label)}
        for node in self.tip_node_list:
            if node in self.tip_list:
                continue
            else:
                self.tip_node_label_dict[node] = "n-" + str(node)

        # The tools to connect the C-P association, and their length.
        self.C_to_P_dict   = {c_node: p_node for c_node, p_node in zip(self.child_list, self.parent_list)}
        self.P_to_paired_C_dict = {num: list(edge_df[edge_df["prox"] == num]["dist"]) for num in self.parent_list}
        self.C_P_len_dict  = {node: edge_len for node, edge_len in zip(self.child_list, self.edge_len_list)}

    def node_x(self, node):
        C_to_P_dict = self.C_to_P_dict
        child_list  = self.child_list
        C_P_len_dict = self.C_P_len_dict
        if node in child_list:
            p_node = C_to_P_dict[node]
            p_node_x = self.node_x(p_node)
            c_to_p_len = C_P_len_dict[node]
            x = p_node_x + c_to_p_len
            return x
        else:
            # in this case, node_num is the toppest node in the tree.
            x = root_len
            return x

    def node_y(self, node):
        tip_list = self.tip_list
        P_to_paired_C_dict = self.P_to_paired_C_dict
        #pass
        if node in tip_list:
            y = node
            return y
        else:
            child_pair = P_to_paired_C_dict[node]
            c_node_y_list = [self.node_y(c_node) for c_node in child_pair]
            y = sum(c_node_y_list) / len(c_node_y_list)
            return y

    # TIP LABEL
    def write_tip_label(self, fontsize = 14, ax = plt):
        tip_list = self.tip_list
        tip_node_label_dict = self.tip_node_label_dict
        x_max = max([self.node_x(node) for node in tip_list])

        for node in tip_list:
            tiplabel=tip_node_label_dict[node]
            ax.text(x = x_max * 1.05, y = node - 0.2, s = tiplabel)
    
    # NODE LABEL
    def write_node_label(self, fontsize = 14, ax = plt):
        node_list = self.tip_node_list
        tip_node_label_dict = self.tip_node_label_dict

        for node in node_list:
            nodelabel=tip_node_label_dict[node]
            ax.text(x = self.node_x(node), y = self.node_y(node) - 0.5, s = nodelabel, fontsize = fontsize)

    # TIP mark
    def plot_tiplabel_mark(self, tiplabel, ax = plt, c = "black", pos_coeff = 1.05, s = 1, y_adj=0):
        tip_list = self.tip_list
        tip_node_label_dict = self.tip_node_label_dict
        nodelist = [node for node in tip_list if tip_node_label_dict[node] == tiplabel]
        if len(nodelist) != 1:
            print("Error in the method 'plot_tiplabel_mark'. tiplable: ", tiplabel, "node_list: ", nodelist)
            print("not unique node for tiplabel.")
            sys.exit()
        node = nodelist[0]
        x_max = max([self.node_x(node) for node in tip_list])

        ax.scatter(x = x_max * pos_coeff, y = node - y_adj, marker = "s", c = c, s = s)

    # upstream from node
    def upstream_line(self, node):
        child_list = self.child_list
        C_to_P_dict = self.C_to_P_dict
        if node in child_list:
            p_node = C_to_P_dict[node]
            upper_line = self.upstream_line(p_node) + [node]
            return upper_line
        else:
            return [node]
    
    #check the association between d_node and a_node
    #d_node: descendant node, a_node: ancestor node
    def check_direct_line(self, d_node, a_node):
        upper_line = self.upstream_line(d_node)
        if a_node not in upper_line:
            #print(a_node, "is not an direct ancestor of", d_node)
            raise ValueError(a_node, "is not an direct ancestor of", d_node)
        else:
            #print(a_node, "is an direct ancestor of", d_node)
            pass
    
    #Plot function from c_node to p_node
    def plot_tree(self, ax = plt, 
                  color_horiz = "k", color_vert = "k", width_horiz = 1, width_vert = 1):
        child_list = self.child_list
        for node in child_list:
            self.plot_C_to_P(node, ax = ax, 
                             color_horiz = color_horiz, color_vert = color_vert, 
                             width_horiz = width_horiz, width_vert = width_vert)

    def plot_C_to_P(self, node, ax = plt,
                    color_horiz = "k", color_vert = "k", width_horiz = 1, width_vert = 1):
        self.plot_C_to_P_horiz(node, ax = ax, color = color_horiz, width = width_horiz)
        self.plot_C_to_P_vert(node, ax = ax, color = color_vert, width = width_vert)

    def plot_C_to_P_horiz(self, node, ax = plt, color = "k", width = 1):
        C_to_P_dict = self.C_to_P_dict
        p_node = C_to_P_dict[node] 
        n_x, np_x = self.node_x(node), self.node_x(p_node)
        n_y, np_y = self.node_y(node), self.node_y(p_node)
        ax.plot([n_x, np_x], [n_y, n_y], c = color, linewidth = width)
        #if ax == False:
        #    plt.plot([n_x, np_x], [n_y, n_y], c = color, linewidth = width)
        #else:
        #    ax.plot([n_x, np_x], [n_y, n_y], c = color, linewidth = width)

    def plot_C_to_P_vert(self, node, ax = plt, color = "k", width = 1):
        C_to_P_dict = self.C_to_P_dict
        p_node = C_to_P_dict[node]
        n_x, np_x = self.node_x(node), self.node_x(p_node)
        n_y, np_y = self.node_y(node), self.node_y(p_node)
        ax.plot([np_x, np_x], [n_y, np_y], c = color, linewidth = width)
        #if ax == False:
        #    plt.plot([np_x, np_x], [n_y, np_y], c = color, linewidth = width)
        #else:
        #    ax.plot([np_x, np_x], [n_y, np_y], c = color, linewidth = width)

    def REFSEQ_set_as_root(self):
        edge_df = self.edge_df
        tip_list = self.tip_list
        tip_node_list       = self.tip_node_list
        tip_node_label_dict = self.tip_node_label_dict
        tip_label = [tip_node_label_dict[node] for node in tip_list]
        C_to_P_dict         = self.C_to_P_dict
        C_P_len_dict        = self.C_P_len_dict

        refseq_node = [node for node in tip_node_list if tip_node_label_dict[node].upper() == "REFSEQ"][0]
        P_refseq_node = C_to_P_dict[refseq_node]
        refseq_to_topnode_len = C_P_len_dict[refseq_node]
        top_refseq_raw = edge_df[(edge_df["prox"] == P_refseq_node) & (edge_df["dist"] == refseq_node)]
        to_be_dropped_index = top_refseq_raw.index
        top_to_refseq_drop = edge_df.drop(to_be_dropped_index).copy()

        root_to_top_raw = pd.DataFrame([[refseq_node, P_refseq_node, refseq_to_topnode_len]], 
                                       columns=["prox", "dist", "length"])

        refseq_root_edge_df = pd.concat([top_to_refseq_drop, root_to_top_raw], ignore_index=True).copy()
        
        if "REFSEQ" in tip_label:
            renewed_tip_label = [elm for elm in tip_label if elm != "REFSEQ"]
        return tree_plot_test.PhyloTree(refseq_root_edge_df, renewed_tip_label)
    #     n_x = node_x(node)
    #     n_y = node_y(node)
    #     np_x = node_x(p_node)
    #     np_y = node_y(p_node)
    #     print(node, p_node)
    #     print(node, p_node)
    #     plt.plot([n_x, np_x], [n_y, np_y])
    #     plt.show()

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
