#!/risapps/rhel7/R/4.1.0/lib64/R/bin/Rscript --vanilla

suppressMessages(library(phytools))
suppressMessages(library(phangorn))

args = commandArgs(trailingOnly=TRUE)

prefix = args[1]
# BOTH file should not contain dummy ROOT or REFSEQ.
InputTree = args[2]
#InputTree = "ROOTED_QueryMTX_SampleMerged_MTXs_phylo_seq77917_REFSEQaddedAsRoot.newick"

shared_var_frac_inOrder = args[3]
#shared_var_frac_inOrder = "shared_var_frac_in_order_PID0010.tsv"

#x_max = as.double(args[3])
x_max = strtoi(args[4])
#y_max = strtoi(args[4])
#pie_size = as.double(args[5])
pie_size = as.double(args[5])
#output_pdf=paste("tree_plot_on_topNodes","xmax",args[3],"ymax",args[4],"piesize",args[5],"pdf", sep =".")
#output_pdf=paste("tree_plot_on_topNodes","xmax",args[3],"piesize",args[4],"pdf", sep =".")
output_pdf=paste(prefix,"tree_plot_on_topNodes","xmax",args[4],"piesize",args[5],"pdf", sep =".")

#print(InputTree)
#print(shared_var_frac_inOrder)
#print(x_max)
#print(y_max)
#print(output_pdf)
#quit()

tree = read.tree(InputTree)
shared_var_frac_df = read.csv(shared_var_frac_inOrder, header = FALSE, sep ="\t")

frac_vec = shared_var_frac_df[,2]

tip_len  = length(tree$tip.label)
node_len = length(tree$node.label)

tip_frac_vec = frac_vec[1:tip_len]
node_frac_vec = frac_vec[(tip_len + 1):(tip_len + node_len)]

# remove the REFSEQ tip from tree, then convert the REFSEQ edgelength to root_edge_length
REFSEQ_ind = which(tree$tip.label == "REFSEQ")
refseq_egde_index = which(tree$edge[,2] == REFSEQ_ind)
root_edge_len = tree$edge.length[refseq_egde_index]
ref_drop_tree = drop.tip(tree, "REFSEQ")
ref_drop_tree$root.edge = root_edge_len

# Remove the index of REFSEQ from tip_frac_vector
tip_frac_vec_after_ref_drop = tip_frac_vec[tree$tip.label != "REFSEQ"]

col_vec = c("red", "white")

pdf(file = output_pdf  , width = 60,height = 60)
#plot(ref_drop_tree, edge.width = 10, x.lim = (c(0, x_max)), y.lim = (c(0, y_max)), no.margin = TRUE, root.edge = TRUE)
plot(ref_drop_tree, edge.width = 10, x.lim = (c(0, x_max)), no.margin = TRUE, root.edge = TRUE)
nodelabels(pie = node_frac_vec, cex = pie_size, piecol = col_vec)
tiplabels(pie = tip_frac_vec_after_ref_drop, cex = pie_size, piecol = col_vec)
dev.off()

