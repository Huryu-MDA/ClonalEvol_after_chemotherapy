#!/risapps/rhel7/R/4.1.0/lib64/R/bin/Rscript --vanilla

suppressMessages(library(phytools))
suppressMessages(library(phangorn))

args = commandArgs(trailingOnly=TRUE)

# BOTH file should not contain dummy ROOT or REFSEQ.
InputTree = args[1]
# InputTree = "PID0002_2010_phyml_trsf_box/ROOTED_QueryMTX_SampleMerged_MTXs_phylo_seq96637_REFSEQaddedAsRoot.newick"

shared_var_frac_inOrder = args[2]
sample_prefix = args[3]

output_tree_plot = paste(sample_prefix, "ref_drop_tree_with_sharedVarsWithPie", "pdf", sep = ".")
#print(output_tree_plot)
#quit()

#shared_var_frac_inOrder = "shared_var_frac_in_order.tsv"

tree = read.tree(InputTree)
#phydat = read.phyDat(InputPhylip)

shared_var_frac_df = read.csv(shared_var_frac_inOrder, header = FALSE, sep ="\t")

frac_vec = shared_var_frac_df[,2]

tip_len  = length(tree$tip.label)
node_len = length(tree$node.label)

#node_len

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

pdf(file = output_tree_plot, width = 60,height = 60)
#dev.new(width = 10, height = 10)
plot(ref_drop_tree, show.tip.label = FALSE, edge.width = 10, no.margin = TRUE, root.edge = TRUE)
#plotTree(tree)
nodelabels(pie = node_frac_vec, cex = 0.5, piecol = col_vec)
tiplabels(pie = tip_frac_vec_after_ref_drop, cex = 0.5, piecol = col_vec)
#nodelabels(thermo = node_frac_vec[1:116], cex = 0.3, horiz = TRUE, piecol = col_vec);
#tiplabels(thermo = tip_frac_vec, cex = 0.3, horiz = TRUE, piecol = col_vec);

dev.off()

quit()

#plotTree(tree)

#pdf(file = "node_label.pdf",width = 30,height = 30)
#plotTree(tree)
#nodelabels()
#dev.off()

#pdf(file = "test_tree.pdf",   # The directory you want to save the file in
#    width = 24, # The width of the plot in inches
#    height = 24) # The height of the plot in inches

#plotTree(tree)
#nodelabels(pie = node_frac_vec, cex = 0.2)
#tiplabels(pie = tip_frac_vec, cex = 0.2)

#dev.off()
#tiplabels(pch = 21, bg = gray(1:23/23), cex = 2, adj = 1.4)
#tiplabels(pch = 19, col = c("yellow", "red", "blue"), adj = 2.5, cex = 2)


