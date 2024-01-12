#!/risapps/rhel7/R/4.1.0/lib64/R/bin/Rscript --vanilla

suppressMessages(library(ape))
suppressMessages(library(adephylo))

# library(ape)
# library(adephylo)

args = commandArgs(trailingOnly=TRUE)
tree_file = args[1]

# print(tree_file)
#quit()

#tr <- read.tree("../data/rodent.tre")
tr <- read.tree(tree_file)
#plot(tr);
#nodelabels();
#tiplabels();

edge_length = tr$edge.length
df_edge <- data.frame(edge_length)
write.table(df_edge, file = "edge_length.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
#df_edge

# Convert the vector to a data frame
tip_lable = tr$tip.label 
df_tip <- data.frame(letter = tip_lable)
write.table(df_tip, file = "tip_label.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

tr_edge_mtx = tr$edge
write.table(tr_edge_mtx, file = "tr_edge_mtx.tsv", row.names = FALSE, col.names = c("prox","dist"), sep = "\t")
