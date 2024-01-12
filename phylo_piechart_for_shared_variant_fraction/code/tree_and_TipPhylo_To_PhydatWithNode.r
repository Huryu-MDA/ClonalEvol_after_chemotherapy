#!/risapps/rhel7/R/4.1.0/lib64/R/bin/Rscript --vanilla

suppressMessages(library(phytools))
suppressMessages(library(phangorn))

args = commandArgs(trailingOnly=TRUE)

# BOTH file should not contain dummy ROOT or REFSEQ.
InputTree = args[1]
#InputTree = "PID0002_2010_phyml_trsf_box/ROOTED_QueryMTX_SampleMerged_MTXs_phylo_seq96637_REFSEQaddedAsRoot.newick"
InputPhylip = args[2]
#InputPhylip = "PID0002_2010_phyml_trsf_box/QueryMTX_SampleMerged_MTXs_phylo_seq96637_REFSEQaddedAsRoot.phy"

#FileTag = strsplit(InputTree, "\\.")
#FileTag = FileTag[[1]][1]
out_phydat = args[3]
#out_phydat = paste("_test_for_seq_anc_tips", "fasta", sep = ".")
#print(InputTree)
#print(InputPhylip)
#quit()
#out_phydat

tree = read.tree(InputTree)
phydat = read.phyDat(InputPhylip)

#is.rooted(tree)

#tree$edge

tree = unroot(tree)
#tree$node.label = as.character(c(1:(tree$Nnode)) + length(tree$tip.label))
#tree$node.label = paste("n-", tree$node.label, sep = "")

fit <- pml(tree, phydat)
fit <- optim.pml(fit, model="GTR", control = pml.control(trace=0))
anc.ml <- ancestral.pml(fit, "ml", return = "phyDat")

write.phyDat(anc.ml, file = out_phydat, format = "fasta")


