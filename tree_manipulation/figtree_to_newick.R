library(ape)
library(treeio)
library(optparse)

op <- OptionParser()
op <- add_option(op, "--figtree", type='character')
op <- add_option(op, "--newicktree", type='character')
args <- parse_args(op)

figtree.file <- args$figtree
new.tree.file <- args$newicktree

figtree <- read.beast(figtree.file)

tree <- figtree@phylo
data <- figtree@data

tree$node.label[
	as.numeric(data$node) - length(tree$tip.label)
] <- data$label

write.tree(tree, new.tree.file)