# View phylogeny
#  Plot, rename tips, root....

# Library ----
library(ape)

# Vars ----
outgroup_pattern <- '^castor_'
wd <- file.path(getwd(), 'pipeline')
input_file <- file.path(wd, '6_supertree', 'supertree_2.tre')

# Read in ----
tree <- read.tree(file = input_file)
tree$edge.length <- NULL
# root tree
outgroup_tips <- tree$tip.label[grepl(pattern = outgroup_pattern,
                                      x = tree$tip.label, ignore.case = TRUE)]
tree <- root(unroot(tree), outgroup = outgroup_tips)
minpx_width <- 768L
png(filename = 'supertree.png', width = minpx_width, height = minpx_width * 3,
    units = 'px')
plot(x = tree, no.margin = TRUE, cex = 0.75)
dev.off()
