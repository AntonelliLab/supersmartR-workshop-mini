# Build a phylogenetic tree Caviomorpha

The analysis is broken up into species-level and backbone phylogenetic trees.

## Process

* `0_restez.R`: Download sections of GenBank
* `1_phylotaR.R`: Fetch sequence clusters
* `2_clusters.R`: Select multiple clusters
* `3_align.R`: Align sequences
* `4_supermatrix.R`: Assemble supermatrices (mono. + backbone)
* `5_phylogeny.R`: Generate trees
* `6_supertree.R`: Assemble the supertree
* `7_view.R`: Plot the resulting tree!

**Notes.**

* You can download the completed first step `1_phylotaR`, see [data/](https://github.com/AntonelliLab/supersmartR-workshop-mini/data)
* To save time, we are restricting sequence retrieval to the `restez`
database (`db_only = TRUE`).

## Target taxon

*A large clade < 500 spp.*

In the default pipeline, we are using [**Hystricomorpha**](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10015),
guinea pigs and allies.

![](https://upload.wikimedia.org/wikipedia/commons/f/fc/Two_Adult_Guinea_Pigs_%28cropped%29.jpg)
