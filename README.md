# Mini `supersmartR` Workshop <img src="https://raw.githubusercontent.com/AntonelliLab/supersmartR/master/logo.png" height="300" align="right"/>

> Originally put together for a mini-workshop on 8 Nov. 2019 in  Gothenburg,
Sweden.

[`supersmartR`](https://github.com/AntonelliLab/supersmartR) is a series of
R packages that form a phylogenetic pipeline.

## <img src="https://raw.githubusercontent.com/AntonelliLab/supersmartR/master/supersmart%20vs%20supersmartr.png" height="450" align="middle"/>

This workshop we will introduce the packages and provide code to run a simple
pipeline to create a phylogenetic tree.

Workshop duration ~1 hr.

* * *

# Prerequisites

* Software
    * [R](https://cran.r-project.org/) (> 3.5)
    * [RStudio](https://www.rstudio.com/)
    * [Desktop Docker](https://docs.docker.com/install/) (Linux containers)
* Basic knowledge
    * R
    * Phylogenetics

(Windows users may struggle installing Docker Desktop, in which case Docker
Toolbox will also work.)

### R packages

Install dependant packages.

```r
# remotes allows installation via GitHub
if (!'remotes' %in% installed.packages()) {
  install.packages('remotes')
}

library(remotes)
# install latest outsider
install_github("antonellilab/outsider.base")
install_github("antonellilab/outsider")
# install latest restez
install_github("hannesmuehleisen/MonetDBLite-R")
install_github("ropensci/restez")
# install latest phylotaR
install_github("ropensci/phylotaR")
# install latest gaius
install_github("antonellilab/gaius")
```

If you're computer is set-up correctly with the R packages and Docker, you
should be able to run the below script without errors.

```r
library(outsider)
repo <- "dombennett/om..hello.world"
module_install(repo = repo, force = TRUE)
hello_world <- module_import(fname = "hello_world", repo = repo)
hello_world()
```

* * *

# Tutorials

## `supersmartR` packages

### `phylotaR`



### `restez`

### `outsider`

### `gaius`

## Pipeline

## Extra

### `ssh`

***

**Learn more: [`supersmartR` GitHub Page](https://github.com/AntonelliLab/supersmartR)**