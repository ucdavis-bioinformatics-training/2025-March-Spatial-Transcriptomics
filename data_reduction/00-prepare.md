---
title: "Prepare R environment for analysis"
author: "UC Davis Bioinformatics Core"
date: "2025-03-23"
output: 
  html_document:
    keep_md: TRUE
    toc: TRUE
---


### Create a new RStudio project

Open RStudio and create a new project, for more info see [Using-Projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects):

*File > New Project > New Directory > New Project*

Name the new directory (e.g. Spatial_transcriptomics), and check "use renv with this project" if present.

Learn more about [renv](https://rstudio.github.io/renv/articles/renv.html).

### Install packages
One of R's many benefits is the large, active user community, which produces and maintains many packages that extend the functionality of base R and provide functions that enable bioinformatic analyses without completely custom code.

The following package installation commands should be run individually, **in the R console**. Many of them will require your input to determine which, if any, dependencies should be updated; for the quickest result, attempt 'n' (none) first.

#### R-universe for arm64 installations

r-universe is a new umbrella project by __rOpenSci__. It uses cross-compiling for arm64 binaries.

<span style="color:blue">For those who are using Macs that have M1/M2/M3 chips, if you have trouble installing the packages and get error that is similar to "ld: warning: ignoring file '/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libR.dylib': found architecture 'arm64', required architecture 'x86_64'", please go to https://r-universe.dev/search/ and search for the packages and use the installation instructions provided there.</span>

#### BiocManager
BiocManager is an interface for the bioinformatics-specific R package repository. We will be using BiocManager to install other packages when possible, rather than the base R function install.packages.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}
```

#### rmarkdown
The rmarkdown package, when used with others like tinytex and knitr, allows you to knit your Rmd document to nicely-formatted reports.

``` r
if (!any(rownames(installed.packages()) == "rmarkdown")){
  BiocManager::install("rmarkdown")
}
library(rmarkdown)
```

#### tinytex
TinyTeX is a small LaTeX distribution for use with R.

``` r
if (!any(rownames(installed.packages()) == "tinytex")){
  BiocManager::install("tinytex")
}
library(tinytex)
```

#### knitr

``` r
if (!any(rownames(installed.packages()) == "knitr")){
  BiocManager::install("knitr")
}
library(knitr)
```

#### kableExtra
The kableExtra package gives the user fine-grained control over table formats. This is useful for knit reports.

``` r
if (!any(rownames(installed.packages()) == "kableExtra")){
  BiocManager::install("kableExtra")
}
library(kableExtra)
```

#### ggplot2
An extremely popular package by the authors of RStudio, ggplot2 produces highly customizable plots.

``` r
if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}
library(ggplot2)
```

#### dplyr
Like ggplot2 and tidyr, dplyr is part of the "tidyverse" by the RStudio authors: a group of packages designed for data analysis and visualization.

``` r
if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}
library(dplyr)
```

#### tidyr

``` r
if (!any(rownames(installed.packages()) == "tidyr")){
  BiocManager::install("tidyr")
}
library(tidyr)
```

#### viridis
viridis produces accessible color palettes.

``` r
if (!any(rownames(installed.packages()) == "viridis")){
  BiocManager::install("viridis")
}
library(viridis)
```

#### hdf5r
HDF5 (heirarchical data format version five) files can be used to store single cell expression data (including output from Cell Ranger). The hdf5r package provides utilities for interacting with the format.

``` r
if (!any(rownames(installed.packages()) == "hdf5r")){
  BiocManager::install("hdf5r")
}
library(hdf5r)
```

#### Seurat
Seurat is an extensive package for the analysis of single cell experiments, from normalization to visualization.

``` r
if (!any(rownames(installed.packages()) == "Seurat")){
  BiocManager::install("Seurat")
}
library(Seurat)
```

#### sf
sf is a package that supports simple features, a standardized way to encode spatial vector data.

``` r
if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  install.packages("sf")
}
```

#### ComplexHeatmap
ComplexHeatmap produces beautiful, highly-customizable heat maps.

``` r
if (!any(rownames(installed.packages()) == "ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)
```

#### biomaRt
This package provides an interface to Ensembl databases.

``` r
if (!any(rownames(installed.packages()) == "biomaRt")){
  BiocManager::install("biomaRt")
}
library(biomaRt)
```

#### limma
Originally developed for microarray data, limma provides functions for linear modeling and differential expression.

``` r
if (!any(rownames(installed.packages()) == "limma")){
  BiocManager::install("limma")
}
library(limma)
```

#### remotes
Some packages (or versions of packages) cannot be installed through Bioconductor. The remotes package contains tools for installing packages from a number of repositories, including GitHub.

``` r
if (!any(rownames(installed.packages()) == "remotes")){
  utils::install.packages("remotes")
}
library(remotes)
```
#### NicheDE
NicheDE performs niche-differential gene expression analysis, i.e. tests if the expression of a given cell type is influenced by being spatially close to another celltype.

``` r
if (!any(rownames(installed.packages()) == "NicheDE")){
  options(timeout=9999999)
  remotes::install_github("kaishumason/NicheDE")
}
library(nicheDE)
```

#### spacexr
spacexr implements computations methods for cell type identification and differential expression on spatial transcriptomics datasets.

``` r
## may need to modify the download time out limit for the installation to be successful
#options(timeout = 600000000)
if (!any(rownames(installed.packages()) == "spacexr")){
remotes::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
```

#### scCustomize
scCustome is a package with collection of functions created and/or curated to aid in visualization and analysis of single-cell data in R

``` r
if (!any(rownames(installed.packages()) == "scCustomize")){
  utils::install.packages("scCustomize")
}
library(scCustomize)
```

#### ape
Analysis of Phylogenetics and Evolution (ape) is used to generate and manipulate phylogenetic trees. In this workshop, we will be using ape to investigate the relationships between clusters.

``` r
if (!any(rownames(installed.packages()) == "ape")){
  utils::install.packages("ape")
}
library(ape)
```

#### openxlsx
The openxlsx package is a suite of tools for reading and writing .xlsx files.

``` r
if (!any(rownames(installed.packages()) == "openxlsx")){
  BiocManager::install("openxlsx")
}
library(openxlsx)
```


#### Verfiy installation
Finally, we can get the session info to ensure that all of the packages were installed and loaded correctly.

``` r
sessionInfo()
```


