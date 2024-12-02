---
title: "Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression"
author: "Bioinformatics Core"
date: "2024-12-02"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Introduction to Single Cell RNA-Seq Part 6: Enrichment and model-based differential expression


## Set up workspace

``` r
library(Seurat)
library(limma)
library(topGO)
library(dplyr)
library(kableExtra)
set.seed(12345)
experiment.aggregate <- readRDS("scRNA_workshop-05.rds")
Idents(experiment.aggregate) <- "finalcluster"
```

## 1. Gene Ontology (GO) Enrichment of Genes Expressed in a Cluster
[Gene Ontology](http://geneontology.org/docs/ontology-documentation/) provides a controlled vocabulary for describing gene products.  Here we use enrichment analysis to identify GO terms that are over-represented among the gene expressed in cells in a given cluster. 

``` r
cluster10 <- subset(experiment.aggregate, idents = '10')
expr <- as.matrix(GetAssayData(cluster10))

# Select genes that are expressed > 0 in at least half of cells
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

# Create topGOdata object
	GOdata <- new("topGOdata",
		ontology = "BP", # use biological process ontology
		allGenes = geneList,
		geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
# Test for enrichment using Fisher's Exact Test
	resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
	GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> GO.ID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Term </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Annotated </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Significant </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Expected </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Fisher </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> GO:0006338 </td>
   <td style="text-align:left;"> chromatin remodeling </td>
   <td style="text-align:right;"> 464 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 5.61 </td>
   <td style="text-align:left;"> 2.6e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0043484 </td>
   <td style="text-align:left;"> regulation of RNA splicing </td>
   <td style="text-align:right;"> 148 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 1.79 </td>
   <td style="text-align:left;"> 0.00041 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903613 </td>
   <td style="text-align:left;"> regulation of protein tyrosine phosphatase activity </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.04 </td>
   <td style="text-align:left;"> 0.00043 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006306 </td>
   <td style="text-align:left;"> DNA methylation </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 0.40 </td>
   <td style="text-align:left;"> 0.00063 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0061061 </td>
   <td style="text-align:left;"> muscle structure development </td>
   <td style="text-align:right;"> 381 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 4.61 </td>
   <td style="text-align:left;"> 0.00067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0031580 </td>
   <td style="text-align:left;"> membrane raft distribution </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.05 </td>
   <td style="text-align:left;"> 0.00086 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0009791 </td>
   <td style="text-align:left;"> post-embryonic development </td>
   <td style="text-align:right;"> 65 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.79 </td>
   <td style="text-align:left;"> 0.00110 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0000512 </td>
   <td style="text-align:left;"> lncRNA-mediated post-transcriptional gene silencing </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00142 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:1903978 </td>
   <td style="text-align:left;"> regulation of microglial cell activation </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00142 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0045060 </td>
   <td style="text-align:left;"> negative thymic T cell selection </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.06 </td>
   <td style="text-align:left;"> 0.00142 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:2000627 </td>
   <td style="text-align:left;"> positive regulation of miRNA catabolic process </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.07 </td>
   <td style="text-align:left;"> 0.00211 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0072711 </td>
   <td style="text-align:left;"> cellular response to hydroxyurea </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.08 </td>
   <td style="text-align:left;"> 0.00292 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0071320 </td>
   <td style="text-align:left;"> cellular response to cAMP </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.31 </td>
   <td style="text-align:left;"> 0.00366 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:2000045 </td>
   <td style="text-align:left;"> regulation of G1/S transition of mitotic cell cycle </td>
   <td style="text-align:right;"> 123 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1.49 </td>
   <td style="text-align:left;"> 0.00370 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0000398 </td>
   <td style="text-align:left;"> mRNA splicing, via spliceosome </td>
   <td style="text-align:right;"> 241 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 2.91 </td>
   <td style="text-align:left;"> 0.00377 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0030889 </td>
   <td style="text-align:left;"> negative regulation of B cell proliferation </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.10 </td>
   <td style="text-align:left;"> 0.00387 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0006376 </td>
   <td style="text-align:left;"> mRNA splice site recognition </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.33 </td>
   <td style="text-align:left;"> 0.00408 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0010720 </td>
   <td style="text-align:left;"> positive regulation of cell development </td>
   <td style="text-align:right;"> 263 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 3.18 </td>
   <td style="text-align:left;"> 0.00452 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0060218 </td>
   <td style="text-align:left;"> hematopoietic stem cell differentiation </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.34 </td>
   <td style="text-align:left;"> 0.00453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GO:0010558 </td>
   <td style="text-align:left;"> negative regulation of macromolecule biosynthetic process </td>
   <td style="text-align:right;"> 1370 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> 16.56 </td>
   <td style="text-align:left;"> 0.00463 </td>
  </tr>
</tbody>
</table>
* Annotated: number of genes (out of all.genes) that are annotated with that GO term
* Significant: number of genes that are annotated with that GO term and meet our criteria for "expressed"
* Expected: Under random chance, number of genes that would be expected to be annotated with that GO term and meeting our criteria for "expressed"
* Fisher: (Raw) p-value from Fisher's Exact Test

## 2. Model-based DE analysis in limma
[limma](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R package for differential expression analysis of bulk RNASeq and microarray data.  We apply it here to single cell data.

Limma can be used to fit any linear model to expression data and is useful for analyses that go beyond two-group comparisons.  A detailed tutorial of model specification in limma is available [here](https://ucdavis-bioinformatics-training.github.io/2021-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes) and in the [limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).


``` r
# filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
cluster10$proper.group <- make.names(cluster10$group)
mm <- model.matrix(~0 + proper.group + S.Score + G2M.Score + percent_MT + nFeature_RNA, data = cluster10[[]])
head(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AAACGCTTCTCTGCTG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.2751272 </td>
   <td style="text-align:right;"> 0.8966284 </td>
   <td style="text-align:right;"> 1.1409396 </td>
   <td style="text-align:right;"> 1038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AACAGGGGTCCCTGAG_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0068710 </td>
   <td style="text-align:right;"> 0.0132416 </td>
   <td style="text-align:right;"> 0.7835455 </td>
   <td style="text-align:right;"> 719 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAGCATCCATCCCACT_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0518141 </td>
   <td style="text-align:right;"> 0.0255571 </td>
   <td style="text-align:right;"> 0.5510810 </td>
   <td style="text-align:right;"> 1391 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAGCGAGCACGAGAAC_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0403860 </td>
   <td style="text-align:right;"> -0.0059914 </td>
   <td style="text-align:right;"> 0.4604758 </td>
   <td style="text-align:right;"> 980 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AAGCGTTCAGCCTATA_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0316891 </td>
   <td style="text-align:right;"> 0.0644262 </td>
   <td style="text-align:right;"> 0.7960199 </td>
   <td style="text-align:right;"> 761 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ACGGTTAGTCTCACAA_A001-C-007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0445603 </td>
   <td style="text-align:right;"> 0.7918382 </td>
   <td style="text-align:right;"> 0.5617978 </td>
   <td style="text-align:right;"> 1199 </td>
  </tr>
</tbody>
</table>

``` r
tail(mm) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> TGGGTTACAAGAATGT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0913972 </td>
   <td style="text-align:right;"> -0.1247840 </td>
   <td style="text-align:right;"> 0.2367798 </td>
   <td style="text-align:right;"> 894 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TGTTTGTTCACTACTT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.1148748 </td>
   <td style="text-align:right;"> -0.1289778 </td>
   <td style="text-align:right;"> 0.4276115 </td>
   <td style="text-align:right;"> 1075 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCCGTGTCCGCTGTT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.1008870 </td>
   <td style="text-align:right;"> -0.0119225 </td>
   <td style="text-align:right;"> 0.6284916 </td>
   <td style="text-align:right;"> 1009 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTTCCAGTCCCAAT_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.0325015 </td>
   <td style="text-align:right;"> -0.1043837 </td>
   <td style="text-align:right;"> 0.3673095 </td>
   <td style="text-align:right;"> 822 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTCTTCCCAGCGTAGA_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0226923 </td>
   <td style="text-align:right;"> 0.0085482 </td>
   <td style="text-align:right;"> 0.5474453 </td>
   <td style="text-align:right;"> 850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TTTACGTGTGTCTTAG_B001-A-301 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> -0.0945735 </td>
   <td style="text-align:right;"> -0.1263372 </td>
   <td style="text-align:right;"> 0.5700326 </td>
   <td style="text-align:right;"> 891 </td>
  </tr>
</tbody>
</table>

``` r
# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit)) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupColorectal.Cancer </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupPolyp </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> G2M.Score </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> percent_MT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> nFeature_RNA </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> CCNL2 </td>
   <td style="text-align:right;"> 0.2242247 </td>
   <td style="text-align:right;"> 0.1990672 </td>
   <td style="text-align:right;"> 0.2470341 </td>
   <td style="text-align:right;"> -0.4213528 </td>
   <td style="text-align:right;"> 0.4293688 </td>
   <td style="text-align:right;"> 0.1358057 </td>
   <td style="text-align:right;"> 0.0002361 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11B </td>
   <td style="text-align:right;"> 0.6242932 </td>
   <td style="text-align:right;"> 0.1148185 </td>
   <td style="text-align:right;"> 0.5999028 </td>
   <td style="text-align:right;"> -0.4532588 </td>
   <td style="text-align:right;"> 0.6436532 </td>
   <td style="text-align:right;"> -0.2417947 </td>
   <td style="text-align:right;"> 0.0000213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC35E2B </td>
   <td style="text-align:right;"> 0.0037741 </td>
   <td style="text-align:right;"> 0.0495442 </td>
   <td style="text-align:right;"> 0.0914102 </td>
   <td style="text-align:right;"> 0.5887303 </td>
   <td style="text-align:right;"> 0.0651147 </td>
   <td style="text-align:right;"> 0.1952115 </td>
   <td style="text-align:right;"> 0.0001180 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CDK11A </td>
   <td style="text-align:right;"> 0.5310057 </td>
   <td style="text-align:right;"> -0.1306109 </td>
   <td style="text-align:right;"> -0.0306059 </td>
   <td style="text-align:right;"> -1.0572770 </td>
   <td style="text-align:right;"> 0.9785408 </td>
   <td style="text-align:right;"> -0.0026322 </td>
   <td style="text-align:right;"> 0.0003089 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NADK </td>
   <td style="text-align:right;"> 0.9047053 </td>
   <td style="text-align:right;"> 0.5465035 </td>
   <td style="text-align:right;"> 0.5874302 </td>
   <td style="text-align:right;"> -0.2504437 </td>
   <td style="text-align:right;"> -0.0059736 </td>
   <td style="text-align:right;"> -0.2493937 </td>
   <td style="text-align:right;"> -0.0001804 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GNB1 </td>
   <td style="text-align:right;"> 0.5492680 </td>
   <td style="text-align:right;"> 0.7200203 </td>
   <td style="text-align:right;"> 0.6733315 </td>
   <td style="text-align:right;"> -0.5749704 </td>
   <td style="text-align:right;"> 0.9643106 </td>
   <td style="text-align:right;"> 0.2342349 </td>
   <td style="text-align:right;"> 0.0001687 </td>
  </tr>
</tbody>
</table>

``` r
# Test 'Normal' - 'Colorectal.Cancer'
contr <- makeContrasts(proper.groupNormal - proper.groupColorectal.Cancer, levels = colnames(coef(fit)))
contr %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> proper.groupNormal - proper.groupColorectal.Cancer </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> proper.groupColorectal.Cancer </td>
   <td style="text-align:right;"> -1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupNormal </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> proper.groupPolyp </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> G2M.Score </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> percent_MT </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nFeature_RNA </td>
   <td style="text-align:right;"> 0 </td>
  </tr>
</tbody>
</table>

``` r
fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30) %>%
	  kable() %>%
	  kable_styling("striped", fixed_thead = TRUE)
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">  </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> logFC </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> AveExpr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> t </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> P.Value </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> adj.P.Val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> B </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> XIST </td>
   <td style="text-align:right;"> 2.1155304 </td>
   <td style="text-align:right;"> 0.8383589 </td>
   <td style="text-align:right;"> 13.661092 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 54.324855 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A2 </td>
   <td style="text-align:right;"> 2.2940612 </td>
   <td style="text-align:right;"> 1.2102978 </td>
   <td style="text-align:right;"> 11.450240 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 40.637829 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PIGR </td>
   <td style="text-align:right;"> 2.1206204 </td>
   <td style="text-align:right;"> 1.3048963 </td>
   <td style="text-align:right;"> 9.653033 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 29.658194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SLC26A3 </td>
   <td style="text-align:right;"> 1.6660558 </td>
   <td style="text-align:right;"> 0.7474037 </td>
   <td style="text-align:right;"> 9.247113 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 27.235093 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCND3 </td>
   <td style="text-align:right;"> 2.1487014 </td>
   <td style="text-align:right;"> 1.1948400 </td>
   <td style="text-align:right;"> 8.593783 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 23.403324 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> GUCA2A </td>
   <td style="text-align:right;"> 1.3966234 </td>
   <td style="text-align:right;"> 0.5377794 </td>
   <td style="text-align:right;"> 8.235661 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 21.345829 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PHGR1 </td>
   <td style="text-align:right;"> 1.5875349 </td>
   <td style="text-align:right;"> 0.7790538 </td>
   <td style="text-align:right;"> 8.094118 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 20.542213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MMP12 </td>
   <td style="text-align:right;"> -1.5478178 </td>
   <td style="text-align:right;"> 0.4498059 </td>
   <td style="text-align:right;"> -8.042518 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 20.250687 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CKB </td>
   <td style="text-align:right;"> 1.8461063 </td>
   <td style="text-align:right;"> 1.2010000 </td>
   <td style="text-align:right;"> 7.809511 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 18.944244 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TYMP </td>
   <td style="text-align:right;"> -1.5086803 </td>
   <td style="text-align:right;"> 0.5986687 </td>
   <td style="text-align:right;"> -7.449051 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 16.957641 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PDE3A </td>
   <td style="text-align:right;"> 1.2645878 </td>
   <td style="text-align:right;"> 0.5870357 </td>
   <td style="text-align:right;"> 7.227138 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 15.757169 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MUC12 </td>
   <td style="text-align:right;"> 1.2465258 </td>
   <td style="text-align:right;"> 0.7228761 </td>
   <td style="text-align:right;"> 6.611899 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 2.00e-07 </td>
   <td style="text-align:right;"> 12.529424 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> RNF213 </td>
   <td style="text-align:right;"> -1.5582052 </td>
   <td style="text-align:right;"> 1.5052670 </td>
   <td style="text-align:right;"> -6.586467 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 2.00e-07 </td>
   <td style="text-align:right;"> 12.399425 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MBNL1 </td>
   <td style="text-align:right;"> 1.6531241 </td>
   <td style="text-align:right;"> 2.4693744 </td>
   <td style="text-align:right;"> 6.362969 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 5.00e-07 </td>
   <td style="text-align:right;"> 11.269562 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NXPE1 </td>
   <td style="text-align:right;"> 1.3190223 </td>
   <td style="text-align:right;"> 0.7592336 </td>
   <td style="text-align:right;"> 6.322599 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 6.00e-07 </td>
   <td style="text-align:right;"> 11.067939 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> FKBP5 </td>
   <td style="text-align:right;"> 1.7002914 </td>
   <td style="text-align:right;"> 1.3458515 </td>
   <td style="text-align:right;"> 6.307107 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 6.00e-07 </td>
   <td style="text-align:right;"> 10.990775 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ZG16 </td>
   <td style="text-align:right;"> 0.9881296 </td>
   <td style="text-align:right;"> 0.4012025 </td>
   <td style="text-align:right;"> 6.191486 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.00e-06 </td>
   <td style="text-align:right;"> 10.418473 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PARP14 </td>
   <td style="text-align:right;"> -1.4952975 </td>
   <td style="text-align:right;"> 1.0742512 </td>
   <td style="text-align:right;"> -6.145782 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.20e-06 </td>
   <td style="text-align:right;"> 10.194040 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MT-CO2 </td>
   <td style="text-align:right;"> -1.1821944 </td>
   <td style="text-align:right;"> 2.3699957 </td>
   <td style="text-align:right;"> -6.134752 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 1.20e-06 </td>
   <td style="text-align:right;"> 10.140028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> UTY </td>
   <td style="text-align:right;"> -1.1717486 </td>
   <td style="text-align:right;"> 0.6887650 </td>
   <td style="text-align:right;"> -6.006045 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 2.20e-06 </td>
   <td style="text-align:right;"> 9.514278 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MAML2 </td>
   <td style="text-align:right;"> 1.5838731 </td>
   <td style="text-align:right;"> 1.5289834 </td>
   <td style="text-align:right;"> 5.912160 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 3.30e-06 </td>
   <td style="text-align:right;"> 9.063133 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SON </td>
   <td style="text-align:right;"> -1.4616397 </td>
   <td style="text-align:right;"> 1.2857812 </td>
   <td style="text-align:right;"> -5.900829 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 3.30e-06 </td>
   <td style="text-align:right;"> 9.008987 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CLCA4 </td>
   <td style="text-align:right;"> 1.0303576 </td>
   <td style="text-align:right;"> 0.4594471 </td>
   <td style="text-align:right;"> 5.763777 </td>
   <td style="text-align:right;"> 0e+00 </td>
   <td style="text-align:right;"> 6.30e-06 </td>
   <td style="text-align:right;"> 8.359471 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SAMD9L </td>
   <td style="text-align:right;"> -1.0383851 </td>
   <td style="text-align:right;"> 0.5141782 </td>
   <td style="text-align:right;"> -5.583241 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.46e-05 </td>
   <td style="text-align:right;"> 7.519295 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> TNFAIP2 </td>
   <td style="text-align:right;"> -1.1456072 </td>
   <td style="text-align:right;"> 0.5712894 </td>
   <td style="text-align:right;"> -5.560157 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.56e-05 </td>
   <td style="text-align:right;"> 7.413159 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> CCDC88A </td>
   <td style="text-align:right;"> -1.2538645 </td>
   <td style="text-align:right;"> 0.8548085 </td>
   <td style="text-align:right;"> -5.553025 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.56e-05 </td>
   <td style="text-align:right;"> 7.380426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PLA2G7 </td>
   <td style="text-align:right;"> -0.9520283 </td>
   <td style="text-align:right;"> 0.3656429 </td>
   <td style="text-align:right;"> -5.536909 </td>
   <td style="text-align:right;"> 1e-07 </td>
   <td style="text-align:right;"> 1.62e-05 </td>
   <td style="text-align:right;"> 7.306572 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ARHGAP15 </td>
   <td style="text-align:right;"> 1.4385488 </td>
   <td style="text-align:right;"> 2.0487528 </td>
   <td style="text-align:right;"> 5.406200 </td>
   <td style="text-align:right;"> 2e-07 </td>
   <td style="text-align:right;"> 2.92e-05 </td>
   <td style="text-align:right;"> 6.712990 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> LINC00996 </td>
   <td style="text-align:right;"> -1.0889424 </td>
   <td style="text-align:right;"> 0.5261653 </td>
   <td style="text-align:right;"> -5.356665 </td>
   <td style="text-align:right;"> 3e-07 </td>
   <td style="text-align:right;"> 3.56e-05 </td>
   <td style="text-align:right;"> 6.490603 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ATP1A1 </td>
   <td style="text-align:right;"> 1.1449796 </td>
   <td style="text-align:right;"> 0.7971706 </td>
   <td style="text-align:right;"> 5.325406 </td>
   <td style="text-align:right;"> 3e-07 </td>
   <td style="text-align:right;"> 3.99e-05 </td>
   <td style="text-align:right;"> 6.350999 </td>
  </tr>
</tbody>
</table>

**Output columns:**

* logFC: log fold change (since we are working with Seurat's natural log transformed data, will be natural log fold change)
* AveExpr: Average expression across all cells in expr2
* t: logFC divided by its standard error
* P.Value: Raw p-value (based on t) from test that logFC differs from 0
* adj.P.Val: Benjamini-Hochberg false discovery rate adjusted p-value
* B: log-odds that gene is DE 

## Save files

``` r
write.csv(GenTable(GOdata, Fisher = resultFisher), file = "cluster10_GOdata.csv")
write.csv(out, file = "cluster10_Normal-Colorectal.Cancer_topTable.csv")
```

## A note on pseudobulk DE

Pseudobulk differential expression uses count data summed across all cells in each sample (typically within each cell type or cluster).  Unlike cell-level DE, pseudobulk DE *requires biological replicates* so we won't perform it on this dataset.

Once counts are summed, pseudobulk data are analyzed like bulk RNASeq data.

Pseudobulk DE may result in better false discovery rate control than cell-level DE, as shown [here](https://www.nature.com/articles/s41467-021-25960-2).

The Seurat function `AggregateExpression()` can be used to sum counts as described [here](https://satijalab.org/seurat/articles/de_vignette).

A tutorial on using limma for bulk RNASeq is available [here](https://ucdavis-bioinformatics-training.github.io/2023-June-RNA-Seq-Analysis/data_analysis/DE_Analysis_mm_with_quizzes).

## Prepare for the next section

#### Download Rmd document

``` r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2024-December-Single-Cell-RNA-Seq-Analysis/main/data_analysis/07-doublet_detection.Rmd", "07-doublet_detection.Rmd")
```

#### Session Information

``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS Ventura 13.7.1
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] org.Hs.eg.db_3.19.1  kableExtra_1.4.0     dplyr_1.1.4         
##  [4] topGO_2.56.0         SparseM_1.82         GO.db_3.19.1        
##  [7] AnnotationDbi_1.66.0 IRanges_2.38.0       S4Vectors_0.42.0    
## [10] Biobase_2.64.0       graph_1.82.0         BiocGenerics_0.50.0 
## [13] limma_3.60.2         Seurat_5.1.0         SeuratObject_5.0.2  
## [16] sp_2.1-4            
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8         
##   [4] magrittr_2.0.3          spatstat.utils_3.0-4    rmarkdown_2.27         
##   [7] zlibbioc_1.50.0         vctrs_0.6.5             ROCR_1.0-11            
##  [10] memoise_2.0.1           spatstat.explore_3.2-7  htmltools_0.5.8.1      
##  [13] sass_0.4.9              sctransform_0.4.1       parallelly_1.37.1      
##  [16] KernSmooth_2.23-22      bslib_0.7.0             htmlwidgets_1.6.4      
##  [19] ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
##  [22] zoo_1.8-12              cachem_1.1.0            igraph_2.0.3           
##  [25] mime_0.12               lifecycle_1.0.4         pkgconfig_2.0.3        
##  [28] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
##  [31] GenomeInfoDbData_1.2.12 fitdistrplus_1.1-11     future_1.33.2          
##  [34] shiny_1.8.1.1           digest_0.6.35           colorspace_2.1-0       
##  [37] patchwork_1.2.0         tensor_1.5              RSpectra_0.16-1        
##  [40] irlba_2.3.5.1           RSQLite_2.3.7           progressr_0.14.0       
##  [43] fansi_1.0.6             spatstat.sparse_3.0-3   httr_1.4.7             
##  [46] polyclip_1.10-6         abind_1.4-5             compiler_4.4.0         
##  [49] bit64_4.0.5             DBI_1.2.3               fastDummies_1.7.3      
##  [52] highr_0.11              MASS_7.3-60.2           tools_4.4.0            
##  [55] lmtest_0.9-40           httpuv_1.6.15           future.apply_1.11.2    
##  [58] goftest_1.2-3           glue_1.7.0              nlme_3.1-164           
##  [61] promises_1.3.0          grid_4.4.0              Rtsne_0.17             
##  [64] cluster_2.1.6           reshape2_1.4.4          generics_0.1.3         
##  [67] gtable_0.3.5            spatstat.data_3.0-4     tidyr_1.3.1            
##  [70] data.table_1.15.4       xml2_1.3.6              XVector_0.44.0         
##  [73] utf8_1.2.4              spatstat.geom_3.2-9     RcppAnnoy_0.0.22       
##  [76] ggrepel_0.9.5           RANN_2.6.1              pillar_1.9.0           
##  [79] stringr_1.5.1           spam_2.10-0             RcppHNSW_0.6.0         
##  [82] later_1.3.2             splines_4.4.0           lattice_0.22-6         
##  [85] bit_4.0.5               survival_3.5-8          deldir_2.0-4           
##  [88] tidyselect_1.2.1        Biostrings_2.72.0       miniUI_0.1.1.1         
##  [91] pbapply_1.7-2           knitr_1.47              gridExtra_2.3          
##  [94] svglite_2.1.3           scattermore_1.2         xfun_0.44              
##  [97] statmod_1.5.0           matrixStats_1.3.0       UCSC.utils_1.0.0       
## [100] stringi_1.8.4           lazyeval_0.2.2          yaml_2.3.8             
## [103] evaluate_0.23           codetools_0.2-20        tibble_3.2.1           
## [106] cli_3.6.2               uwot_0.2.2              systemfonts_1.1.0      
## [109] xtable_1.8-4            reticulate_1.39.0       munsell_0.5.1          
## [112] jquerylib_0.1.4         GenomeInfoDb_1.40.1     Rcpp_1.0.12            
## [115] globals_0.16.3          spatstat.random_3.2-3   png_0.1-8              
## [118] parallel_4.4.0          blob_1.2.4              ggplot2_3.5.1          
## [121] dotCall64_1.1-1         listenv_0.9.1           viridisLite_0.4.2      
## [124] scales_1.3.0            ggridges_0.5.6          crayon_1.5.2           
## [127] leiden_0.4.3.1          purrr_1.0.2             rlang_1.1.3            
## [130] KEGGREST_1.44.0         cowplot_1.1.3
```
