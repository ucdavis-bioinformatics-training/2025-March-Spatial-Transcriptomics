---
title: "Spatial Transcriptomics Part 7: Cell Segmentation"
author: "UCD Bioinformatics Core"
date: "`r Sys.Date()`"
output:
    html_document:
      keep_md: TRUE
      toc: TRUE
---

# Spatial Transcriptomics Part 7: Cell segmentation

Cell segmentation is an essential task that allows us to study cellular biology at a normal and pathological conditions, since cell is the smallest unit of life. It allows analysis such as cell count, cell type, cell shape, etc. and in turn to study how these festures change over time and in response to different pathological conditions.

Generating accurate cell segmentation is more difficult than one may think. There are many factors that pose challenge to the task. For example, cell shapes differ significantly that no mathematical model has been successful in representing all possibilities. Cell sizes have a wide range as well. In addition, the quality of imaging at the cell boundary is more likely to have weak gradients when cells are in close proximity that makes it difficult to separate them. The methods developed for this task all have their own limitations that makes it difficult to know which method one should choose.

Even though all image-based spatial transcriptomics platform try to develop the best approach for their assays, there is still cases where some portion of the cells are not segmented accurately under close inspection. As an example, 10X clearly states that their most advanced segmentation algorithm, Xenium Multimodal Cell Segmentation algorithm, is not performing as well for large cells with a diameter greater than 100 microns, such as muscle, adipocytes and signet ring cells. In addition, when nuclei stain is weak in some cell types, such as dorsal root ganglion, it doesn't perform optimally either.

There are many methods that have been developed to tackle this challenge. There are two main sources of information that are used in cell segmentation. One is imaging information and the other is the spatial transcriptomics information. [Cellpose](https://www.nature.com/articles/s41592-025-02595-5), [Baysor](https://www.nature.com/articles/s41587-021-01044-w) and [StarDist](https://arxiv.org/abs/1806.03535) use primarily image information. [ST-CellSeg](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012254), [ClusterMap](https://pubmed.ncbi.nlm.nih.gov/34625546/) leverage the spatial transcriptomics to segment cells without imaging information. [Baysor](https://www.nature.com/articles/s41587-021-01044-w) can take into account both sources of information and makes it the most versatile tool.

First, there are two softwares that I would like to introduce to you. One is [QuPath](https://qupath.github.io/). It is a very good tool for image exploration and it has many functions that can aide in annotating images. There is a [great video](https://www.youtube.com/watch?v=dZYjIOI76WM) from an export. The other is [napari](https://napari.org/stable/index.html), which has a [github repo](https://github.com/napari/napari) for installation. Both are great tools for spatial transcriptomics.

Second, Xenium Explorer is also a great place to start by visualizing your tissues, cells, transcripts, especially when the multimodal imaging is used for an assay. Let take a look at [one example](https://www.10xgenomics.com/datasets/preview-data-xenium-prime-gene-expression).

---
