
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="vignettes/figures/home.png" width="100%" style="margin-top: 20px;" />

<!-- badges: start -->

[![R
\>3.5](https://img.shields.io/badge/R-%3E3.5-success.svg)](https://www.r-project.org/)
<a href='#devtools'>![installed with
devtools](https://img.shields.io/badge/installed%20with-devtools-blueviolet.svg)</a>
[![GitHub
codesize](https://img.shields.io/github/languages/code-size/idrblab/ANPELA)](https://github.com/idrblab/ANPELA/releases)
[![WebServer](https://img.shields.io/badge/Web_Server-ANPELA-blue)](https://idrblab.org/anpela2025/)
<!-- badges: end -->

# How to Use `ANPELA`?

## Contents

- [Overview](#overview)
- [Citations](#citations)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Example Data](#example-data)
- [Usage and Examples](#usage-and-examples)

## Overview

Single-cell proteomics (SCP) has emerged as a powerful technique that
significantly advances our understanding of complex biological systems
with new level of granularity. Because of the extreme difficulty in
processing SCP data, ANPELA was developed for identifying the optimal
workflow based on a well-designed assessment strategy.<br><br>
**ANPELA** is a **user-centric** and **application-oriented** tool which
is capable of navigating the data processing for SCP. **ANPELA 3.0** has
significantly improved its practicality, focusing primarily on the
following points: **1**. <u>multi-scenarios deployment</u> (versatile
choices meet diverse user needs); **2**. <u>data security</u> (local
execution ensures data confidentiality); **3**. <u>open source</u>
(modular codes facilitate readers’ free editing); and **4**.
<u>user-friendly interface</u> (the visual interface enhances user
application).<br><br> The local software and webserver of ANPELA are
available at <a href="https://idrblab.org/anpela2025">
https://idrblab.org/anpela2025</a>.

## Citations

You can cite the `ANPELA` publication as follows:

> **ANPELA: Significantly Enhanced Quantification Tool for
> Cytometry-Based Single-Cell Proteomics.**
>
> Ying Zhang, Huaicheng Sun, Xichen Lian, Jing Tang, Feng Zhu\*.
>
> *Advanced Science*. 2023 May;10(15):e2207061. <br>doi:
> [10.1002/advs.202207061](https://onlinelibrary.wiley.com/doi/10.1002/advs.202207061);
> PubMed ID: 36950745.

## Repo Contents

- R: R package code.
- man: package manual for help in R session.
- vignettes: R vignettes for R session html help pages.

## System Requirements

### Hardware

- **Minimum**:<br> OS: Windows® 32/64-bit, macOS®, Linux<br> Processor:
  2.0 GHz<br> Memory: 8 GB RAM<br> Storage: 2 GB available hard-disk
  space<br><br>
- **Recommended**: <br> OS: Windows® 32/64-bit, macOS®, Linux<br>
  Processor: 4.0 GHz<br> Memory: 32 GB RAM<br> Storage: 8 GB available
  hard-disk space<br><br>
- An internet connection is necessary to download the ANPELA desktop
  software or R packages (including the ANPELA R package and other
  prerequisite R packages).

### Software

- *R v4.0 or later*. The installation file of *R* language, named
  ‘*R*-x.y.z.tar.gz’ (where x, y, and z represent the version numbers of
  the software release), retrieved from CRAN R-project website
  (<https://cran.r-project.org>), which is compatible with user’s
  operating system (the latest version of this protocol was developed
  and tested using ‘*R* 4.4.1’).

- *RStudio*. The RStudio installation file ‘RStudio-x.y.z.zip’ (where x,
  y, and z represent the version numbers of the software release) from
  RStudio website (<https://www.rstudio.com/>), which is compatible with
  user’s operating system.

- *RTools*. The RTools installation file ‘rtoolsx-y-z.zip’ (where x, y,
  and z denote the specific version numbers of the software release) can
  be obtained from the official RTools website
  (<https://cran.r-project.org/bin/windows/Rtools> /);

- *Prerequisite R packages*. This protocol necessitates a range of *R*
  packages (as described in [Installation Guide](#installation-guide)),
  with certain packages obtained from CRAN
  (<https://cran.r-project.org>), others installed from Bioconductor
  (<https://bioconductor.org/>), and the remaining accessible from
  GitHub (<https://github.com/>).

## Installation Guide

This section provides a comprehensive guide to the installation and
configuration of `ANPELA` *R* package. It outlines the necessary tools
required to ensure a seamless and accurate installation of the protocol.

- **Tool 1.** ***R*** **language, RStudio and RTools**

Install the R language, RStudio and RTools using their installation
files, which are compatible with user’s operating system.

- **Tool 2. The required** ***R*** **Packages**

Install a variety of *R* packages imported in this protocol.

Installed from ***CRAN*** (can also from other repositories): dplyr,
doParallel, rstan, Rtsne, pastecs, cowplot, ggpubr, gridExtra,
MLmetrics, fossil, clusterCrit, VennDiagram, stringr, bbmle, mc2d,
parallel, doSNOW, foreach, igraph, mclust, pheatmap, magrittr, withr,
tree, sampling and ellipse.

Installation commands:

``` r
CRAN_packages <- c("dplyr", "doParallel", "rstan", "Rtsne", "pastecs", "cowplot",
                   "ggpubr", "gridExtra", "MLmetrics", "fossil", "clusterCrit",
                   "VennDiagram","stringr", "bbmle", "mc2d", "doSNOW", "foreach",
                   "igraph","mclust", "pheatmap", "magrittr", "withr", "tree",
                   "sampling", "ellipse")
install.packages(CRAN_packages, dependencies = TRUE)
```

Installed from ***Bioconductor*** (can also from other repositories):
flowCore, limma, SCORPIUS, slingshot, destiny, flowStats, flowAI,
flowCut, flowClean, spillR, PeacoQC, systemPipeR, flowTrans and flowVS.

Installation commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

Bioconductor_packages <- c("flowCore", "limma", "SCORPIUS", "slingshot", "destiny", 
                           "flowStats", "flowAI", "flowCut", "flowClean", "spillR",
                           "PeacoQC", "systemPipeR", "flowTrans", "flowVS")

BiocManager::install(Bioconductor_packages, ask = FALSE)
```

Installed from ***GitHub*** (can also from other repositories):
Rtsne.multicore, Rphenograph, CytoSpill and autospill.

Installation commands:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("RGLab/Rtsne.multicore")
devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("qmiao19/CytoSpill")
devtools::install_github("carlosproca/autospill")
```

- **Tool 3. ANPELA** ***R*** **Package**

Install the `ANPELA` package by running the following command in
RStudio:

``` r
devtools::install_github("idrblab/ANPELA", build_vignettes = TRUE)
```

During the installation of `ANPELA`, the appearance of the error message
“ERROR: dependency ‘package_name’ is not available for package
‘`ANPELA`’” indicates that the required imported package,
‘package_name’, has not been successfully installed. Users should refer
to the detailed reinstallation instructions described in **Tool 2** to
resolve this issue and ensure the proper installation of the missing
package before proceeding with the `ANPELA` installation.

## Example Data

You can use the example data provided in `ANPELA` to try it out.

**Example data acquired from flow cytometry for CSI analysis**
([Download_287
MB](https://idrblab.cn/anpela2025/example_data/FC_CSI.zip)): A
single-cell proteomic dataset involving 23 markers of fresh thymus
tissue obtained from nine patients undergoing elective thymectomy,
including three myasthenia gravis (MG) patients and six healthy controls
(non-MG).

**Example data acquired from mass cytometry for CSI analysis**
([Download_174
MB](https://idrblab.cn/anpela2025/example_data/MC_CSI.zip)): A
single-cell proteomic dataset consisting 35 surface markers of two
antigen specific T cell lines generated by naïve CD4+ T cells stimulated
with tolerogenic dendritic cells (tolDCs) or mature inflammatory myeloid
dendritic cells (mDCs) pulsed with proinsulin peptide.

**Example data acquired from flow cytometry for PTI analysis**
([Download_3.70
MB](https://idrblab.cn/anpela2025/example_data/FC_PTI.zip)): A
single-cell proteomic temporal dataset using 10 cell surface markers in
vitro hematopoietic differentiation system from human embryonic stem
cells (HUES9) at 6 sequential timepoints (0, 2, 4, 6, 8, 10 day) to
capture cells at different developmental stages.

**Example data acquired from mass cytometry for PTI analysis**
([Download_882
KB](https://idrblab.cn/anpela2025/example_data/MC_PTI.zip)): A gated
single-cell proteomic temporal dataset encoding the activation dynamics
of 14 CD8+ cell intracellular markers perturbed by tetradecanoylphorbol
acetate (PMA)/ionomycin at 8 sequential timepoints (0, 1, 5, 15, 30, 60,
120 and 240 min).

## Usage and Examples

For the usage and examples of `ANPELA`, users can refer to the vignette
“How to Use `ANPELA`” built in the package.

``` r
vignette("ANPELA")
```

For details of each method provided by `ANPELA`, users can refer to the
vignette “Methods_Introduction” built in the package.

``` r
vignette("Methods_Introduction")
```

For the functions provided by `ANPELA`, users can refer to other
vignettes.

``` r
vignette("Getmarker")
vignette("Process")
vignette("FCprocess")
vignette("MCprocess")
vignette("Assess")
vignette("CSIassess")
vignette("PTIassess")
vignette("Ranking")
vignette("FC_CSI")
vignette("MC_CSI")
vignette("FC_PTI")
vignette("MC_PTI")
```
