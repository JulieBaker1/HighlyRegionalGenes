Introduction
============

Select genes that have biological significance is useful for downstream
analysis.Here we bring up a feature selection method called
`HighlyRegionalGenes`which identify genes with regional patterns.The
method mainly contains three steps:

-   construct SNN
-   calculate score and rank genes
-   select the regional distributed genes

We use iteration to make the SNN more precise.Below is the workflow.

![workflow.png](https://github.com/JulieBaker1/HighlyRegionalGenes/blob/master/images/61e0650f339ebd15ded6f3ea569ef67.png)

System Requirements
-------------------

-   R version: &gt;= 3.6.1
-   Seurat version: &gt;3.1.4

If you would like to run the codes that compare HRG with other feature
selection methods, M3drop and parallel packages need to be installes. It
is necessary to notify you that those codes will cost much time.

Installation
------------

the package can be installed directly from the github.

    install.packages("devtools")
    library(devtools)
    install_github("JulieBaker1/HighlyRegionalGenes")

usage
-----

The vignette of HighlyRegionalGenes can be found in the project
[wiki](https://github.com/JulieBaker1/HighlyRegionalGenes/wiki).
