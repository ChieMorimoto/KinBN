# KinBN

KinBN is a free software (GNU General Public License v3.0) for kinship analysis based on the Bayesian network . It can be applied to short tandem repeat (STR) markers commonly used in forensic genetics to calculate the likelihood ratio (LR) considering the linkage between loci, mutation, and drop-out. The software is a graphical user interface written in R language. The software has been validated under various conditions.

## Changes in v2.0.0
Added functionalities:
* The software can be used as the R-package KinBN.
* The software can calculate the LR value while accounting for allele drop-out.
* Project data can be saved and loaded as required.

Minor changes:
* Layout of the software has changed.
* The user can set the calculation conditions, such as the estimation method of allele frequency, in the Setting tab.

## Getting started
1. Install the R language  (v.4.3 or v.4.2) available on the R Development Core Team website (http://www.R project.org).

2. Enter the code below the R console to install KinBN and other necessary packages.
```r
install.packages('https://github.com/ChieMorimoto/KinBN/releases/download/v2.0.0/KinBN_2.0.0.zip',repos=NULL,type='win.binary')
install.packages(c(“tcltk2”, “bnlearn”, “gRbase”, “gRain”, “tkrplot”, “kinship2”))
```

3. Enter the codes below to start the graphical user interface (GUI).
```r
library(KinBN)
KinBN()
```
