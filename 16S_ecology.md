## Packages and dependencies

```r
###Install phyloseq###

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")


###Install ape###

install.packages("ape")


###Install qiime2R###

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")


# Installing from GitHub requires the remotes package
install.packages("remotes")
# Windows users will also need to have RTools installed! http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/

# To install the latest version:
remotes::install_github("david-barnett/microViz")

# To install a specific "release" version of this package, e.g. an old version 
remotes::install_github("david-barnett/microViz@0.12.0") 

###Install Microbiome###

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiome")


###Install vegan###
install.packages("vegan")
```

### Import packages
```r
library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(ape)
library(microbiome)
library(dplyr)
```

