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

Import articacts to initialize a phyloseq object

```{r}
physeq <- qza_to_phyloseq(
  features="/Users/bensakdinan/Desktop/MIMM_617/00_data/feature-sample-filtered-table.qza",
  tree="/Users/bensakdinan/Desktop/MIMM_617/00_data/rooted-tree.qza",
  taxonomy="/Users/bensakdinan/Desktop/MIMM_617/00_data/filtered-taxonomy.qza",
  metadata="/Users/bensakdinan/Desktop/MIMM_617/00_data/sample-metadata.txt"
)
```
## Section 1
Now that we have initialized a phyloseq object, we can explore the data. 

Q1A) How many taxa (ASVs) and samples do we have in the dataset?
 - 1840 taxa in 152 samples
```{r}
physeq
```

Now let's look at the metadata in the physeq object

Q1B) How many time points are in the experiment
  - 4 (0, 5, 13, 31)
Q1C) What are the three diet treatments the mice are given?
  - Inulin, resistant starch, and cellulose
Q1C) Where are the mouse vendors from?
 - Guangdong, Beijing, Shanghai
```{r}
sample_data(physeq)
metadata <- as.data.frame(sample_data(physeq))
metadata
uniq_metadata <- unique(metadata$Day)
uniq_metadata
```

## Section 2: Exploring Microbiome Composition

Let's make a simple bar plot to get a general visualization of our data
```{r}
# If a row in physeq is unknown, it will be assigned the name "uncultured" 
fixed_physeq <- tax_fix(physeq)
fixed_physeq <- fixed_physeq %>% 
  tax_fix(unknowns = c("uncultured"))

# Assigning the taxa level as "Phylum"
comp_barplot(fixed_physeq, tax_level = "Phylum")

# comp_barplot by default displays the top 8 most abundant taxa, but this can be increased to 12
comp_barplot(fixed_physeq, tax_level = "Phylum", n_taxa = 12)

# By family 
comp_barplot(fixed_physeq, tax_level = "Family", n_taxa = 12)

# By genus 
comp_barplot(fixed_physeq, tax_level = "Genus", n_taxa = 12)

# By species
fixed_physeq %>% tax_fix(unknowns = c("uncultured Genus"))
comp_barplot(fixed_physeq, tax_level = "Species")
```

Question 2A) Which bacterial phyla tend to be the most abundant in the dataset?
 - Bacteroidota
Question 2B) What about at the family and genus levels?
 - Some samples are abundant in Bacteroidaceae and others in Muribaculaceae
Question 2C) Try doing this at species level, anything strange to report of the species names?
 - Species-levle resolution is not available for this dataset 
 
## Temporarily skip section 3

## Section 4: 

```{r}
beijing_physeq <- subset_samples(fixed_physeq, vendor == "Beijing")
shanghai_physeq <- subset_samples(fixed_physeq, vendor == "Shanghai")
guangdong_physeq <- subset_samples(fixed_physeq, vendor == "Guangdong")

comp_barplot(beijing_physeq, tax_level = "Family")
comp_barplot(shanghai_physeq, tax_level = "Family")
comp_barplot(guangdong_physeq, tax_level = "Family")

# comp_barplot(beijing_physeq, tax_level = "Family", group_by = "Group", facet_by = "Day", n_taxa = 12)
# comp_barplot(shanghai_physeq, tax_level = "Family", group_by = "Group", facet_by = "Day", n_taxa = 12)
# comp_barplot(guangdong_physeq, tax_level = "Family", group_by = "Group", facet_by = "Day", n_taxa = 12)
```
Question 4A) Look through different taxonomic levels - What differences in the relative abundance do you observe based on mouse vendor?
 - At the family level, there are differences in the abudances of different families between vendor location
Question 4B) In the published paper they describe differences specifically in Muribaculaceae (polysaccharide-degrading). Describe the differences in this family between vendors
 - The relative abundance of Muribaculaceae is comparable between Beijing and Shanghai, but it seems that Muribaculaceae is significantly increased in Guangdong.
Question 4C) In general do you see variability between replicates?
 - Yes, there is some small level of variance between replicates, but this variance initally appears to be small
Question 4D) In these studies, cellulose is the control diet (all mice are kept on cellulose before day 0, and then either switch to inulin, resistant starch or are kept on cellulose). Is microbiome composition stable over time in this control? Is this different across taxonomic levels?
 - I can see significant changes in the composition in resistant starch, in inulin, and at day 31 in cellulose
Question 4E) What effects do inulin and resistant starch have on the microbiome? Are these constant over time? Do they differ based on mouse vendor?
```{r}
comp_barplot(shanghai_physeq, tax_level = "Family", group_by = "Group", facet_by = "Day", n_taxa = 12)
comp_barplot(guangdong_physeq, tax_level = "Family", group_by = "Group", facet_by = "Day", n_taxa = 12)
```

