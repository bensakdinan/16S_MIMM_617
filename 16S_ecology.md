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

## Section 5: Normalizing depth
- Why is it important to normalize read depth?

Question 5A) What would be a good sequencing depth to subsample at? 
 - 37605 sequencing. I reordered my fixed_physeq data by number of sequences, 
   from lowest to highest. Therefore, the number of seqs at the first index
   is our minimum sequencing depth
Question 5B) In general, what would the lowest acceptable sampling depth where we could have adequate detection of all ASVs?
 - Looking at my rarefaction curve, it seems that the sequencing depth begins to 
   saturate at around 15,000 reads
```{r}
sequencing_depth <- sample_sums(fixed_physeq)
sequencing_depth_df <- as.data.frame(sequencing_depth)

# Reorder this dataset from lowest to highest
sequencing_depth_ordered  <- sequencing_depth_df[order(sequencing_depth_df$sequencing_depth), ]
print(sequencing_depth_ordered)

# Rarefaction
matrix_phyloseq <- as(t(otu_table(physeq)), "matrix")

## Calculate row sums (sequencing depth of each sample)

row_sums <- rowSums(matrix_phyloseq)

# Order rows by row sums and take the 50 rows with the lowest sums
subset_matrix_phyloseq <- matrix_phyloseq[order(row_sums)[1:20], ]



#now we have a matrix that we can put into the rarecurve function! 
#step = 1000  means that the function calculates ASV richness every 1000 sequences until the total number of species for that sample is reached

rarecurve(subset_matrix_phyloseq, step = 1000, label = FALSE)

# Rarefy by the minimum sequencing depth
rarefied_phyloseq <- rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)))

```

## Section 7: Between Sample Diversity

Question 7A) Imagine two points on the PCoA plot above. In case one, these two points are distributed vertically and in case 2, the points are distributed horizontally (equidistant to case 1). In which case is there greater distance between points?
 - Horizontal PCoA. The horizontal axis explains 31.3% of variance 
Question 7B) Which metadata has the largest effect on beta-diversity based on Bray-Curtis ordinations?
 - vendor
Question 7C) Qualitatively, which distance metric yields higher differences between the groups. Why would this be the case?
 - Weighted UniFrac
Question 7D) In which vendor is the diet effects most pronounced? Does resistant start or inulin generally have a larger effect compared to the control cellulose diet?
 - Beijing, Inulin
Question 7E) In the shanghai vendor, try to determine if the effects of inulin are static, or change over time (hint: you’ll need to play around with filtering)
 -  Changes over time
Question 7F) Which metadata groups have a significant effect on beta-diversity? Which has the largest effect?
 - Diet (Group 2) has a significant effect on beta diversity.
Question 7G) In which vendor is the effect of diet the largest? Does this match your observations with the corresponding PCoA plot?
 - Beijing

```{r}
# This returns a distance 
bray_distance <- distance(rarefied_phyloseq, method="bray")


# Using this distance matrix, we can build a ordination plot
ordination_pcoa_bray <- ordinate(rarefied_phyloseq, "PCoA", "bray")
plot_ordination <- plot_ordination(rarefied_phyloseq, ordination_pcoa_bray)
plot_ordination(rarefied_phyloseq, ordination_pcoa_bray, color = "" )
print(plot_ordination)

#subset for each vendor 
# can use the != operator to select samples that ARE NOT day 0 (baseline)
bejing_physeq_rarefied <- subset_samples(rarefied_phyloseq, vendor == "Beijing" & Day != 0 )
shanghai_physeq_rarefied <- subset_samples(rarefied_phyloseq, vendor == "Shanghai" & Day!= 0)
guangdong_physeq_rarefied <- subset_samples(rarefied_phyloseq, vendor == "Guangdong" & Day!= 0)


#create ordinations for each vendor 
ordination_pcoa_bray_beijing <- ordinate(bejing_physeq_rarefied, "PCoA", "bray")
ordination_pcoa_bray_shanghai <- ordinate(shanghai_physeq_rarefied, "PCoA", "bray")
ordination_pcoa_bray_guangdong <- ordinate(guangdong_physeq_rarefied, "PCoA", "bray")



#plot 
plot_ordination(bejing_physeq_rarefied, ordination_pcoa_bray_beijing, color = "Group")
plot_ordination(shanghai_physeq_rarefied, ordination_pcoa_bray_shanghai, color = "Group")
plot_ordination(guangdong_physeq_rarefied, ordination_pcoa_bray_guangdong, color = "Group")


# Calculate bray curtis distance matrix
bray_distance <- distance(rarefied_phyloseq, method = "bray")

# make a data frame from the sample_data sow we can feed it into 
metadata_df <- data.frame(sample_data(rarefied_phyloseq))

# Adonis PERMANOVA test
#Requires,the bray_distance matrix, the group of interest to compare and the metadata dataframe 

adonis2(bray_distance ~ Group, data = metadata_df)

# Calculate bray curtis distance matrix
bray_distance <- distance(rarefied_phyloseq, method = "bray")

# make a data frame from the sample_data sow we can feed it into 
metadata_df <- data.frame(sample_data(rarefied_phyloseq))

# Adonis PERMANOVA test
#Requires,the bray_distance matrix, the group of interest to compare and the metadata dataframe 

adonis2(bray_distance ~ Group, data = metadata_df)

bejing_physeq_rarefied <- subset_samples(fixed_physeq, vendor == "Beijing" & Day != 0 )
shanghai_physeq_rarefied <- subset_samples(fixed_physeq, vendor == "Shanghai" & Day!= 0)
guangdong_physeq_rarefied <- subset_samples(fixed_physeq, vendor == "Guangdong" & Day!= 0)


# Calculate bray curtis distance matrix
bray_distance_beijing <- distance(bejing_physeq_rarefied, method = "bray")
bray_distance_shanghai <- distance(shanghai_physeq_rarefied, method = "bray")
bray_distance_guangdong <- distance(guangdong_physeq_rarefied, method = "bray")

# make a data frame from the sample_data sow we can feed it into 
metadata_df_b <- data.frame(sample_data(bejing_physeq_rarefied))
metadata_df_s <- data.frame(sample_data(shanghai_physeq_rarefied))
metadata_df_g <- data.frame(sample_data(guangdong_physeq_rarefied))


adonis(bray_distance_beijing ~ Group, data = metadata_df_b)
adonis(bray_distance_shanghai ~ Group, data = metadata_df_s)
adonis(bray_distance_guangdong ~ Group, data = metadata_df_g)

```
