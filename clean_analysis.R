# Clean analysis for Loutet et al. 2024
# Last updated: 20 November 2024

# # Set the path correctly to where we want to work ---------------------------

setwd('/Users/samloutet/Desktop/HudsonBay2024')
# All data should be here


# Load the packages necessary -----------------------------------------------

#if (!requireNamespace("BiocManager", quietly = TRUE)) 
#install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome"))
library(ggpubr) # plotting
library(tibble)
library(picante) # alpha diversity
library(microbiome) # alpha diversity
library(phyloseq) # functions for manipulating amplicon data
library(tidyverse) # general data manipulation and plotting
library (readr) # importing sequence data into r
library(seqinr) #lets you work with fasta files
library(ape) # manipulating data sets
library(vegan) # statistical tools
library(RColorBrewer) #useful for plots
#install.packages("viridis") #useful for plots
library(viridis) #useful for plots
library(ggplot2) #useful for plots
library(dplyr)
#install.packages("devtools")
#devtools::install_github("jfq3/QsRutils", build_vignettes = TRUE) #for Good's coverage
library(QsRutils)
library(cowplot)
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")
library(metagMisc)


# Import all the files from QIIME -------------------------------------------

# Sequences
count_table <- read_tsv(file="Export/table/table.tsv", skip = 1)

# And specify that the first column of data are rownames
count_table <- column_to_rownames(count_table, var = colnames(count_table)[1])

# Rarecurve expects the dataset as a dataframe so we need to use as.data.frame
count_table_df <- as.data.frame(count_table)

# Metadata
sample_info_tab<-read_tsv('sample-metadata.tsv')

#remove first row with character designations
sample_info_tab<-sample_info_tab[-1,]

# Convert metadata rows which are numeric to numeric
sample_info_tab$'lat-decimal' <- as.numeric(sample_info_tab$'lat-decimal')
sample_info_tab$'long-decimal' <- as.numeric(sample_info_tab$'long-decimal')
sample_info_tab$'temperature-C' <- as.numeric(sample_info_tab$'temperature-C')
sample_info_tab$'salinity-PSU' <- as.numeric(sample_info_tab$'salinity-PSU')
sample_info_tab$'coastline-km' <- as.numeric(sample_info_tab$'coastline-km')
sample_info_tab$'ice-distance-km' <- as.numeric(sample_info_tab$'ice-distance-km')
sample_info_tab$'days-since-ice' <- as.numeric(sample_info_tab$'days-since-ice')

# Order the metadata rows similarly to in count_table
sample_info_tab$'sample-id' <- factor(sample_info_tab$'sample-id', levels = c("NaMah-V4-16", "NaMah-V4-18", "NaMah-V4-19", "NaMah-V4-24", "NaMah-V4-26", "NaMah-V4-28", "NaMah-V4-29", "NaMah-V4-3", "NaMah-V4-33", "NaMah-V4-4", "NaMah-V4-5", "NaMah-V4-6", "NaMah-V4-7", "NaMah-V4-8", "NaMah-V4-9", "NaMah-V4rep-27", "NaMah-V4rep-30", "NaMah-V4rep-31", "NaMah-V4rep-32"))

#Taxonomy
taxonomy <- read_tsv(file="Export/taxonomy/taxonomy.tsv")
taxonomy_df<-data.frame(taxonomy)

#Tree
tree = read_tree("Export/exported-tree/tree.nwk")


# Rarefaction ---------------------------------------------------------------

# Use rarecurve, from the Vegan package.

# Plot the rarefaction curves
# Add a veritical line to the plot indicating the fewest # of sequences from 
# any sample
# Export for use in ggplot
otu_curve = rarecurve(t(count_table_df), step=100, cex=0.5, lwd=2, ylab="ASVs", xlab="Sequence Reads", label = TRUE, tidy = TRUE)
otu_curve$Site <- factor(otu_curve$Site) 

# Make it easier to see the real site names
sites_key <- data.frame(sample_info_tab$`sample-id`, sample_info_tab$`site-name`, row.names = NULL)
names(sites_key)[1] <- "Site"
names(sites_key)[2] <- "Name"
otu_curve <- inner_join(otu_curve, sites_key, by = "Site")

# Make colour bars
darks1 <- brewer.pal(6, "Dark2")
darks2 <- brewer.pal(7, "Set2")
darks3 <- brewer.pal(6, "Accent")
darks4 <- brewer.pal(12, "Paired")
darks5 <- brewer.pal(8, "Set1")
darks6 <- brewer.pal(11, "Set3")
darks <- c(darks4, darks2, darks3, darks1, darks5, darks6, darks4, darks2, darks3, darks1, darks5, darks6, darks4)

# Ensure the levels of the 'Name' variable match the order in the 'darks' palette
otu_curve$Name <- factor(otu_curve$Name, levels = unique(otu_curve$Name))  # Set factor levels for Name

# Extract only the last point for each site
label_points <- otu_curve %>%
  group_by(Site) %>%
  slice_max(Sample) %>% # Select the row with the maximum Sample value for each Site
  ungroup()

# Identify overlapping labels and nudge them manually
label_points <- label_points %>%
  mutate(
    # Adjust positions for specific labels to avoid overlap
    nudge_x_17 = ifelse(Name == "17", 100, 0),  # Nudge right for label 17
    nudge_y_11 = ifelse(Name == "11", -10, 0)    # Nudge down for label 11
  )

# Plot
rare_plot <- ggplot(otu_curve, aes(x=Sample, y=Species, color=Name)) +
  geom_line(linewidth=0.75) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.75, fill = NA), legend.position = "none", legend.key = element_blank(), legend.justification = "top") +
  theme(text = element_text(size = 14)) +
  scale_color_manual(values=darks, name = "site") +
  geom_label(
    data = label_points,
    aes(label = Name, fill = Name, color=Name), # Add fill matching line color
    color = "black", # Label text color
    fill = darks[1:19], 
    size = 3,
    nudge_x = label_points$nudge_x_17,  # Apply nudges
    nudge_y = label_points$nudge_y_11   # Apply nudges
  ) +
  xlab("number of sequences") + ylab("number of ASVs")

print(rare_plot)

# identify the number of reads in the sample with the fewest reads.
min(rowSums(t(count_table_df)))

#create phyloseq object for OTU table
count_table_phyloseq_df = otu_table(count_table_df, taxa_are_rows = TRUE)

#set rarefaction depth for samples equal the number of reads in the sample with
# the fewest reads
sums = sample_sums(count_table_phyloseq_df)
raredepth = min(sums[sums > 0])

# rarefy without replacement
# set seed and record this for your records so others can reproduce
# the seed I'm starting with is 1
count_table_phyloseq_df=rarefy_even_depth(count_table_phyloseq_df, rngseed = 1,sample.size = raredepth, replace = F)


# Move Each Taxonomic Level to Its Own Column -------------------------------

head(taxonomy_df)

taxonomy_organized <- taxonomy_df %>%
  #get rid of the letters in front of each taxon name
  mutate(taxonomy=str_replace_all(string=Taxon, pattern=c("d__|p__|c__|o__|f__|g__|s__" ), replacement="")) %>%
  #separate into columns
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species"), sep=";") %>%
  select(-Taxon, -Confidence) %>%
  column_to_rownames(var = 'Feature.ID')

head(taxonomy_organized)


# Create Phyloseq Objects ---------------------------------------------------

# Ordinations are done using the phyloseq package, which first requires making
# phyloseq objects of data
TAX<-tax_table(as.matrix(taxonomy_organized))
sample_info_tab_df<-data.frame(sample_info_tab, row.names=sample_info_tab$`sample-id`)

# Order the metadata rows similarly to in count_table
sample_info_tab_df <- sample_info_tab_df %>% arrange(factor(sample.id, levels = c("NaMah-V4-16", "NaMah-V4-18", "NaMah-V4-19", "NaMah-V4-24", "NaMah-V4-26", "NaMah-V4-28", "NaMah-V4-29", "NaMah-V4-3", "NaMah-V4-33", "NaMah-V4-4", "NaMah-V4-5", "NaMah-V4-6", "NaMah-V4-7", "NaMah-V4-8", "NaMah-V4-9", "NaMah-V4rep-27", "NaMah-V4rep-30", "NaMah-V4rep-31", "NaMah-V4rep-32")))

META<-sample_data(sample_info_tab_df)
ASV<-as.matrix(count_table_phyloseq_df)

#Make one phyloseq object, which contains all 4 objects
ps <- phyloseq(ASV, META, TAX, tree)


# Clean up phyloseq object taxonomy -----------------------------------------

#Filter out those ambigious Kingdom annotations
ps <- subset_taxa(ps, Kingdom %in% c("Bacteria", "Archaea"))
table(tax_table(ps)[, "Kingdom"], exclude = NULL)
tax.clean <- data.frame(tax_table(ps))

# Get a count of the number of taxa in each phylum.
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Rename unknown or NA taxa to lowest known rank
# rename NA to last known taxa

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste( tax.clean[i,1], "(k)", sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste(tax.clean[i,2], "(p)", sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste(tax.clean[i,3], "(c)", sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste(tax.clean[i,4], "(o)", sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste(tax.clean[i,5], "(f)", sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], "(g)", sep = "_")
  }
}
tax_table(ps) <- as.matrix(tax.clean)
#rename uncultured to last known taxa
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == " uncultured"){
    kingdom <- paste( tax.clean[i,1], "(k)", sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == " uncultured"){
    phylum <- paste(tax.clean[i,2], "(p)", sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == " uncultured"){
    class <- paste(tax.clean[i,3], "(c)", sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == " uncultured"){
    order <- paste(tax.clean[i,4], "(o)", sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == " uncultured"){
    family <- paste(tax.clean[i,5], "(f)", sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == " uncultured"){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], "(g)", sep = "_")
  }
}
tax_table(ps) <- as.matrix(tax.clean)

#rename unknown to last known taxa
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == " Unknown_Phylum"){
    kingdom <- paste( tax.clean[i,1], "(k)", sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == " Unknown_Class"){
    phylum <- paste(tax.clean[i,2], "(p)", sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == " Unknown_Order"){
    class <- paste(tax.clean[i,3], "(c)", sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,4] == " Gammaproteobacteria_Incertae_Sedis"){
    class <- paste(tax.clean[i,3], "(c)", sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == " Unknown_Family"){
    order <- paste(tax.clean[i,4], "(o)", sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == " Unknown_Genus"){
    family <- paste(tax.clean[i,5], "(f)", sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == " Unknown_Species"){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], "(g)", sep = "_")
  }
}
tax_table(ps) <- as.matrix(tax.clean)


# Reroot tree ---------------------------------------------------------------

# first define function from link above to find furthest outgroup
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  
  # tablify parts of tree that we need.
  treeDT <-
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>%
    cbind(data.table(id = tree.unrooted$tip.label))
  
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup) }

# then run on my phyloseq tree
my.tree <- phy_tree(ps)
out.group <- pick_new_outgroup(my.tree)
out.group

# Then use this outgroup to root the tree
new.tree <- ape::root(my.tree, outgroup=out.group, resolve.root=TRUE)

# and convert to dichotomy tree
new.tree2 <- ape::multi2di(new.tree)
phy_tree(ps) <- new.tree2
phy_tree(ps)


# Calculate microbial diversity ---------------------------------------------

# Change the Phylum, Class, Order, Family as necessary for calculations

# Get the relative abundances
diversity_data <-
  ps %>%
  tax_glom("Species") %>%
  transform_sample_counts(function(x)100* x / sum(x)) %>%
  psmelt() %>%
  as_tibble()

# Highest abundance: all samples pooled together
diversity_data %>%
  group_by(OTU) %>%
  summarise(Abundance = mean(Abundance)) %>%
  arrange(-Abundance)

# Sanity check: is total abundance of each sample 100%?
diversity_data %>%
  group_by(Sample) %>%
  summarise(Abundance = sum(Abundance))

# Get most abundant family for each sample individually and see how many samples
# they are the most abundant in
diversity_data %>%
  group_by(site.name) %>%
  arrange(-Abundance) %>%
  slice(1) %>%
  select(Order) %>%
  ungroup() %>%
  dplyr::count(Order, name = "n_samples") %>%
  arrange(-n_samples)

# Beta diversity ------------------------------------------------------------

# Include the control for the ordinations

# Get the metadata out as separate object
ps.meta <- meta(ps)

# Add the row names as a new column for easy integration later.
ps.meta$sam_name <- rownames(ps.meta)

# Do a NMDS with weighted Unifrac distance for samples
out.wunifrac.nmds <- ordinate(ps, "NMDS", "wunifrac")

# Get the metadata stuff
env <- select(ps.meta, salinity.PSU, temperature.C)
env$salinity.PSU <- as.numeric(env$salinity.PSU)
env$temperature.C <- as.numeric(env$temperature.C)

# Run the envfit function -- removes the control when doing this, don't worry
env_nmds = envfit(out.wunifrac.nmds, env, permutations = 999,na.rm = TRUE)
env_nmds

# Save the environmental vectors and scale them properly -- changed the scale from 0.1 to 0.05 for better fitting
env_nmds_vectors = as.data.frame(scores(env_nmds, "vectors")) * ordiArrowMul(env_nmds, fill = 0.05, at = c(0,0))

# Make the plot by salinity
wunifrac_nmds_plot = plot_ordination(ps, out.wunifrac.nmds, color="salinity.PSU") +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  theme(text = element_text(size = 14)) +  
  geom_hline(yintercept=0, linetype = 'dashed', linewidth = 0.50, color = 'lightgrey') +
  geom_vline(xintercept=0, linetype = 'dashed', linewidth = 0.50, color = 'lightgrey') +
  geom_point(size = 6, alpha= 0.75) +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c() +
  labs(color = "salinity") +
  geom_segment(aes(x = 0, y = 0, xend = env_nmds_vectors[1,1], yend = env_nmds_vectors[1,2]), linewidth =0.5, alpha = 0.5, colour = "grey30", arrow = arrow(length=unit(0.30,"cm"), type = "closed")) +
  geom_segment(aes(x = 0, y = 0, xend = env_nmds_vectors[2,1], yend = env_nmds_vectors[2,2]), linewidth =0.5, alpha = 0.5, colour = "grey30", arrow = arrow(length=unit(0.30,"cm"), type = "closed"))

print(wunifrac_nmds_plot)

# Make the plot by temperature
wunifrac_nmds_plot_temp = plot_ordination(ps, out.wunifrac.nmds, color="temperature.C") +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  theme(text = element_text(size = 14)) +  
  geom_hline(yintercept=0, linetype = 'dashed', linewidth = 0.50, color = 'lightgrey') +
  geom_vline(xintercept=0, linetype = 'dashed', linewidth = 0.50, color = 'lightgrey') +
  geom_point(size = 6, alpha= 0.75) +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c(option="plasma") +
  labs(color = "temperature (C)") +
  geom_segment(aes(x = 0, y = 0, xend = env_nmds_vectors[1,1], yend = env_nmds_vectors[1,2]), linewidth =0.5, alpha = 0.5, colour = "grey30", arrow = arrow(length=unit(0.30,"cm"), type = "closed")) +
  geom_segment(aes(x = 0, y = 0, xend = env_nmds_vectors[2,1], yend = env_nmds_vectors[2,2]), linewidth =0.5, alpha = 0.5, colour = "grey30", arrow = arrow(length=unit(0.30,"cm"), type = "closed"))

print(wunifrac_nmds_plot_temp)

# Arrange them
plot_grid(wunifrac_nmds_plot + theme(legend.justification = c(0,1)), wunifrac_nmds_plot_temp + theme(legend.justification = c(0,1)), ncol=1, align='hv')

# PERMANOVA -----------------------------------------------------------------

# Permutational multivariate analysis of variance (PERMANOVA) tells us if there is a statistical
# difference between metadata groups. PERMANOVA is done by applying the adonis function to a
# distance matrix.

# Normally would check for a sufficient level of homogeneity of dispersion between groups if you
# had a discrete variables, but all of our variables are continuous and it doesn't make sense to
# group them. 

# If you do have discrete variables, you would use the following command:
# anova(betadisper(dist.wuf, sample_info_tab_df$discrete_variable))
# and if the p value is significant (below 0.05) then adonis could be impacted.
# You could disclose this caveat that adonis could be unreliable if your p is 
# significant. 

# We also won't use a strata because there is no obvious reason for clustering co-variables
# to be making a difference in our data. There is only one sample at each site. 

# Must remove the control from the phyloseq object and the metadata for adonis to work because 
# there is no metadata associated with it. 
surface.ps = subset_samples(ps, site.name != "Control")
dist.wuf.surface <- phyloseq::distance(surface.ps, "wunifrac")
info_surface <- subset(sample_info_tab_df, site.name != "Control")

# Use adonis to check salinity, temperature and coastline
adonis2(dist.wuf.surface ~ info_surface$salinity.PSU, perm=999)
adonis2(dist.wuf.surface ~ info_surface$temperature.C, perm=999)

# Alpha diversity -----------------------------------------------------------

# Find the number of unique ASVs in each sample
unique_count = colSums(count_table_df != 0)
unique_count = as.data.frame(unique_count)

# Add this to the metadata 
unique_meta = merge(sample_info_tab_df, unique_count, by = "row.names", all = TRUE)

# Remove control for further analysis, it has 164 unique ASVs (average is 163.15 including control, 163.11 without)
unique_surface_meta <- subset(unique_meta, site.name != "Control")

# Add function to round values
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

# Calculate ASVs vs metadata groups

# Surface ASVs vs salinity
test_asv_sal_surface = cor.test(unique_surface_meta$unique_count, unique_surface_meta$salinity.PSU, method="spearman", exact=FALSE)
spearman_label_asvsal = paste("rho =", round_df(test_asv_sal_surface$estimate, 4), "\n", "p =", round_df(test_asv_sal_surface$p.value,4), " ")

# Surface ASVs vs temperature
test_asv_temp_surface = cor.test(unique_surface_meta$unique_count, unique_surface_meta$temperature.C, method="spearman", exact=FALSE)
spearman_label_asvtemp = paste("rho =", round_df(test_asv_temp_surface$estimate, 4), "\n", "p =", round_df(test_asv_temp_surface$p.value,4), " ")

# Surface ASVs vs distance to coast
test_asv_coast_surface = cor.test(unique_surface_meta$unique_count, unique_surface_meta$coastline.km, method="spearman", exact=FALSE)
spearman_label_asvcoast = paste("rho =", round_df(test_asv_coast_surface$estimate, 4), "\n", "p =", round_df(test_asv_coast_surface$p.value,4), " ")

#Non-phylogenetic Diversities
ps.div <- microbiome::alpha(ps, index = "all")
ps.div

# get the metadata out as separate object
ps.meta <- meta(ps)

# Add the row names as a new column for easy integration later.
ps.meta$sam_name <- rownames(ps.meta)

# Add the row names to diversity table
ps.div$sam_name <- rownames(ps.div)

#Phylogenetic Diversities
ps.asvtab <- as.data.frame(ps@otu_table)

# We first need to check if the tree is rooted or not
ps@phy_tree

#Yes, it is rooted
ps.tree<-ps@phy_tree

# t(ou_table) transposes the table for use in picante and the tree file comes from the first 
# code chunck we used to read tree file
df.pd <- pd(t(ps.asvtab), ps.tree,include.root=T)

# Add the results of PD to this file
ps.meta$Phylogenetic_Diversity <- df.pd$PD

# merge these two data frames into one
div.df <- merge(ps.div,ps.meta, by = "row.names")

# check the tables
#each column name is a different alpha diversity metric
colnames(div.df)

# Exclude the control
div_surface <- subset(div.df, site.name != "Control")

# Calculate PD vs metadata groups

# PD vs salinity
test_pd_sal_surface = cor.test(div_surface$Phylogenetic_Diversity, div_surface$salinity.PSU, method = "spearman", exact = FALSE)
spearman_label_pdsal = paste("rho =", round_df(test_pd_sal_surface$estimate, 4), "\n", "p =", round_df(test_pd_sal_surface$p.value,4), " ")

# PD vs temperature
test_pd_temp_surface = cor.test(div_surface$Phylogenetic_Diversity, div_surface$temperature.C, method = "spearman", exact = FALSE)
spearman_label_pdtemp = paste("rho =", round_df(test_pd_temp_surface$estimate, 4), "\n", "p =", round_df(test_pd_temp_surface$p.value,4), " ")

# PD vs distance to coast
test_pd_coast_surface = cor.test(div_surface$Phylogenetic_Diversity, div_surface$coastline.km, method = "spearman", exact = FALSE)
spearman_label_pdcoast = paste("rho =", round_df(test_pd_coast_surface$estimate, 4), "\n", "p =", round_df(test_pd_coast_surface$p.value,4), " ")

# Test Chao1 vs environmental vairbales
test_chao1_sal_surface = cor.test(div_surface$chao1, div_surface$salinity.PSU, method = "spearman", exact = FALSE)
spearman_label_chao1sal =  paste("rho =", round_df(test_chao1_sal_surface$estimate, 4), "\n", "p =", round_df(test_chao1_sal_surface$p.value,4), " ")
test_chao1_temp_surface = cor.test(div_surface$chao1, div_surface$temperature.C, method = "spearman", exact = FALSE)
spearman_label_chao1temp =  paste("rho =", round_df(test_chao1_temp_surface$estimate, 4), "\n", "p =", round_df(test_chao1_temp_surface$p.value,4), " ")
test_chao1_coast_surface = cor.test(div_surface$chao1, div_surface$coastline.km, method = "spearman", exact = FALSE)
spearman_label_chao1coast =  paste("rho =", round_df(test_chao1_coast_surface$estimate, 4), "\n", "p =", round_df(test_chao1_coast_surface$p.value,4), " ")

# Test Shannon vs environmental variables
test_shan_sal_surface = cor.test(div_surface$diversity_shannon, div_surface$salinity.PSU, method = "spearman", exact = FALSE)
spearman_label_shansal = paste("rho =", round_df(test_shan_sal_surface$estimate, 4), "\n", "p =", round_df(test_shan_sal_surface$p.value,4), " ")
test_shan_temp_surface = cor.test(div_surface$diversity_shannon, div_surface$temperature.C, method = "spearman", exact = FALSE)
spearman_label_shantemp =  paste("rho =", round_df(test_shan_temp_surface$estimate, 4), "\n", "p =", round_df(test_shan_temp_surface$p.value,4), " ")
test_shan_coast_surface = cor.test(div_surface$diversity_shannon, div_surface$coastline.km, method = "spearman", exact = FALSE)
spearman_label_shancoast =  paste("rho =", round_df(test_shan_coast_surface$estimate, 4), "\n", "p =", round_df(test_shan_coast_surface$p.value,4), " ")

# Make alpha diversity heatmap

# Name the metadata groups
env_variables <- c('salinity', 'temperature')

# Get the rho values for each test
ASVs_surface <- c(test_asv_sal_surface$estimate, test_asv_temp_surface$estimate)
faith_surface <- c(test_pd_sal_surface$estimate, test_pd_temp_surface$estimate)
chao1_surface <- c(test_chao1_sal_surface$estimate, test_chao1_temp_surface$estimate)
shannon_surface <- c(test_shan_sal_surface$estimate, test_shan_temp_surface$estimate)

# Get the p values for each test
ASVs_surface_p <- c(test_asv_sal_surface$p.value, test_asv_temp_surface$p.value)
faith_surface_p <- c(test_pd_sal_surface$p.value, test_pd_temp_surface$p.value)
chao1_surface_p <- c(test_chao1_sal_surface$p.value, test_chao1_temp_surface$p.value)
shannon_surface_p <- c(test_shan_sal_surface$p.value, test_shan_temp_surface$p.value)

# Assemble the data frame
alpha_test <- c("# of ASVs","# of ASVs","Faith's PD","Faith's PD","Chao1","Chao1","Shannon Index","Shannon Index")
alpha_env <- c('salinity', 'temperature','salinity', 'temperature','salinity', 'temperature','salinity', 'temperature')
alpha_rho <- c(ASVs_surface, faith_surface, chao1_surface, shannon_surface)
alpha_p <- c(ASVs_surface_p, faith_surface_p, chao1_surface_p, shannon_surface_p)
surface_alpha_diversity <- data.frame(alpha_test, alpha_env, alpha_rho, alpha_p)

# Make the environment variables a factor to order them how you like
surface_alpha_diversity$alpha_env <- as.factor(surface_alpha_diversity$alpha_env)
surface_alpha_diversity$alpha_env <- factor(surface_alpha_diversity$alpha_env, levels = levels(surface_alpha_diversity$alpha_env)[c(2,1)])
alpha_env <- as.factor(alpha_env)
alpha_env <- factor(alpha_env, levels = levels(alpha_env)[c(2,1)])

# Order the diversity tests as you like
surface_alpha_diversity$alpha_test <- as.factor(surface_alpha_diversity$alpha_test)
surface_alpha_diversity$alpha_test <- factor(surface_alpha_diversity$alpha_test, levels = levels(surface_alpha_diversity$alpha_test)[c(1,3,2,4)])

# Round the values to 4 digits
surface_alpha_diversity$alpha_rho <- round_df(surface_alpha_diversity$alpha_rho, 4)
surface_alpha_diversity$alpha_p <- round_df(surface_alpha_diversity$alpha_p, 4)

# Set the limits of the rho colourway
limit <- max(abs(surface_alpha_diversity$alpha_rho)) * c(-1, 1)

# remove the Chao1 from the heatmap because we used DADA2 which removes singletons
surface_alpha_diversity <- surface_alpha_diversity[-c(5,6),]

# Plot graphs

# ASVs vs salinity
asvs_surface_salinity = ggplot(unique_surface_meta, aes(unique_count, salinity.PSU, color=salinity.PSU)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c() +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3, label = spearman_label_asvsal, color="grey30", size=5) +
  labs(color = "salinity") +
  xlab("number of ASVs") + ylab("salinity") +
  theme(legend.justification = c(0,1))

print(asvs_surface_salinity)

# ASVs vs temperature
asvs_surface_temp = ggplot(unique_surface_meta, aes(unique_count, temperature.C, color=temperature.C)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c(option="plasma") +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_asvtemp, color="grey30", size=5) +
  labs(color = "temperature (C)") +
  xlab("number of ASVs") + ylab("temperature (C)") +
  theme(legend.justification = c(0,1))

print(asvs_surface_temp)

# PD vs salinity
pd_surface_salinity = ggplot(div_surface, aes(Phylogenetic_Diversity, salinity.PSU, color=salinity.PSU, label=site.name)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c() +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_pdsal, color="grey30", size=5) +
  labs(color = "salinity") +
  xlab("Faith's PD") + ylab("salinity") +
  theme(legend.justification = c(0,1))

print(pd_surface_salinity)

# PD vs temperature
pd_surface_temp = ggplot(div_surface, aes(Phylogenetic_Diversity, temperature.C, color=temperature.C)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c(option="plasma") +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_pdtemp, color="grey30",size=5) +
  labs(color = "temperature (C)") +
  xlab("Faith's PD") + ylab("temperature (C)") +
  theme(legend.justification = c(0,1))

print(pd_surface_temp)

# Chao1 vs salinity
chao1_surface_salinity = ggplot(div_surface, aes(chao1, salinity.PSU, color=salinity.PSU, label=site.name)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c() +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_chao1sal, color="grey30", size=5) +
  labs(color = "salinity") +
  xlab("Chao1") + ylab("salinity") +
  theme(legend.justification = c(0,1))

print(chao1_surface_salinity)

# Chao1 vs temperature
chao1_surface_temp = ggplot(div_surface, aes(chao1, temperature.C, color=temperature.C)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c(option="plasma") +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_chao1temp, color="grey30", size=5) +
  labs(color = "temperature (C)") +
  xlab("Chao1") + ylab("temperature (C)") +
  theme(legend.justification = c(0,1))

print(chao1_surface_temp)

# Shannon vs salinity
shan_surface_salinity = ggplot(div_surface, aes(diversity_shannon, salinity.PSU, color=salinity.PSU, label=site.name)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c() +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_shansal, color="grey30",size=5) +
  labs(color = "salinity") +
  xlab("Shannon Index") + ylab("salinity") +
  theme(legend.justification = c(0,1))

print(shan_surface_salinity)

# Shannon vs temperature
shan_surface_temp = ggplot(div_surface, aes(diversity_shannon, temperature.C, color=temperature.C)) +
  theme(text = element_text(size = 14)) +
  theme(panel.background=element_blank() ,panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', size = 0.75, fill = NA), legend.key = element_blank(), legend.justification = "top") +
  scale_shape_manual(values = c(16)) +
  scale_color_viridis_c(option="plasma") +
  geom_point(size = 6, alpha= 0.75) +
  annotate("text", x = Inf, y = -Inf, hjust=1.1, vjust=-0.3,  label = spearman_label_shantemp, color="grey30", size=5) +
  labs(color = "temperature (C)") +
  xlab("Shannon Index") + ylab("temperature (C)") +
  theme(legend.justification = c(0,1))

print(shan_surface_temp)

# Arrange them
asv_pd_shan <- plot_grid(asvs_surface_salinity, asvs_surface_temp, pd_surface_salinity, pd_surface_temp, shan_surface_salinity, shan_surface_temp, ncol=2, align='hv')
print(asv_pd_shan)

# Heatmaps -----------------------------------------------------------------

# Transform counts to relative abundance for the top 20 ASVs
relative_ASVs <- transform_sample_counts(ps, function(x) x/sum(x)*100)
relative_ASVs <- prune_taxa(names(sort(taxa_sums(relative_ASVs),TRUE)[1:20]), relative_ASVs)

# Get the data from the plot_heatmap function for the surface samples
raw_heatmap_relative = plot_heatmap(relative_ASVs, method = "NMDS", distance = "wunifrac", taxa.order = "Family", sample.label = "site.name", sample.order = "site.name", trans = NULL) +
  theme(panel.background=element_blank(),panel.grid = element_blank()) +
  labs(x="sample", y="species", fill="abundance") +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(raw_heatmap_relative)

# Extract the data we want
relative_heatmap_data <- raw_heatmap_relative$data
relative_heatmap_data = subset(relative_heatmap_data, select=c(OTU,Abundance,site.name,salinity.level,salinity.PSU,temperature.C, temperature.level, coastline.km, coastline.level, Species, Order, Family, Phylum, Class))

# Make it look pretty

# Get the species correctly sorted 
sample_names(relative_ASVs) <- gsub("-","_", sample_names(relative_ASVs))
sorted_taxa_relative <- phyloseq_to_df(relative_ASVs, sorting="taxonomy")

# Replace the unique species names and edit them accordingly
relative_species_old = unique(sorted_taxa_relative$Species)
relative_species_new = gsub("_", " ", relative_species_old)
relative_species_new = gsub(")", "", relative_species_new)
relative_species_new = gsub("(f", " sp.", relative_species_new, fixed = TRUE)
relative_species_new = gsub(" (g", " sp.", relative_species_new, fixed = TRUE)
relative_species_new = gsub("(g", " sp.", relative_species_new, fixed = TRUE)
relative_species_new = gsub("(c", " sp.", relative_species_new, fixed = TRUE)
relative_species_new = gsub("Clade Ia", "SAR11 Clade Ia", relative_species_new, fixed = TRUE)
relative_species_new = gsub("uncultured Cytophagia", "Cytophagia sp.", relative_species_new, fixed = TRUE)
relative_species_new = gsub("uncultured Pelagibacterales", "Pelagibacterales sp.", relative_species_new, fixed = TRUE)

# Add new column for the index value
new_index = c(1:380)
relative_heatmap_data$new_index = new_index

# Put the edited names in the data we are using for the heatmap - this works
for (i in 1:380){
  if (relative_heatmap_data$Species[i] %in% relative_species_old){
    name_id = which(relative_species_old == relative_heatmap_data$Species[i])
    relative_heatmap_data$Species <- replace(relative_heatmap_data$Species, i, relative_species_new[name_id])
    relative_heatmap_data$new_index <- replace(relative_heatmap_data$new_index, i, name_id)
  }
}

# Make a new column to replace new OTU code in 
relative_ASVcode = c(relative_heatmap_data$Species)
relative_heatmap_data$ASVcode = relative_ASVcode

# Get the last four digits of the OTU
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
code = substrRight(as.character(relative_heatmap_data$OTU), 4)

# Make the new column so that it has the species + the last four digits of the OTU code
for (i in 1:380){
  new_code = paste(relative_heatmap_data$Species[i], " (ASV #", code[i], ")", sep="")
  relative_heatmap_data$ASVcode <- replace(relative_heatmap_data$ASVcode, i, new_code)
}

# Make the column a factor so that it does not sort it
relative_heatmap_data$ASVcode <- as.factor(relative_heatmap_data$ASVcode)
relative_heatmap_data$OTU <- as.factor(relative_heatmap_data$OTU)

# Make a factor of 50 with the unique ASV codes in the correct order - have to do this by hand, in reverse order
# Print the ugly heatmap with no species labels (should be OTU) then search ASVcode_ordered for the last four
# digits in each OTU code on that heatmap, starting at the bottom, and put the index of the ASV code in for ordering
ASVcode_ordered = unique(relative_heatmap_data$ASVcode)
ASVcode_ordered = ASVcode_ordered[c(13,15,18,9,3,19,14,20,5,12,4,2,6,7,16,10,1,8,17,11)]

# Get colours
heatmap_col_relative = c("#D90166", "#7B0393","#0B3AED", "#077399", "#3D9152", "#AFC01F", "#E6CB2B", "#FFAC1F", "#FF7A0A", "#80461B", "#626363")

# Edit the order so that it doesn't have underscores
relative_heatmap_data$Order = gsub("_", " ", relative_heatmap_data$Order, fixed=TRUE)

# Edit the class so that it doesn't have underscores
relative_heatmap_data$Class = gsub("_", " ", relative_heatmap_data$Class, fixed=TRUE)

# Edit the family so that it doesn't have underscores
relative_heatmap_data$Family = gsub("_", " ", relative_heatmap_data$Family, fixed=TRUE)
relative_heatmap_data$Family = gsub("Clade I", "SAR11 Clade I", relative_heatmap_data$Family, fixed=TRUE)

# Edit the phylum so that it doesn't have underscores
relative_heatmap_data$Phylum = gsub("_", " ", relative_heatmap_data$Phylum, fixed=TRUE)

# Make the site name ordered properly
relative_heatmap_data$site.name <- factor(relative_heatmap_data$site.name, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "14", "15", "16", "17", "18", "19", "Control"))

# Use ggplot to plot the new heatmap
relative_heatmap_new <- ggplot(data=relative_heatmap_data, aes(x=reorder(site.name,salinity.PSU), y=OTU)) +
  theme(text = element_text(size = 14)) +
  geom_point(aes(color=Family, size=ifelse(Abundance==0, NA,Abundance))) +
  facet_grid(~salinity.level,scales = "free_x", space = "free_x") +
  scale_color_manual(values = heatmap_col_relative, name="family") +
  scale_size(range=c(0,6), name="relative abundance (%)", breaks=c(0,10,20,30,40,50,60), labels=c("0","10", "20", "30", "40", "50", "60")) +
  labs(x="sample",y="", subtitle = 'salinity') +
  scale_y_discrete(labels=ASVcode_ordered) +
  theme(panel.background=element_blank(), panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.50, fill = NA), legend.justification = "top", legend.key=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.subtitle = element_text(hjust = 0.5, vjust=-2), strip.background = element_rect(fill = "white")) +
  guides(color = guide_legend(override.aes = list(size = 5)))

print(relative_heatmap_new)

# Use ggplot to plot the new heatmap
relative_heatmap_new_temp <- ggplot(data=relative_heatmap_data, aes(x=reorder(site.name,temperature.C), y=OTU)) +
  theme(text = element_text(size = 14)) +
  geom_point(aes(color=Family, size=ifelse(Abundance==0, NA,Abundance))) +
  facet_grid(~temperature.level,scales = "free_x", space = "free_x") +
  scale_color_manual(values = heatmap_col_relative, name="family") +
  scale_size(range=c(0,6), name="relative abundance (%)", breaks=c(0,10,20,30,40,50,60), labels=c("0","10", "20", "30", "40", "50", "60")) +
  labs(x="sample",y="", subtitle = 'temperature (C)') +
  scale_y_discrete(labels=ASVcode_ordered) +
  theme(panel.background=element_blank(), panel.grid = element_line(color = 'white', linewidth = 0.50), panel.border = element_rect(colour = 'black', linewidth = 0.50, fill = NA), legend.justification = "top", legend.key=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.subtitle = element_text(hjust = 0.5, vjust=-2), strip.background = element_rect(fill = "white")) +
  guides(color = guide_legend(override.aes = list(size = 5)))

print(relative_heatmap_new_temp)

## Calculate the abundances for specific ASVs

ASV_e623 <- relative_heatmap_data %>% select(Abundance, ASVcode, site.name) %>% 
  filter(ASVcode == " Polaribacter sp. (ASV #e623)") #%>%
#filter(site.name == "14" | site.name == "15")

mean(ASV_e623$Abundance)
min(ASV_e623$Abundance)
max(ASV_e623$Abundance)

ASV_397e <- relative_heatmap_data %>% select(Abundance, ASVcode, site.name) %>% 
  filter(ASVcode == " SAR11 Clade Ia sp. (ASV #397e)")

mean(ASV_397e$Abundance)

ASV_29af <- relative_heatmap_data %>% select(Abundance, ASVcode, site.name) %>% 
  filter(ASVcode == " Comamonadaceae sp. (ASV #29af)")

mean(ASV_29af$Abundance)

ASV_1af0 <- relative_heatmap_data %>% select(Abundance, ASVcode, site.name) %>% 
  filter(ASVcode == " Nitrincolaceae sp. (ASV #1af0)")

mean(ASV_1af0$Abundance)




