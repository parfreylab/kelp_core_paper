## load libraries
library(phyloseq)
library(tidyverse)
library(indicspecies)
library(data.table)

dephyloseq = function(phylo_obj){
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  ## how many metadata columns you have
  metacols = ncol(meta)+1
  ## get out the otu table
  ## if your metadta is empty after running this, you need to use
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  #otu = as.data.frame(as.matrix(phylo_obj@otu_table))
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  ## get out the taxonomy file
  tax = as.data.frame(phylo_obj@tax_table)
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_name")
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_name", values_to="asv_abundance")
  ## Join the metadata and otu table with the taoxnomy table
  mot = full_join(mo, tax)
  ## Specify the output for the dephyloseq funciton
  output = mot
}

## read in non-rarefied data
saccharinaData = readRDS("/Users/evankohn/Desktop/NSERC work/kelpcore_phyloseq_new.RDS")
saccdata_core = subset_samples(saccharinaData,
                          kelp_type == "wild")
saccdata_core = tax_glom(saccdata_core, taxrank = "Genus")
## caulcuate sample sums of  dataset
saccdata_core@sam_data$read_depth = sample_sums(saccdata_core)

## getting data needed
metadata = as.data.frame(saccdata_core@sam_data)
#rotating data set and then reformating as dataframe
otu = as.data.frame(saccdata_core@otu_table) %>% t() %>% as.data.frame()#rotating data set and then reformating as dataframe

## metadata_core = subset(metadata,
##             metadata$nocore == "core")
## otu_core = 

## create a value with the comparison stored
comparison = metadata$sample_type2

## run indval to compare the host species
indval <- multipatt(otu, comparison, duleg = TRUE, control = how(nperm=999))

# Get indval statistic
indval.str <- as.data.frame(indval$str)
indval.str$rn <- rownames(indval.str)

# get p-value
indval.stat <- as.data.frame(indval$sign) #get dataframe of indval statistic
indval.stat$rn <- rownames(indval.stat) # make column of ASVs

# Prevalence as dataframe
indval.prev <- as.data.frame(indval$A)
# extract rownames into column
setDT(indval.prev, keep.rownames = TRUE)[]

## rename columns
colnames(indval.prev) <- paste0("prev.", colnames(indval.prev))
names(indval.prev)[names(indval.prev) == 'prev.rn'] <- 'rn'
# Specificity as dataframe
indval.sp <- as.data.frame(indval$B)
# extract rownames into column
setDT(indval.sp, keep.rownames = TRUE)[]

## rename columns
colnames(indval.sp) <- paste0("sp.", colnames(indval.sp))
names(indval.sp)[names(indval.sp) == 'sp.rn'] <- 'rn'
# Join statistics together
str.and.stat = full_join(indval.str, indval.stat,
                         by="rn")
prev.and.fid = full_join(indval.prev, indval.sp,
                         by="rn")
indval_table = full_join(str.and.stat, prev.and.fid,
                         by="rn")

## get taxonomy table form phyloseq object
tax = as.data.frame(saccdata_core@tax_table)
## get the ASV ID into a new column
tax = tax %>% rownames_to_column(var="asv_name")
## rename columns fromt he IndVal output to join with taxonomy
names(indval_table)[names(indval_table) == 'rn'] <- 'asv_name'
## merge with taxonomy
indval_table= inner_join(indval_table, tax)

## Make indVal relative to sacc
indval_table_new <- indval_table %>%
  mutate(stat = if_else(s.env == 1, 1 - stat, stat))

## save a .csv with the indval output
write.csv(indval_table_new, "/Users/evankohn/Desktop/NSERC work/indval_output_genus_UPDATED.csv")



