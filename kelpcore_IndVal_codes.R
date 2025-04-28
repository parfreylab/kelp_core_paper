
####IMPORTANT: first, make sure you have data in phyloseq format####
project_data <- readRDS("kelpcore_phyloseq_new.RDS") #to start from here, later, run this command

### INDVAL Indicator Species Analysis ###

library(dplyr)
library(vegan)
library(labdsv)
library(indicspecies)
library(phyloseq)


### from HERE after revising metadata####
metadata <- read.csv(file="genus_level_final.sam.csv")  ### after adding columns manually
View (metadata)
otu_table<- read.csv(file="otu_rev.csv")
#otu_table<- read.csv(file="otu_rev_genus.csv")
View (otu_table)


### subsamples from metadata ### 
metadata <- subset(metadata, sample_type2 %in% c("meristem", "env"))
metadata <- subset(metadata, location != "LAB")
metadata <- subset(metadata, location != "BAM")
#metadata <- subset(metadata, kelp_type != "outplanted")
metadata <- subset(metadata, practitioner != "katy")
metadata <- subset(metadata, kelp_type = "hatchery")
metadata <- subset(metadata, sample_type %in% c("kelp", "water", "rock"))

##OPTIONAL: subsetting samples for evenness
library(dplyr)

# Define a function to sample up to 10 rows for each time point
sample_max_10 <- function(df) {
  if (nrow(df) <= 10) {
    return(df)
  } else {
    return(df[sample(1:nrow(df), 10), ])
  }
}

# Separate the data into two parts: one for subsampling and one to keep intact
metadata_to_sample <- metadata %>%
  filter(location == "GWS", sample_type2 == "meristem")

metadata_others <- metadata %>%
  filter(!(location == "GWS" & sample_type2 == "meristem"))

# Apply the sampling function only to the data that needs subsampling
metadata_sampled <- metadata_to_sample %>%
  group_by(time) %>%
  group_map(~ sample_max_10(.x)) %>%
  bind_rows()

# Combine the sampled data with the rest of the data
metadata <- bind_rows(metadata_sampled, metadata_others)


#View(metadata)
levels(metadata$sample_type2)
metadata$sample_type2 = factor(metadata$sample_type2)
levels(metadata$sample_type2)
table(metadata$sample_type2)

master_table <- left_join(metadata, otu_table, by = "sample_id")
#View(as.data.frame(master_table))
#write.csv(master_table, file="mater_stable_genus.csv", row.names=T)
## remove NA##
#na.omit(master_table)
#na.omit(master_table)
## replace NA with 0##
master_table[is.na(master_table)] <- 0

write.csv(master_table, file="master_table.csv", row.names=T)

### Loading your data (this should contain metadata + taxa abundances)
#fucus <- read.table("your_sample(rows)_taxa(columns)_table.txt", header = T)
kelp <- master_table
head(kelp)

### Set factors you are interested in i.e. want to get core for fucus but using seawater as control, then use the column with info on sample_type

## use one of the interesting info
sample_type <-kelp$sample_type2

class(sample_type)
levels(sample_type)

### Creating an object to store abundances only so you can run the analysis (i.e. remove 7 first columns of metadata with dplyr)
kelp_abund <- kelp %>% dplyr::select(-(1:10)) 
head(kelp_abund)

#### Multipatt analysis: indval with fucus and water ####
multipatt.kelp <- multipatt(kelp_abund, sample_type, control = how(nperm=999))
summary(multipatt.kelp)

## get output and save it
indaval_output <- capture.output(summary(multipatt.kelp, indvalcomp=TRUE))
write.table(as.data.frame(indaval_output), file = "~/desktop/test2.txt", quote=F, row.names=F, col.names=T, sep="n")



