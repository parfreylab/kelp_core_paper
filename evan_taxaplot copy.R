##### set up #####
library(tidyverse)
install.packages("phyloseq")
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")
library(phyloseq)
library(plyr)
library(qualpalr)
library(ggh4x)
#library(vegan)
#library(ggpattern)
library(ggpubr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12),
                                   axis.text.x = element_text(colour = "black", size = 10),
                                   legend.text = element_text(size = 9, colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   # axis.ticks = element_blank(),
                                   panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank()
                             ))

## load functions
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
  tax = tax %>% rownames_to_column(var="ASVid")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "ASVid", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}


## read in data
all = readRDS("/Users/evankohn/Desktop/NSERC work/kelpcore_phyloseq_new.RDS")
#View(all@sam_data)

print(c(unique(all@sam_data$kelp_type)))


##### SUBSET SAMPLES #####

kelp = subset_samples(all, sample_type=="kelp" &
                        kelp_type %in% c("wild", "outplanted"))

print(c(unique(kelp@sam_data$kelp_type)))

#### SET UP GROUPS FOR NURSERY KELP ####

## group at genus level
kelp.gen = tax_glom(kelp, taxrank = "Genus")

## calculate read depth
kelp.gen@sam_data$read_depth = sample_sums(kelp.gen)

## get data out of phyloseq
kelpdf = dephyloseq(kelp.gen)

## calculate relative abundance
kelpdf$ra = as.numeric(kelpdf$asv_abundance)/as.numeric(kelpdf$read_depth)

## make plotnames
kelpdf$plotnames = paste0(kelpdf$Order,"; ", kelpdf$Genus)

## summarize data by taxaplot group type. 
kelp.sum = ddply(kelpdf, c("plotnames", "kelp_type"),
                summarise,
                sum = sum(ra))

## sort data by relative abundance. This is how the loop will pick the mos tabundant taxa
sorted = kelp.sum[order(-kelp.sum$sum),]


#### PICK TOP TAXA FOR EACH GROUP ####
## make empty dataframe to store output from the loop
top.df = NULL

## start loop
for(i in c(unique(kelp.sum$kelp_type))) {
  for(j in i) {
    
    ## subset dataframe by samples
    #!# Remeber to change te substrate to your group! 
    sample = subset(sorted, sorted$kelp_type %in% c(j))
    
    ## get top 15 genera
    top = sample[c(1:10),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)
    
    ## close loop 
  }
}


##### SUBSET AND JOIN DATAFRAMES #####
toplist = c(unique(top.df$plotnames))

## only keep target taxa and change others to Others
kelpdf.tops = subset(kelpdf, kelpdf$plotnames %in% c(toplist))

## calculate sample number for each group
samplegroups = ddply(kelpdf.tops, c("kelp_type", "location", "time", "Row.names"),
           summarise,
           nrows=length(ra))
samplegroups = ddply(samplegroups, c("kelp_type", "location", "time"),
                     summarise,
                     nsamples=length(nrows))


## add the sampel counts to the data
kelpdf.tops = left_join(kelpdf.tops, samplegroups)

## group the dataframe to get mean RA for each group
tp = ddply(kelpdf.tops, c("kelp_type", "location", "time", "plotnames", "nsamples"),
           summarise,
           meanra = mean(ra))

## calculate others
tpothers = ddply(tp,  c("kelp_type", "location", "time", "nsamples"),
                 summarise,
                 sumra = sum(meanra))
tpothers$meanra = 1-tpothers$sumra
tpothers$plotnames = "Others"

## combine mean per-top taxa relative abundance with the others caluclation
alldata = full_join(tp, tpothers)

#### GET COLORS FOR TAXAPLOT #####

# 1. find out how many colors you need
numcol <- length(unique(alldata$plotnames))

# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(3)

# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. Extract hex colors
hex = as.data.frame(newpal$hex)
colnames(hex) <- c("taxa_color")

# 5. Get list of taxa
tops = as.data.frame(c(unique(alldata$plotnames)))
colnames(tops) <- c("plotnames")

# 6. Join color list and taxa names
topcolors = cbind(tops, hex)

# 7. for the "others" plot name, replace that with grey 90 (this is just an astetic thing)
topcolors[topcolors$plotnames == "Others",]$taxa_color <- "grey90"


## freshet colors
topcolors[topcolors$plotnames == "Verrucomicrobiales; Persicirhabdus",]$taxa_color = "#3C2266"
topcolors[topcolors$plotnames == "Caulobacterales; Litorimonas",]$taxa_color ="#989433"
  topcolors[topcolors$plotnames == "Flavobacteriales; Maribacter",]$taxa_color = "#378E90"
  topcolors[topcolors$plotnames == "Thiotrichales; Cocleimonas",]$taxa_color ="#8F4234"
  #topcolors[topcolors$plotnames == "Granulosicoccales; Granulosicoccus",]$taxa_color ="#1C523A"
  #topcolors[topcolors$plotnames == "Caulobacterales; Robiginitomaculum",]$taxa_color = "#8A3E75"
  #topcolors[topcolors$plotnames == "Alteromonadales; Pseudoalteromonas",]$taxa_color = "#534A1C"
  #topcolors[topcolors$plotnames == "Gammaproteobacteria; Gammaproteobacteria",]$taxa_color ="#339448"
    topcolors[topcolors$plotnames == "Verrucomicrobiales; Rubritalea",]$taxa_color = "#44272C"
    #topcolors[topcolors$plotnames == "Rhodobacterales; Yoonia-Loktanella",]$taxa_color ="#273542"


# 8. Make an object R can pull form for the colors
plotcolors <- topcolors$taxa_color
names(plotcolors) <- topcolors$plotnames

## ORDER THE TAXA SO OTHERS ARE AT THE BOTTOM #####
## order by decreasing relative abundance
alldata = alldata[order(-alldata$meanra),]

## get list of factors in order
natural.genus.order = as.list(c(unique(alldata$plotnames)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")


# create new column with date and sample for x-axis
## x=paste0(time," (",nsamples,")")

alldata$sampletime = paste0(alldata$time," (",as.character(alldata$nsamples),")")

print(c(unique(alldata$sampletime)))

## order dates
alldata$sampletime = factor(alldata$sampletime, levels=c("JAN_2019 (3)","MAR_2019 (8)", "JULY_2020 (9)"   ,"AUG_2020 (9)","SEP_2020 (8)"  , "NOV_2020 (10)","JAN_2021 (10)","FEB_2021 (5)", "MAR_2021 (10)", "APR_2021 (19)", "APR_2021 (9)", "APR_2021 (12)", "APR_2021 (8)", "APR_2021 (5)", "MAY_2021 (10)", "MAY_2021 (9)", "JUN_2021 (3)", "JUN_2021 (11)", "JUN_2021 (4)", "JUN_2021 (10)", "JUN_2021 (5)", "JUN_2021 (7)", "JUN_2021 (9)", "JULY_2021 (2)","JULY_2021 (3)", "NOV_2021 (16)", "NOV_2021 (10)"))


## set plot_names levels and order by realtive abundance (do this last!)
plot.order = unlist(plot.order)
alldata$plotnames = factor(alldata$plotnames, levels=c(plot.order))

##### MAKE PLOT #####
## remove emppty data
alldata = subset(alldata, alldata$meanra !="NaN")

## abbreviate locations
alldata <- alldata %>%
  mutate(location = case_when(
    location == "Talbot" ~ "TB",
    location == "East West Bay" ~ "EWB",
    location == "site 1 Loughboro Inlet" ~ "LI 1",
    location == "site 2 Loughboro Inlet" ~ "LI 2",
    location == "Hood Head Farm" ~ "HHF",
    location == "Interfor" ~ "INT",
    TRUE ~ location 
  ))

## reorder plot
alldata$kelp_type <- factor(alldata$kelp_type, levels = c("wild", "outplanted"))


ggplot(alldata, aes(x = sampletime,
               y=meanra,
               fill=plotnames))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=plotcolors)+
   facet_nested(.~kelp_type+location, scales="free", space="free")+
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  labs(x= NULL, y="Mean Relative Abundance", fill= "Taxa (Order;Genus)")
