library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the data.table package for better metadata manipulation
library(data.table)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)
#load the vegan library
library(vegan)
library(dplyr)
# filter and reformat data frames
library(tibble)
# Needed for converting column to row names

# Read RDS after DADA2 pipelines
project_data <- readRDS("kelpcore_project_month_added2_1500.RDS")


# subset data (entire datasets but not lab, and other hosts)
project_data <- project_data %>%
  subset_samples(location != "LAB") %>%
  subset_samples(sample_type != "Alaria") %>%
  subset_samples(sample_type != "Alaria_line") %>%
  subset_samples(sample_type != "Alaria_seedstring") %>%
  subset_samples(sample_type != "kelp_line") %>%
  subset_samples(kelp_type != "hatchery_disease") %>%
  #subset_samples(location != "Cascadia") %>%
  #subset_samples(kelp_type != "hatchery") %>%
  #subset_samples(sample_type2 == "meristem")
  subset_samples(sample_type != "nereocystis") %>%
  #subset_samples(sample_type2 != "env")
  subset_samples(sample_type != "nereocystis_line")


# subset data (only kelp samples)
project_data <- project_data %>%
  subset_samples(location != "LAB") %>%
  #subset_samples(practitioner != "katy") %>%
  #subset_samples(practitioner == "katy") %>%
  #subset_samples(sample_type %in% c("kelp"))
  subset_samples(sample_type %in% c("kelp", "kelp_seedstring"))


# subset data (-nocore; katy's ouplanted but initial months)
project_data <- project_data %>%
  subset_samples(nocore != "nocore")

# subset data (katy's ouplanted but initial months)
project_data <- project_data %>%
  subset_samples(sample_type3 != "tissue")
# subset data (katy's ouplanted but initial months)
project_data <- project_data %>%
  subset_samples(kelp_type == "outplanted")

# subset data (only hatchery)
project_data <- project_data %>%
  subset_samples(kelp_type == "hatchery") %>%
  subset_samples(sample_type == "kelp_seedstring")

# subset data (only overlapping months between wild and outplanted)
project_data <- project_data %>%
  subset_samples(time %in% c("APR_2021","JUN_2021", "NOV_2021"))

# subset data (only kelp at GWS for timeseries)
project_data <- project_data %>%
  subset_samples(location == "GWS") %>%
  subset_samples(sample_type %in% c("kelp", "rock"))

# subset data (only wild-kelp)
project_data <- project_data %>%
  subset_samples(kelp_type %in% "wild") %>%
  subset_samples(sample_type %in% "kelp")

# subset data (only compariable month/age)
project_data <- project_data %>%
  subset_samples(comparison %in% "yes")

#### rarefy data ####
set.seed(24)
project_data.rarefied <- rarefy_even_depth(project_data, sample.size = 1500)

#### Create Plotting Objects ####
# 1. reshape data based on taxonomic level you are interested in, and select the top N taxa to show in the plot
taxonomy_plot_obj <- project_data %>%
  tax_glom(taxrank = "asv") # agglomerate at your rank of interest (Rank5 is roughly equal to "family" in our example)

#OPTIONAL, OFTEN RECOMMENDED: select top taxa
# recommend roughly 20 taxa maximum, since it becomes more difficult to distinguish colours with more taxa than that
topOTUs <- names(sort(taxa_sums(taxonomy_plot_obj), TRUE)[1:15]) #where N is the number of taxa you want to retain for plotting purposes
# filter taxa present to just those in the topOTUs list
taxonomy_plot_obj <- prune_taxa(topOTUs, taxonomy_plot_obj) #note that this method of selecting the top OTUs will not show you the proportion of "non-top" OTUs in the plot. for a method that does this, see the last section of this guide.

# REQUIRED: transform to relative abundance, melt to long format for plotting
taxonomy_plot_obj <- taxonomy_plot_obj %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
  psmelt() %>%                                          # Melt to long format
  arrange(Genus)                                        # Arrange by the rank you are going to use in the plot

taxonomy_plot_obj$Genus <- as.character(taxonomy_plot_obj$Genus) #must assign this as character in order to introduce a new string "other"
taxonomy_plot_obj$Genus[-which(taxonomy_plot_obj$OTU %in% topOTUs)] <- "other"
taxonomy_plot_obj$Genus <- factor(taxonomy_plot_obj$Genus)

# OPTIONAL: order levels you are interested in the way you would like them plotted. many ways to do this.
taxonomy_plot_obj$time <- factor(taxonomy_plot_obj$time, levels = c("0DAY", "4DAY", "14DAY", "21DAY", "28DAY")) #method 1

taxonomy_plot_obj$sample_type <- factor(taxonomy_plot_obj$sample_type, levels = c("none", "Maribacter", "Litorimonas", "Sulfitobacter1", "PSD", "Vibrio")) #method 1

taxonomy_plot_obj$Genus <- factor(taxonomy_plot_obj$Genus, levels = c("ASV2.Maribacter", "ASV14.Maribacter", "ASV717.Maribacter", "ASV12.Sulfitobacter", "ASV55.Sulfitobacter", "ASV67.Sulfitobacter", "ASV103.Sulfitobacter", "ASV142.Sulfitobacter", "ASV376.Sulfitobacter", "ASV3.Litorimonas","ASV18.Litorimonas" ,"ASV49.Litorimonas", "ASV5.Pseudoalteromonas", "ASV34.Vibrio", "Nonlabens", "Octadecabacter", "Cellulophaga", "Cocleimonas", "Colwellia", "Ectothiorhodospiraceae", "Granulosicoccus", "Hyphomonadaceae", "Neptuniibacter", "Persicirhabdus", "Roseobacter", "Rubritalea")) #method 1

taxonomy_plot_obj$sample_type <- factor(taxonomy_plot_obj$sample_type, levels = c("none", "Maribacter", "Litorimonas", "Sulfitobacter1", "Sulfitobacter2", "PSD", "Nonlabens", "Octadecabacter", "Vibrio"))

##abundance table
write.csv(taxonomy_plot_obj, file="taxonomy_plot_summary_data2.csv", row.names=T)
##read after the modification
taxonomy_plot_obj<- read.csv(file="taxonomy_plot_summary_data.csv")


#### Make Plots ####
#example plot showing relative abundance of taxa, with panels dividing samples according to two factors (FACTOR_1 and FACTOR_2)

pdf("Taxonomy_plot_kelpcore_Lab_modified_ASVs_Genus.pdf"
    , width = 13 # Default is 7
    , height = 10 # Change to 10; make it taller
)
ggplot(taxonomy_plot_obj, aes(x = time, y = Abundance, fill = Genus)) + 
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  facet_wrap(~ sample_type, drop=TRUE, scales="free_x") +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) + #geom_bar controls the bar plots themselves
  scale_y_continuous(expand = c(0.005,0.005)) + #this controls the y axis scale, for bigger margins above and below, increase the numbers provided
  scale_fill_manual(values=c("#e28200","#fa8351", "#e4c098","#009f9a",
                             "#14dec3",
                             "#05865a",
                             "#2bd2e7",
                             "#68caca",
                             "#68d2e2", "#cc728f",
                             "#ee88f3",
                             "#ff60c2", "#6681ff", "#e04d36", "#00106e", "#35618c", "#a6d6ea", 
                             "#dcc6a7",
                             "#88aee1",
                             "#aec7a3",
                             "#95bbef", 
                             "#d3eed4",
                             "#d6bee2",
                             "#8ecdc7",
                             "#e7b8b7",
                             "#abc5bf",
                             "#acb5d1", "#2d2032")) +
  #scale_fill_manual(values = cbPalette) + #example with colourblind palette
  #theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  ggtitle("My Plot Title") + #plot title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")) #another "theme" command. a lot of extra control can be established here. this line ensures that there is no padding in the multi-plot grid
dev.off()

##order designation##
project_data$time <- factor(project_data$time, level = c("JAN_2019", "JULY_2020",	"AUG_2020",	"SEP_2020",	"NOV_2020",	"JAN_2021", "FEB_2021",	"MAR_2021",	"APR_2021",	"MAY_2021",	"JUN_2021", "JULY_2021", "OCT_2021",	"NOV_2021"))


###########Beta Diversity (NMDS plots) ##############
set.seed(6)
NMDS.bray <- ordinate(
  physeq = project_data,
  method = "NMDS",
  distance = "bray"
) # you can choose different methods and distance metrics, see the ordinate help page for details. this function works with "phyloseq" class objects.

#### making beta div plots ####
#we get more plotting control if we don't use the phyloseq plotting functions for ordination plots, and instead add the results of the ordination to our existing metadata
NMDS <- as.data.frame(unclass(sample_data(project_data)))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2

###############NMDS PLOT: Sample Type###############
pdf("NMDS_additional analysis6_NOV2024.pdf"
    , width = 7 # Default is 7
    , height = 6 # Change to 10; make it taller
)

ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color=month, shape=kelp_type)) +
  geom_point(size=5, alpha = 0.7) + 
  
  scale_colour_manual(values=c("#00BFFF","#32CD32","#DC143C", "#FF4500", "#FFA500", "#FFD700","#800080", "#483D8B", "#4682B4", "#5F9EA0", "#39FF14", "#000000", "#DC143C", "#8B0000", "#FF4500", "#FFA500", "#800080", "#ADFF2F")) +
  #scale_colour_manual(values=c("#FFA500", "#800080","#3CB371", "#FFA500", "#800080", "#ADFF2F")) +
  #scale_colour_manual(values=c("#238B45", "#00539D")) +
  
  #facet_wrap(~ time, drop=TRUE, scales="free") +
  scale_shape_manual(values=c(16, 8, 15)) +
  theme_bw() +
  labs(title="NMDS_kelpcore_sampletype_time_coloured") + 
  scale_fill_manual(values = c("2020" = "solid", "2021" = "open")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#xlim(0.8, -0.8) + ylim(-0.4, 0.5)
dev.off()

###############NMDS PLOT: filled color###############
pdf("NMDS_additional analysis9_Oct2024.pdf"
    , width = 9 # Default is 7
    , height = 6 # Change to 10; make it taller
)
p <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, color=time, shape=location, fill=sample_type)) +
  geom_point(size=5, alpha=0.85, stroke=1.5) + 
  scale_colour_manual(values=cbPalette) +  # Use cbPalette for time color
  
  # Shapes for 'location' and fill for 'sample_type'
  scale_shape_manual(values=c(21, 24, 25, 22)) +  
  scale_fill_manual(values=c("#006400", "grey", "blue")) +  # Ensure fill in the legend reflects correct colors
  
  theme_bw() +
  labs(title="NMDS_kelpcore_sampletype_time_coloured") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  
  guides(fill=guide_legend(override.aes=list(shape=21), ncol=2),  # Make the fill legend two columns
         color=guide_legend(ncol=2))  # Make the color legend two columns

p

dev.off()


##PERMANOVA TEST##

#subsampling for targetting comparisons and match with NMDS results
project_data <- project_data %>%
  #subset_samples(sample_type2 %in% c("meristem", "env")) %>%
  subset_samples(kelp_type != "outplanted") %>%
  subset_samples(sample_type %in% c("kelp", "kelp_seedstring"))

project_data <- project_data %>%
  subset_samples(sample_type %in% c("kelp_seedstring", "kelp")) %>%
  subset_samples(kelp_type != "outplanted")

project_bray <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))


#adonis(project_bray ~ Type, data=sample_df, method="bray")
adonis2(project_bray ~ sample_type3,
        data= sample_df , permutations=9999, by = "margin")


## betadispersion test
#need a distance matrix and the sample data (metadata) to start
project_bray.rarefied <- phyloseq::distance(project_data, method = "bray")
sample_df <- data.frame(sample_data(project_data))
beta.FACTOR1 <- betadisper(project_bray.rarefied, sample_df$sample_type) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
b3.1=permutest(beta.FACTOR1)
capture.output(b3.1,file="betadispersion_wild_nuresery.doc")
beta.FACTOR2 <- betadisper(project_bray.rarefied, sample_df$blade_health) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
b3.2=permutest(beta.FACTOR2)
capture.output(b3.2,file="betadispersion_field2020-bladehealth.doc")



### bubble plot with relative abundance and prevalence####
###upload data
Long_np<- read.csv("Averages_for_Dot_plot_all_location.csv")
Long_np<- read.csv("Averages_for_Dot_plot_outplant_time.csv")
Long_np<- read.csv("Averages_for_Dot_plot_wild_seasonality_expanded.csv")
head(Long_np)
Long_np<- read.csv("Averages_for_Dot_plot_all_wild_core_specificity.csv")
##Genus
Long_np<- read.csv("Averages_for_Dot_plot_wild_seasonality_expanded_genus.csv")
Long_np<- read.csv("Averages_for_Dot_plot_wild_seasonality_expanded_ASV.csv")
Long_np<- read.csv("Averages_for_Dot_plot_outplant_time_2groups.csv")
Long_np<- read.csv("Averages_for_Dot_plot_outplant_time_ASV.csv")


#### Go to wide/fat format
Long_np$taxa <- factor(Long_np$taxa, level = c("ASV1.Persicirhabdus","ASV2.Maribacter","ASV3.Litorimonas","ASV4.Cocleimonas","ASV6.Granulosicoccus","ASV7.Hellea","ASV10.Rubritalea","ASV11.Ectothiorhodospiraceae","ASV15.Granulosicoccus", "ASV19.Blastopirellula","ASV27.Loktanella","ASV29.Marixanthomonas"))

Long_np$taxa <- factor(Long_np$taxa, level = c("ASV1.Persicirhabdus","Others.Persicirhabdus","ASV2.Maribacter","Others.Maribacter","Litorimonas.ASV3","Others.Litorimonas","ASV6.Granulosicoccus","ASV15.Granulosicoccus","Others.Granulosicoccus","ASV7.Hellea","Others.Hellea","ASV4.Cocleimonas","Others.Cocleimonas","Ectothiorhodospiraceae.ASV11","ASV29.Marixanthomonas"))

Long_np$taxa <- factor(Long_np$taxa, level = c("ASV10.Rubritalea","ASV19.Blastopirellula","ASV27.Loktanella"))

Long_np$sample <- factor(Long_np$sample, level = c("GWS_JULY_2020","GWS_AUG_2020","GWS_SEP_2020", "GWS_NOV_2020","GWS_JAN_2021","GWS_MAR_2021","GWS_APR_2021","GWS_MAY_2021","GWS_JUN_2021","GWS_JULY_2021","GWS_NOV_2021","GWS_rock","GWS_water","TB_APR_2021","TB_MAY_2021", "TB_JUN_2021","TB_rock","TB_water","LHP_APR_2021","LHP_MAY_2021","LHP_JUN_2021", "LHP_JULY_2021","LHP_rock","LHP_water","SCP_APR_2021","SCP_MAY_2021","SCP_JUN_2021","SCP_JULY_2021","SCP_rock","SCP_water"))

Long_np$sample <- factor(Long_np$sample, level = c("LI1_FAB_2021","LI1_APR_2021","LI1_JUN_2021","LI2_FAB_2021","LI2_APR_2021","LI2_JUN_2021","EWB_FAB_2021","EWB_APR_2021","EWB_JUN_2021","INT_FAB_2021","INT_APR_2021","INT_JUN_2021","TAL_FAB_2021","TAL_APR_2021","TAL_JUN_2021","HHF_JAN_2019","HHF_MAR_2019","BAM_DEC_2021", "NIC_DEC_2020","NOAA_JAN_2019","Cascadia_OCT_2021"))

###Make dot plot 
pdf("bubble_plot_seasonality_industrial_genus.pdf"
    , width = 12 # Default is 7
    , height = 8 # Change to 10; make it taller
)
#png("bubble_plot_all_new_prevalance_test.png", width = 11, height = 10, units = 'in', res = 600)

Long_np$Relative_Abundance <- as.numeric(as.character(Long_np$abundance))
Dot_plot_Survey_core<- ggplot(Long_np, aes(y= taxa, x = sample)) + 
  geom_point(aes(size = abundance, color = prevalence))  +
  scale_colour_gradient(low = "white", high = "black") +
  #facet_wrap(~ location, drop=TRUE, scales="free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  #scale_size(range=c(2,8)) + theme_bw() +
  scale_size(range=c(1,9), breaks=c(0.001,0.01,0.1,0.5)) + theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black"))  +
  theme(axis.text.x=element_text(size=12, angle= 45, hjust= 1), axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=8, face="bold"), axis.title.y=element_text(size=9, face="bold")) + ylab("Core Taxa") + xlab("Time Series") + ggtitle("")

Dot_plot_Survey_core

dev.off()


#################################
### kelp core box plot with scatter dots####
theme_set(theme_bw())
setwd("~/desktop/Manipulation for Kelp and Bacteria/R/box_plot")
data <- read.csv("categorical_core_effects.csv")
data <- read.csv("categorical_core_effects_focual_species.csv")

# Plot
pdf("boxplot_gameto_scaled_corenoncore.pdf"
    , width = 8 # Default is 7
    , height = 11 # Change to 10; make it taller
)
data %>%
  ggplot( aes(x=taxa, y=sporo, fill=taxa)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="grey20", size=1, alpha=0.9) +
  theme_bw() +
  scale_y_continuous(limit = c(0,30)) +
  facet_wrap(~ round, drop=TRUE, scales="free") +
  theme(
    legend.position="none",
    plot.title = element_text(size=13)
  ) +
  #stat_compare_means(method = "anova") +
  stat_compare_means(method = "t.test") + 
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 1) +
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  ggtitle("A boxplot with jitter") +
  xlab("")

dev.off()

#subset each round
df_group1 <- subset(df, round == "1")
df_group2 <- subset(df, round == "2")
df_group3 <- subset(df, round == "3")
df_group4 <- subset(df, round == "4")
# Perform the ANOVA
aov_result <- aov(gameto ~ taxa, data = df_group4)
summary(aov_result)
# a,b,c,d group assignment for each categorical group
library(agricolae)
tukey_result <- HSD.test(aov_result, "taxa", group = TRUE)
print(tukey_result)

df <- data.frame(data)
compare_means(gameto ~ taxa, data = df, method = "pairwise.t.test")
# Pairwise comparison against reference
stat<-compare_means(gameto ~ taxa,  data = df, ref.taxa = "control",
                    method = "anova")


my_comparisons <- list( c("control", "core"), c("core", "non_core"), c("control", "non_core") )


# Calculate the mean and standard deviation of the mpg variable by cylinder
summary <- data %>%
  group_by(taxa, round) %>%
  summarize(mean = mean(sporo), sd = sd(sporo))
#subset adjusted/non-adjusted
data <- subset(data, group == "adjusted")
data <- subset(data, taxa != "control")


data$s_log2fold <- as.numeric(data$s_log2fold)
n <- as.numeric(data$s_log2fold)

pdf("barplot_sporo_log2fold.pdf"
    , width = 8 # Default is 7
    , height = 9 # Change to 10; make it taller
)

ggplot(data, aes(x = taxa, y = sporo, fill=taxa)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  #scale_y_continuous(breaks = seq(-1,1,0.2),limits = c(-1,1),
  #labels = scales::percent)+
  theme_bw() +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(.9)) +
  facet_wrap(~ round, drop=TRUE, scales="free") +
  theme(axis.text.x = element_text(size = 15, angle = -45, hjust = 0)) +
  xlab("")

dev.off()

## Volcano plot (log fold change)

tmp <- read.csv(file="log2fold_final_rawdata.csv")

tmp <- tmp %>% filter(treatment != "Control")
tmp <- tmp %>% filter(s_pvalue != "1")
#tmp <- tmp %>% filter(trial == "1")

View(tmp)
# remove rows that contain NA values
de <- tmp[complete.cases(tmp), ]
# The basic scatter plot: x is "log2FoldChange", y is "pvalue"
p <- ggplot(data=de, aes(x=s_log2FoldChange, y=-log10(s_BH_P_value))) + geom_point() + theme_minimal()
p

p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="grey80") +
  geom_hline(yintercept=-log10(0.05), col="grey80")
p2
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$s_log2FoldChange > 0.6 & de$s_BH_P_value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$s_log2FoldChange < -0.6 & de$s_BH_P_value < 0.05] <- "DOWN"


# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de, aes(x=s_log2FoldChange, y=-log10(s_BH_P_value), col=diffexpressed)) + geom_point() + theme_minimal()
p
# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="grey80") +
  geom_hline(yintercept=-log10(0.05), col="grey80")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("#FF4500", "black", "royalblue"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("#FF4500", "black", "royalblue")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=s_log2FoldChange, y=-log10(s_BH_P_value), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()



library(ggrepel)
# plot adding up all layers we have seen so far
pdf("logfoldchange_coculture_sporo_all_nov.21.23_label.pdf"
    , width = 9 # Default is 7
    , height = 5 # Change to 10; make it taller
)

ggplot(data=de, aes(x=s_log2FoldChange, y=-log10(s_BH_P_value), col=diffexpressed, label=sample)) +
  geom_point(alpha=0.8, size=5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~ trial, drop=TRUE, scale="free_y") +
  geom_text_repel(data=subset(de, s_BH_P_value < 0.05), size = 3, max.overlaps = 20) +
  theme(axis.text.x=element_text(size=20)) +
  theme(axis.text.y=element_text(size=20)) +
  scale_color_manual(values=c("#FF4500", "black", "royalblue")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="grey80") +
  theme_bw() +
  geom_hline(yintercept=-log10(0.05), col="grey80") +
  scale_y_log10() 

dev.off()



