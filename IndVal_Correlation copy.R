install.packages("lme4")
install.packages("Matrix")
library(ggplot2)
library(tidyverse)
library(scales)

## read in data, clean it, merge it ####
CoCultureData <- read.csv("/Users/evankohn/Desktop/NSERC work/Chapter3_co-culture_results_Appendix3.Table4_lwp.csv")
IndValData <- read.csv("/Users/evankohn/Desktop/NSERC work/indval_output_genus_UPDATED.csv")
AuxData <- read.csv("/Users/evankohn/Desktop/Honors/data and scripts/Bacteria_IAA_Production_Data+Genus.csv")
FrequencyData <- read.csv("/Users/evankohn/Desktop/NSERC work/Genus_Frequency.csv")

## Filtering and Merging ####
## Remove unidentified and control treatments 
CoCultureData %>% 
  filter(Genus != "Unidentified",
         Genus != "Control") -> CoCultureData_Clean
IndValData %>% 
  filter(stat != "NA") ->IndValData_Clean
AuxData %>% 
  filter(Genus != "NA") %>% 
  mutate(across(IAA, ~ replace(., . < 0, 0))) %>% 
  subset(select = -c(Genus, Family)) -> AuxDataClean
colnames(AuxDataClean)[1] <- "Isolate.Identifier"

## Merge coculture data with frequency data for mixed model
Merged_data <- merge(CoCultureData_Clean, IndValData_Clean, by = "Genus", all = TRUE)

## Filter out Genera that weren't tested
Merged_data %>% 
  filter(Isolate.Identifier != "NA") -> Merged_data_clean

## Check Correlations ####
## filter to t1
Merged_data_clean %>% 
  filter(Co.culture.Trial == "1") -> Trial1_data

## correlation double check....
cor.test(Trial1_data$Log2.Fold.Change, Trial1_data$stat)
## p = .002457, r = 0.40
cor.test(Trial1_data$Gametophyte_fold_change, Trial1_data$stat)
## p = 0.026, r = 0.30 

## regression checks
t1_regression <- lm(Log2.Fold.Change. ~ stat, data = Trial1_data)
summary(t1_regression)
## p = 0.00246, R^2 = 0.16

t1_regression_gameto <- lm(Gametophyte_fold_change ~ stat, data = Trial1_data)
summary(t1_regression_gameto)
## p = 0.0266, R^2 = 0.09

## Linear Mixed effects Model ####
library(Matrix)
library(lme4)
install.packages("lmerTest")
library(lmerTest)

model <- lmer(Log2.Fold.Change. ~ stat + (1 | Co.culture.Trial), data = Merged_data_clean)
summary(model)
anova(model)

## OlD: p = 0.0281, F = 4.9792 , indval stat is a useful predictor of log2.fold.change in sporo
## NEW: p = 0.0005158, F = 12.957

model2 <- lmer(Gametophyte_fold_change ~ stat + (1 | Co.culture.Trial), data = Merged_data_clean)
summary(model2)
anova(model2)
## p = 0.01673, F = 5.937

## multiple linear regression ####
# make dataset
Merged_data_Aux <- merge(Merged_data_clean, AuxDataClean, by = "Isolate.Identifier")

Merged_data_Aux %>% 
  filter(Isolate.Identifier != "NA") -> Merged_data_clean_Aux

multipleReg <- lm(Log2.Fold.Change. ~ stat + IAA, data = Merged_data_clean_Aux)
summary(multipleReg)

## Merge coculture data with frequency data
Merged_data_freq <- merge(CoCultureData_Clean, FrequencyData, by = "Genus", all = TRUE)

## Filter out Genera that weren't tested
Merged_data_freq %>% 
  filter(Isolate.Identifier != "NA") -> Merged_data_clean_freq
Merged_data_clean_freq %>% 
  filter(Co.culture.Trial == "1") -> Regression_data_freq


## Aux, IndVal Relationship ####

## load and filter IAA yes/no data
IAA_yn_Data <- read.csv("/Users/evankohn/Desktop/Honors/data and scripts/IAA_Jungsoo_data_YN.csv")
IAA_yn_Data %>% 
  drop_na(IAA) -> Trial_1_yn_data

## merge with stat data
Merged_data_Auxyn <- merge(IndValData_Clean, IAA_yn_Data, by = "Genus")
Merged_data_Auxyn$IAA_Yes_No <- as.factor(Merged_data_Auxyn$IAA_Yes_No)

Merged_data_Auxyn %>% 
  filter(IAA != "NA") -> Aux_yesno_clean

## Check t-test assumptions:
library(car)
ggplot(Aux_yesno_clean, aes(x = stat)) +   
  geom_histogram(binwidth = 0.25)
leveneTest(data = Aux_yesno_clean, stat ~ IAA_Yes_No, center = mean)
## p = 0.095, variance equal

t.test(stat ~ IAA_Yes_No, data = Aux_yesno_clean, var.equal = FALSE)
# p = 0.08415


## plots
Aux1 <- ggplot(Aux_yesno_clean, aes(x=factor(IAA_Yes_No, levels = c("Yes", "No")), y=stat, fill=IAA_Yes_No)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.6) +  # Add jittered points
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, fill="black") +
  scale_fill_manual(values = c("#FF9999", "#9999FF")) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 8),     # Adjust Y-axis title size
    legend.position = "none",   
    axis.text = element_text(size = 6) ) +
  xlab("")+
  ylab("Strength of Association with S. latissima")+
  theme(legend.position = "none") +
  annotate("text", x = 1, y = 1, label = "a", size=4, color="black")+
  annotate("text", x = 2, y = 0.8, label = "b", size=4, color="black")
Aux1

t.test(Log2.Fold.Change.sporo ~ IAA_Yes_No, data = Trial_1_yn_data, var.equal = TRUE)
# p = 0.5825

t.test(Log2.Fold.Change.gameto ~ IAA_Yes_No, data = Trial_1_yn_data, var.equal = TRUE)
# p = 0.3052

Aux2 <- ggplot(Trial_1_yn_data, aes(x=factor(IAA_Yes_No, levels = c("Yes", "No")), y=Log2.Fold.Change.sporo, fill=IAA_Yes_No)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.6) +  # Add jittered points
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, fill="black") +
  scale_fill_manual(values = c("#FF9999", "#9999FF")) +
  theme_minimal() +
  ylim(-4, 3)+
  theme(
    axis.title.y = element_text(size = 8), 
    axis.title.x = element_text(size = 8),
    legend.position = "none",   
    axis.text = element_text(size = 6)) +
  xlab("Auxin Production?")+
  ylab("Log2FoldChange in Sporophyte Number")+
  theme(legend.position = "none")

Aux3 <- ggplot(Trial_1_yn_data, aes(x=factor(IAA_Yes_No, levels = c("Yes", "No")), y=Log2.Fold.Change.gameto, fill=IAA_Yes_No)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.6) +  # Add jittered points
  stat_summary(fun=mean, geom="point", shape=21, size=2.5, fill="black") +
  scale_fill_manual(values = c("#FF9999", "#9999FF")) +
  theme_minimal() +
  ylim(-4, 3)+
  theme(
    axis.title.y = element_text(size = 8),
    legend.position = "none",   
    axis.text = element_text(size = 6)) +
  xlab("")+
  ylab("Log2FoldChange in Gametophyte Coverage")+
  theme(legend.position = "none")
Aux3

## Plot ####

# Custom colors with alpha for each trial
trial_colors <- c("1" = alpha("orange", 0.5), 
                  "2" = alpha("blue", 0.5), 
                  "3" = alpha("#4DAF4A", 0.5), 
                  "4" = alpha("red", 0.5))

## Sporo
plot1 <-ggplot(Merged_data_clean, aes(x = stat, y = Log2.Fold.Change., color = as.factor(Co.culture.Trial))) +
  geom_point() +
  theme_minimal() +
  theme(aspect.ratio = 1/1)+
  labs(color = "Trial", x = NULL) +
  ylab("Log2FoldChange in Number of Sporophytes") +
  geom_smooth(aes(group = Co.culture.Trial), method = lm, se = FALSE, size = .75) +
  geom_smooth(method = lm, color = alpha("black", 0.9), se = FALSE, size = 1.5) +
  scale_color_manual(values = trial_colors)  
plot1

## Gameto
plot2 <- ggplot(Merged_data_clean, aes(x = stat, y = Gametophyte_fold_change, color = as.factor(Co.culture.Trial))) +
  geom_point() +
  theme_minimal() +
  theme(aspect.ratio = 1/1)+
  labs(color = NULL) +
  xlab("Strength of Association with S. latissima (IndVal)") +
  ylab("Log2FoldChange in Gametophyte Coverage") +
  geom_smooth(aes(group = Co.culture.Trial), method = lm, se = FALSE, size = .75) +
  geom_smooth(method = lm, color = alpha("black", 0.9), se = FALSE, size = 1.5) +
  scale_color_manual(values = trial_colors, guide = "none")  
plot2

## combining plots
library(cowplot)

legend <- get_legend(plot1)

# Combine plots without their legends
combined_plots <- plot_grid(plot1 + theme(legend.position = "none"), 
                            plot2 + theme(legend.position = "none"), 
                            ncol = 1, align = "v", labels = c("A", "B"), label_size = 10)


# Add the legend to the right of the combined plots
final_plot <- plot_grid(combined_plots, legend, rel_widths = c(3, 0.5))
final_plot
ggsave(file = "/Users/evankohn/Desktop/stat_cocultureFig_NEW.JPEG", plot = final_plot, dpi = 800, units = "in", width = 6, height = 8)

## Aux plots
combined_Aux_plots <- plot_grid(Aux1, Aux2, Aux3, ncol = 3, labels = c("A", "B", "C"), label_size = 10)
combined_Aux_plots

ggsave(file = "/Users/evankohn/Desktop/IAAplot1.JPEG", plot = combined_Aux_plots, dpi = 800, units = "in", width = 4, height = 3)

