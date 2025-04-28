library(ggplot2)
library(tidyverse)


## Sporophyte plot #####
CoCultureData <- read.csv("/Users/evankohn/Desktop/NSERC work/Chapter3_co-culture_results_Appendix3.Table4_lwp.csv")

# Adjust P-values 
CoCultureData <- CoCultureData %>%
  mutate(Trial1Pval = ifelse(Co.culture.Trial == "1", T.TEST.P.value, NA)) %>% 
  mutate(Trial2Pval = ifelse(Co.culture.Trial == "2", T.TEST.P.value, NA)) %>% 
  mutate(Trial3Pval = ifelse(Co.culture.Trial == "3", T.TEST.P.value, NA)) %>% 
  mutate(Trial4Pval = ifelse(Co.culture.Trial == "4", T.TEST.P.value, NA))

CoCultureData <- CoCultureData %>%
  mutate(AdjustedPval1 = ifelse(!is.na(Trial1Pval), p.adjust(Trial1Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval2 = ifelse(!is.na(Trial2Pval), p.adjust(Trial2Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval3 = ifelse(!is.na(Trial3Pval), p.adjust(Trial3Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval4 = ifelse(!is.na(Trial4Pval), p.adjust(Trial4Pval, method = "BH"), NA))

CoCultureDataNEWP <- CoCultureData %>%
  mutate(NewPVal = coalesce(AdjustedPval1, AdjustedPval2, AdjustedPval3, AdjustedPval4))

## plot
trial_colors <- c("1" = alpha("orange", 0.5), 
                  "2" = alpha("blue", 0.5), 
                  "3" = alpha("#4DAF4A", 0.5), 
                  "4" = alpha("red", 0.5))

sporophyte_plot <- ggplot(data = CoCultureDataNEWP, aes(x = Log2.Fold.Change., y = -log10(NewPVal), color = as.factor(Co.culture.Trial))) + 
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +  
  geom_point() +
  scale_y_log10() +
  xlab("Log2FoldChange in Sporophyte Number") +  
  ylab(NULL)+
  geom_hline(yintercept = -log10(0.05), col = "BLACK") +
  scale_color_manual(values = trial_colors) +
  geom_point(aes(color = ifelse(Log2.Fold.Change. >= -0.6 & Log2.Fold.Change. <= 0.6, 
                                "grey", as.factor(Co.culture.Trial)))) +
  labs(color = "Co-Culture Trial") +
  theme_minimal()
  


## Gametophyte Plot #### 

CoCultureData_G <- read.csv("/Users/evankohn/Desktop/NSERC work/GametoP_Chapter3_co-culture_results_Appendix3.Table4_lwp.csv")

## Adjust p-values
CoCultureData_G <- CoCultureData_G %>%
  mutate(Trial1Pval = ifelse(Co.culture.Trial == "1", Gameto_pVal, NA)) %>% 
  mutate(Trial2Pval = ifelse(Co.culture.Trial == "2", Gameto_pVal, NA)) %>% 
  mutate(Trial3Pval = ifelse(Co.culture.Trial == "3", Gameto_pVal, NA)) %>% 
  mutate(Trial4Pval = ifelse(Co.culture.Trial == "4", Gameto_pVal, NA))

CoCultureData_G <- CoCultureData_G %>%
  mutate(AdjustedPval1 = ifelse(!is.na(Trial1Pval), p.adjust(Trial1Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval2 = ifelse(!is.na(Trial2Pval), p.adjust(Trial2Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval3 = ifelse(!is.na(Trial3Pval), p.adjust(Trial3Pval, method = "BH"), NA)) %>% 
  mutate(AdjustedPval4 = ifelse(!is.na(Trial4Pval), p.adjust(Trial4Pval, method = "BH"), NA))

CoCultureDataNEWP_G <- CoCultureData_G %>%
  mutate(NewPVal = coalesce(AdjustedPval1, AdjustedPval2, AdjustedPval3, AdjustedPval4))

## plot

gametophyte_plot <- ggplot(data= CoCultureDataNEWP_G, aes(x=Gametophyte_fold_change, y=-log10(NewPVal),, color = as.factor(Co.culture.Trial))) + 
  annotate("rect", xmin = -0.6, xmax = 0.6, ymin = 0, ymax = Inf, fill = "grey", alpha = 0.5) +
  geom_point() +
  scale_y_log10()+
  theme_minimal() +
  xlab("Log2FoldChange in Gametophyte Coverage") +
  geom_hline(yintercept=-log10(0.05), col="BLACK")+
  scale_color_manual(values = trial_colors)  +
  geom_point(aes(color = ifelse(Gametophyte_fold_change >= -0.6 & Gametophyte_fold_change <= 0.6, 
                                "grey", as.factor(Co.culture.Trial)))) +
  theme(legend.position = "none")
gameto

## combine
library(cowplot)
volcanoplot <- plot_grid(gameto, sporo, ncol = 2, labels = c("A", "B"))
volcanoplot

ggsave(file = "/Users/evankohn/Desktop/volcano.JPEG", plot = volcanoplot, dpi = 800, units = "in", width = 8, height = 4)



