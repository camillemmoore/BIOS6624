#############################################################
## Program Name: Descriptives.R						                 ##
## Purpose: To read in the clean dataset and create graphs ##
##     		 and descriptive statistics				               ##
## Created by: Nichole Carlson	                           ##
## (R code by: Kevin Josey and updated by Camille Moore)	 ##
#############################################################

# Dependencies
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(gt) 
library(gtsummary) # this is a package to make tables of descriptive stats
library(flextable) # this can help print tables to docx (MS Word)
library(tableone) # this is another fairly easy to use descriptive stats package
library(psych) # used for corr.test
library(knitr) # for writing tables
library(kableExtra) # for writing tables

###You will have to set your working directory####
setwd("/Users/mooreca/Documents/BIOS6624/PSAExample/")

### read in the first cleaned data set (psaclean.csv)

## read in the clean data from our data directory
psaclean <- read.csv("DataProcessed/psaclean.csv")
  
## create scatter plots with PSA as the outcome
p1 <- ggplot(psaclean, aes(x=age, y=psa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Age (years)')+ylab('PSA (mg/ml)')+
  ggtitle('PSA vs. Age')+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(psaclean, aes(x=wt, y=psa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Prostate Weight (g)')+ylab('PSA (mg/ml)')+
  ggtitle('PSA vs. Prostate Weight')+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(psaclean, aes(x=cavol, y=psa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Cancer Volume (cc)')+ylab('PSA (mg/ml)')+
  ggtitle('PSA vs. Cancer Volume')+
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(psaclean, aes(x=bph, y=psa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Benign prostatic hyperplasia/cm2')+ylab('PSA (mg/ml)')+
  ggtitle('PSA vs. Benign prostatic hyperplasia')+
  theme(plot.title = element_text(hjust = 0.5))

p5 <- ggplot(psaclean, aes(x=cappen, y=psa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Capsular penetration (cm)')+ylab('PSA (mg/ml)')+
  ggtitle('PSA vs. Capsular penetration')+
  theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(psaclean, aes(x=factor(gleason), y=psa))+
  geom_boxplot(outlier.size = 0)+
  geom_beeswarm()+
  theme_classic()+
  xlab('Gleason Score')+ylab('PSA (mg/ml)')+
  ggtitle('PSA by Gleason Score')+
  theme(plot.title = element_text(hjust = 0.5))

p7 <- ggplot(psaclean, aes(x=svi_factor, y=psa))+
  geom_boxplot(outlier.size = 0)+
  geom_beeswarm()+
  theme_classic()+
  xlab('Seminal Vesicle Invasion')+ylab('PSA (mg/ml)')+
  ggtitle('PSA by SVI')+
  theme(plot.title = element_text(hjust = 0.5))

pdf('Output/WS3_plots_rawscale.pdf',
    height=10, width=7)
ggarrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2, nrow=4)
dev.off()


## Table of descriptive statistics
### There are other Table 1 packages besides gtsummary.  Use one you like.
### Can print to screen to see things #####
### Write to a file - options for word(docx) using flextable package, 
### html using gt package and csv.  

### This table has a fairly comprehensive list of summary stats
### May be overkill for a final report, but good for data exploration
tbl <- tbl_summary(
  data = psaclean,
  
  include = c(psa, lpsa, cavol, wt, age, bph, cappen, svi, gleason),
  
  # Specify types
  type = list(
    c(svi, gleason) ~ "categorical",
    all_continuous() ~ "continuous2"
  ),
  
  # Multi-line statistics for continuous, n (%) for categorical
  statistic = list(
    all_continuous() ~ c(
      "{mean} ({sd})",             # Mean (SD)
      "{median} ({p25}, {p75})",  # Median (Q1, Q3)
      "{min}-{max}"                      # Min, Max
    ),
    all_categorical() ~ "{n} ({p}%)"  # Count (%) for categorical
  ),
  
  # Show missing values if any
  missing = "ifany"
)

tbl

### save to word using flextable package
ft <- as_flex_table(tbl)
flextable::save_as_docx(ft, path = "Output/psaclean_descriptive_stats.docx")

### save to .csv
tbl_df <- as_tibble(tbl)
write.csv(tbl_df, "Output/psaclean_descriptive_stats.csv", row.names = FALSE)

### save to html using gt package
gt_tbl <- as_gt(tbl)
gt::gtsave(gt_tbl, "Output/psaclean_descriptive_stats.html")


## read in the final clean data from our data directory after addressing findings in initial descriptive statistics
psaclean2 <- read.csv("DataProcessed/psaclean2.csv")
  
## create scatter plots with log(PSA) as the outcome (log transforming the outcome => approximately normal)

p1 <- ggplot(psaclean2, aes(x=age, y=lpsa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Age (years)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Age')+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(psaclean2, aes(x=wt, y=lpsa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Prostate Weight (g)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Prostate Weight')+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(psaclean2, aes(x=cavol, y=lpsa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Cancer Volume (cc)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Cancer Volume')+
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(psaclean2, aes(x=bph, y=lpsa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Benign prostatic hyperplasia/cm2')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Benign prostatic hyperplasia')+
  theme(plot.title = element_text(hjust = 0.5))

p5 <- ggplot(psaclean2, aes(x=cappen, y=lpsa))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Capsular penetration (cm)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Capsular penetration')+
  theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(psaclean2, aes(x=factor(gleason), y=lpsa))+
  geom_boxplot(outlier.size = 0)+
  geom_beeswarm()+
  theme_classic()+
  xlab('Gleason Score')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) by Gleason Score')+
  theme(plot.title = element_text(hjust = 0.5))


p7 <- ggplot(psaclean2, aes(x=svi_factor, y=lpsa))+
  geom_boxplot(outlier.size = 0)+
  geom_beeswarm()+
  theme_classic()+
  xlab('Seminal Vesicle Invasion')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) by SVI')+
  theme(plot.title = element_text(hjust = 0.5))


pdf('Output/WS3_plots_logscale.pdf',
    height=10, width=7)
ggarrange(p1, p2, p3, p4, p5, p6, p7, ncol = 2, nrow=4)
dev.off()


## Re-do descriptive stats tables now that outlier is removed
tbl2 <- tbl_summary(
  data = psaclean2,
  
  include = c(psa, lpsa, cavol, wt, age, bph, cappen, svi, gleason),
  
  # Specify types
  type = list(
    c(svi, gleason) ~ "categorical",
    all_continuous() ~ "continuous2"
  ),
  
  # Multi-line statistics for continuous, n (%) for categorical
  statistic = list(
    all_continuous() ~ c(
      "{mean} ({sd})",             # Mean (SD)
      "{median} ({p25}, {p75})",  # Median (Q1, Q3)
      "{min}-{max}"                      # Min, Max
    ),
    all_categorical() ~ "{n} ({p}%)"  # Count (%) for categorical
  ),
  
  # Show missing values if any
  missing = "ifany"
)

tbl2

### save to word using flextable package
ft2 <- as_flex_table(tbl2)
flextable::save_as_docx(ft2, path = "Output/psaclean2_descriptive_stats.docx")

### save to .csv
tbl_df2 <- as_tibble(tbl2)
write.csv(tbl_df2, "Output/psaclean2_descriptive_stats.csv", row.names = FALSE)

### save to html using gt package
gt_tbl2 <- as_gt(tbl2)
gt::gtsave(gt_tbl2, "Output/psaclean2_descriptive_stats.html")

# Make pearson correlation coefficients
cor_mat <- psaclean2 %>%
  dplyr::select(psa, lpsa, cavol, wt, age, bph, cappen) %>%
  corr.test() # depends on library(psych)

#Write to a file for future use: Not pretty but useful.
sink("Output/CorrelationTablesClean.txt")
print(cor_mat)
sink()

## create scatter plots with log(PSA) as the outcome and an interaction for svi
p1 <- ggplot(psaclean2, aes(x=age, y=lpsa, group=svi_factor, color=svi_factor))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Age (years)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Age')+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(psaclean2, aes(x=wt, y=lpsa, group=svi_factor, color=svi_factor))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Prostate Weight (g)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Prostate Weight')+
  theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(psaclean2, aes(x=cavol, y=lpsa, group=svi_factor, color=svi_factor))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Cancer Volume (cc)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Cancer Volume')+
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(psaclean2, aes(x=bph, y=lpsa, group=svi_factor, color=svi_factor))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Benign prostatic hyperplasia/cm2')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Benign prostatic hyperplasia')+
  theme(plot.title = element_text(hjust = 0.5))

p5 <- ggplot(psaclean2, aes(x=cappen, y=lpsa, group=svi_factor, color=svi_factor))+
  geom_point()+
  geom_smooth(method='lm', se = F)+
  theme_classic()+
  xlab('Capsular penetration (cm)')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) vs. Capsular penetration')+
  theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(psaclean2, aes(x=factor(gleason), y=lpsa, fill=svi_factor))+
  geom_boxplot(outlier.size = 0)+
  geom_beeswarm()+
  theme_classic()+
  xlab('Gleason Score')+ylab('ln(PSA)')+
  ggtitle('ln(PSA) by Gleason Score')+
  theme(plot.title = element_text(hjust = 0.5))+facet_grid(.~svi_factor)


pdf('Output/WS5_plots_svi_interaction.pdf',
    height=10, width=7)
ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow=3, common.legend = T)
dev.off()


###Maybe make a table 1 by SVI.  We would need to put gleason as a facctor to add it.
###Uses the tableone R package
Table1bySVI <- CreateTableOne(vars = c("psa", "lpsa", "cavol", "wt", "age", "bph", "cappen","gleason"),
                              factorVars = 'gleason',
                              strata=c("svi_factor"),
                              #addOverall = T, # can be used to add the overall
                              data=psaclean2)
### Using Mean SD
tab_mean <- print(Table1bySVI, missing=T, showAllLevels=T, nonnormal = F, 
                  quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

### Using Median IQR
tab_median <- print(Table1bySVI, missing=T, showAllLevels=T, nonnormal = T, 
                    quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

### Save ouput as csv, word or html
write.csv(tab_mean, "Output/psaclean2_tableone_by_svi.csv", row.names = TRUE)
write.csv(tab_median, "Output/psaclean2_tableone_medians_by_svi.csv", row.names = TRUE)

### Example word
ft <- flextable(data.frame(Variable=rownames(tab_mean), as.data.frame(tab_mean)))
save_as_docx(ft, path = "Output/psaclean2_tableone_by_svi.docx")

### Example html
kbl(tab_mean, format = "html") %>%
  kableExtra::save_kable("Output/psaclean2_tableone_by_svi.html")
