library(dplyr)
library(tidyverse)
library(readxl)
library('haven')
library(ggplot2)
library(tibble)
library(magrittr)
library(knitr)
library(ggpubr)
library(ggsignif)
library(ggpmisc)
library(interactions)
library(effects)
library(jtools)
library(patchwork)
library(lme4)
library(lmerTest)
library(sjPlot)

#######

## Load and organize data
dataset = read_excel('/Users/noahschweitzer/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Microglia_project/final_combined_results20250722.xlsx')


FileName = unique(dataset$FileName)
files = c("Mouse0507","Mouse0506","Mouse0917","Mouse0924","Mouse0926","Mouse1003",
          "Mouse1031","Mouse1029","Mouse1024","Mouse1022","Mouse1017","Mouse1105",
          "Mouse1107","Mouse1119", "Mouse1203", "Mouse1205", "Mouse20250325","Mouse20250327")

dataset$Age <- ifelse(dataset$FileName %in% files, "old", "young")

dataset <- dataset %>%
  mutate(APOE_Sex = factor(case_when(
    APOE == "APOE3" & Sex== "Female" ~ 'APOE3 Female',
    APOE == "APOE3" & Sex== "Male" ~ 'APOE3 Male',
    APOE == "APOE4" & Sex== "Female" ~ 'APOE4 Female',
    APOE == "APOE4" & Sex== "Male" ~ 'APOE4 Male',
  )))




dataset = subset(dataset, MinBranchLength>0 & MaxBranchLength>1)
dataset$t_index = 1+((dataset$MinBranchLength - dataset$MaxBranchLength)/(dataset$MinBranchLength+ dataset$MaxBranchLength))


dataset <- dataset %>%
  mutate(log_MaxBranchLength = log(MaxBranchLength),
         log_AvgBranchLength = log(AvgBranchLength))  # You can add 1 if needed to avoid log(0)

# Compute means grouped by FileName and Time
dataset_mean <- dataset %>%
  group_by(FileName, Time) %>%
  summarise(mean_FullCellComplexity = mean(FullCellComplexity, na.rm = TRUE),
            mean_log_MaxBranchLength = mean(log_MaxBranchLength, na.rm = TRUE),
            mean_log_AvgBranchLength = mean(log_AvgBranchLength, na.rm = TRUE),
            mean_t_index = mean(t_index, na.rm = TRUE),
            .groups = "drop")

# Reattach categorical information
summary_vars <- dataset %>%
  select(FileName, Time, APOE, Sex, APOE_Sex, Injection, Age) %>%
  distinct()

# Join to create full summary dataframe
mean_values_by_time <- left_join(summary_vars, dataset_mean, by = c("FileName", "Time"))

attach(dataset)

mean_values_by_time <- mean_values_by_time %>%
  mutate(APOE_Sex_Age = factor(case_when(
    APOE_Sex == 'APOE3 Female' & Age == "young" ~ "APOE3F\nYoung",
    APOE_Sex == 'APOE3 Male' & Age == "young" ~ "APOE3M\nYoung",
    APOE_Sex == 'APOE4 Female' & Age == "young" ~ "APOE4F\nYoung",
    APOE_Sex == 'APOE4 Male' & Age == "young" ~ "APOE4M\nYoung",
    APOE_Sex == 'APOE3 Female' & Age == "old" ~ "APOE3F\nOld",
    APOE_Sex == 'APOE3 Male' & Age == "old" ~ "APOE3M\nOld",
    APOE_Sex == 'APOE4 Female' & Age == "old" ~ "APOE4F\nOld",
    APOE_Sex == 'APOE4 Male' & Age == "old" ~ "APOE4M\nOld",
  )))


###########
### Pre injection quantifications

ggplot(subset(mean_values_by_time,Injection=="PreInjection" & Age=="young"),aes(x=Time,y=mean_log_MaxBranchLength,color=Sex,group=FileName))+
  facet_grid(~APOE)+geom_point()+geom_smooth(linewidth=2,span=2)+
  ylab("log Max Branch Length (pixels)")+xlab("Time (Frame #)")+  scale_color_manual(values=c("#F8766D","#619CFF"))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20),
        strip.text = element_text(size = 18))+ggtitle("Young Mice")

ggplot(subset(mean_values_by_time,Injection=="PreInjection" & Age=="old"),aes(x=Time,y=mean_log_MaxBranchLength,color=Sex,group=FileName))+
  facet_grid(~APOE)+geom_point()+geom_smooth(linewidth=2,span=2)+
  ylab("log Max Branch Length (pixels)")+xlab("Time (Frame #)")+  scale_color_manual(values=c("#F8766D","#619CFF"))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20),
        strip.text = element_text(size = 18))+ggtitle("Old Mice")

ggplot(subset(mean_values_by_time,Injection=="PreInjection" & Age=="young"),aes(x=Time,y=mean_FullCellComplexity,color=Sex,group=FileName))+
  facet_grid(~APOE)+geom_point()+geom_smooth(linewidth=2,span=2)+
  ylab("Full Cell Complexity")+xlab("Time (Frame #)")+  scale_color_manual(values=c("#F8766D","#619CFF"))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20),
        strip.text = element_text(size = 18))+ggtitle("Young Mice")

ggplot(subset(mean_values_by_time,Injection=="PreInjection" & Age=="old"),aes(x=Time,y=mean_FullCellComplexity,color=Sex,group=FileName))+
  facet_grid(~APOE)+geom_point()+geom_smooth(linewidth=2,span=2)+
  ylab("Full Cell Complexity")+xlab("Time (Frame #)")+  scale_color_manual(values=c("#F8766D","#619CFF"))+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20),
        strip.text = element_text(size = 18))+ggtitle("Old Mice")




my_comparisons <- list( c("APOE3 Female","APOE3 Male"),
                        c("APOE3 Female","APOE4 Female"))

ggplot(subset(mean_values_by_time,  Age=="young" & Injection=="PreInjection"), aes(x=APOE_Sex,y=mean_FullCellComplexity,fill=Sex))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.38, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=6)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.03)+  scale_fill_manual(values=c("#F8766D","#619CFF"))+
  ylab("Full Cell Complexity")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20))+ggtitle("Young Mice")


my_comparisons <- list( c("APOE3 Female","APOE3 Male"),c("APOE4 Female","APOE4 Male"),
                        c("APOE3 Female","APOE4 Female"),c("APOE3 Male","APOE4 Male") )
ggplot(subset(mean_values_by_time,  Age=="old" & Injection=="PreInjection"), aes(x=APOE_Sex,y=mean_FullCellComplexity,fill=Sex))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.38, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=6)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.03)+  scale_fill_manual(values=c("#F8766D","#619CFF"))+
  ylab("Full Cell Complexity")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20))+ggtitle("Old Mice")



ggplot(subset(mean_values_by_time,  Age=="old" & Injection=="PreInjection"), aes(x=APOE_Sex,y=mean_log_MaxBranchLength,fill=Sex))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.38, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=6)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.03)+scale_fill_manual(values=c("#F8766D","#619CFF"))+
  ylab("log Max Branch Length")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20))




ggplot(subset(mean_values_by_time,  Age=="old" & Injection=="PreInjection"), aes(x=APOE_Sex,y=mean_t_index,fill=Sex))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.38, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=6)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.003)+scale_fill_manual(values=c("#F8766D","#619CFF"))+
  ylab("T-Index")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20))




my_comparisons <- list( c("APOE3F\nYoung","APOE3F\nOld"),c("APOE3M\nYoung","APOE3M\nOld"),
                        c("APOE4F\nYoung","APOE4F\nOld"))


mean_values_by_time <- mean_values_by_time %>% mutate(APOE_Sex_Age = factor(APOE_Sex_Age,levels=c("APOE3F\nYoung","APOE3F\nOld",
                                                                      "APOE3M\nYoung","APOE3M\nOld","APOE4F\nYoung","APOE4F\nOld","APOE4M\nYoung","APOE4M\nOld")))


ggplot(subset(mean_values_by_time, Injection=="PreInjection"), aes(x=APOE_Sex_Age,y=mean_FullCellComplexity,fill=Sex))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.38, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=6)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.003)+scale_fill_manual(values=c("#F8766D","#619CFF"))+
  ylab("Full Cell Complexity")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_blank(),
        title = element_text(size = 20))

##############

### Post injection quantifications


x = subset(mean_values_by_time, Injection=="PostInjection")
uniquelist = unique(x$FileName)
experiment <- mean_values_by_time %>% filter(FileName %in% uniquelist)


file_counts <- subset(experiment,Age=="old") %>%
  group_by(APOE_Sex) %>%
  summarise(NumUniqueFiles = n_distinct(FileName), .groups = "drop")
print(file_counts)

# Create a new column AdjustedTime initialized with NA
experiment$AdjustedTime <- NA

# Loop through each unique FileName
for (file in unique(experiment$FileName)) {
  # Subset data for the current FileName
  file_data <- experiment[experiment$FileName == file, ]
  
  # Find the minimum Time where Injection is "PostInjection"
  min_time <- min(file_data$Time[file_data$Injection == "PostInjection"], na.rm = TRUE)
  
  # Adjust Time and store it in AdjustedTime
  experiment$AdjustedTime[experiment$FileName == file] <- 
    experiment$Time[experiment$FileName == file] - min_time
}


experiment <- experiment %>%
  mutate(AdjustedTime_binned = factor(case_when(
    AdjustedTime<=0 ~ 'Baseline',
    AdjustedTime>0 & AdjustedTime<=20 ~  "0-20 min",
    AdjustedTime>20 & AdjustedTime<=40 ~  "20-40 min",
    AdjustedTime>40 ~ "40-80 min",
  )))


experiment <- experiment %>% mutate(AdjustedTime_binned = factor(AdjustedTime_binned,levels=c("Baseline","0-20 min",
                                                                                              "20-40 min","40-80 min")))

my_comparisons <- list( c("Baseline","0-20 min"), c("Baseline","20-40 min"), c("Baseline","40-80 min"))



experiment <- experiment %>%
  mutate(APOE_Sex = factor(case_when(
    APOE == "APOE3" & Sex== "Female" ~ 'APOE3\nFemale',
    APOE == "APOE3" & Sex== "Male" ~ 'APOE3\nMale',
    APOE == "APOE4" & Sex== "Female" ~ 'APOE4\nFemale',
    APOE == "APOE4" & Sex== "Male" ~ 'APOE4\nMale',
  )))


# Step 1: Get baseline value per FileName (summary-level)
baseline_df <- experiment %>%
  filter(AdjustedTime_binned == "Baseline") %>%
  distinct(FileName, .keep_all = TRUE) %>%
  select(FileName, baseline_t_index = mean_t_index)

# Step 2: Join and compute percent change
experiment <- experiment %>%
  left_join(baseline_df, by = "FileName") %>%
  mutate(t_index_pct_change = 100 * (mean_t_index - baseline_t_index) / baseline_t_index)


# Step 1: Get baseline value per FileName (summary-level)
baseline_df <- experiment %>%
  filter(AdjustedTime_binned == "Baseline") %>%
  distinct(FileName, .keep_all = TRUE) %>%
  select(FileName, baseline_FullCellComplexity = mean_FullCellComplexity)

# Step 2: Join and compute percent change
experiment <- experiment %>%
  left_join(baseline_df, by = "FileName") %>%
  mutate(FullCellComplexity_pct_change = 100 * (mean_FullCellComplexity - baseline_FullCellComplexity) / baseline_FullCellComplexity)



# Step 1: Get baseline value per FileName (summary-level)
baseline_df <- experiment %>%
  filter(AdjustedTime_binned == "Baseline") %>%
  distinct(FileName, .keep_all = TRUE) %>%
  select(FileName, baseline_MaxBranchLength = mean_log_MaxBranchLength)

# Step 2: Join and compute percent change
experiment <- experiment %>%
  left_join(baseline_df, by = "FileName") %>%
  mutate(MaxBranchLength_pct_change = 100 * (mean_log_MaxBranchLength - baseline_MaxBranchLength) / baseline_MaxBranchLength)


ggplot(subset(experiment, Age=="old"), aes(x=AdjustedTime_binned,y=t_index_pct_change))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.4, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 1.5)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )+ylab("T-Index % change")

ggplot(subset(experiment, Age=="old"), aes(x=AdjustedTime_binned,y=FullCellComplexity_pct_change))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.2, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 1.5)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )+ylab("Full Cell Complexity % change")

ggplot(subset(experiment, Age=="old"), aes(x=AdjustedTime_binned,y=MaxBranchLength_pct_change))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.2, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.75)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )+ylab("Log Max Branch Length % change")


ggplot(subset(experiment, Age=="old"), aes(x=AdjustedTime_binned,y=log_MaxBranchLength))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.4, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 0.005)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )



ggplot(subset(experiment, Age=="young"), aes(x=AdjustedTime_binned,y=FullCellComplexity_pct_change))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.2, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 1.5)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )+ylab("Full Cell Complexity % change")+ggtitle("Young Mice")


ggplot(subset(experiment, Age=="old"), aes(x=AdjustedTime_binned,y=FullCellComplexity_pct_change))+geom_violin(position='dodge')+
  stat_summary(fun = "mean", geom = "crossbar",width= 0.2, position = position_dodge(width = 0.9),show.legend = FALSE)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",size=5)+
  geom_dotplot(binaxis='y',position="dodge", stackdir='center', binwidth = 1.5)+facet_grid(APOE_Sex ~ ., switch = "y") +  # <- Moves strip text to right
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title = element_text(size = 22),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_blank(),
    title = element_text(size = 20),
    strip.text.y.left = element_text(size = 18,angle=0),  # optional: hide default left strip
    strip.placement = "outside" # <- Ensures strip is placed on right
  )+ylab("Full Cell Complexity % change")+ggtitle("Old Mice")

