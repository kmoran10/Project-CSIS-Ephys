
# M-current

library(Matrix)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(car)

source("functions/geom_boxjitter.R")



rie <- read.csv("rawdata/rmp_and_epsc.csv")
rie$subslice <- paste(rie$subject, rie$slice, sep = "")
rie <- rie %>% relocate(subslice)

mcd <- read.csv("rawdata/mcurrent.csv")

rie.id <- rie %>% select(subslice,group,sex)

tmcd <- as.data.frame(t(mcd))

names(tmcd) <- lapply(tmcd[1, ], as.character)
tmcd <- tmcd[-1,]
tmcd1 <- tibble::rownames_to_column(tmcd, "ID")



tmcd2 <- tmcd1 %>%
  separate(ID, into = c("mc.tp", "subject"), sep = "_", remove = TRUE) %>%
  select(subject, everything()) %>%
  separate(mc.tp, into = c("drug", "tp"), sep = "\\.", remove = TRUE) %>%
  select(subject, everything()) %>% 
  rename(subslice = subject)

mcd3 <- left_join(rie.id,tmcd2, by="subslice")

mcd4 <- mcd3[!is.na(mcd3$drug),]



mc.pro <- mcd4 %>%
  filter(tp %in% c("Y1", "Y2", "Y3", "Y4")) %>%
  mutate(
    drug_group = case_when(
      drug %in% c("TTX1", "TTX2") ~ "TTX",
      drug %in% c("XE1", "XE2") ~ "XE",
      TRUE ~ NA_character_
    ),
    group = group
  ) %>%
  pivot_longer(cols = 6:16) %>%
  # Remove NA values before calculations
  filter(!is.na(value)) %>%
  mutate(period = if_else(tp %in% c("Y1", "Y2"), "Y1Y2", "Y3Y4")) %>%
  group_by(subslice, group, sex, drug_group, period, name) %>%
  summarise(
    avg = mean(value, na.rm = TRUE),  # na.rm for additional protection
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = period, values_from = avg) %>%
  mutate(
    diff = if_else(
      is.na(Y1Y2) | is.na(Y3Y4),
      NA_real_,  # Return NA if either period is completely missing
      Y1Y2 - Y3Y4
    )
  ) %>%
  group_by(subslice, group, sex, drug_group, name) %>%
  summarise(
    avg_diff = mean(diff, na.rm = TRUE),  # Will ignore NA diffs
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = name, values_from = avg_diff) %>%
  rename(drug = drug_group)



mcpro.long <- mc.pro %>% 
  pivot_longer(col = c(5:15), names_to = "mV", values_to = "pA")
mcpro.long$mV <- as.numeric(mcpro.long$mV)
mcpro.long$group <- as.factor(mcpro.long$group)



mc.diff <- mc.pro %>%
  select(subslice, group, sex, drug, 5:15) %>%  # Select relevant columns
  pivot_longer(cols = -(1:4), names_to = "variable") %>%
  pivot_wider(names_from = drug, values_from = value) %>%
  mutate(MC = TTX - XE) %>%
  select(subslice, group, sex, variable, MC) %>%
  pivot_wider(names_from = variable, values_from = MC)


mc.diff.long <- mc.diff %>% 
  pivot_longer(col = c(4:14), names_to = "mV", values_to = "pA")
mc.diff.long$mV <- as.numeric(mc.diff.long$mV)
mc.diff.long$group <- as.factor(mc.diff.long$group)


# First, identify which subslice values to exclude
subslices_to_exclude <- rie %>%
  filter(RMP.xe > -40) %>%
  pull(subslice) %>%
  unique()

# Then filter both dataframes to exclude those subslices
mcpro.long_filtered <- mcpro.long %>%
  filter(!subslice %in% subslices_to_exclude)

mc.diff.long_filtered <- mc.diff.long %>%
  filter(!subslice %in% subslices_to_exclude)



####graphs

mcpro.long %>% 
  filter(drug == "TTX") %>% 
  ggplot(., aes(mV, pA, color = group, linetype = sex)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = T) +
  theme_classic() +
  ggtitle("Outward Current (TTX)")+ 
  geom_hline(yintercept=0)


mcpro.long %>% 
  filter(drug == "TTX") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", size=1, ylab = "pA (+/- SEM)")+
  stat_compare_means(aes(group = group, label=..p.adj..), method="t.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c( "***", "**", "*", "")),
                     label.y= 40 , size=8) +
  ggtitle("Outward Current (TTX)")
## not grouping by sex

mcpro.long %>% 
  filter(drug == "TTX") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", size=1, ylab = "pA (+/- SEM)")+
  # stat_compare_means(aes(group = group, label=..p.adj..), method="t.test", label = "p.signif", 
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
  #                                       symbols = c( "***", "**", "*", "")),
  #                    label.y= 40 , size=8) +
  ggtitle("Outward Current (TTX)")



