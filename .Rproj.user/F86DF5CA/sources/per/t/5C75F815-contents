
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


#group*sex col
#mcpro.long$gxs <- paste0(mcpro.long$sex, "_",mcpro.long$group)

mcpro.long$sex <- as.factor(mcpro.long$sex)

#####################################################  Graphs

# mcpro.long %>% 
#   filter(drug == "TTX") %>% 
#   ggplot(., aes(mV, pA, color = group, linetype = sex)) +
#   geom_point() +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = T) +
#   theme_classic() +
#   ggtitle("Outward Current (TTX)")+ 
#   geom_hline(yintercept=0)


## TTX
mcpro.long %>% 
  filter(drug == "TTX") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", linewidth=1, ylab = "pA (+/- SEM)")+
  stat_compare_means(aes(group = group, label=..p.adj..), method="t.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c( "***", "**", "*", "")),
                     label.y= 40 , size=8) +
  ggtitle("Outward Current (TTX)")

mcpro.long %>% 
  filter(drug == "TTX") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="sex", color="sex", linewidth=1, ylab = "pA (+/- SEM)")+
  stat_compare_means(aes(group = group, label=..p.adj..), method="t.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c( "***", "**", "*", "")),
                     label.y= 40 , size=8) +
  ggtitle("Outward Current (TTX)")


mcpro.long %>%
  filter(drug == "TTX") %>%
  ggplot(aes(x = mV, y = pA, color = sex, linetype = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  labs(title = "Outward Current (TTX)", y = "pA (+/- SEM)") +
  theme_classic()



## XE
mcpro.long %>% 
  filter(drug == "XE") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (XE)")

mcpro.long %>% 
  filter(drug == "XE") %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="sex", color="sex", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (XE)")


mcpro.long %>%
  filter(drug == "XE") %>%
  ggplot(aes(x = mV, y = pA, color = sex, linetype = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  labs(title = "Outward Current (XE)", y = "pA (+/- SEM)") +
  theme_classic()



#### TTX-XE
mcpro.long_filtered %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (TTX-XE)")

mcpro.long_filtered %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="sex", color="sex", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (TTX-XE)")


mcpro.long_filtered %>%
  ggplot(aes(x = mV, y = pA, color = sex, linetype = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "Outward Current (TTX-XE)", y = "pA (+/- SEM)") +
  theme_classic()+
  theme(plot.title = element_text(size = 20)) 



### actual analysis -- no diffs, but still very low Ns
mc.diff.long_filtered2 <- mc.diff.long_filtered

mc.diff.long_filtered2$mvsq <- mc.diff.long_filtered2$mV^2



mc.dif.mm.fullint <- lmer(pA~group * sex * (mV+mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.fullint)


mc.dif.mm.int <- lmer(pA~group * sex + mV + mvsq + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.int)


mc.dif.mm.main <- lmer(pA~group + sex + mV + mvsq + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.main)

## Check for linearity and homoscedasticity
#plot(mc.dif.mm.main, type = c("p", "smooth"))
#
## Check for normality of residuals
#qqnorm(resid(mc.dif.mm.main))
#qqline(resid(mc.dif.mm.main))
##theyre *fine*

###### correlation time. extract -30mV M-current pA and attach to rie

rie2 <- mc.diff.long_filtered %>%
  filter(mV == -30) %>%
  rename(max_MC_pA = pA) %>%
  select(1,2,5) %>%
  left_join(rie,.)

write.csv(rie2, "rawdata/rie_and_maxMC.csv")



rie2 %>% 
  filter(RMP.xe < -40) %>% 
  filter(epsc.events <1000) %>% 
  ggplot(aes(group,max_MC_pA, fill = sex, color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="Max M-Current @ -30 mV",x="Group", y = "pA") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 
### i need to double check this is correct......

