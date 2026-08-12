

## first grubbs testing each metric to remove outliers before running tests 



library(tidyverse)
library(ggpubr)
library(broom)
library(effectsize)
library(outliers)


rie <- read.csv("rawdata/rmp_and_epsc.csv")
rie$subslice <- paste(rie$subject, rie$slice, sep = "")
rie <- rie %>% relocate(subslice)  %>% 
  mutate(sex = recode(sex, "M" = "Male", "F" = "Female"))

#### RMP

rmp.tp <- rie %>% 
  select(1,4,7:11) %>% 
  pivot_longer(col = c(4:7), names_to = "RMP.timepoint", values_to = "RMP")

rmp.tp %>%
  ggplot(aes(RMP.timepoint,RMP, fill = sex, color = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="Resting Membrane Potential",x="Group", y = "RMP (mV)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



in.re.tp <- rie %>% 
  select(1,4, 7:11, 32:39) %>% 
  mutate(
    in.re.start = (inp.start.10 - inp.start.n10) / 20,
    in.re.pre = (inp.pre.10 - inp.pre.n10) / 20,
    in.re.ttx = (inp.ttx.10 - inp.ttx.n10) / 20,
    in.re.xe = (inp.xe.10 - inp.xe.n10) / 20
  ) %>% 
  select(1:7, 16:19) %>%
  pivot_longer(col = c(8:11), names_to = "Inp.re.timepoint", values_to = "Inp.re")

in.re.tp2 <- in.re.tp %>%
  filter(!(
    (Inp.re.timepoint == "in.re.start" & RMP.start > -40) |
      (Inp.re.timepoint == "in.re.pre" & RMP.pre > -40) |
      (Inp.re.timepoint == "in.re.ttx" & RMP.ttx > -40) |
      (Inp.re.timepoint == "in.re.xe" & RMP.xe > -40)
  ))


inre.wide <- in.re.tp2 %>% select(1,8,9) %>% pivot_wider(names_from = Inp.re.timepoint, values_from = Inp.re)

riei <- merge(rie, inre.wide, by = "subslice")



riei2 <- riei %>% 
  select(subslice, subject, slice, sex, group, RMP.pre, in.re.pre, epsc.events, epsc.amp, epsc.auc.pAms, epsc.flag)


grubbs.test(type = 11, riei2$RMP.pre) #no outliers
grubbs.test(type = 11, riei2$in.re.pre) #no outliers
grubbs.test(type = 11, riei2$epsc.events) # remove the >2000s?
grubbs.test(type = 11, riei2$epsc.amp) #no outliers
grubbs.test(type = 11, riei2$epsc.auc.pAms) #remove c0850A?

gt.events.riei2 <- riei2 %>% filter(epsc.events < 1000)
grubbs.test(type = 11, gt.events.riei2$epsc.events) # remove the >1000s works fine
grubbs.test(type = 10, gt.events.riei2$epsc.events) # remove the >1000s works fine (there was still outliers in the one tailed if removing only >2000s)

gt.auc.riei2 <- riei2 %>% filter(subslice != "c0850A")
grubbs.test(type = 11, gt.auc.riei2$epsc.auc.pAms) #remove c0850A works fine


### RMP 
head(riei2)

rmp.pre <- riei2 %>% select(subslice, sex, group, RMP.pre)

rmp.lm.int <- lm(RMP.pre ~ sex*group, data = rmp.pre)
summary(rmp.lm.int)

rmp.lm.main <- lm(RMP.pre ~ sex + group, data = rmp.pre)
summary(rmp.lm.main)

rmp.pre %>% 
  group_by(sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(RMP.pre[group == unique(group)[1]]),
    sd1 = sd(RMP.pre[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(RMP.pre[group == unique(group)[2]]),
    sd2 = sd(RMP.pre[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(RMP.pre ~ group)),
    t_stat = t_test[[1]]$statistic,
    df = t_test[[1]]$parameter,
    p_value = t_test[[1]]$p.value,
    diff = mean2 - mean1,
    ci_low = t_test[[1]]$conf.int[1],
    ci_high = t_test[[1]]$conf.int[2],
    
    # Cohen's d (manual calculation)
    pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1 + n2 - 2)),
    cohens_d = diff/pooled_sd,
    
    .groups = "drop"
  ) %>%
  select(-t_test)

rmp.pre %>%
  ggplot(aes(sex,RMP.pre, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="A. Resting Membrane Potential",x="Sex", y = "RMP (mV)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) + 
  annotate("text", x = 1, y = -55, label = "**", size = 13) + 
  annotate("segment", x = 0.75, xend = 2.25, y = -51, color = "black", size = 1) + 
  annotate("text", x = 1.5, y = -52, label = "Interaction Effect", size = 5)+ 
  annotate("text", x = 1.5, y = -50.7, label = "**", size = 12)+
  scale_x_discrete(limits = c("Male", "Female"))


### INPUT RES
head(riei2)

inpre.pre <- riei2 %>% select(subslice, sex, group, in.re.pre)

inpre.pre.lm.int <- lm(in.re.pre ~ sex*group, data = inpre.pre)
summary(inpre.pre.lm.int)

inpre.pre.lm.main <- lm(in.re.pre ~ sex + group, data = inpre.pre)
summary(inpre.pre.lm.main)

inpre.pre %>% 
  group_by(sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(in.re.pre[group == unique(group)[1]]),
    sd1 = sd(in.re.pre[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(in.re.pre[group == unique(group)[2]]),
    sd2 = sd(in.re.pre[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(in.re.pre ~ group)),
    t_stat = t_test[[1]]$statistic,
    df = t_test[[1]]$parameter,
    p_value = t_test[[1]]$p.value,
    diff = mean2 - mean1,
    ci_low = t_test[[1]]$conf.int[1],
    ci_high = t_test[[1]]$conf.int[2],
    
    # Cohen's d (manual calculation)
    pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1 + n2 - 2)),
    cohens_d = diff/pooled_sd,
    
    .groups = "drop"
  ) %>%
  select(-t_test)

inpre.pre %>%
  ggplot(aes(sex,in.re.pre, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  labs(title="B. Input Resistance",x="Sex", y = "Input Resistance (MΩ)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) + 
  annotate("text", x = 1, y = 6, label = "**", size = 13)+ 
  annotate("segment", x = 0.75, xend = 2.25, y = 8, color = "black", size = 1) + 
  annotate("text", x = 1.5, y = 7.8, label = "Interaction Effect", size = 5)+ 
  annotate("text", x = 1.5, y = 8.05, label = "**", size = 12) +
  scale_x_discrete(limits = c("Male", "Female"))



### group diffs and figs for EPSC measures. 
## per grubbs tests - for events, filtering >1000. for amp, no outliers. for AUC removing c0850A

## events - 
riei2 %>% 
  select(subslice, sex, group, epsc.events) %>% 
  filter(epsc.events < 1000) %>% 
  group_by(sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(epsc.events[group == unique(group)[1]]),
    sd1 = sd(epsc.events[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(epsc.events[group == unique(group)[2]]),
    sd2 = sd(epsc.events[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(epsc.events ~ group)),
    t_stat = t_test[[1]]$statistic,
    df = t_test[[1]]$parameter,
    p_value = t_test[[1]]$p.value,
    diff = mean2 - mean1,
    ci_low = t_test[[1]]$conf.int[1],
    ci_high = t_test[[1]]$conf.int[2],
    
    # Cohen's d (manual calculation)
    pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1 + n2 - 2)),
    cohens_d = diff/pooled_sd,
    
    .groups = "drop"
  ) %>%
  select(-t_test)

events.riei2 <- riei2 %>% 
  select(subslice, sex, group, epsc.events) %>% 
  filter(epsc.events < 1000)

events.lm.int <- lm(epsc.events ~ sex*group, data = events.riei2)
summary(events.lm.int)

events.lm.main <- lm(epsc.events ~ sex + group, data = events.riei2)
summary(events.lm.main) # NO MAIN EFFECT EVEN IN FILTERED DATA

rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <1000) %>% 
  mutate(sex = recode(sex, "F" = "Female", "M" = "Male")) %>% 
  ggplot(aes(sex,epsc.events, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  labs(title="C. EPSC Events",x="Sex", y = "Number of Events") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_x_discrete(limits = c("Male", "Female"))



## amp - 
riei2 %>% 
  select(subslice, sex, group, epsc.amp) %>% 
  group_by(sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(epsc.amp[group == unique(group)[1]]),
    sd1 = sd(epsc.amp[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(epsc.amp[group == unique(group)[2]]),
    sd2 = sd(epsc.amp[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(epsc.amp ~ group)),
    t_stat = t_test[[1]]$statistic,
    df = t_test[[1]]$parameter,
    p_value = t_test[[1]]$p.value,
    diff = mean2 - mean1,
    ci_low = t_test[[1]]$conf.int[1],
    ci_high = t_test[[1]]$conf.int[2],
    
    # Cohen's d (manual calculation)
    pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1 + n2 - 2)),
    cohens_d = diff/pooled_sd,
    
    .groups = "drop"
  ) %>%
  select(-t_test)

riei2 %>% 
  ggplot(aes(sex,epsc.amp, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  labs(title="D. EPSC Amplitude",x="Sex", y = "Current (pA)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



## auc - 
riei2 %>% 
  select(subslice, sex, group, epsc.auc.pAms) %>% 
  filter(subslice != "c0850A") %>% 
  group_by(sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(epsc.auc.pAms[group == unique(group)[1]]),
    sd1 = sd(epsc.auc.pAms[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(epsc.auc.pAms[group == unique(group)[2]]),
    sd2 = sd(epsc.auc.pAms[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(epsc.auc.pAms ~ group)),
    t_stat = t_test[[1]]$statistic,
    df = t_test[[1]]$parameter,
    p_value = t_test[[1]]$p.value,
    diff = mean2 - mean1,
    ci_low = t_test[[1]]$conf.int[1],
    ci_high = t_test[[1]]$conf.int[2],
    
    # Cohen's d (manual calculation)
    pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2)/(n1 + n2 - 2)),
    cohens_d = diff/pooled_sd,
    
    .groups = "drop"
  ) %>%
  select(-t_test)



riei2 %>% 
  filter(subslice != "c0850A") %>% 
  ggplot(aes(sex,epsc.auc.pAms, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  labs(title="E. EPSC AUC",x="Sex", y = "AUC (pA-ms)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 
