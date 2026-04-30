
# prelim group diffs

library(tidyverse)
library(ggpubr)
library(broom)
library(effectsize)


rie <- read.csv("rawdata/rmp_and_epsc.csv")
rie$subslice <- paste(rie$subject, rie$slice, sep = "")
rie <- rie %>% relocate(subslice)

#### RMP

rmp.tp <- rie %>% 
  select(1,4,7:11) %>% 
  pivot_longer(col = c(4:7), names_to = "RMP.timepoint", values_to = "RMP")

rmp.tp %>%
  ggplot(aes(RMP.timepoint,RMP, fill = sex, color = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="Resting Membrane Potential",x="Group", y = "RMP (mV)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

## NEED TO DO LM AND NORMAL T TEST -- ESP @ JUST START

rmp.start <- rmp.tp %>% filter(RMP.timepoint == "RMP.start") %>% select(!RMP.timepoint)

rmp.lm.int <- lm(RMP ~ sex*group, data = rmp.start)
summary(rmp.lm.int)

rmp.lm.main <- lm(RMP ~ sex + group, data = rmp.start)
summary(rmp.lm.main)


rmp.tp %>%
  group_by(RMP.timepoint, sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(RMP[group == unique(group)[1]]),
    sd1 = sd(RMP[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(RMP[group == unique(group)[2]]),
    sd2 = sd(RMP[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(RMP ~ group)),
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
## ok - will go back and add marker to graph?




#### input res (10mV minus -10mV div by 20)

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


in.re.tp2 %>%
  ggplot(aes(Inp.re.timepoint,Inp.re, fill = sex, color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="Input Resistance",x="Group", y = "Input Resistance (MΩ)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

# all well under 20, looks good



in.re.tp2.start <- in.re.tp2 %>% filter(Inp.re.timepoint == "in.re.start") %>% select(!Inp.re.timepoint)

inre.lm.int <- lm(Inp.re ~ sex*group, data = in.re.tp2.start)
summary(inre.lm.int)

inre.lm.main <- lm(Inp.re ~ sex + group, data = in.re.tp2.start)
summary(inre.lm.main)



in.re.tp2 %>%
  na.omit() %>% 
  group_by(Inp.re.timepoint, sex) %>%
  summarise(
    # Group stats
    group1 = unique(group)[1],
    mean1 = mean(Inp.re[group == unique(group)[1]]),
    sd1 = sd(Inp.re[group == unique(group)[1]]),
    n1 = sum(group == unique(group)[1]),
    
    group2 = unique(group)[2],
    mean2 = mean(Inp.re[group == unique(group)[2]]),
    sd2 = sd(Inp.re[group == unique(group)[2]]),
    n2 = sum(group == unique(group)[2]),
    
    # T-test
    t_test = list(t.test(Inp.re ~ group)),
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
## ok - will go back and add marker to graph?








##### EPSCs (save as 350x450)

##event number
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <1000) %>% 
  ggplot(aes(group,epsc.events, fill = sex, color = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC Events",x="Group", y = "Number of Events") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 
# excluding the 3 cells that had 1000+ events

epsc.events.lm.int <- lm(epsc.events ~ sex*group, data = rie)
summary(epsc.events.lm.int)

epsc.events.lm.main <- lm(epsc.events ~ sex + group, data = rie)
summary(epsc.events.lm.main)
# no diffs, may have to come back later to do piecewise t.tests




##epsc amplitude
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <1000) %>% 
  ggplot(aes(group,epsc.amp, fill = sex, color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC Amplitude",x="Group", y = "Current (pA)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

epsc.amp.lm.int <- lm(epsc.amp ~ sex*group, data = rie)
summary(epsc.amp.lm.int)

epsc.amp.lm.main <- lm(epsc.amp ~ sex + group, data = rie)
summary(epsc.amp.lm.main)
# no diffs, may have to come back later to do piecewise t.tests




##epsc AUC pAms
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <1000) %>% 
  ggplot(aes(group,epsc.auc.pAms, fill = sex, color=group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC AUC",x="Group", y = "AUC (pA-ms)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

epsc.auc.lm.int <- lm(epsc.auc.pAms ~ sex*group, data = rie)
summary(epsc.auc.lm.int)

epsc.auc.lm.main <- lm(epsc.auc.pAms ~ sex + group, data = rie)
summary(epsc.auc.lm.main)
# no diffs, may have to come back later to do piecewise t.tests




