

# Male and D Female select behavior comparisons 


library(tidyverse)
library(ggpubr)
library(broom)
library(effectsize)


oft <- read.csv("rawdata/oft.csv")
epm <- read.csv("rawdata/epm.csv")
nsf <- read.csv("rawdata/nsf.csv")

a <- merge(oft, epm, by = c("subject","sex","group"), all = TRUE)
behs <- merge(a, nsf, by = c("subject","sex","group"), all = TRUE)


behs %>% 
  ggplot(aes(sex,oft_10cm_center_pct_time, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#55a0fd", "#f64ed5")) +
  labs(title="OFT 10cm center",x="Sex", y = "Seconds") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) +
  scale_x_discrete(limits = c("Male", "D. Female"))

center.int <- lm(oft_10cm_center_pct_time ~ sex*group, data = behs)
summary(center.int)

center.main <- lm(oft_10cm_center_pct_time ~ sex + group, data = behs)
summary(center.main)
## NOTHING



behs %>% 
  ggplot(aes(sex,oft_corners_pct_time, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#55a0fd", "#f64ed5")) +
  labs(title="OFT corners & time ",x="Sex", y = "& time") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) +
  scale_x_discrete(limits = c("Male", "D. Female"))

corner.int <- lm(oft_corners_pct_time ~ sex*group, data = behs)
summary(corner.int) ########## DEFINITE INTERACTION HERE

corner.main <- lm(oft_corners_pct_time ~ sex + group, data = behs)
summary(corner.main)




behs %>% 
  ggplot(aes(sex,epm_open_time, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#55a0fd", "#f64ed5")) +
  labs(title="EPM open time",x="Sex", y = "Seconds") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) +
  scale_x_discrete(limits = c("Male", "D. Female"))

epm.open.int <- lm(epm_open_time ~ sex*group, data = behs)
summary(epm.open.int)

epm.open.main <- lm(epm_open_time ~ sex + group, data = behs)
summary(epm.open.main)
## NOTHING



behs %>% 
  ggplot(aes(sex,novel_eat_latency, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#55a0fd", "#f64ed5")) +
  labs(title="latency to eat in novel arena",x="Sex", y = "Seconds") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) +
  scale_x_discrete(limits = c("Male", "D. Female"))

novel.lat.int <- lm(novel_eat_latency ~ sex*group, data = behs)
summary(novel.lat.int)

novel.lat.main <- lm(novel_eat_latency ~ sex + group, data = behs)
summary(novel.lat.main)
## BOTH A SEX AND STRESS MAIN EFFECT - NO INTERACTION. 




