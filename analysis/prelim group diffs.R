
# prelim group diffs

library(tidyverse)
library(ggpubr)

source("functions/geom_boxjitter.R")


rie <- read.csv("rawdata/rmp_and_epsc.csv")
rie$subslice <- paste(rie$subject, rie$slice, sep = "")
rie <- rie %>% relocate(subslice)

#### RMP

rmp.tp <- rie %>% 
  select(1,4,7:11) %>% 
  pivot_longer(col = c(4:7), names_to = "RMP.timepoint", values_to = "RMP")

rmp.tp %>%
  ggplot(aes(RMP.timepoint,RMP, fill = sex, color = group))+
  #stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
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
  #stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
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


##### EPSCs (save as 350x450)

##event number
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <3000) %>% 
  ggplot(aes(group,epsc.events, fill = sex))+
  #stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  #scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC Events",x="Group", y = "Number of Events") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 
# excluding the 1 cell that had 3000+ events


##epsc amplitude
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <3000) %>% 
  ggplot(aes(group,epsc.amp, fill = sex))+
  #stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  #scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC Amplitude",x="Group", y = "Current (pA)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



##epsc AUC pAms
rie %>% 
  filter(RMP.pre < -40) %>% 
  filter(epsc.events <3000) %>% 
  ggplot(aes(group,epsc.auc.pAms, fill = sex))+
  #stat_compare_means(method = "t.test", size = 6, label.x = 1.5) +
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, 
                 alpha = 1,
                 width = 0.5,
                 jitter.height = 0.02, jitter.width = 0.02, errorbar.draw = TRUE,
                 position = position_dodge(0.8)) +
  scale_fill_manual(values=c("hotpink", "skyblue")) +
  #scale_color_manual(values=c("gray", "black"))+
  labs(title="EPSC AUC",x="Group", y = "AUC (pA-ms)") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 



### still need to do all actual analysis on RMP, EPSC events, EPSC amp, and EPSC AUP
