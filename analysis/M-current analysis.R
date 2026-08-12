
# M-current

library(tidyverse)



rie <- read.csv("rawdata/rmp_and_epsc.csv")
rie$subslice <- paste(rie$subject, rie$slice, sep = "")
rie <- rie %>% relocate(subslice)

rie <- rie %>% mutate(sex = recode(sex, "F" = "Female", "M" = "Male"))  
  

library(Matrix)
library(lme4)
library(lmerTest)
library(ggpubr)
library(car)
library(emmeans)

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


# i have a gut feeling i did something wrong here maybe swapping some values so i am going to redo this M-current calculation code
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
mc.diff.long_filtered %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="group", color="group", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (TTX-XE)")

mc.diff.long_filtered %>% 
  ggline(., x = "mV", y = "pA", add = "mean_se", group="sex", color="sex", linewidth=1, ylab = "pA (+/- SEM)")+
  ggtitle("Outward Current (TTX-XE)")

# 800x600
mc.diff.long_filtered %>%
  ggplot(aes(x = mV, y = pA, color = sex, linetype = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "A. Outward Current (TTX-XE)", y = "pA (+/- SEM)") +
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) 


######################################################################### inserting grubbs testing and grubbs tested filtered analysis here, as well as removing -25mV timepoint. 

grubbs.test(type = 11, mc.diff$'-25') #none
grubbs.test(type = 11, mc.diff$'-30') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-35') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-40') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-45') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-50') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-55') #some - fix - c0880B
grubbs.test(type = 11, mc.diff$'-60') #none
grubbs.test(type = 11, mc.diff$'-65') #none
grubbs.test(type = 11, mc.diff$'-70') #none
grubbs.test(type = 11, mc.diff$'-75') #some - fix - c0824A or c0818C

# first just removing c0880B -- this IS the cell that I thought might be too posterior anyway

mc.diff.f1 <- mc.diff %>% filter(subslice != "c0880B")
mcpro.long.f1 <- mcpro.long %>% filter(subslice != "c0880B")

grubbs.test(type = 11, mc.diff.f1$'-25') #some now? gonna remove entire pool anyway so dw
grubbs.test(type = 10, mc.diff.f1$'-30') #some - barely - c0880A
grubbs.test(type = 10, mc.diff.f1$'-35') #some - barely - c0830B
grubbs.test(type = 11, mc.diff.f1$'-40') #none
grubbs.test(type = 11, mc.diff.f1$'-45') #none
grubbs.test(type = 11, mc.diff.f1$'-50') #some - barely - c0896B
grubbs.test(type = 11, mc.diff.f1$'-55') #none
grubbs.test(type = 11, mc.diff.f1$'-60') #none
grubbs.test(type = 11, mc.diff.f1$'-65') #none
grubbs.test(type = 11, mc.diff.f1$'-70') #none
grubbs.test(type = 11, mc.diff.f1$'-75') #some - some? - c0824A or c0818C

# ok - none of these occur in multiple places. Therefore, I'm just gonna keep all of them. 

mc.diff.f1.wide <- mc.diff.f1 %>% select(!("-25"))

mc.diff.f1.long <- mc.diff.f1 %>% 
  pivot_longer(col = c(4:14), names_to = "mV", values_to = "pA") %>% filter(mV != "-25")
mc.diff.f1.long$mV <- as.numeric(mc.diff.f1.long$mV)
mc.diff.f1.long$group <- as.factor(mc.diff.f1.long$group)


mc.diff.f1.long%>%
  ggplot(aes(x = mV, y = pA, color = group, linetype = sex)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "A. XE991 Sensitive Current (M-current)", y = "pA (+/- SEM)") +
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#55a0fd", "#f64ed5")) 


mc.diff.f1.long2 <- mc.diff.f1.long

mc.diff.f1.long2$mvsq <- mc.diff.f1.long2$mV^2

# ****
mc.dif.f1.mm.fullint <- lmer(pA~group * sex * poly(mV, 2) + (1 + poly(mV, 2)| subslice), data=mc.diff.f1.long2, na.action=na.exclude)
summary(mc.dif.f1.mm.fullint)
AIC(mc.dif.f1.mm.fullint) # lower AIC -good

mc.dif.f1.mm.fullint.nrs <- lmer(pA~group * sex * poly(mV, 2) + (1 | subslice), data=mc.diff.f1.long2, na.action=na.exclude)
summary(mc.dif.f1.mm.fullint.nrs)
AIC(mc.dif.f1.mm.fullint.nrs)


mc.dif.f1.mm.int1 <- lmer(pA~ group*sex + sex*poly(mV, 2) + group*poly(mV, 2) + (1 + mV| subslice), data=mc.diff.f1.long2, na.action=na.exclude)
summary(mc.dif.f1.mm.int1)


mc.dif.f1.mm.main <- lmer(pA~group + sex + poly(mV, 2) + (1 + poly(mV, 2)| subslice), data=mc.diff.f1.long2, na.action=na.exclude)
summary(mc.dif.f1.mm.main)



# males only
mc.diff.f1.long2.m <- mc.diff.f1.long2 %>% filter(sex=="Male")
mc.diff.f1.long2.f <- mc.diff.f1.long2 %>% filter(sex=="Female")

mc.dif.f1.mm.fullint.m <- lmer(pA~group * poly(mV, 2) + (1 + poly(mV, 2) | subslice), data=mc.diff.f1.long2.m, na.action=na.exclude)
summary(mc.dif.f1.mm.fullint.m)
AIC(mc.dif.f1.mm.fullint.m) # lower AIC -good

mc.dif.f1.mm.fullint.m.nrs <- lmer(pA~group * poly(mV, 2) + (1 | subslice), data=mc.diff.f1.long2.m, na.action=na.exclude)
summary(mc.dif.f1.mm.fullint.m.nrs)
AIC(mc.dif.f1.mm.fullint.m.nrs)


#*****
mc.dif.f1.mm.main.m <- lmer(pA~group + poly(mV, 2) + (1 + poly(mV, 2) | subslice), data=mc.diff.f1.long2.m, na.action=na.exclude)
summary(mc.dif.f1.mm.main.m) # main effect of group = [b = -4.45 ± 2.37, n = 37, p = 0.069]
AIC(mc.dif.f1.mm.main.m) # lower AIC -good
anova(mc.dif.f1.mm.main.m, type = 2)
emmeans(mc.dif.f1.mm.main.m, ~ group)

mc.dif.f1.mm.main.m.nrs <- lmer(pA~group + poly(mV, 2) + (1 | subslice), data=mc.diff.f1.long2.m, na.action=na.exclude)
summary(mc.dif.f1.mm.main.m.nrs)
AIC(mc.dif.f1.mm.main.m.nrs)



mc.diff.f1.long2.m %>%
  ggplot(aes(x = mV, y = pA, color = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "B. Males: XE991 Sensitive Current (M-current)", y = "pA (+/- SEM)") +
  theme_classic()+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#55a0fd", "#f64ed5")) 


## asking tommy and troy about how they have analyzed these data before... 
## saving m-c data to send to them 
write.csv(mc.diff.f1.wide, "rawdata/mc.diff.f1.wide.csv")
write.csv(mc.diff.f1.long, "rawdata/mc.diff.f1.long.csv")

wideforprism <- mcpro.long.f1 %>% pivot_wider(names_from = subslice, values_from = pA)

wideforprism.ttx.f.ctrl <- wideforprism %>% filter(drug == "TTX", sex == "Female", group == "ctrl") %>% select(where(~ !all(is.na(.))))
wideforprism.xe.f.ctrl <- wideforprism %>% filter(drug == "XE", sex == "Female", group == "ctrl") %>% select(where(~ !all(is.na(.))))
wideforprism.ttx.m.ctrl <- wideforprism %>% filter(drug == "TTX", sex == "Male", group == "ctrl") %>% select(where(~ !all(is.na(.))))
wideforprism.xe.m.ctrl <- wideforprism %>% filter(drug == "XE", sex == "Male", group == "ctrl") %>% select(where(~ !all(is.na(.))))

wideforprism.ttx.f.expt <- wideforprism %>% filter(drug == "TTX", sex == "Female", group == "expt") %>% select(where(~ !all(is.na(.))))
wideforprism.xe.f.expt <- wideforprism %>% filter(drug == "XE", sex == "Female", group == "expt") %>% select(where(~ !all(is.na(.))))
wideforprism.ttx.m.expt <- wideforprism %>% filter(drug == "TTX", sex == "Male", group == "expt") %>% select(where(~ !all(is.na(.))))
wideforprism.xe.m.expt <- wideforprism %>% filter(drug == "XE", sex == "Male", group == "expt") %>% select(where(~ !all(is.na(.))))


write.csv(wideforprism.ttx.f.ctrl, "rawdata/data for prism/wideforprism.ttx.f.ctrl.csv")
write.csv(wideforprism.xe.f.ctrl, "rawdata/data for prism/wideforprism.xe.f.ctrl.csv")
write.csv(wideforprism.ttx.m.ctrl, "rawdata/data for prism/wideforprism.ttx.m.ctrl.csv")
write.csv(wideforprism.xe.m.ctrl, "rawdata/data for prism/wideforprism.xe.m.ctrl.csv")
write.csv(wideforprism.ttx.f.expt, "rawdata/data for prism/wideforprism.ttx.f.expt.csv")
write.csv(wideforprism.xe.f.expt, "rawdata/data for prism/wideforprism.xe.f.expt.csv")
write.csv(wideforprism.ttx.m.expt, "rawdata/data for prism/wideforprism.ttx.m.expt.csv")
write.csv(wideforprism.xe.m.expt, "rawdata/data for prism/wideforprism.xe.m.expt.csv")


# ## testing just an anova....
# # in all
# aov1 <- lm(pA ~ group*sex*mV,
#          data = mc.diff.f1.long2)
# Anova(aov1, type = 2)
# 
# aov2 <- lm(pA ~ group+sex+mV,
#             data = mc.diff.f1.long2)
# Anova(aov2, type = 2)
# 
# #in males
# aov3 <- lm(pA ~ group*mV,
#             data = mc.diff.f1.long2.m)
# Anova(aov3, type = 2) #no int
# 
# aov4 <- lm(pA ~ group+mV,
#             data = mc.diff.f1.long2.m)
# Anova(aov4, type = 2) ###### SIG IN MALES
# 
# #in females
# aov5 <- lm(pA ~ group*mV,
#             data = mc.diff.f1.long2.f)
# Anova(aov5, type = 2)
# 
# aov6 <- lm(pA ~ group+mV,
#             data = mc.diff.f1.long2.f)
# Anova(aov6, type = 2)
### LOL JK FORGOT THIS IS REPEATED MEASURES I GOTTA USE MY LMER



######################################################################### END inserting grubbs testing and grubbs tested filtered analysis here, as well as removing -25mV timepoint. 

### actual analysis --  
mc.diff.long_filtered2 <- mc.diff.long_filtered

mc.diff.long_filtered2$mvsq <- mc.diff.long_filtered2$mV^2



mc.dif.mm.fullint <- lmer(pA~group * sex * (mV+mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.fullint)


mc.dif.mm.int1 <- lmer(pA~ group*sex + sex*(mV + mvsq) + group*(mV+mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.int1)


mc.dif.mm.main <- lmer(pA~group + sex + (mV + mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
summary(mc.dif.mm.main)



## analysis JUST MALES
# mc.diff.long_filtered2_males <- mc.diff.long_filtered2 %>% filter(sex=="M")
# 
# mc.dif.mm.fullint.MALE <- lmer(pA~group * (mV+mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
# summary(mc.dif.mm.fullint.MALE)
# 
# mc.dif.mm.main.MALE <- lmer(pA~group + (mV*mvsq) + (1 | subslice), data=mc.diff.long_filtered2, na.action=na.exclude)
# summary(mc.dif.mm.main.MALE)

## Check for linearity and homoscedasticity
#plot(mc.dif.mm.main, type = c("p", "smooth"))
#
## Check for normality of residuals
#qqnorm(resid(mc.dif.mm.main))
#qqline(resid(mc.dif.mm.main))
##theyre *fine*

###### correlation time. extract -40mV M-current pA and attach to rie

rie2 <- mc.diff.long_filtered %>%
  filter(mV == -40) %>%
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
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="Max M-Current @ -40 mV",x="Group", y = "pA") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=24, hjust = 0.4)) 

max.mc.pa.lm.int <- lm(max_MC_pA ~ sex*group, data = rie2)
summary(max.mc.pa.lm.int) #no effect here -- BUT should still report given full interaction model had definite effect. 

max.mc.pa.lm.main <- lm(max_MC_pA ~ sex+group, data = rie2)
summary(max.mc.pa.lm.main)



# ## max MC in just males t.test
# rie2.males <- rie2 %>% filter(sex == "Male")
# 
# t.test(max_MC_pA ~ group, data = rie2.males)
# 
# ## max MC in just STRESS (compare sex) t.test
# rie2.males <- rie2 %>% filter(sex == "Male")
# 
# t.test(max_MC_pA ~ group, data = rie2.males)

# CSIS_XE_sensitive_current_4Troy <- mc.diff %>% arrange(group,sex)
# write.csv(CSIS_XE_sensitive_current_4Troy, "results/CSIS_XE_sensitive_current_4Troy.csv")



##### POSTER FIGS

#750 x 600
mc.diff.long_filtered %>%
  filter(sex == "Male") %>% 
  ggplot(aes(x = mV, y = pA, color = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "A. Outward Current (TTX-XE): Males", y = "pA (+/- SEM)") +
  theme_classic()+
  coord_cartesian(ylim = c(-10,30)) +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) 


mc.diff.long_filtered %>%
  filter(sex == "Female") %>% 
  ggplot(aes(x = mV, y = pA, color = group)) +
  # Line connecting means
  stat_summary(geom = "line", fun = mean, size = 1) +
  # Points at means
  stat_summary(geom = "point", fun = mean, size = 2) +
  # Error bars (mean ± SEM)
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
  labs(title = "B. Outward Current (TTX-XE): Females", y = "pA (+/- SEM)") +
  theme_classic()+
  coord_cartesian(ylim = c(-10,30))+
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) 





# 600 x 750 or 750 x 600?
rie2 %>% 
  ggplot(aes(sex,max_MC_pA, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 1.5) +
  scale_fill_manual(values=c("#f64ed5", "#55a0fd")) +
  scale_color_manual(values=c("gray", "black"))+
  labs(title="C. Max M-Current @ -40 mV",x="Group", y = "pA") +
  theme_classic() +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) + 
  annotate("text", x = 1, y = 55, label = "p=0.13", size = 6) +
  scale_x_discrete(limits = c("Male", "Female"))






# rie2 %>%
#   ggplot(aes(x = RMP.pre, y = max_MC_pA, color = group, linetype = sex)) +
#   geom_smooth(method = "lm", se=F)+
#   labs(title = "D. M-Current * RMP Interaction", y = "pA (+/- SEM)", x = "RMP") +
#   theme_classic()+
#   coord_cartesian(ylim = c(0,35)) +
#   theme(axis.line = element_line(colour = 'black', size = 1),
#         axis.ticks = element_line(colour = "black", size = 1),
#         #legend.position="none",
#         axis.text=element_text(size=16),
#         axis.title=element_text(size=18,face="bold"),
#         plot.title = element_text(size=28, hjust = 0.4)) +
#   scale_color_manual(values=c("#f64ed5", "#55a0fd")) 
# 
# 
# 
# max.mc.by.RMP.lm.int <- lm(max_MC_pA ~ RMP.pre*sex*group, data = rie2)
# summary(max.mc.by.RMP.lm.int)


#750 x 600
rie2 %>%
  filter(sex == "Male") %>% 
  ggplot(aes(x = RMP.pre, y = max_MC_pA, color = group)) +
  geom_smooth(method = "lm", se=F)+
  labs(title = "D. M-Current * RMP Regression: Males", y = "Maximum M-current (pA @ -40 mV)", x = "RMP") +
  theme_classic()+
  geom_point(size = 3) +
  coord_cartesian(ylim = c(5,60)) +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) 


rie2.M <- rie2 %>% filter(sex == "Male") 

max.mc.by.RMP.lm.int.MALE <- lm(max_MC_pA ~ RMP.pre*group, data = rie2.M)
summary(max.mc.by.RMP.lm.int.MALE)

max.mc.by.RMP.lm.main.MALE <- lm(max_MC_pA ~ RMP.pre+group, data = rie2.M)
summary(max.mc.by.RMP.lm.main.MALE)




# 750 x 600
rie2 %>%
  filter(sex == "Female") %>% 
  ggplot(aes(x = RMP.pre, y = max_MC_pA, color = group)) +
  geom_smooth(method = "lm", se=F)+
  labs(title = "E. M-Current * RMP Regression: Females", y = "Maximum M-current (pA @ -40 mV)", x = "RMP") +
  theme_classic()+
  geom_point(size = 3) +
  coord_cartesian(ylim = c(5,60)) +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) +
  annotate("text", x = -60, y = 45, label = "Interaction p=0.03", size = 6)


rie2.F <- rie2 %>% filter(sex == "Female") 

max.mc.by.RMP.lm.int.FEM <- lm(max_MC_pA ~ RMP.pre*group, data = rie2.F)
summary(max.mc.by.RMP.lm.int.FEM)





















#### INPUT RESISTANCE * MAX MC --- SAME INTERACTION EFFECT AS MC * RMP THO
in.re.tp <- rie2 %>% 
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

in.re.tp3 <- in.re.tp2 %>% select(1,8,9) %>% pivot_wider(names_from = Inp.re.timepoint, values_from = Inp.re)

rie3 <- merge(rie2,in.re.tp3, by = "subslice")


rie3 %>%
  ggplot(aes(x = in.re.pre, y = max_MC_pA, color = group, linetype = sex)) +
  geom_smooth(method = "lm", se=F)+
  labs(title = "M-Current * Input Resistance", y = "pA (+/- SEM)", x = "Input Resistance (MΩ)") +
  theme_classic()+
  coord_cartesian(ylim = c(0,35)) +
  theme(axis.line = element_line(colour = 'black', size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        #legend.position="none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        plot.title = element_text(size=28, hjust = 0.4)) +
  scale_color_manual(values=c("#f64ed5", "#55a0fd")) 



max.mc.by.RMP.lm.int <- lm(max_MC_pA ~ RMP.pre*sex*group, data = rie3)
summary(max.mc.by.RMP.lm.int)





### i need to double check this is all correct and do runs with outliers and flags removed......
