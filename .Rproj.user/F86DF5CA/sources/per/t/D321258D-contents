
# M-current

library(Matrix)
library(lme4)
library(lmerTest)
library(tidyverse)
library(ggpubr)
library(car)




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


time_cols <- paste0("-", c(75, 70, 65, 60, 55, 50, 45, 40, 35, 30, 25))


df_long <- mcd4 %>%
  pivot_longer(cols = all_of(time_cols),
               names_to = "time",
               values_to = "reading")

# Create drug category and tp category, then average per subslice, group, sex, time, drug_cat, tp_cat
averaged <- df_long %>%
  mutate(
    drug_cat = case_when(
      drug %in% c("TTX1", "TTX2") ~ "TTX",
      drug %in% c("XE1", "XE2")  ~ "XE"
    ),
    tp_cat = case_when(
      tp %in% c("Y1", "Y2") ~ "Y1Y2",
      tp %in% c("Y3", "Y4") ~ "Y3Y4"
    )
  ) %>%
  group_by(subslice, group, sex, time, drug_cat, tp_cat) %>%
  summarise(avg_reading = mean(reading, na.rm = TRUE), .groups = "drop")

# Pivot wider so that each combination of drug_cat and tp_cat becomes a column
result <- averaged %>%
  unite("col_name", drug_cat, tp_cat, sep = "_avg_") %>%
  pivot_wider(id_cols = c(subslice, group, sex, time),
              names_from = col_name,
              values_from = avg_reading)

# Arrange by subslice and time (so time order is preserved)
result <- result %>%
  mutate(time = factor(time, levels = time_cols)) %>%
  arrange(subslice, group, sex, time)