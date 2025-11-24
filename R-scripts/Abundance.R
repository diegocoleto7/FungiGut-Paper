library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
ps_mico <- read_rds("Batch/MiCoP_corrected.rds")
ps_fun <- read_rds("Batch/phyloseq_object_corrected.rds")
#View(otu_table(ps_mico))
#View(otu_table(ps_fun))

# Filtrado
thr <- 0.0001

# FungiGut
df_fun <- psmelt(ps_fun) %>%
  mutate(
    Species   = sub(" ", "_", Species),
    OTU       = sub(" ", "_", OTU),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(Abundance >= thr) %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

# MiCoP
df_mico <- psmelt(ps_mico) %>%
  mutate(
    Species   = sub(" ", "_", Species),
    OTU       = sub(" ", "_", OTU),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(Abundance >= thr) %>%
  group_by(Sample) %>%
  mutate(Abundance = 100 * Abundance / sum(Abundance, na.rm = TRUE)) %>%
  ungroup()


df_fun$Group  <- as.factor(df_fun$Group)
df_mico$Group <- as.factor(df_mico$Group)

tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")


# FungiGut
meta_fun <- df_fun %>%
  distinct(Sample, Group)

tax_fun <- df_fun %>%
  select(Species, all_of(tax_cols)) %>%
  distinct(Species, .keep_all = TRUE)

tab_fun_by_group <- df_fun %>%
  select(Sample, Species, Abundance) %>%
  complete(Sample, Species, fill = list(Abundance = 0)) %>%
  left_join(tax_fun, by = "Species") %>%
  left_join(meta_fun, by = "Sample") %>%
  group_by(Group, Species) %>%
  summarise(
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Abundance_mean = mean(Abundance, na.rm = TRUE),
    .groups = "drop") %>% arrange(Group, desc(Abundance_mean))

tab_fun_celiac <- tab_fun_by_group %>% filter(Group == "Celiac") %>% select(-Group)
tab_fun_metacardis <- tab_fun_by_group %>% filter(Group == "Non_Celiac") %>% select(-Group)

# MiCoP
meta_mico <- df_mico %>%
  distinct(Sample, Group)

tax_mico <- df_mico %>%
  select(Species, all_of(tax_cols)) %>%
  distinct(Species, .keep_all = TRUE)

tab_mico_by_group <- df_mico %>%
  select(Sample, Species, Abundance) %>%
  complete(Sample, Species, fill = list(Abundance = 0)) %>%
  left_join(tax_mico, by = "Species") %>%
  left_join(meta_mico, by = "Sample") %>%
  group_by(Group, Species) %>%
  summarise(
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Abundance_mean = mean(Abundance, na.rm = TRUE),
    .groups = "drop") %>% arrange(Group, desc(Abundance_mean))

tab_mico_celiac     <- tab_mico_by_group %>% filter(Group == "Celiac") %>% select(-Group)
tab_mico_metacardis <- tab_mico_by_group %>% filter(Group == "Non_Celiac") %>% select(-Group)

# Save csv
write_csv(tab_fun_celiac, file = "Tables//FungiGut_Celiac_abundance_mean.csv")

write_csv(tab_fun_metacardis, file = "Tables/FungiGut_Metacardis_abundance_mean.csv")

write_csv(tab_mico_celiac, file = "Tables/MiCoP_Celiac_abundance_mean.csv")

write_csv(tab_mico_metacardis, file = "Tables/MiCoP_Metacardis_abundance_mean.csv")



# FungiGut
prev_fun <- df_fun %>%
  group_by(Group, Species) %>%
  summarise(
    n_samples_present = n_distinct(Sample[Abundance > 0]),
    .groups = "drop")

# MiCoP
prev_mico <- df_mico %>%
  group_by(Group, Species) %>%
  summarise(
    n_samples_present = n_distinct(Sample[Abundance > 0]),
    .groups = "drop")

View(prev_mico)
View(prev_fun)
