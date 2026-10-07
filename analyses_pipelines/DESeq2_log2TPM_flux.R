#dodatek log2FoldChange in padj iz DESeq2 v log2TPM
library(dplyr)
library(readxl)
library(stringr)
library(tidyr)


#branje rezultatov posameznih pipeline-ov in preimenovanje stolpcev za boljšo preglednost
log2TPM <- read.csv2("v2_vsi_podsistemi_gene_auc_log2.csv")
colnames(log2TPM)[colnames(log2TPM) == "wilcox_padj"] <- "wilcox_padj_log2TPM"
colnames(log2TPM)[colnames(log2TPM) == "AUC"] <- "AUC_log2TPM"

DESeq2 <- read.csv2("anotirani_kandidati_human1.csv")
colnames(DESeq2)[colnames(DESeq2) == "SUBSYSTEM"] <- "subsystem"
colnames(DESeq2)[colnames(DESeq2) == "padj"] <- "padj_DESeq2"
colnames(DESeq2)[colnames(DESeq2) == "baseMean"] <- "mean_norm_counts"

flux <- read.csv2("all_podsistemi_reaction_auc_flux.csv")
colnames(flux)[colnames(flux) == "wilcox_padj"] <- "wilcox_padj_flux"
colnames(flux)[colnames(flux) == "AUC"] <- "AUC_flux"

human1 <- read_excel("Human-GEM.xlsx")
library(dplyr)

#___________________________________________
#dodatek log2TPM rezultatov k DESeq2 rezultatom
#__________________________________________

# iz DESeq2 vzemi samo želene stolpce
DESeq2_selected <- DESeq2 %>%
  select(gene_id, mean_norm_counts, log2FoldChange, padj_DESeq2, subsystem) %>%
  distinct(gene_id, subsystem, .keep_all = TRUE)

# iz log2TPM vzemi samo želene stolpce
log2TPM_selected <- log2TPM %>%
  select(gene_id, subsystem, AUC_log2TPM, Rank_log2TPM_norm_distance, wilcox_padj_log2TPM) %>%
  distinct(gene_id, subsystem, .keep_all = TRUE)

# osnova je DESeq2, da se ne izgubijo geni, ki niso v human1 / log2TPM
DESeq2_z_log2TPM <- DESeq2_selected %>%
  left_join(
    log2TPM_selected,
    by = c(
      "gene_id" = "gene_id",
      "subsystem" = "subsystem"
    )
  ) %>%
  distinct()

# shrani rezultat
write.csv2(DESeq2_z_log2TPM, "DESeq2_z_log2TPM.csv", row.names = FALSE)

#___________________________________________
#dodatek flux rezultatov k DESeq2_log2TPM rezultatom
#__________________________________________

# priprava human1: reakcija -> gen -> podsistem

human1_gene_reaction <- human1 %>%
  select(ID, SUBSYSTEM, `GENE ASSOCIATION`) %>%
  filter(!is.na(ID), !is.na(SUBSYSTEM), !is.na(`GENE ASSOCIATION`)) %>%
  mutate(
    `GENE ASSOCIATION` = str_replace_all(`GENE ASSOCIATION`, "\\(|\\)", ""),
    `GENE ASSOCIATION` = str_replace_all(`GENE ASSOCIATION`, "\\band\\b|\\bor\\b", ";")
  ) %>%
  separate_rows(`GENE ASSOCIATION`, sep = ";") %>%
  mutate(`GENE ASSOCIATION` = str_trim(`GENE ASSOCIATION`)) %>%
  filter(`GENE ASSOCIATION` != "") %>%
  distinct(ID, SUBSYSTEM, `GENE ASSOCIATION`)

# optional sanity check
cat("Število unikatnih reaction-subsystem-gene povezav v Human1:\n")
print(nrow(human1_gene_reaction))

# join flux rezultatov na Human1
# pogoji: reaction_id == ID subsystem == SUBSYSTEM
flux_gene_level <- flux %>%
  left_join(
    human1_gene_reaction,
    by = c(
      "reaction_id" = "ID",
      "subsystem" = "SUBSYSTEM"
    )
  ) %>%
  distinct()

# preveri koliko flux reakcij se je mapiralo na gene
cat("\nŠtevilo vrstic po joinu flux -> Human1:\n")
print(nrow(flux_gene_level))

cat("\nŠtevilo vrstic brez gena po joinu:\n")
print(sum(is.na(flux_gene_level$`GENE ASSOCIATION`)))

# preimenuj GENE ASSOCIATION v gene_id za lažji končni join
flux_gene_level <- flux_gene_level %>%
  rename(gene_id = `GENE ASSOCIATION`)

# dodatek flux rezultatov k DESeq2_z_log2TPM, pogoji: gene_id == gene_id subsystem == subsystem

DESeq2_z_log2TPM_z_flux <- DESeq2_z_log2TPM %>%
  left_join(
    flux_gene_level %>%
      select(
        gene_id,
        subsystem,
        reaction_id,
        AUC_flux,
        n_samples,
        wilcox_p,
        wilcox_padj_flux,
        Rank_flux_norm_distance
      ) %>%
      distinct(),
    by = c(
      "gene_id" = "gene_id",
      "subsystem" = "subsystem"
    )
  ) %>%
  distinct()

#odstranitev nepotrebnih stolpcev
DESeq2_z_log2TPM_z_flux <- DESeq2_z_log2TPM_z_flux %>%
  select(-wilcox_p)
DESeq2_z_log2TPM_z_flux <- DESeq2_z_log2TPM_z_flux %>%
  select(-n_samples)
#premestitev stolpcev za boljšo preglednost
DESeq2_z_log2TPM_z_flux <- DESeq2_z_log2TPM_z_flux %>%
  relocate(reaction_id, .after = gene_id)

DESeq2_z_log2TPM_z_flux <- DESeq2_z_log2TPM_z_flux %>%
  relocate(subsystem, .after = reaction_id)

write.csv2(
  DESeq2_z_log2TPM_z_flux,
  "DESeq2_z_log2TPM_z_flux.csv",
  row.names = FALSE
)
#___________________________________________
#analiza rezultatov
#__________________________________________
#strog kriterij- upoštevanje fluxov in tpm

# naloži glavno tabelo
rezultati <- read.csv2("DESeq2_z_log2TPM_z_flux.csv", stringsAsFactors = FALSE)

# STROG FILTER- tako na log2TPM in flux morata veljati pogoja
najboljsi_kandidati <- rezultati %>%
  filter(
    !is.na(wilcox_padj_log2TPM),
    !is.na(AUC_log2TPM),
    !is.na(wilcox_padj_flux),
    !is.na(AUC_flux),
    wilcox_padj_log2TPM < 0.01,
    AUC_log2TPM > 0.8,
    wilcox_padj_flux < 0.01,
    AUC_flux > 0.8
  ) %>%
  distinct()

# hiter pregled
cat("Število kandidatov:", nrow(najboljsi_kandidati), "\n")
print(head(najboljsi_kandidati, 20))

# shrani v csv
write.csv2(
  najboljsi_kandidati,
  "najboljsi_kandidati_log2TPM_in_flux.csv",
  row.names = FALSE
)

#MILEJŠI FILTER- izpolnjeni pogoji na log2TPM ali fluxi
kandidati_milejsi <- rezultati %>%
  filter(
    (
      !is.na(wilcox_padj_log2TPM) &
        !is.na(AUC_log2TPM) &
        wilcox_padj_log2TPM < 0.01 &
        AUC_log2TPM > 0.9
    )
    |
      (
        !is.na(wilcox_padj_flux) &
          !is.na(AUC_flux) &
          wilcox_padj_flux < 0.01 &
          AUC_flux > 0.9
      )
  ) %>%
  distinct(gene_id, reaction_id, subsystem, .keep_all = TRUE)

# pregled
cat("Število kandidatov:", nrow(kandidati_milejsi), "\n")

# shrani
write.csv2(
  kandidati_milejsi,
  "kandidati_milejsi_log2TPM_ali_flux.csv",
  row.names = FALSE
)
#grafični prikaz najboljših kandidatov

library(dplyr)
library(ggplot2)
library(pROC)
library(readxl)

#__________________________________________
# 1. BRANJE PODATKOV
#__________________________________________

kandidati <- read.csv2("DESeq2_z_log2TPM_z_flux.csv", stringsAsFactors = FALSE)
CPTAC_log2 <- read.csv("CPTAC_clean_log2.csv", check.names = FALSE, stringsAsFactors = FALSE)
flux_clean <- read.csv("flux_clean.csv", check.names = FALSE, stringsAsFactors = FALSE)
human1 <- read_excel("Human-GEM.xlsx")

#__________________________________________
# 2. IZBIRA POSAMEZNEGA KANDIDATA
#__________________________________________

izbran_gen <- "ENSG00000176153"
izbrana_reakcija <- "MAR04116"

# poišči točno ta par v master tabeli
candidate_row <- kandidati %>%
  filter(gene_id == izbran_gen, reaction_id == izbrana_reakcija)

if (nrow(candidate_row) == 0) {
  stop("Izbrani gene_id + reaction_id par ni bil najden v master tabeli.")
}

if (nrow(candidate_row) > 1) {
  warning("Za ta gene_id + reaction_id obstaja več vrstic. Uporabljam prvo.")
  candidate_row <- candidate_row[1, ]
}

print(candidate_row)

izbran_podsistem <- candidate_row$subsystem[1]

# dodatne info iz Human1
gene_info <- human1 %>%
  select(`GENE ASSOCIATION`, SUBSYSTEM) %>%
  filter(SUBSYSTEM == izbran_podsistem)

reaction_info <- human1 %>%
  select(ID, NAME, SUBSYSTEM, EQUATION) %>%
  filter(ID == izbrana_reakcija)

print(reaction_info)

#__________________________________________
# 3. PRIPRAVA GENE PODATKOV (log2TPM)
#__________________________________________

if (!(izbran_gen %in% colnames(CPTAC_log2))) {
  stop(paste("Gen", izbran_gen, "ne obstaja v CPTAC_clean_log2.csv"))
}

gene_data <- CPTAC_log2 %>%
  select(id_pacienta, tip_vzorca, all_of(izbran_gen)) %>%
  rename(expr = all_of(izbran_gen)) %>%
  mutate(tip_vzorca = factor(tip_vzorca, levels = c("normal", "tumor")))

summary(gene_data$expr)
table(gene_data$tip_vzorca)

# ROC za gen
roc_gene_data <- gene_data %>%
  mutate(label = ifelse(tip_vzorca == "tumor", 1, 0))

roc_gene <- roc(
  response = roc_gene_data$label,
  predictor = roc_gene_data$expr,
  quiet = TRUE
)

auc_gene <- as.numeric(auc(roc_gene))
auc_gene_display <- ifelse(auc_gene < 0.5, 1 - auc_gene, auc_gene)

#__________________________________________
# 4. PRIPRAVA REACTION PODATKOV (flux)
#__________________________________________

if (!(izbrana_reakcija %in% colnames(flux_clean))) {
  stop(paste("Reakcija", izbrana_reakcija, "ne obstaja v flux_clean.csv"))
}

reaction_data <- flux_clean %>%
  select(patient_id, type, all_of(izbrana_reakcija)) %>%
  rename(flux = all_of(izbrana_reakcija)) %>%
  mutate(type = factor(type, levels = c("normal", "tumor")))

summary(reaction_data$flux)
table(reaction_data$type)

# ROC za reakcijo
roc_reaction_data <- reaction_data %>%
  mutate(label = ifelse(type == "tumor", 1, 0))

roc_reaction <- roc(
  response = roc_reaction_data$label,
  predictor = roc_reaction_data$flux,
  quiet = TRUE
)

auc_reaction <- as.numeric(auc(roc_reaction))
auc_reaction_display <- ifelse(auc_reaction < 0.5, 1 - auc_reaction, auc_reaction)

#__________________________________________
# 5. GRAFI ZA GEN
#__________________________________________

# box + jitter
ggplot(gene_data, aes(x = tip_vzorca, y = expr, color = tip_vzorca)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  labs(
    title = paste("log2TPM gen", izbran_gen),
    subtitle = paste("Subsystem:", izbran_podsistem),
    x = "Tip vzorca",
    y = "log2(TPM + 1)"
  ) +
  theme_minimal()

# density
ggplot(gene_data, aes(x = expr, fill = tip_vzorca)) +
  geom_density(alpha = 0.4) +
  labs(
    title = paste("Gostota porazdelitve log2TPM za", izbran_gen),
    subtitle = paste("Subsystem:", izbran_podsistem),
    x = "log2(TPM + 1)",
    y = "Gostota"
  ) +
  theme_minimal()

# histogram
ggplot(gene_data, aes(x = expr, fill = tip_vzorca)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  labs(
    title = paste("Histogram log2TPM za", izbran_gen),
    subtitle = paste("Subsystem:", izbran_podsistem),
    x = "log2(TPM + 1)",
    y = "Število vzorcev"
  ) +
  theme_minimal()

# ROC
plot(
  roc_gene,
  main = paste("ROC za gen", izbran_gen,
               "- AUC =", round(auc_gene_display, 4))
)
abline(a = 0, b = 1, lty = 2, col = "grey")

# paired plot
ggplot(gene_data, aes(x = tip_vzorca, y = expr, group = id_pacienta)) +
  geom_line(alpha = 0.35) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = paste("Paired log2TPM za", izbran_gen),
    subtitle = paste("Subsystem:", izbran_podsistem),
    x = "Tip vzorca",
    y = "log2(TPM + 1)"
  ) +
  theme_minimal()

#__________________________________________
# 6. GRAFI ZA REAKCIJO
#__________________________________________

reaction_name <- ifelse(nrow(reaction_info) > 0, reaction_info$NAME[1], NA)
reaction_subsystem <- ifelse(nrow(reaction_info) > 0, reaction_info$SUBSYSTEM[1], NA)

# box + jitter
ggplot(reaction_data, aes(x = type, y = flux, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  labs(
    title = paste("Flux reakcije", izbrana_reakcija),
    subtitle = paste(reaction_name, "|", reaction_subsystem),
    x = "Tip vzorca",
    y = "Flux"
  ) +
  theme_minimal()

# density
ggplot(reaction_data, aes(x = flux, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(
    title = paste("Gostota porazdelitve fluxa za", izbrana_reakcija),
    subtitle = paste(reaction_name, "|", reaction_subsystem),
    x = "Flux",
    y = "Gostota"
  ) +
  theme_minimal()

# histogram
ggplot(reaction_data, aes(x = flux, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  labs(
    title = paste("Histogram fluxa za", izbrana_reakcija),
    subtitle = paste(reaction_name, "|", reaction_subsystem),
    x = "Flux",
    y = "Število vzorcev"
  ) +
  theme_minimal()

# ROC
plot(
  roc_reaction,
  main = paste("ROC za reakcijo", izbrana_reakcija,
               "- AUC =", round(auc_reaction_display, 4))
)
abline(a = 0, b = 1, lty = 2, col = "grey")

# paired plot
ggplot(reaction_data, aes(x = type, y = flux, group = patient_id)) +
  geom_line(alpha = 0.35) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = paste("Paired flux za", izbrana_reakcija),
    subtitle = paste(reaction_name, "|", reaction_subsystem),
    x = "Tip vzorca",
    y = "Flux"
  ) +
  theme_minimal()




rezultati <- read.csv2("DESeq2_z_log2TPM_z_flux.csv")