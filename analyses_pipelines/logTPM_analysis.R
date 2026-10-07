# ============================================================
# PIPELINE V2: log2TPM + AUC + paired Wilcoxon + pravilna BH korekcija
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(pROC)

# ============================================================
# 1) VHODNI PODATKI
# ============================================================

CPTAC <- read.csv(
  "CPTAC-3_final.csv",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
human1 <- read_excel("Human-GEM.xlsx")

# vzemi gene iz CPTAC (stolpci od 3 naprej)
cptac_genes <- colnames(CPTAC)[3:ncol(CPTAC)]

# metadata
meta <- CPTAC[, 1:2, drop = FALSE]

# ekspresijski del
expr <- CPTAC[, 3:ncol(CPTAC), drop = FALSE]

# originalna imena genov
gene_names <- colnames(expr)

# odstrani verzije ENSGxxxx.xx -> ENSGxxxx
gene_names_clean <- sub("\\..*", "", gene_names)

# preimenuj stolpce
colnames(expr) <- gene_names_clean

# vse pretvori v numeric
expr[] <- lapply(expr, as.numeric)

# združi gene z istim imenom kot povprečje TPM
expr_collapsed <- sapply(split(seq_along(colnames(expr)), colnames(expr)), function(idx) {
  if (length(idx) == 1) {
    expr[[idx]]
  } else {
    rowMeans(expr[, idx, drop = FALSE], na.rm = TRUE)
  }
})

expr_collapsed <- as.data.frame(expr_collapsed, check.names = FALSE)
# sestavi začasni CPTAC po združevanju genov
CPTAC_gene <- cbind(meta, expr_collapsed)
#========================================================
# 2) ČIŠČENJE VZORCEV:
# združi več normal ali več tumor vzorcev istega pacienta
# kot povprečno TPM vrednost, ampak samo za paciente,
# ki imajo prisoten oba tipa vzorca (normal in tumor)
#========================================================

# preštej vzorce na pacienta
patient_counts <- CPTAC_gene %>%
  count(id_pacienta, tip_vzorca) %>%
  tidyr::pivot_wider(names_from = tip_vzorca,
                     values_from = n,
                     values_fill = 0)

# obdrži samo paciente, ki imajo vsaj en normal in vsaj en tumor
valid_patients <- patient_counts %>%
  filter(normal >= 1, tumor >= 1) %>%
  pull(id_pacienta)

CPTAC_gene <- CPTAC_gene %>%
  filter(id_pacienta %in% valid_patients)

# združi po pacientu in tipu vzorca kot povprečje TPM
CPTAC <- CPTAC_gene %>%
  group_by(id_pacienta, tip_vzorca) %>%
  summarise(across(everything(), ~ mean(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

#========================================================
# 3) priprava cptac_genes po čiščenju
#========================================================

cptac_genes <- colnames(CPTAC)[3:ncol(CPTAC)]
write.csv(CPTAC, "CPTAC_clean.csv", row.names = FALSE)

CPTAC_log2 <- CPTAC
CPTAC_log2[, 3:ncol(CPTAC_log2)] <- log2(CPTAC_log2[, 3:ncol(CPTAC_log2)] + 1)
write.csv(CPTAC_log2, "CPTAC_clean_log2.csv", row.names = FALSE)
CPTAC <-  CPTAC_log2


gene_subsystem <- human1 %>%
  select(`GENE ASSOCIATION`, SUBSYSTEM) %>% #izbereš le stolpca GENE ASSOCIATION in SUBSYSTEM
  filter(!is.na(`GENE ASSOCIATION`), !is.na(SUBSYSTEM)) %>% #obdržijo se tiste vrstice ki imajo oba podatka
  mutate(
    gene_association = str_replace_all(`GENE ASSOCIATION`, "\\(|\\)", ""), #odstraniš ( ali )
    gene_association = str_replace_all(gene_association, "\\band\\b|\\bor\\b", ";") #zamenjaš and ali or z ;
  ) %>%
  separate_rows(gene_association, sep = ";") %>% #ločitev gene glede na ; -- vsak gen dobi svojo vrstico in asociacijo
  mutate(gene_association = str_trim(gene_association)) %>% #odstrani odvečmne presledke na začetku in koncu
  filter(gene_association != "") #odstrani prazne vrstice


gene_subsystem <- gene_subsystem %>% rename(gene_id = gene_association) #preimenuje gene_association v gene_id

write.csv2(gene_subsystem, "TPM_human1_annotation.csv", row.names = FALSE)

gene_subsystem_cptac <- gene_subsystem %>% filter(gene_id %in% cptac_genes) #vzame le tiste gene, ki se nahajajo v CPTAC datoteki in v TPM_human1_annotation.csv
gene_subsystem_cptac <- gene_subsystem_cptac %>% distinct(gene_id, SUBSYSTEM) #odstrani duplikate vrstic
write.csv2(gene_subsystem_cptac, "TPM_CPTAC_human1_annotation.csv", row.names = FALSE)
#_________________________________________________
#izračun evklidove distance za vsak podsistem
#________________________________________________
gene_matrix <- CPTAC[, 3:ncol(CPTAC)]
subsystem_list <- split(gene_subsystem_cptac$gene_id,gene_subsystem_cptac$SUBSYSTEM)
patients <- unique(CPTAC$id_pacienta)

results <- data.frame()

for(p in patients){ #for zanka 
  normal_row <- which(CPTAC$id_pacienta == p & CPTAC$tip_vzorca == "normal") #poiščeš vrstico normalnega vzorca tega pacienta
  tumor_row  <- which(CPTAC$id_pacienta == p & CPTAC$tip_vzorca == "tumor") #poiščeš vrstico tumorskega vzorca tega pacienta
  
  if(length(normal_row)==1 & length(tumor_row)==1){ #če ima le eno vrstico normalnega in eno vrstico tumorksega vzorca
    for(sub in names(subsystem_list)){ #druga for zanka za pregled podsistemov
      genes <- subsystem_list[[sub]] #vzame gene tega podsistema
      genes <- genes[genes %in% colnames(gene_matrix)] #preveri če se ti geni podsistema nahajajo v CPTAC 
      if(length(genes) > 1){ #preveri da sta vsaj dva gena v podsistemu
        normal_vec <- as.numeric(gene_matrix[normal_row, genes]) #vzame TPM vrednosti normalnega vzorca za gene tega podsistema
        tumor_vec  <- as.numeric(gene_matrix[tumor_row, genes]) 
        d <- sqrt(sum((tumor_vec - normal_vec)^2)) #izračuna evklidova distanca, vsi geni so v enakem vrsten redu, zato se lahko matriki normal_vec in tumor_vec preprosto odštejeta
        results <- rbind(results,
                         data.frame(patient=p,
                                    subsystem=sub,
                                    distance=d,
                                    n_genes=length(genes)))
      }
    }
  }
}
results <- results[order(-results$distance),] 


results$distance_norm <- results$distance / sqrt(results$n_genes)
summary_results <- results %>%
  group_by(subsystem) %>%
  summarise(
    mean_distance = mean(distance_norm),
    median_distance = median(distance_norm)
  ) %>%
  arrange(desc(median_distance))
install.packages("ggplot2")
library(ggplot2)

ggplot(summary_results[1:15,],
       aes(x=reorder(subsystem, median_distance),
           y=median_distance)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(x="Subsystem",
       y="Median tumor-normal distance")

results <- results[order(-results$distance_norm),] 
write.csv2(results, "distanca_pacientov_log2.csv", row.names = FALSE)
write.csv2(summary_results, "mediana_distanc_podsistemov_log2.csv", row.names = FALSE)
mediana_podsistemov <- summary_results
#__________________________________________
#analiza podatkov
#________________________________________

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(pROC)

# load data
CPTAC <- read.csv("CPTAC_clean_log2.csv")
results <- read.csv2("distanca_pacientov_log2.csv")
gene_subsystem_cptac <- read.csv2("TPM_CPTAC_human1_annotation.csv")
mediana_podsistemov <- read.csv2("mediana_distanc_podsistemov_log2.csv")
kandidati <- read.csv2("top10_podsistemi_gene_auc_log2.csv")

# PARAMETRI
izbran_pacient <- "C3L-02888"
izbran_podsistem <- "Retinol metabolism"
outlier_gene <- "ENSG00000025423"

# geni podsistema
genes_sub <- gene_subsystem_cptac %>%
  filter(SUBSYSTEM == izbran_podsistem) %>%
  pull(gene_id) %>%
  unique()

genes_sub <- genes_sub[genes_sub %in% colnames(CPTAC)]

# podatki
patient_data <- CPTAC %>%
  filter(id_pacienta == izbran_pacient) %>%
  dplyr::select(id_pacienta, tip_vzorca, all_of(genes_sub))

# long format
patient_long <- patient_data %>%
  pivot_longer(
    cols = -c(id_pacienta, tip_vzorca),
    names_to = "gene_id",
    values_to = "TPM"
  ) %>%
  filter(!is.na(TPM)) %>%
  mutate(color_group = case_when(
    gene_id == outlier_gene ~ "outlier",
    tip_vzorca == "normal" ~ "normal",
    tip_vzorca == "tumor" ~ "tumor"
  ))

# jitter plot
ggplot(patient_long,
       aes(x = tip_vzorca, y = log10(TPM))) +
  geom_jitter(
    aes(color = color_group),
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(values = c(
    "outlier" = "green",
    "normal" = "red",
    "tumor" = "blue"
  )) +  geom_text(
    data = subset(patient_long, gene_id == outlier_gene),
    aes(label = gene_id),
    vjust = -1,
    color = "green",
    size = 3
  ) +
  labs(
    title = paste("TPM genov v podsistemu:", izbran_podsistem),
    subtitle = paste("Pacient:", izbran_pacient),
    x = "Tip vzorca",
    y = "log10(TPM+1)",
    color = "Skupina"
  ) +
  theme_minimal()

CPTAC %>%
  select(id_pacienta, tip_vzorca, ENSG00000170835) %>%
  ggplot(aes(x = tip_vzorca, y = ENSG00000170835)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

ggplot(CPTAC,
       aes(x = tip_vzorca,
           y = log10(ENSG00000170835 + 1))) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

wilcox.test(
  ENSG00000170835 ~ tip_vzorca,
  data = CPTAC
)

# jitter graph
ggplot(patient_long,
       aes(x = tip_vzorca, y = TPM)) +
  geom_jitter(
    aes(color = color_group),
    width = 0.12,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(values = c(
    "outlier" = "green",
    "normal" = "red",
    "tumor" = "blue"
  )) +  geom_text(
    data = subset(patient_long, gene_id == outlier_gene),
    aes(label = gene_id),
    vjust = -1,
    color = "green",
    size = 3
  ) +
  labs(
    title = paste("TPM genov v podsistemu:", izbran_podsistem),
    subtitle = paste("Pacient:", izbran_pacient),
    x = "Tip vzorca",
    y = "TPM",
    color = "Skupina"
  ) +
  theme_minimal()


# ROC
roc_data <- CPTAC %>%
  mutate(label = ifelse(tip_vzorca == "tumor", 1, 0))

roc_obj <- roc(roc_data$label, roc_data$ENSG00000203805)

plot(roc_obj)
auc(roc_obj)


CPTAC %>%
  select(id_pacienta, tip_vzorca, ENSG00000203805) %>%
  ggplot(aes(x = tip_vzorca, y = ENSG00000203805)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

ggplot(CPTAC,
       aes(x = tip_vzorca,
           y = log10(ENSG00000170835 + 1))) +
  geom_boxplot() +
  geom_jitter(width = 0.1)

wilcox.test(
  ENSG00000170835 ~ tip_vzorca,
  data = CPTAC
)
ggplot(CPTAC,
       aes(x = log2(ENSG00000025423 + 1),
           fill = tip_vzorca)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = "Gostota porazdelitve",
    x = "log2(TPM + 1)"
  ) +
  theme_minimal()

ggplot(CPTAC,
       aes(x = log10(ENSG00000100994 + 1),
           color = tip_vzorca)) +
  stat_ecdf(size = 1) +
  scale_color_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = "ECDF – kumulativna porazdelitev",
    x = "log10(TPM + 1)",
    y = "F(x)"
  ) +
  theme_minimal()

#ROC KRIVULJA
install.packages("pROC")
library(pROC)
roc_data <- CPTAC %>%
  mutate(
    label = ifelse(tip_vzorca == "tumor", 1, 0)
  )
roc_obj <- roc(
  response = roc_data$label,
  predictor = roc_data$ENSG00000100994
)

plot(roc_obj, col = "blue", lwd = 3,
     main = "ROC curve – ENSG00000100994")
abline(a = 0, b = 1, lty = 2, col = "grey")
auc(roc_obj)
#vzamem top 10 podsistemov glede na mediano normalizirane distance iz mediana_podsistemov. Nato pa vsakemu genu v teh podsistemih pripišem AUC
library(dplyr)
library(pROC)

# 1) namesto n() vstavimo koliko top podsistemov želimo nadaljno analizirati, če se pusti n, se vzamejo vsi podsistemi
top10_subsystems <- mediana_podsistemov %>%
  arrange(desc(median_distance)) %>%
  slice(1:n()) %>%
  pull(subsystem)

# 2) geni iz teh podsistemov
top_genes_df <- gene_subsystem_cptac %>%
  filter(SUBSYSTEM %in% top10_subsystems) %>%
  distinct(SUBSYSTEM, gene_id)

# 3) priprava labela
CPTAC_auc <- CPTAC %>%
  filter(tip_vzorca %in% c("normal", "tumor")) %>%
  mutate(label = ifelse(tip_vzorca == "tumor", 1, 0))

# 4) izračun AUC za vsak gen
gene_auc_results <- data.frame()

for (i in seq_len(nrow(top_genes_df))) {
  
  sub_i <- top_genes_df$SUBSYSTEM[i]
  gene_i <- top_genes_df$gene_id[i]
  
  if (gene_i %in% colnames(CPTAC_auc)) {
    
    tmp <- CPTAC_auc %>%
      select(id_pacienta, tip_vzorca, label, all_of(gene_i)) %>%
      filter(!is.na(.data[[gene_i]]))
    
    # AUC rabi vsaj 2 razreda
    if (length(unique(tmp$label)) == 2) {
      
      roc_obj <- roc(
        response = tmp$label,
        predictor = tmp[[gene_i]],
        quiet = TRUE
      )
      
      auc_val <- as.numeric(auc(roc_obj))
      
      # obrni smer, da bo AUC vedno >= 0.5
      if (auc_val < 0.5) {
        auc_val <- 1 - auc_val
      }
      
      gene_auc_results <- rbind(
        gene_auc_results,
        data.frame(
          subsystem = sub_i,
          gene_id = gene_i,
          AUC = auc_val,
          n_samples = nrow(tmp)
        )
      )
    }
  }
}

# 5) uredi po podsistemu in AUC
gene_auc_results <- gene_auc_results %>% arrange(desc(AUC))

# 6) shrani
write.csv2(gene_auc_results, "vsi_podsistemi_gene_auc_log2.csv", row.names = FALSE)
# dodal sem Rank_log2TPM_norm_distance v datoteko mediana_distanc_podsistemov_log. Zdaj pa želim ta rank dodati tudi v  gene_auc_results
mediana_podsistemov <- read.csv2(
  "mediana_distanc_podsistemov_log2.csv",
  stringsAsFactors = FALSE
)

# uredi po median_distance in dodaj rank
mediana_podsistemov <- mediana_podsistemov %>%
  arrange(desc(median_distance)) %>%
  mutate(Rank_log2TPM_norm_distance = row_number())

write.csv2(
  mediana_podsistemov,
  "mediana_distanc_podsistemov_log2.csv",
  row.names = FALSE
)
gene_auc_results <- gene_auc_results %>%
  left_join(
    mediana_podsistemov %>%
      select(subsystem, Rank_log2TPM_norm_distance),
    by = "subsystem"
  )
write.csv2(gene_auc_results, "vsi_podsistemi_gene_auc_log2.csv", row.names = FALSE)


#vizualizacija tistih genov s top AUC
library(dplyr)
library(ggplot2)

izbran_podsistem <- "Acylglycerides metabolism"
izbran_gen <- "ENSG00000166391"

plot_data <- CPTAC %>%
  select(id_pacienta, tip_vzorca, all_of(izbran_gen)) %>%
  filter(!is.na(.data[[izbran_gen]])) %>%
  rename(TPM = all_of(izbran_gen))

ggplot(plot_data, aes(x = tip_vzorca, y = TPM, color = tip_vzorca)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = paste("TPM gena", izbran_gen),
    subtitle = paste("Podsistem:", izbran_podsistem),
    x = "Tip vzorca",
    y = "TPM"
  ) +
  theme_minimal()

ggplot(plot_data, aes(x = tip_vzorca, y = log10(TPM + 1), color = tip_vzorca)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = paste("Izražanje gena", izbran_gen),
    subtitle = paste("Podsistem:", izbran_podsistem),
    x = "Tip vzorca",
    y = "log10(TPM + 1)"
  ) +
  theme_minimal()

install.packages("ggpubr")
library(ggpubr)

plot_data <- CPTAC %>%
  select(id_pacienta, tip_vzorca, all_of(izbran_gen)) %>%
  filter(tip_vzorca %in% c("normal", "tumor")) %>%
  filter(!is.na(.data[[izbran_gen]])) %>%
  rename(TPM = all_of(izbran_gen)) %>%
  mutate(
    logTPM = log10(TPM + 1),
    label = ifelse(tip_vzorca == "tumor", 1, 0)
  )

# 1) density plot
ggplot(plot_data,
       aes(x = logTPM, fill = tip_vzorca)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = "Gostota porazdelitve",
    subtitle = paste("Gen:", izbran_gen, "| Podsistem:", izbran_podsistem),
    x = "log10(TPM + 1)",
    y = "density"
  ) +
  theme_minimal()

# 2) ROC curve
roc_obj <- roc(
  response = plot_data$label,
  predictor = plot_data$TPM,
  quiet = TRUE
)

# če je AUC < 0.5, obrnemo smer, da bo interpretacija lepša
auc_val <- as.numeric(auc(roc_obj))

if (auc_val < 0.5) {
  roc_obj <- roc(
    response = plot_data$label,
    predictor = -plot_data$TPM,
    quiet = TRUE
  )
  auc_val <- as.numeric(auc(roc_obj))
}

plot(
  roc_obj,
  col = "blue",
  lwd = 3,
  main = paste("ROC curve –", izbran_gen)
)
abline(a = 0, b = 1, lty = 2, col = "grey")

auc_val
write.csv2(gene_auc_results, "top10_podsistemi_gene_auc_log2.csv", row.names = FALSE)



#primerjava top 10 log2TPm in TPM kandidatov
log2TPM <- read.csv2("top10_podsistemi_gene_auc_log2.csv",check.names = FALSE, stringsAsFactors = FALSE)
TPM <- read.csv2("top10_podsistemi_gene_auc.csv", check.names = FALSE, stringsAsFactors = FALSE)
CPTAC <- read.csv("CPTAC_clean.csv", check.names = FALSE, stringsAsFactors = FALSE)
CPTAC_log2 <- read.csv("CPTAC_clean_log2.csv", check.names = FALSE, stringsAsFactors = FALSE)
#__________________________________________
# PRIPRAVA PODATKOV ZA WILCOXON TEST
# VERZIJA 2: BH korekcija na unikatnih genih
#__________________________________________
# ============================================================
# GLOBALNA ANALIZA log2TPM NA RAVNI UNIKATNIH GENOV
# AUC + paired Wilcoxon + globalna BH-korekcija
# ============================================================

library(dplyr)
library(pROC)

# ------------------------------------------------------------
# 1) VHODNI PODATKI
# ------------------------------------------------------------

# očiščene in log2-transformirane vrednosti TPM
CPTAC_log2 <- read.csv(
  "CPTAC_clean_log2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# povezave med geni Human1 in CPTAC
# uporabimo samo za določitev vseh metabolnih genov,
# podsistemov v končni tabeli ne ohranimo
gene_annotation <- read.csv2(
  "TPM_CPTAC_human1_annotation.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2) SEZNAM VSEH UNIKATNIH METABOLNIH GENOV
# ------------------------------------------------------------

unique_genes <- gene_annotation %>%
  distinct(gene_id) %>%
  filter(gene_id %in% colnames(CPTAC_log2))

cat(
  "Število unikatnih metabolnih genov za analizo:",
  nrow(unique_genes),
  "\n"
)


# ------------------------------------------------------------
# 3) PRIPRAVA PARNIH NORMALNIH IN TUMORSKIH VZORCEV
# ------------------------------------------------------------

normal_log <- CPTAC_log2 %>%
  filter(tip_vzorca == "normal") %>%
  arrange(id_pacienta)

tumor_log <- CPTAC_log2 %>%
  filter(tip_vzorca == "tumor") %>%
  arrange(id_pacienta)

# preveri, da imata skupini enako število pacientov
stopifnot(nrow(normal_log) == nrow(tumor_log))

# preveri, da so pacienti v popolnoma enakem vrstnem redu
stopifnot(
  identical(
    as.character(normal_log$id_pacienta),
    as.character(tumor_log$id_pacienta)
  )
)

cat(
  "Število parnih pacientov:",
  nrow(normal_log),
  "\n"
)


# ------------------------------------------------------------
# 4) PRAZNA TABELA REZULTATOV
# ena vrstica = en unikaten gen
# ------------------------------------------------------------

log2TPM_global <- unique_genes %>%
  mutate(
    AUC = NA_real_,
    n_pairs = NA_integer_,
    n_samples = NA_integer_,
    wilcox_p = NA_real_
  )


# ------------------------------------------------------------
# 5) AUC IN PARNI WILCOXONOV TEST
# ------------------------------------------------------------

for (i in seq_len(nrow(log2TPM_global))) {
  
  gene <- log2TPM_global$gene_id[i]
  
  normal_vals <- normal_log[[gene]]
  tumor_vals  <- tumor_log[[gene]]
  
  # ohranimo samo paciente, pri katerih sta prisotni obe vrednosti
  valid_pairs <- !is.na(normal_vals) & !is.na(tumor_vals)
  
  n_pairs <- sum(valid_pairs)
  
  log2TPM_global$n_pairs[i] <- n_pairs
  log2TPM_global$n_samples[i] <- 2 * n_pairs
  
  # če v metodah piše "vsaj 5 parov", mora biti >= 5
  if (n_pairs >= 5) {
    
    normal_valid <- normal_vals[valid_pairs]
    tumor_valid  <- tumor_vals[valid_pairs]
    
    # --------------------------------------------------------
    # PARNI WILCOXONOV TEST
    # --------------------------------------------------------
    
    wilcox_result <- tryCatch(
      wilcox.test(
        tumor_valid,
        normal_valid,
        paired = TRUE,
        exact = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(wilcox_result)) {
      log2TPM_global$wilcox_p[i] <- wilcox_result$p.value
    }
    
    
    # --------------------------------------------------------
    # AUC
    # uporabimo ista veljavna para kot pri Wilcoxonovem testu
    # --------------------------------------------------------
    
    auc_values <- c(
      normal_valid,
      tumor_valid
    )
    
    auc_labels <- c(
      rep(0, n_pairs),   # normal
      rep(1, n_pairs)    # tumor
    )
    
    roc_result <- tryCatch(
      roc(
        response = auc_labels,
        predictor = auc_values,
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(roc_result)) {
      
      auc_value <- as.numeric(auc(roc_result))
      
      # AUC naj predstavlja moč ločevanja ne glede na smer
      if (auc_value < 0.5) {
        auc_value <- 1 - auc_value
      }
      
      log2TPM_global$AUC[i] <- auc_value
    }
  }
}


# ------------------------------------------------------------
# 6) GLOBALNA BH-KOREKCIJA
# čez vse unikatne analizirane metabolne gene
# ------------------------------------------------------------

log2TPM_global$wilcox_padj <- p.adjust(
  log2TPM_global$wilcox_p,
  method = "BH"
)


# ------------------------------------------------------------
# 7) UREJANJE REZULTATOV
# ------------------------------------------------------------

log2TPM_global <- log2TPM_global %>%
  arrange(desc(AUC))


# ------------------------------------------------------------
# 8) SHRANJEVANJE
# ------------------------------------------------------------

write.csv2(
  log2TPM_global,
  "globalni_unikatni_geni_auc_wilcox_log2TPM.csv",
  row.names = FALSE
)
# koliko genov ima p < 0.01
sum(log2TPM_unique$wilcox_p < 0.01, na.rm = TRUE)

# koliko genov ima padj < 0.01
sum(log2TPM_unique$wilcox_padj < 0.01, na.rm = TRUE)

# koliko genov ima padj < 0.01 in AUC > 0.9
# tukaj uporabimo tabelo na ravni gen-podsistem, zato lahko isti gen šteje večkrat
sum(log2TPM_v2$wilcox_padj < 0.01 & log2TPM_v2$AUC > 0.9, na.rm = TRUE)

# če želiš šteti unikatne gene z AUC > 0.9 in padj < 0.01:
log2TPM_v2 %>%
  distinct(gene_id, AUC, wilcox_padj) %>%
  summarise(
    n_unique_genes_padj_001_AUC_09 =
      sum(wilcox_padj < 0.01 & AUC > 0.9, na.rm = TRUE)
  )