#______________________________________________________________________________
# PORAZDELITEV NORMALIZIRANIH EVKLIDOVIH RAZDALJ
# PRI TREH IZBRANIH METABOLNIH PODSISTEMIH
#______________________________________________________________________________

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)


#______________________________________________________________________________
# 1) BRANJE PODATKOV
#______________________________________________________________________________

df_raw <- read.csv2(
  "distanca_pacientov_log2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

df <- df_raw %>%
  select(
    patient,
    subsystem,
    distance,
    n_genes,
    distance_norm
  ) %>%
  mutate(
    patient = as.character(patient),
    subsystem = as.character(subsystem),
    distance_norm = as.numeric(
      str_replace(
        as.character(distance_norm),
        ",",
        "."
      )
    )
  ) %>%
  filter(
    !is.na(patient),
    !is.na(subsystem),
    !is.na(distance_norm)
  )


#______________________________________________________________________________
# 2) IZBRANI METABOLNI PODSISTEMI
#______________________________________________________________________________

selected_subsystems <- tribble(
  ~subsystem, ~subsystem_si,
  
  "Formation and hydrolysis of cholesterol esters",
  "Sinteza in hidroliza holesterol estrov",
  
  "Protein assembly",
  "Sestava proteinov",
  
  "Protein degradation",
  "Degradacija proteinov"
)


#______________________________________________________________________________
# 3) PRIPRAVA PODATKOV
#______________________________________________________________________________

df_hist <- df %>%
  inner_join(
    selected_subsystems,
    by = "subsystem"
  ) %>%
  mutate(
    subsystem_si = factor(
      subsystem_si,
      levels = c(
        "Sinteza in hidroliza holesterol estrov",
        "Sestava proteinov",
        "Degradacija proteinov"
      )
    )
  )


#______________________________________________________________________________
# 4) GRAF
#______________________________________________________________________________

p_hist <- ggplot(
  df_hist,
  aes(x = distance_norm)
) +
  
  geom_histogram(
    binwidth = 0.5,
    boundary = 0,
    fill = "grey80",
    color = "grey35",
    linewidth = 0.35
  ) +
  
  facet_wrap(
    ~ subsystem_si,
    nrow = 1,
    scales = "fixed"
  ) +
  
  labs(
    x = "Normalizirana Evklidska razdalja log2TPM",
    y = "Število pacientov"
  ) +
  
  theme_minimal(
    base_size = 12,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0,
      margin = margin(b = 10)
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold",
      margin = margin(r = 8)
    ),
    
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    
    strip.text = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )


# Prikaz grafa
p_hist


#______________________________________________________________________________
# 5) DIMENZIJE KONČNEGA PNG
#______________________________________________________________________________

final_width  <- 11
final_height <- 5.5
final_dpi    <- 600


#______________________________________________________________________________
# 6) SHRANJEVANJE
#______________________________________________________________________________

ggsave(
  filename = "porazdelitev_izbrani_podsistemi_log2TPM.png",
  plot = p_hist,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = final_dpi,
  bg = "white",
  limitsize = FALSE
)
#========================================================
# HEATMAPS: log2TPM in PRETOKI
# Enak vrstni red pacientov v obeh grafih
#========================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(forcats)
library(scales)
library(viridis)

#--------------------------------------------------------
# 1) LOG2TPM
#--------------------------------------------------------

df_raw <- read.csv2(
  "distanca_pacientov_log2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

df <- df_raw %>%
  select(patient, subsystem, distance, n_genes, distance_norm) %>%
  mutate(
    patient = as.character(patient),
    subsystem = as.character(subsystem),
    distance_norm = as.numeric(str_replace(as.character(distance_norm), ",", "."))
  ) %>%
  filter(
    !is.na(patient),
    !is.na(subsystem),
    !is.na(distance_norm)
  )

# Število podsistemov za prikaz
top_n <- 10

# Prevodi za log2TPM - enaki kot v diplomi
subsystem_translation_log2 <- tribble(
  ~subsystem, ~subsystem_si,
  "Formation and hydrolysis of cholesterol esters", "Sinteza in hidroliza holesterol estrov",
  "Acylglycerides metabolism", "Metabolizem acilgliceridov",
  "Retinol metabolism", "Metabolizem retinola",
  "Protein assembly", "Sestava proteinov",
  "Protein modification", "Modifikacija proteinov",
  "Vitamin A metabolism", "Metabolizem vitamina A",
  "Starch and sucrose metabolism", "Metabolizem škroba in saharoze",
  "Protein degradation", "Degradacija proteinov",
  "Ascorbate and aldarate metabolism", "Metabolizem askorbata in aldarata",
  "Miscellaneous", "Ostale reakcije"
)

# Top 10 podsistemov po mediani normalizirane distance
top_subsystems_log2 <- df %>%
  group_by(subsystem) %>%
  summarise(
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance)) %>%
  slice_head(n = top_n) %>%
  pull(subsystem)

# Podatki za top 10 log2TPM
df_top10_log2 <- df %>%
  filter(subsystem %in% top_subsystems_log2) %>%
  left_join(subsystem_translation_log2, by = "subsystem") %>%
  mutate(
    subsystem_si = if_else(is.na(subsystem_si), subsystem, subsystem_si)
  )

# Vrstni red podsistemov za log2TPM
subsystem_order_log2 <- df_top10_log2 %>%
  group_by(subsystem_si) %>%
  summarise(
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance)) %>%
  pull(subsystem_si)

#--------------------------------------------------------
# 2) SKUPNI VRSTNI RED PACIENTOV
#    Tukaj ga določimo iz log2TPM top 10
#--------------------------------------------------------

patient_order_common <- df_top10_log2 %>%
  group_by(patient) %>%
  summarise(
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance)) %>%
  pull(patient)

#--------------------------------------------------------
# 3) HEATMAP ZA LOG2TPM
#--------------------------------------------------------

p_top10_log2 <- df_top10_log2 %>%
  mutate(
    subsystem_si = factor(
      subsystem_si,
      levels = subsystem_order_log2
    ),
    patient = factor(
      patient,
      levels = rev(patient_order_common)
    )
  ) %>%
  ggplot(
    aes(
      x = subsystem_si,
      y = patient,
      fill = distance_norm
    )
  ) +
  geom_tile(color = "black", linewidth = 0.15) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",
    name = "Normalizirana\nEvklidska\nrazdalja"
  ) +
  scale_x_discrete() +
  scale_y_discrete() +
  labs(
    x = "Metabolni podsistem",
    y = "Pacient",
  ) +
  theme_minimal(
    base_size = 12,
    base_family = "Helvetica"
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0
    ),
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 7
    ),
    axis.text.x = element_text(
      size = 9,
      face = "plain",
      colour = "black",
      angle = 55,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 5)
    ),
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

p_top10_log2

ggsave(
  filename = "heatmap_top10_podsistemi_log2_slovensko.jpg",
  plot = p_top10_log2,
  width = 11,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

#--------------------------------------------------------
# 4) PRETOKI
#--------------------------------------------------------

df_flux_raw <- read.csv2(
  "distanca_pacientov_flux.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

df_flux <- df_flux_raw %>%
  select(patient, subsystem, distance, distance_norm) %>%
  mutate(
    patient = as.character(patient),
    subsystem = as.character(subsystem),
    distance_norm = as.numeric(str_replace(as.character(distance_norm), ",", "."))
  ) %>%
  filter(
    !is.na(patient),
    !is.na(subsystem),
    !is.na(distance_norm)
  )

# Prevodi za pretoke - enaki kot v diplomi
subsystem_translation_flux <- tribble(
  ~subsystem, ~subsystem_si,
  "Pentose Phosphate Pathway", "Pentoza-fosfatna pot",
  "Pentose phosphate pathway", "Pentoza-fosfatna pot",
  "Glycolysis / Gluconeogenesis", "Glikoliza/glukoneogeneza",
  "Pyrimidine metabolism", "Metabolizem pirimidinov",
  "Pentose and glucuronate interconversions", "Pretvorbe pentoze in glukuronata",
  "Arginine and proline metabolism", "Metabolizem arginina in prolina",
  "Exchange/demand reactions", "Izmenjevalne reakcije",
  "Exchange reactions", "Izmenjevalne reakcije",
  "Exchange, demand and sink reactions", "Izmenjevalne reakcije",
  "Alanine, aspartate and glutamate metabolism", "Metabolizem alanina, aspartata in glutamata",
  "Miscellaneous", "Ostale reakcije",
  "Glycine, serine and threonine metabolism", "Metabolizem glicina, serina in treonina",
  "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
  "Cikel citronske kisline ter \nmetabolizem glioksilata in dikarboksilatov"
)

#--------------------------------------------------------
# 4) TOP 10 PODSISTEMOV ZA PRETOKE
#--------------------------------------------------------

# Top 10 podsistemov po mediani normalizirane distance

top_subsystems_flux <- df_flux %>%
  group_by(subsystem) %>%
  summarise(
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance)) %>%
  slice_head(n = top_n) %>%
  pull(subsystem)


# Podatki za top 10 pretokov

df_top10_flux <- df_flux %>%
  filter(subsystem %in% top_subsystems_flux) %>%
  left_join(subsystem_translation_flux, by = "subsystem") %>%
  mutate(
    subsystem_si = if_else(
      is.na(subsystem_si),
      subsystem,
      subsystem_si
    )
  )


# Vrstni red podsistemov za pretoke
# Najbolj spremenjen podsistem bo prvi na levi strani

subsystem_order_flux <- df_top10_flux %>%
  group_by(subsystem_si) %>%
  summarise(
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance)) %>%
  pull(subsystem_si)


#--------------------------------------------------------
# 5) KONTROLA: ali so pacienti enaki v obeh datasetih?
#--------------------------------------------------------

setdiff(
  unique(df_top10_log2$patient),
  unique(df_top10_flux$patient)
)

setdiff(
  unique(df_top10_flux$patient),
  unique(df_top10_log2$patient)
)

# Če oba izpisa vrneta character(0), je vse OK.


#--------------------------------------------------------
# 6) HEATMAP ZA PRETOKE
# Uporabimo ISTI patient_order_common kot pri log2TPM
#--------------------------------------------------------

p_top10_flux <- df_top10_flux %>%
  mutate(
    subsystem_si = factor(
      subsystem_si,
      levels = subsystem_order_flux
    ),
    patient = factor(
      patient,
      levels = rev(patient_order_common)
    )
  ) %>%
  ggplot(
    aes(
      x = subsystem_si,
      y = patient,
      fill = distance_norm
    )
  ) +
  geom_tile(
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",
    name = "Normalizirana\nEvklidska\nrazdalja"
  ) +
  scale_x_discrete() +
  scale_y_discrete() +
  labs(
    x = "Metabolni podsistem",
    y = "Pacient",
  ) +
  theme_minimal(
    base_size = 12,
    base_family = "Helvetica"
  ) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 7
    ),
    
    axis.text.x = element_text(
      size = 9,
      face = "plain",
      colour = "black",
      angle = 55,
      hjust = 1,
      vjust = 1,
      margin = margin(t = 5)
    ),
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      10, 20, 10, 10
    )
  )

p_top10_flux


#--------------------------------------------------------
# 7) SHRANJEVANJE JPG
#--------------------------------------------------------

ggsave(
  filename = "heatmap_top10_podsistemi_flux_slovensko.jpg",
  plot = p_top10_flux,
  width = 11,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
#združitev heatmapa tpm in pretokov
library(patchwork) #install.packages("patchwork")

p_combined <- (p_top10_log2 / p_top10_flux) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      size = 14,
      face = "bold",
      family = "Helvetica",
      colour = "black",
      hjust = 0,
      vjust = 0.5
    ),
    plot.tag.position = c(-0.01, 0.995)
  )

ggsave(
  filename = "heatmap_top10_log2_flux_skupaj.jpg",
  plot = p_combined,
  width = 12,
  height = 19,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)
p_combined
# =========================
# MA plot in volcano plot brez označenih genov
# =========================

library(dplyr)
library(ggplot2)
library(scales)

# -------------------------
# 1) Nastavitve pragov
# -------------------------

padj_cutoff <- 0.01
lfc_cutoff <- 1

# -------------------------
# 2) Branje celotnih DESeq2 rezultatov
# -------------------------

res_all <- read.csv2(
  "DESeq2_results_shrunk_all.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -------------------------
# 3) Priprava podatkov
# -------------------------

res_all <- res_all %>%
  mutate(
    baseMean = as.numeric(baseMean),
    log2FoldChange = as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    pvalue = as.numeric(pvalue),
    
    padj_plot = ifelse(is.na(padj), NA, pmax(padj, .Machine$double.xmin)),
    neg_log10_padj = -log10(padj_plot),
    
    significance = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "Povišano v tumorju",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Znižano v tumorju",
      TRUE ~ "Ni značilno"
    )
  ) %>%
  filter(!is.na(baseMean), baseMean > 0)


p_ma <- ggplot(res_all, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(
    aes(color = significance),
    alpha = 0.55,
    size = 1.1
  ) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  scale_x_log10(
    labels = label_scientific()
  ) +
  scale_color_manual(
    values = c(
      "Povišano v tumorju" = "firebrick3",
      "Znižano v tumorju" = "steelblue3",
      "Ni značilno" = "grey70"
    )
  ) +
  labs(
    title = "MA prikaz diferencialnega izražanja genov",
    subtitle = paste0(
      "Ob upoštevanju pogojev padj < ", padj_cutoff,
      " in |log2FC| ≥ ", lfc_cutoff
    ),
    x = "Povprečje normaliziranih odčitkov",
    y = "log2FC tumor / normalno",
    color = "Kategorija"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_ma

ggsave(
  "DESeq2_MA_plot.png",
  p_ma,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "DESeq2_MA_plot.pdf",
  p_ma,
  width = 8,
  height = 6
)



p_volcano <- ggplot(res_all, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(
    aes(color = significance),
    alpha = 0.55,
    size = 1.1
  ) +
  geom_vline(
    xintercept = c(-lfc_cutoff, lfc_cutoff),
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "black",
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = c(
      "Povišano v tumorju" = "firebrick3",
      "Znižano v tumorju" = "steelblue3",
      "Ni značilno" = "grey70"
    )
  ) +
  labs(
    x = "LFC",
    y = expression(-log[10]("padj_DESeq2")),
    color = "Kategorija"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p_volcano

ggsave(
  "DESeq2_volcano_plot.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "DESeq2_volcano_plot.pdf",
  p_volcano,
  width = 8,
  height = 6
)

DESeq2_rezultati <- read.csv2(
  "DESeq2_significant_genes_all_padj_log2fc.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#GRAF PORAZDELITVE VITAMINA A IN ŠKROBA/SAHAROZE
#========================================================
# PORAZDELITVE NORMALIZIRANIH DISTANC: prikaz IQR
# Y-os = število pacientov
# X-os = jasno določeni intervali normalizirane distance
#========================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)
library(scales)
library(tidyr)

#--------------------------------------------------------
# 1) Branje podatkov
#--------------------------------------------------------

df_raw <- read.csv2(
  "distanca_pacientov_log2.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

df <- df_raw %>%
  select(patient, subsystem, distance, n_genes, distance_norm) %>%
  mutate(
    patient = as.character(patient),
    subsystem = as.character(subsystem),
    distance_norm = as.numeric(str_replace(as.character(distance_norm), ",", "."))
  ) %>%
  filter(
    !is.na(patient),
    !is.na(subsystem),
    !is.na(distance_norm)
  )

#--------------------------------------------------------
# 2) Izbrana podsistema in slovenski prevodi
#--------------------------------------------------------

selected_subsystems <- tribble(
  ~subsystem, ~subsystem_si,
  "Vitamin A metabolism", "Metabolizem vitamina A",
  "Starch and sucrose metabolism", "Metabolizem škroba in saharoze"
)

df_iqr <- df %>%
  inner_join(selected_subsystems, by = "subsystem") %>%
  mutate(
    subsystem_si = factor(
      subsystem_si,
      levels = c(
        "Metabolizem vitamina A",
        "Metabolizem škroba in saharoze"
      )
    )
  )

#--------------------------------------------------------
# 3) Nastavitev širine intervala histograma
#--------------------------------------------------------

bin_width <- 0.25

x_max <- ceiling(max(df_iqr$distance_norm, na.rm = TRUE) / bin_width) * bin_width

x_breaks <- seq(
  from = 0,
  to = x_max,
  by = bin_width
)

#--------------------------------------------------------
# 4) Izračun Q1, mediane, Q3 in IQR
# type = 7 je primerljiv z Excelovo QUARTILE.INC funkcijo
#--------------------------------------------------------

iqr_stats <- df_iqr %>%
  group_by(subsystem_si) %>%
  summarise(
    n = n(),
    q1 = quantile(distance_norm, probs = 0.25, na.rm = TRUE, type = 7),
    mediana = median(distance_norm, na.rm = TRUE),
    q3 = quantile(distance_norm, probs = 0.75, na.rm = TRUE, type = 7),
    povprecje = mean(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    iqr = q3 - q1,
    label = paste0(
      "Q1 = ", comma(q1, decimal.mark = ",", accuracy = 0.001), "\n",
      "Mediana = ", comma(mediana, decimal.mark = ",", accuracy = 0.001), "\n",
      "Q3 = ", comma(q3, decimal.mark = ",", accuracy = 0.001), "\n",
      "IQR = ", comma(iqr, decimal.mark = ",", accuracy = 0.001)
    )
  )

iqr_lines <- iqr_stats %>%
  select(subsystem_si, q1, mediana, q3) %>%
  pivot_longer(
    cols = c(q1, mediana, q3),
    names_to = "statistika",
    values_to = "xintercept"
  ) %>%
  mutate(
    statistika = recode(
      statistika,
      q1 = "Q1",
      mediana = "Mediana",
      q3 = "Q3"
    ),
    statistika = factor(statistika, levels = c("Q1", "Mediana", "Q3"))
  )

#--------------------------------------------------------
# 5) Graf: histogram s številom pacientov
#--------------------------------------------------------

p_iqr_count <- ggplot(df_iqr, aes(x = distance_norm)) +
  
  # Zasenčeno območje IQR: od Q1 do Q3
  geom_rect(
    data = iqr_stats,
    aes(xmin = q1, xmax = q3, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = "lightblue",
    alpha = 0.25
  ) +
  
  # Histogram: y-os je število pacientov v intervalu
  geom_histogram(
    binwidth = bin_width,
    boundary = 0,
    closed = "left",
    fill = "grey85",
    color = "grey35",
    linewidth = 0.25
  ) +
  
  # Navpične črte za Q1, mediano in Q3
  geom_vline(
    data = iqr_lines,
    aes(xintercept = xintercept, linetype = statistika),
    linewidth = 0.65,
    color = "black"
  ) +
  
  # Oznake z vrednostmi Q1, mediana, Q3 in IQR
  geom_label(
    data = iqr_stats,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.05,
    size = 3,
    label.size = 0.25,
    fill = "white"
  ) +
  
  facet_wrap(
    ~ subsystem_si,
    ncol = 1,
    scales = "free_y"
  ) +
  
  scale_linetype_manual(
    values = c(
      "Q1" = "dashed",
      "Mediana" = "solid",
      "Q3" = "dashed"
    )
  ) +
  
  scale_x_continuous(
    breaks = x_breaks,
    labels = comma_format(decimal.mark = ",", accuracy = 0.01),
    limits = c(0, x_max),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_y_continuous(
    breaks = pretty_breaks(n = 6),
    expand = expansion(mult = c(0, 0.10))
  ) +
  
  labs(
    x = "Normalizirana Evklidska distanca",
    y = "Število pacientov",
    linetype = "Legenda"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0, margin = margin(b = 8)),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 8)),
    axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 9),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

p_iqr_count

#--------------------------------------------------------
# 6) Shranjevanje
#--------------------------------------------------------

ggsave(
  filename = "porazdelitev_IQR_vitaminA_skrob_sukroza_st_pacientov.pdf",
  plot = p_iqr_count,
  width = 9,
  height = 7,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "porazdelitev_IQR_vitaminA_skrob_sukroza_st_pacientov.png",
  plot = p_iqr_count,
  width = 9,
  height = 7,
  units = "in",
  dpi = 300
)

#___________________________________________
#štetje transkriptov v DESeq2
#________________________________________
library(dplyr)
library(stringr)

# Branje DESeq2 rezultatov
DESeq2_rezultati <- read.csv2(
  "DESeq2_significant_genes_all_padj_log2fc.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Pretvorba log2FoldChange v numeric, če slučajno ni že pravilno prebran
DESeq2_rezultati <- DESeq2_rezultati %>%
  mutate(
    log2FoldChange = as.numeric(str_replace(as.character(log2FoldChange), ",", "."))
  )

# Štetje vrstic/genov glede na smer log2FoldChange
stevilo_po_smeri <- DESeq2_rezultati %>%
  summarise(
    skupaj_genov = sum(!is.na(log2FoldChange)),
    pozitivni_log2FC = sum(log2FoldChange > 0, na.rm = TRUE),
    negativni_log2FC = sum(log2FoldChange < 0, na.rm = TRUE),
    log2FC_enak_0 = sum(log2FoldChange == 0, na.rm = TRUE),
    manjkajoci_log2FC = sum(is.na(log2FoldChange))
  )

cat("Skupaj diferencialno izraženih genov:", stevilo_po_smeri$skupaj_genov, "\n")
cat("Povišano izraženi v tumorju, log2FC > 0:", stevilo_po_smeri$pozitivni_log2FC, "\n")
cat("Znižano izraženi v tumorju, log2FC < 0:", stevilo_po_smeri$negativni_log2FC, "\n")

#______________________________________________________________
# TPM rezultati: združen prikaz dveh AUC grafov v enem PNG
#______________________________________________________________

library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

#--------------------------------------------------------
# 1) Branje podatkov
#--------------------------------------------------------

tpm_rezultati <- read.csv2(
  "vsi_podsistemi_gene_auc_log2.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

log2TPM_v2 <- read.csv2(
  "globalni_unikatni_geni_auc_wilcox_log2TPM.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_annotation <- read.csv2(
  "TPM_CPTAC_human1_annotation.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#--------------------------------------------------------
# 2) Povzetek števila genov
#--------------------------------------------------------

log2TPM_v2 %>%
  summarise(
    n_vseh_genov = n(),
    n_padj_001 = sum(wilcox_padj < 0.01, na.rm = TRUE),
    n_AUC_080 = sum(AUC >= 0.80, na.rm = TRUE),
    n_AUC_090 = sum(AUC >= 0.90, na.rm = TRUE),
    
    n_padj_001_AUC_080 = sum(
      wilcox_padj < 0.01 & AUC >= 0.80,
      na.rm = TRUE
    ),
    
    n_padj_001_AUC_090 = sum(
      wilcox_padj < 0.01 & AUC >= 0.90,
      na.rm = TRUE
    )
  )

#--------------------------------------------------------
# 3) Prvi graf: vsi unikatni geni
#--------------------------------------------------------

plot_data_global <- log2TPM_v2 %>%
  filter(
    !is.na(AUC),
    !is.na(wilcox_padj),
    wilcox_padj > 0
  ) %>%
  distinct(gene_id, AUC, wilcox_padj) %>%
  mutate(
    minus_log10_padj = -log10(wilcox_padj)
  )

p_auc_global <- ggplot(
  plot_data_global,
  aes(x = AUC, y = minus_log10_padj)
) +
  geom_point(alpha = 0.6, size = 1.2) +
  
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0.90,
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed"
  ) +
  
  coord_cartesian(xlim = c(0.5, 1)) +
  
  labs(
    x = "AUC",
    y = expression(-log[10](Wilcoxon_padj_TPM))
  ) +
  
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 20, 10, 10)
  )

#--------------------------------------------------------
# 4) Prevodi podsistemov
#--------------------------------------------------------

subsystem_translation_log2 <- tribble(
  ~subsystem, ~subsystem_si,
  "Formation and hydrolysis of cholesterol esters", "Sinteza in hidroliza holesterol estrov",
  "Acylglycerides metabolism", "Metabolizem acilgliceridov",
  "Retinol metabolism", "Metabolizem retinola",
  "Protein assembly", "Sestava proteinov",
  "Protein modification", "Modifikacija proteinov",
  "Vitamin A metabolism", "Metabolizem vitamina A",
  "Starch and sucrose metabolism", "Metabolizem škroba in saharoze",
  "Protein degradation", "Degradacija proteinov",
  "Ascorbate and aldarate metabolism", "Metabolizem askorbata in aldarata",
  "Miscellaneous", "Ostale reakcije"
)

#--------------------------------------------------------
# 5) Če je stolpec SUBSYSTEM, ga preimenuj
#--------------------------------------------------------

if ("SUBSYSTEM" %in% colnames(gene_annotation)) {
  gene_annotation <- gene_annotation %>%
    rename(subsystem = SUBSYSTEM)
}

#--------------------------------------------------------
# 6) Drugi graf: geni najbolj spremenjenih podsistemov
#--------------------------------------------------------

plot_data_subsystems <- gene_annotation %>%
  distinct(gene_id, subsystem) %>%
  inner_join(
    log2TPM_v2 %>% select(gene_id, AUC, wilcox_padj),
    by = "gene_id"
  ) %>%
  inner_join(
    subsystem_translation_log2,
    by = "subsystem"
  ) %>%
  filter(
    !is.na(AUC),
    !is.na(wilcox_padj),
    wilcox_padj > 0
  ) %>%
  mutate(
    minus_log10_padj = -log10(wilcox_padj),
    subsystem_si = factor(
      subsystem_si,
      levels = subsystem_translation_log2$subsystem_si
    )
  )

p_auc_subsystems <- ggplot(
  plot_data_subsystems,
  aes(x = AUC, y = minus_log10_padj, color = subsystem_si)
) +
  geom_point(alpha = 0.65, size = 1.6) +
  
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed"
  ) +
  
  geom_vline(
    xintercept = 0.90,
    linetype = "dotted"
  ) +
  
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed"
  ) +
  
  coord_cartesian(xlim = c(0.5, 1)) +
  
  labs(
    x = "AUC",
    y = expression(-log[10](Wilcox_padj_TPM)),
    color = "Metabolni podsistem"
  ) +
  
  theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 20, 10, 10)
  )

#--------------------------------------------------------
# 7) Združitev obeh grafov v eno sliko
#    (/ = eden pod drugim)
#--------------------------------------------------------

p_combined_auc <- (p_auc_global / p_auc_subsystems) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(
      size = 14,
      face = "bold",
      family = "Helvetica",
      colour = "black",
      hjust = 0
    ),
    plot.tag.position = c(-0.01, 0.99)
  )

# Prikaz
p_combined_auc

#--------------------------------------------------------
# 8) Nastavitve izvoza
#    Tukaj lahko prosto spreminjate končne dimenzije PNG
#--------------------------------------------------------

final_width  <- 11   # širina v inch
final_height <- 12   # višina v inch
final_dpi    <- 600  # ločljivost

#--------------------------------------------------------
# 9) Shranjevanje v PNG
#--------------------------------------------------------

ggsave(
  filename = "AUC_Wilcox_grafa_skupaj.png",
  plot = p_combined_auc,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = final_dpi,
  bg = "white",
  limitsize = FALSE
)



#prikaz grafa distribucije
library(dplyr)
library(ggplot2)

gen <- "ENSG00000117394"

plot_data <- CPTAC_log2 %>%
  select(id_pacienta, tip_vzorca, all_of(gen)) %>%
  filter(!is.na(.data[[gen]])) %>%
  rename(expr = all_of(gen))

ggplot(plot_data, aes(x = expr, fill = tip_vzorca)) +
  geom_histogram(
    position = "identity",
    alpha = 0.5,
    bins = 20
  ) +
  scale_fill_manual(values = c("normal" = "red", "tumor" = "blue")) +
  labs(
    title = paste("Porazdelitev izražanja gena", gen),
    x = "log2(TPM + 1)",
    y = "Število pacientov",
    fill = "Tip vzorca"
  ) +
  theme_minimal()

#grafi za pretoke
#______________________________________________________________________________
# FLUX REZULTATI: AUC IN GLOBALNO KORIGIRANE P-VREDNOSTI
# Združen prikaz vseh reakcij in reakcij top 10 podsistemov
#______________________________________________________________________________

library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)


#______________________________________________________________________________
# 1) NALOŽI REZULTATE REAKCIJ
#______________________________________________________________________________

reaction_flux <- read.csv2(
  "all_podsistemi_reaction_auc_flux.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


#______________________________________________________________________________
# 2) PREVODI TOP 10 PODSISTEMOV
# Povezava je narejena po ranku podsistema
#______________________________________________________________________________

subsystem_translation_flux <- tribble(
  ~Rank_flux_norm_distance, ~subsystem_si,
  1L,  "Pentoza-fosfatna pot",
  2L,  "Glikoliza/glukoneogeneza",
  3L,  "Metabolizem pirimidinov",
  4L,  "Pretvorbe pentoze in glukuronata",
  5L,  "Metabolizem arginina in prolina",
  6L,  "Izmenjevalne reakcije",
  7L,  "Metabolizem alanina, aspartata in glutamata",
  8L,  "Ostale reakcije",
  9L,  "Metabolizem glicina, serina in treonina",
  10L, "Cikel citronske kisline ter\nmetabolizem glioksilata in dikarboksilatov"
)


#______________________________________________________________________________
# 3) GRAF a): VSE UNIKATNE REAKCIJE
#______________________________________________________________________________

plot_flux_all <- reaction_flux %>%
  distinct(
    reaction_id,
    AUC,
    wilcox_padj
  ) %>%
  filter(
    !is.na(AUC),
    !is.na(wilcox_padj),
    wilcox_padj > 0
  ) %>%
  mutate(
    minus_log10_padj = -log10(wilcox_padj)
  )


graf_flux_all <- ggplot(
  plot_flux_all,
  aes(
    x = AUC,
    y = minus_log10_padj
  )
) +
  
  geom_point(
    alpha = 0.6,
    size = 1.6
  ) +
  
  # meja AUC = 0,80
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed"
  ) +
  
  # meja AUC = 0,90
  geom_vline(
    xintercept = 0.90,
    linetype = "dotted"
  ) +
  
  # meja prilagojene p-vrednosti = 0,01
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed"
  ) +
  
  coord_cartesian(
    xlim = c(0.5, 1)
  ) +
  
  labs(
    x = "AUC",
    y = expression(-log[10](Wilcoxon_padj_pretoki))
  ) +
  
  theme_bw(
    base_size = 11,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      10, 20, 10, 10
    )
  )


#______________________________________________________________________________
# 4) PRIPRAVA PODATKOV ZA GRAF b):
# REAKCIJE DESETIH NAJBOLJ SPREMENJENIH PODSISTEMOV
#______________________________________________________________________________

plot_flux_top10 <- reaction_flux %>%
  
  # obdržimo samo reakcije top 10 podsistemov
  filter(
    Rank_flux_norm_distance %in% 1:10
  ) %>%
  
  # dodamo slovenska imena podsistemov
  inner_join(
    subsystem_translation_flux,
    by = "Rank_flux_norm_distance"
  ) %>%
  
  # odstranimo podvojitve iste reakcije znotraj istega podsistema
  distinct(
    reaction_id,
    Rank_flux_norm_distance,
    subsystem_si,
    AUC,
    wilcox_padj
  ) %>%
  
  filter(
    !is.na(AUC),
    !is.na(wilcox_padj),
    wilcox_padj > 0
  ) %>%
  
  mutate(
    minus_log10_padj = -log10(wilcox_padj),
    
    subsystem_si = factor(
      subsystem_si,
      levels = subsystem_translation_flux$subsystem_si
    )
  )


#______________________________________________________________________________
# 5) BARVE ZA TOP 10 PODSISTEMOV
#______________________________________________________________________________

flux_colors <- c(
  
  "Pentoza-fosfatna pot" =
    "#E41A1C",
  
  "Glikoliza/glukoneogeneza" =
    "#377EB8",
  
  "Metabolizem pirimidinov" =
    "#4DAF4A",
  
  "Pretvorbe pentoze in glukuronata" =
    "#984EA3",
  
  "Metabolizem arginina in prolina" =
    "#FF7F00",
  
  "Izmenjevalne reakcije" =
    "#A65628",
  
  "Metabolizem alanina, aspartata in glutamata" =
    "#F781BF",
  
  "Ostale reakcije" =
    "#1B9E77",
  
  "Metabolizem glicina, serina in treonina" =
    "#D9A400",
  
  "Cikel citronske kisline ter\nmetabolizem glioksilata in dikarboksilatov" =
    "#666666"
)


#______________________________________________________________________________
# 6) GRAF b): REAKCIJE TOP 10 PODSISTEMOV
#______________________________________________________________________________

graf_flux_top10 <- ggplot(
  plot_flux_top10,
  aes(
    x = AUC,
    y = minus_log10_padj,
    color = subsystem_si
  )
) +
  
  geom_point(
    alpha = 0.75,
    size = 2
  ) +
  
  # meja AUC = 0,80
  geom_vline(
    xintercept = 0.80,
    linetype = "dashed"
  ) +
  
  # meja AUC = 0,90
  geom_vline(
    xintercept = 0.90,
    linetype = "dotted"
  ) +
  
  # meja prilagojene p-vrednosti = 0,01
  geom_hline(
    yintercept = -log10(0.01),
    linetype = "dashed"
  ) +
  
  scale_color_manual(
    values = flux_colors,
    drop = FALSE
  ) +
  
  coord_cartesian(
    xlim = c(0.5, 1)
  ) +
  
  labs(
    x = "AUC",
    y = expression(-log[10](Wilcox_padj_pretoki)),
    color = "Metabolni podsistem"
  ) +
  
  theme_bw(
    base_size = 11,
    base_family = "Helvetica"
  ) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 10,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    legend.title = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 11
    ),
    
    plot.margin = margin(
      10, 20, 10, 10
    )
  )


#______________________________________________________________________________
# 7) ZDRUŽITEV GRAFOV
#
# a) vse unikatne reakcije
# b) reakcije desetih najbolj spremenjenih podsistemov
#______________________________________________________________________________

graf_flux_combined <- (
  graf_flux_all /
    graf_flux_top10
) +
  
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  
  theme(
    plot.tag = element_text(
      size = 14,
      face = "bold",
      family = "Helvetica",
      colour = "black",
      hjust = 0
    ),
    
    plot.tag.position = c(
      -0.01,
      0.98
    )
  )


# Prikaz skupnega grafa
graf_flux_combined


#______________________________________________________________________________
# 8) NASTAVITVE KONČNEGA PNG
#
# Tukaj lahko poljubno spreminjate končne dimenzije.
#______________________________________________________________________________

final_width  <- 11
final_height <- 12
final_dpi    <- 600


#______________________________________________________________________________
# 9) SHRANJEVANJE
#______________________________________________________________________________

ggsave(
  filename = "flux_AUC_Wilcox_grafa_skupaj.png",
  plot = graf_flux_combined,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = final_dpi,
  bg = "white",
  limitsize = FALSE
)

#število reakcij v 4.5
reaction_flux %>%
  summarise(
    n_vseh_reakcij = n(),
    
    n_padj_001 = sum(
      wilcox_padj < 0.01,
      na.rm = TRUE
    ),
    
    n_AUC_080 = sum(
      AUC >= 0.80,
      na.rm = TRUE
    ),
    
    n_AUC_090 = sum(
      AUC >= 0.90,
      na.rm = TRUE
    ),
    
    n_padj_001_AUC_080 = sum(
      wilcox_padj < 0.01 &
        AUC >= 0.80,
      na.rm = TRUE
    ),
    
    n_padj_001_AUC_090 = sum(
      wilcox_padj < 0.01 &
        AUC >= 0.90,
      na.rm = TRUE
    )
  )

#INTEGRIRANI REZULTATI
#TABELA raznih parametrov in rezultatov
library(dplyr)

rezultati <- read.csv2(
  "DESeq2_z_log2TPM_z_flux.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# izbrani pragovi
padj_prag <- 0.01
AUC_prag  <- 0.80
LFC_prag  <- 1


pregled_stevila <- bind_rows(
  
  # ----------------------------------------------------------
  # DESeq2
  # ----------------------------------------------------------
  rezultati %>%
    filter(
      !is.na(padj_DESeq2),
      !is.na(log2FoldChange),
      padj_DESeq2 < padj_prag,
      abs(log2FoldChange) >= LFC_prag
    ) %>%
    summarise(
      Analiza = "DESeq2",
      Pogoji = "padj < 0,01; |log2FC| ≥ 1",
      `Število genov` = n_distinct(gene_id),
      `Število reakcij` = NA_integer_
    ),
  
  
  # ----------------------------------------------------------
  # log2TPM
  # ----------------------------------------------------------
  rezultati %>%
    filter(
      !is.na(wilcox_padj_log2TPM),
      !is.na(AUC_log2TPM),
      wilcox_padj_log2TPM < padj_prag,
      AUC_log2TPM >= AUC_prag
    ) %>%
    summarise(
      Analiza = "log2TPM",
      Pogoji = "Wilcox padj < 0,01; AUC ≥ 0,80",
      `Število genov` = n_distinct(gene_id),
      `Število reakcij` = NA_integer_
    ),
  
  
  # ----------------------------------------------------------
  # pretoki
  # ----------------------------------------------------------
  rezultati %>%
    filter(
      !is.na(wilcox_padj_flux),
      !is.na(AUC_flux),
      wilcox_padj_flux < padj_prag,
      AUC_flux >= AUC_prag
    ) %>%
    summarise(
      Analiza = "Pretoki",
      Pogoji = "Wilcox padj < 0,01; AUC ≥ 0,80",
      `Število genov` = NA_integer_,
      `Število reakcij` = n_distinct(
        reaction_id,
        na.rm = TRUE
      )
    ),
  
  
  # ----------------------------------------------------------
  # DESeq2 + log2TPM
  # ----------------------------------------------------------
  rezultati %>%
    filter(
      !is.na(padj_DESeq2),
      !is.na(log2FoldChange),
      !is.na(wilcox_padj_log2TPM),
      !is.na(AUC_log2TPM),
      
      padj_DESeq2 < padj_prag,
      abs(log2FoldChange) >= LFC_prag,
      wilcox_padj_log2TPM < padj_prag,
      AUC_log2TPM >= AUC_prag
    ) %>%
    summarise(
      Analiza = "DESeq2 + log2TPM",
      Pogoji = "DESeq2 + log2TPM pogoji",
      `Število genov` = n_distinct(gene_id),
      `Število reakcij` = NA_integer_
    ),
  
  
  # ----------------------------------------------------------
  # DESeq2 + log2TPM + pretoki
  # ----------------------------------------------------------
  rezultati %>%
    filter(
      !is.na(padj_DESeq2),
      !is.na(log2FoldChange),
      !is.na(wilcox_padj_log2TPM),
      !is.na(AUC_log2TPM),
      !is.na(wilcox_padj_flux),
      !is.na(AUC_flux),
      
      padj_DESeq2 < padj_prag,
      abs(log2FoldChange) >= LFC_prag,
      wilcox_padj_log2TPM < padj_prag,
      AUC_log2TPM >= AUC_prag,
      wilcox_padj_flux < padj_prag,
      AUC_flux >= AUC_prag
    ) %>%
    summarise(
      Analiza = "DESeq2 + log2TPM + pretoki",
      Pogoji = "Vsi trije sklopi pogojev",
      `Število genov` = n_distinct(gene_id),
      `Število reakcij` = n_distinct(
        reaction_id,
        na.rm = TRUE
      )
    )
)

print(pregled_stevila)

write.csv2(
  pregled_stevila,
  "pregled_stevila_kandidatov.csv",
  row.names = FALSE
)

#COUPLED GRAFI ZA PRIKAZ SMERI SPREMEMBE METABOLNIH PRETOKOV PRI TOP KANDIDATIH INTEGRIRANE TABELE
#______________________________________________________________________________
# PAIRED PLOTI METABOLNIH PRETOKOV ZA 4 IZBRANE REAKCIJE
#
# Rdeča = večji pretok v tumorskem vzorcu
# Modra = večji pretok v normalnem vzorcu
# Siva  = enak pretok
#
# Reakcije:
# MAR06628
# MAR06022
# MAR01976
# MAR04116
#______________________________________________________________________________

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)


#______________________________________________________________________________
# 1) BRANJE PODATKOV
#______________________________________________________________________________

flux_clean <- read.csv(
  "flux_clean.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


#______________________________________________________________________________
# 2) IZBRANE REAKCIJE
#______________________________________________________________________________

izbrane_reakcije <- c(
  "MAR06628",
  "MAR06022",
  "MAR01976",
  "MAR04116"
)


# Preveri, ali vse reakcije obstajajo v tabeli

manjkajoce_reakcije <- setdiff(
  izbrane_reakcije,
  colnames(flux_clean)
)

if (length(manjkajoce_reakcije) > 0) {
  
  stop(
    paste(
      "Naslednje reakcije manjkajo v flux_clean.csv:",
      paste(manjkajoce_reakcije, collapse = ", ")
    )
  )
}


#______________________________________________________________________________
# 3) PRIPRAVA PODATKOV
#______________________________________________________________________________

flux_selected <- flux_clean %>%
  
  select(
    patient_id,
    type,
    all_of(izbrane_reakcije)
  ) %>%
  
  pivot_longer(
    cols = all_of(izbrane_reakcije),
    names_to = "reaction_id",
    values_to = "flux"
  ) %>%
  
  filter(
    !is.na(patient_id),
    !is.na(type),
    !is.na(flux)
  )


#______________________________________________________________________________
# 4) OBDRŽIMO LE PACIENTE Z OBEMA VZORCEMA
#______________________________________________________________________________

flux_paired <- flux_selected %>%
  
  group_by(
    reaction_id,
    patient_id
  ) %>%
  
  filter(
    n_distinct(type) == 2
  ) %>%
  
  ungroup()


#______________________________________________________________________________
# 5) KONTROLA ŠTEVILA PACIENTOV
#______________________________________________________________________________

kontrola_parov <- flux_paired %>%
  
  group_by(reaction_id) %>%
  
  summarise(
    n_pacientov = n_distinct(patient_id),
    .groups = "drop"
  )

print(kontrola_parov)


#______________________________________________________________________________
# 6) DOLOČITEV SMERI SPREMEMBE ZA VSAKEGA PACIENTA
#______________________________________________________________________________

flux_smeri <- flux_paired %>%
  
  select(
    reaction_id,
    patient_id,
    type,
    flux
  ) %>%
  
  pivot_wider(
    names_from = type,
    values_from = flux
  ) %>%
  
  mutate(
    smer = case_when(
      
      tumor > normal ~
        "Povečan v tumorskem vzorcu",
      
      normal > tumor ~
        "Povečan v normalnem vzorcu",
      
      TRUE ~
        "Enak pretok"
    )
  )


#______________________________________________________________________________
# 7) POVZETEK SMERI SPREMEMB
#______________________________________________________________________________

povzetek_smeri <- flux_smeri %>%
  
  group_by(reaction_id) %>%
  
  summarise(
    
    n_pacientov = n(),
    
    n_povecan_tumor = sum(
      smer == "Povečan v tumorskem vzorcu"
    ),
    
    n_povecan_normal = sum(
      smer == "Povečan v normalnem vzorcu"
    ),
    
    n_enak = sum(
      smer == "Enak pretok"
    ),
    
    delez_povecan_tumor =
      100 * n_povecan_tumor / n_pacientov,
    
    .groups = "drop"
  )

print(povzetek_smeri)


#______________________________________________________________________________
# 8) SMER SPREMEMBE DODAMO NAZAJ POSAMEZNIM MERITVAM
#______________________________________________________________________________

flux_paired_plot <- flux_paired %>%
  
  left_join(
    
    flux_smeri %>%
      select(
        reaction_id,
        patient_id,
        smer
      ),
    
    by = c(
      "reaction_id",
      "patient_id"
    )
  ) %>%
  
  mutate(
    
    type = factor(
      type,
      levels = c(
        "normal",
        "tumor"
      ),
      labels = c(
        "Normalni",
        "Tumorski"
      )
    ),
    
    reaction_id = factor(
      reaction_id,
      levels = izbrane_reakcije
    ),
    
    smer = factor(
      smer,
      levels = c(
        "Povečan v tumorskem vzorcu",
        "Povečan v normalnem vzorcu",
        "Enak pretok"
      )
    )
  )


#______________________________________________________________________________
# 9) BARVE GLEDE NA SMER SPREMEMBE
#______________________________________________________________________________

barve_smeri <- c(
  
  "Povečan v tumorskem vzorcu" =
    "#D73027",
  
  "Povečan v normalnem vzorcu" =
    "#4575B4",
  
  "Enak pretok" =
    "#777777"
)


#______________________________________________________________________________
# 10) FUNKCIJA ZA POSAMEZEN PAIRED PLOT
#______________________________________________________________________________

naredi_paired_plot <- function(reaction) {
  
  podatki <- flux_paired_plot %>%
    filter(
      reaction_id == reaction
    )
  
  
  stats <- povzetek_smeri %>%
    filter(
      reaction_id == reaction
    )
  
  
  # Podnaslov posameznega panela
  
  podnaslov_reakcije <- paste0(
    "Večji pretok v tumorskem vzorcu pri ",
    stats$n_povecan_tumor,
    "/",
    stats$n_pacientov,
    " pacientih (",
    round(
      stats$delez_povecan_tumor,
      1
    ),
    " %)"
  )
  
  
  ggplot(
    podatki,
    aes(
      x = type,
      y = flux,
      group = patient_id,
      color = smer
    )
  ) +
    
    # povezava normalni -> tumorski
    
    geom_line(
      linewidth = 0.55,
      alpha = 0.65
    ) +
    
    # posamezne vrednosti
    
    geom_point(
      size = 1.8,
      alpha = 0.85
    ) +
    
    # barve glede na smer
    
    scale_color_manual(
      values = barve_smeri,
      drop = FALSE
    ) +
    
    labs(
      title = reaction,
      subtitle = podnaslov_reakcije,
      x = "Tip vzorca",
      y = "Ocenjeni metabolni pretok",
      color = "Smer spremembe"
    ) +
    
    theme_minimal(
      base_size = 12,
      base_family = "Helvetica"
    ) +
    
    theme(
      
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0
      ),
      
      plot.subtitle = element_text(
        size = 10,
        hjust = 0,
        margin = margin(
          b = 8
        )
      ),
      
      axis.title.x = element_text(
        size = 11,
        face = "bold",
        margin = margin(
          t = 7
        )
      ),
      
      axis.title.y = element_text(
        size = 11,
        face = "bold",
        margin = margin(
          r = 7
        )
      ),
      
      axis.text.x = element_text(
        size = 10,
        colour = "black"
      ),
      
      axis.text.y = element_text(
        size = 9,
        colour = "black"
      ),
      
      legend.title = element_text(
        size = 10,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
}


#______________________________________________________________________________
# 11) IZDELAVA POSAMEZNIH GRAFOV
#______________________________________________________________________________

p_MAR06628 <- naredi_paired_plot(
  "MAR06628"
)

p_MAR06022 <- naredi_paired_plot(
  "MAR06022"
)

p_MAR01976 <- naredi_paired_plot(
  "MAR01976"
)

p_MAR04116 <- naredi_paired_plot(
  "MAR04116"
)


#______________________________________________________________________________
# 12) UREDITEV LEGENDE
#
# Legendo odstranimo iz prvih treh grafov.
# Samo MAR04116 jo ohrani.
# Tako se lahko v končni figuri pojavi samo ena legenda.
#______________________________________________________________________________

p_MAR06628 <- p_MAR06628 +
  theme(
    legend.position = "bottom",
    legend.title = element_text(
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      size = 9
    )
  )

p_MAR06022 <- p_MAR06022 +
  theme(
    legend.position = "none"
  )

p_MAR01976 <- p_MAR01976 +
  theme(
    legend.position = "none"
  )

p_MAR04116 <- p_MAR04116 +
  theme(
    legend.position = "none"
  )



#______________________________________________________________________________
# 13) POSTAVITEV KONČNE FIGURE
#
# A B
# C D
# G G
#
# G = prostor za eno samo skupno legendo
#______________________________________________________________________________

layout_flux <- "
AB
CD
GG
"


p_flux_paired_combined <- wrap_plots(
  
  A = p_MAR06628,
  B = p_MAR06022,
  C = p_MAR01976,
  D = p_MAR04116,
  G = guide_area(),
  
  design = layout_flux,
  
  guides = "collect",
  
  heights = c(
    1,
    1,
    0.12
  )
  
) +
  
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) &
  
  theme(
    
    plot.tag = element_text(
      size = 14,
      face = "bold",
      family = "Helvetica",
      colour = "black",
      hjust = 0
    ),
    
    plot.tag.position = c(
      0.01,
      0.98
    )
  )


#______________________________________________________________________________
# 14) PRIKAZ KONČNE FIGURE
#______________________________________________________________________________

p_flux_paired_combined


#______________________________________________________________________________
# 15) DIMENZIJE KONČNEGA PNG
#
# Tukaj lahko spreminjate velikost slike.
#______________________________________________________________________________

final_width  <- 11
final_height <- 9
final_dpi    <- 600


#______________________________________________________________________________
# 16) SHRANJEVANJE
#______________________________________________________________________________

ggsave(
  filename = "paired_flux_4_reakcije_smer.png",
  plot = p_flux_paired_combined,
  width = final_width,
  height = final_height,
  units = "in",
  dpi = final_dpi,
  bg = "white",
  limitsize = FALSE
)