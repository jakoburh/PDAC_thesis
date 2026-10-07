# =========================
# 0) Namestitev paketov
# =========================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("apeglm", quietly = TRUE)) BiocManager::install("apeglm")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("biomaRt", quietly = TRUE)) {BiocManager::install("biomaRt")}
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)
library(DESeq2)
library(apeglm)
library(biomaRt)
library(dplyr)
library(tidyr) 
library(BiocManager) 
# =========================
# 1) Uvoz podatkov in urejanje podatkov
# =========================
#združevanje podatkov iz PAAD in CPTAC
paad <- read.csv("PAAD_raw_counts_transposed.csv", row.names = 1, check.names = FALSE)
cptac <- read.csv2("CPTAC-3_raw_counts_transposed.csv", row.names = 1, check.names = FALSE)

rownames(paad) <- trimws(rownames(paad))
rownames(cptac) <- trimws(rownames(cptac))

common_genes <- intersect(rownames(paad), rownames(cptac))

paad2 <- paad[common_genes, , drop = FALSE]
cptac2 <- cptac[common_genes, , drop = FALSE]

all_counts <- cbind(paad2, cptac2)

all_counts[] <- lapply(all_counts, function(x) as.numeric(as.character(x)))

# pretvori v data.frame in shrani originalne gene ID-je v stolpec
all_counts <- as.data.frame(all_counts)
all_counts$gene_id <- rownames(all_counts)

# odstrani verzije genov v stolpcu, ne v rownames
all_counts$gene_id <- sub("\\..*", "", all_counts$gene_id)

# združi podvojene gene po odstranitvi verzij
all_counts <- all_counts %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum), .groups = "drop")

# gene_id nazaj v rownames
all_counts <- as.data.frame(all_counts)
rownames(all_counts) <- all_counts$gene_id
all_counts$gene_id <- NULL

write.csv2(all_counts, "all_raw_counts.csv", row.names = TRUE)
# BRANJE METAPODATKOV IN COUNT PODATKOV metadata je ločena z ;
colData <- read.csv("all_DESeq2_cisti.csv", sep = ";", stringsAsFactors = FALSE)

# raw counts
counts_data <- read.csv2("all_raw_counts.csv", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) #row.names = 1 zagotovi da se tisti prazen kvadrat, kjer se nahajajo imena stolpcev, ne šteje kot ime vrstice

rownames(colData) <- colData$sample
rownames(colData)
all(colnames(counts_data) %in% rownames(colData)) #preverimo da so vse vrstice v colData v stolpcih counts_data
colData <- colData[colnames(counts_data), ] #naredimo vrstni red teh vrstic enak stolpcem
all(colnames(counts_data) == rownames(colData)) #preverimo če vrstni red štima

# =========================
# 2) Naredimo matriko in izpeljemo izračune
# =========================
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ database + Type)
#prefiltriranje- odstranitev tistih ki imajo manj kot 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds [keep,]

#za določanje factor level- da se primerja zdrav vzorec s tumorskim
dds$Type <- relevel(dds$Type, ref ="normal")
dds$database <- relevel(dds$database, ref = "CPTAC-3")   # ali "TCGA"
#da zaženemo DESeq
dds <- DESeq(dds)
#potreben tudi collapse technical replicates, če jih imamo
#technical replicate so counti izmerjeni na istem vzorcu da preverimo machine nois, pcr bias, sequencing noise...
#tabela vrednosti padj, mean, log2foldchange...
# =========================
# 3)rezultati in analiza
# =========================
res <- results(dds, contrast = c("Type","tumor","normal"))
summary(res)


res001 <- results(dds, alpha = 0.01) #pri p-vrednosti manj kot 0,01
summary(res001)
#table(res$padj < 0.05, useNA="ifany")
#hist(res$pvalue, breaks=50)
#plotMA(res2)
resultsNames(dds)

res_shrunk <- lfcShrink( dds,coef = "Type_tumor_vs_normal",type = "apeglm") #prilagodi log2FoldChange za gene z malo podatki ali z veliko variance, kar povzroči nestabilen log2FoldChange
summary(res_shrunk)
# =========================
# 4) Osnovni pregledi
# =========================
table(res$padj < 0.01, useNA = "ifany")
hist(res$pvalue, breaks = 50, main = "Histogram of p-values", xlab = "p-value")
plotMA(res_shrunk, ylim = c(-5, 5))


vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "Type") #graf povprečja normaliziranih countov in log2FC

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_shrunk_df <- as.data.frame(res_shrunk)
res_shrunk_df$gene_id <- rownames(res_shrunk_df)

sig_res <- subset(
  res_shrunk_df,
  !is.na(padj) & padj < 0.01 & abs(log2FoldChange) >= 1
)
table(sig_res$padj < 0.01, useNA = "ifany")
hist(sig_res$pvalue, breaks = 50, main = "Histogram of p-values", xlab = "p-value")
plotMA(res_shrunk, ylim = c(-5, 5))
write.csv2(res_df, "DESeq2_results_raw_all.csv", row.names = FALSE)
write.csv2(res_shrunk_df, "DESeq2_results_shrunk_all.csv", row.names = FALSE)
write.csv2(sig_res, "DESeq2_significant_genes_all.csv", row.names = FALSE)




# -----------------------------
# 1) Razvrščanje po najnižji padj
# -----------------------------
sig_by_padj <- sig_res[order(sig_res$padj), ]

# -----------------------------
# 2) Razvrščanje po najvišjem |log2FC|
# -----------------------------
sig_by_log2fc <- sig_res[order(-abs(sig_res$log2FoldChange)), ]

# -----------------------------
# 3) Kombinacija: nizek padj + visok |log2FC|
# score = -log10(padj) * |log2FC|
# -----------------------------
sig_res$combined_score <- -log10(sig_res$padj) * abs(sig_res$log2FoldChange)
sig_by_combined <- sig_res[order(-sig_res$combined_score), ]

# -----------------------------
# po želji samo prvih 20
# -----------------------------
top20_padj <- head(sig_by_padj, 20)
top20_log2fc <- head(sig_by_log2fc, 20)
top20_combined <- head(sig_by_combined, 20)

# pregled v R
head(top20_padj, 20)
head(top20_log2fc, 20)
head(top20_combined, 20)

# -----------------------------
# izvoz v csv
# -----------------------------
write.csv2(sig_by_padj, "DESeq2_significant_genes_all_padj.csv", row.names = FALSE)
write.csv2(sig_by_log2fc, "DESeq2_significant_genes_all_log2fc.csv", row.names = FALSE)
write.csv2(sig_by_combined, "DESeq2_significant_genes_all_padj_log2fc.csv", row.names = FALSE)

write.csv2(top20_padj, "top20_genes_by_padj.csv", row.names = FALSE)
write.csv2(top20_log2fc, "top20_genes_by_log2fc.csv", row.names = FALSE)
write.csv2(top20_combined, "top20_genes_by_combined_score.csv", row.names = FALSE)

#dodajanje angleških imen namesto ENSG


# =========================
# 5) Anotacija genov z biomaRt
# =========================

# preberi datoteko s signifikantnimi kandidati
podatki <- read.csv2(
  "DESeq2_significant_genes_all_padj_log2fc.csv",
  stringsAsFactors = FALSE
)

# poskrbi, da je gene_id tekst
podatki$gene_id <- as.character(podatki$gene_id)

# če bi se slučajno še kje pojavile verzije ENSG, jih odstranimo
podatki$gene_id <- sub("\\..*", "", podatki$gene_id)

# unikatni geni za poizvedbo v Ensembl
genes <- unique(podatki$gene_id)

# povezava na Ensembl preko biomaRt
mart <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# pridobi anotacije
annotations <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
  filters = "ensembl_gene_id",
  values = genes,
  mart = mart
)

# združi DESeq2 rezultate z anotacijami
koncni_kandidati_anotirani <- merge(
  podatki,
  annotations,
  by.x = "gene_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)

# razvrsti po combined score
koncni_kandidati_anotirani <- koncni_kandidati_anotirani[
  order(-koncni_kandidati_anotirani$combined_score),
]

# izvozi
write.csv2(
  koncni_kandidati_anotirani,
  "anotirani_kandidati.csv",
  row.names = FALSE
)

# =========================
# 6) Dodajanje podsistemov iz Human1
# =========================
koncni_kandidati_anotirani <- read.csv2("anotirani_kandidati.csv")

# preberi Human-GEM
human <- readxl::read_excel("Human-GEM.xlsx")

# izberi samo pomembna stolpca
human_sub <- human[, c("GENE ASSOCIATION", "SUBSYSTEM")]

# preimenuj stolpca, da se izogneš težavam s presledki
colnames(human_sub) <- c("GENE_ASSOCIATION", "SUBSYSTEM")

# odstrani oklepaje iz GPR zapisov
human_sub$GENE_ASSOCIATION <- gsub("\\(|\\)", "", human_sub$GENE_ASSOCIATION)

# zamenjaj "and" in "or" z ";" da lahko gene razbijemo v posamezne vrstice
human_sub$GENE_ASSOCIATION <- gsub(" and | or ", ";", human_sub$GENE_ASSOCIATION)

# odstrani odvečne presledke
human_sub$GENE_ASSOCIATION <- trimws(human_sub$GENE_ASSOCIATION)

# razbij GPR v posamezne gene
gene_map <- human_sub %>%
  tidyr::separate_rows(GENE_ASSOCIATION, sep = ";") %>%
  dplyr::rename(gene = GENE_ASSOCIATION)

# počisti presledke še enkrat
gene_map$gene <- trimws(gene_map$gene)

# odstrani prazne gene, če obstajajo
gene_map <- gene_map[gene_map$gene != "" & !is.na(gene_map$gene), ]

# združi anotirane DESeq2 kandidate s Human1 mapo
# POZOR: tukaj mapiramo po po ENSG ID
anotirani_kandidati_human1 <- merge(
  koncni_kandidati_anotirani,
  gene_map,
  by.x = "gene_id",
  by.y = "gene",
  all.x = TRUE
)

# razvrsti po combined score
anotirani_kandidati_human1 <- anotirani_kandidati_human1[
  order(-anotirani_kandidati_human1$combined_score),
]

# izvozi
write.csv2(
  anotirani_kandidati_human1,
  "anotirani_kandidati_human1.csv",
  row.names = FALSE
)


#____________________________________________________________________
#dodatna analiza rezultatov
#____________________________________________________________________
library(dplyr)
library(ggplot2)

# preberi podatke
counts_data <- read.csv2("all_raw_counts.csv",
                         row.names = 1,
                         check.names = FALSE,
                         stringsAsFactors = FALSE)

colData <- read.csv("all_DESeq2_cisti.csv",
                    sep = ";",
                    stringsAsFactors = FALSE)

# uredi colData
rownames(colData) <- colData$sample

# preveri usklajenost
all(colnames(counts_data) %in% rownames(colData))
colData <- colData[colnames(counts_data), ]
all(colnames(counts_data) == rownames(colData))

# izberi gen
gene_id <- "ENSG00000118137"

# preveri, če gen obstaja
gene_id %in% rownames(counts_data)

# vzemi counte za ta gen
gene_counts <- as.numeric(counts_data[gene_id, ])

# naredi tabelo za ggplot
plot_df <- data.frame(
  sample = colnames(counts_data),
  count = gene_counts,
  log2_count = log2(gene_counts + 1),
  Type = colData$Type,
  database = colData$database
)

#grafi brez DESeq2

ggplot(plot_df, aes(x = Type, y = count, color = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = paste("Raw counts za gen", gene_id),
    x = "Tip vzorca",
    y = "Raw counts"
  ) +
  theme_minimal()

ggplot(plot_df, aes(x = Type, y = log2_count, color = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = paste("log2(raw count + 1) za gen", gene_id),
    x = "Tip vzorca",
    y = "log2(count + 1)"
  ) +
  theme_minimal()

ggplot(plot_df, aes(x = log2_count, fill = Type)) +
  geom_density(alpha = 0.4) +
  labs(
    title = paste("Gostota porazdelitve za gen", gene_id),
    x = "log2(count + 1)",
    y = "Gostota"
  ) +
  theme_minimal()

#grafi z DESeq2

norm_counts <- counts(dds, normalized = TRUE)

plot_df_norm <- data.frame(
  sample = colnames(norm_counts),
  norm_count = as.numeric(norm_counts[gene_id, ]),
  log2_norm_count = log2(as.numeric(norm_counts[gene_id, ]) + 1),
  Type = colData$Type,
  database = colData$database
)

ggplot(plot_df_norm, aes(x = Type, y = log2_norm_count, color = Type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  labs(
    title = paste("log2(normalized count + 1) za gen", gene_id),
    x = "Tip vzorca",
    y = "log2(normalized count + 1)"
  ) +
  theme_minimal()

ggplot(plot_df_norm, aes(x = log2_norm_count, fill = Type)) +
  geom_density(alpha = 0.4) +
  labs(
    title = paste("Gostota porazdelitve (log2 normalized counts) za", gene_id),
    x = "log2(normalized count + 1)",
    y = "Gostota"
  ) +
  theme_minimal()

#izračun auc
library(pROC)
library(dplyr)

gene_id <- "ENSG00000118137"

# normalized counts iz DESeq2
norm_counts <- counts(dds, normalized = TRUE)

# dataframe
auc_df <- data.frame(
  count = as.numeric(norm_counts[gene_id, ]),
  log2_count = log2(as.numeric(norm_counts[gene_id, ]) + 1),
  Type = colData$Type
)

auc_df$Type <- factor(auc_df$Type, levels = c("normal", "tumor"))

roc_obj <- roc(auc_df$Type, auc_df$log2_count)

auc(roc_obj)

#POGLED DOLOČENEGA GENA, KI NI DIFERENCIALNO IZRAŽEN PO DESeq2
gene_of_interest <- "ENSG00000132840"

gene_of_interest %in% rownames(counts_data) #ali je ta gen sploh v podatkih?

res_shrunk_df[res_shrunk_df$gene_id == gene_of_interest, ]