#__________________________________________
# namestitev knjižnic
#__________________________________________
library(dplyr)
library(tidyr)
library(readxl)
#__________________________________________
# 1. Priprava vhodnih fluxov za obdelavo
#__________________________________________
flux_mapa <- "C:/Users/Nitro 17/Desktop/diplomska naloga/Jakobb/flux/data"

#meta_podatki vsebujejo  podatke o patient_id, sample_id in type. Cilj je vsaki datoteki s fluxi pripisati te podatke, kjer bodo v vrsticah imena reakcij
meta_podatki <- read.csv2("meta_podatki.csv", stringsAsFactors = FALSE)

# ustvari standardiziran sample_id za lažje ujemanje z imeni datotek
meta_podatki$sample_id_std <- gsub("-", "_", meta_podatki$sample_id)

# preveri, ali meta_podatki vsebuje pričakovane stolpce
if (!all(c("sample_id", "type", "patient_id") %in% colnames(meta_podatki))) {stop("meta_podatki nima pričakovanih stolpcev: sample_id, type, patient_id")}

#prebere imena datotek v tej mapi 
flux_datoteke <- list.files(path = flux_mapa,pattern = "\\.csv$", full.names = TRUE)

# če je meta_podatki.csv slučajno v isti mapi, ga odstrani iz seznama
flux_datoteke <- flux_datoteke[basename(flux_datoteke) != "meta_podatki.csv"]

# kontrolni izpis
cat("Najdenih vseh csv datotek v flux mapi:", length(flux_datoteke), "\n")
cat("Primer imen datotek:\n")
print(head(basename(flux_datoteke), 10))

vsi_fluxi <- list()                 # sem bomo shranili uspešno prebrane in ujemajoče se vzorce
manjkajoci_meta <- character()      # datoteke, ki nimajo ujemanja v meta_podatki
napacne_datoteke <- character()     # datoteke, ki nimajo stolpcev reaction in flux ali se ne preberejo
vse_reakcije <- character()         # unija vseh reakcij čez vse veljavne vzorce

for (trenutna_datoteka in flux_datoteke) {
  sample_id <- tools::file_path_sans_ext(basename(trenutna_datoteka))   # iz imena datoteke brez .csv dobi sample_id
  meta_vrstica <- meta_podatki %>% filter(sample_id_std == !!sample_id) #poiščeš ujemajočo vrstico imena datoteke z vrstico v meta podatkih
  if (nrow(meta_vrstica) == 0) { #če ni ujemanja med imenom datoteke in meta podatki, preskoči to datoteko
    manjkajoci_meta <- c(manjkajoci_meta, sample_id)   
    next}
  if (nrow(meta_vrstica) > 1) { # če se isti sample_id pojavi večkrat v meta_podatki, ustavi program
    stop(paste("Sample_id", sample_id, "se v meta_podatki pojavi večkrat.")) }
  flux_df <- try(read.csv(trenutna_datoteka, stringsAsFactors = FALSE), silent = TRUE) #preberi datoteko ki se pojavi v metab podatkih
  if (inherits(flux_df, "try-error")) {  # če datoteke ni bilo mogoče pravilno prebrati, jo zabeleži in preskoči
    napacne_datoteke <- c(napacne_datoteke, sample_id)
    next
  }
  if (!all(c("reaction", "flux") %in% colnames(flux_df))) {   # preveri, ali datoteka vsebuje stolpca reaction in flux
    napacne_datoteke <- c(napacne_datoteke, sample_id)
    next
  }
  if (any(duplicated(flux_df$reaction))) { # preveri, ali so v datoteki podvojene reakcije
    stop(paste("Datoteka", sample_id, "vsebuje podvojene reakcije."))
  }
  
  vsi_fluxi[[sample_id]] <- list( # shrani podatke trenutnega vzorca v seznam
    sample_id = sample_id,
    type = meta_vrstica$type,
    patient_id = meta_vrstica$patient_id,
    flux_df = flux_df
  )
  vse_reakcije <- union(vse_reakcije, flux_df$reaction)  # dodaj reakcije trenutnega vzorca v skupni seznam vseh reakcij
}

#hitro vmesno preverjanje rezultatov
cat("\nŠtevilo uspešno prebranih in ujemajočih se vzorcev:", length(vsi_fluxi), "\n") #cat je le lepši print funkcija, le da je potrebno roćno dodajati nove vrsti z \n 
cat("\nŠtevilo datotek brez ujemanja v meta_podatki:", length(manjkajoci_meta), "\n")
if (length(manjkajoci_meta) > 0) {
  cat("Primer datotek brez ujemanja:\n")
  print(head(manjkajoci_meta, 20))
}
cat("\nŠtevilo napačnih / neberljivih flux datotek:", length(napacne_datoteke), "\n")
if (length(napacne_datoteke) > 0) {
  cat("Primer napačnih datotek:\n")
  print(head(napacne_datoteke, 20))
}
# če ni nobene veljavne datoteke, ustavi program
if (length(vsi_fluxi) == 0) {
  stop("Noben flux vzorec ni bil uspešno prebran in usklajen z meta_podatki.")
}

vse_reakcije <- sort(vse_reakcije) # uredi vse reakcije po abecedi za bolj pregledno matriko
cat("\nSkupno število unikatnih reakcij v veljavnih vzorcih:", length(vse_reakcije), "\n")

seznam_vrstic <- list()
#za vsako datoteko naredimo prazen vector z NA -> dodamo vrednostifluxov na tista mesta, kjer obstajajo -> pretvorimo vector v dataframe -> dodamo meta podatke -> shranimo kot vrstico končne tabele -> ponovimo za vse datoteke
for (ime_vzorca in names(vsi_fluxi)) {
  podatki <- vsi_fluxi[[ime_vzorca]] # vzemi podatke trenutnega vzorca
  flux_df <- podatki$flux_df
  poln_flux <- rep(NA, length(vse_reakcije))# ustvari prazen named vector z NA za vse reakcije (vector= enodimenzionalna tabel. npr poln_flux(NA, NA, NA...) velikosti dolžine vseh reakcij)
  names(poln_flux) <- vse_reakcije #za vsak sample naredimo ta prazen vector
  poln_flux[flux_df$reaction] <- flux_df$flux # vstavi obstoječe flux vrednosti na prava mesta v praznem vectorju
  ena_vrstica <- as.data.frame(t(poln_flux)) # pretvori v eno vrstico
  # dodaj meta stolpce
  ena_vrstica$sample_id <- podatki$sample_id
  ena_vrstica$type <- podatki$type
  ena_vrstica$patient_id <- podatki$patient_id
  ena_vrstica <- ena_vrstica %>% # preuredi stolpce, da so meta podatki na začetku
    select(sample_id, type, patient_id, everything())
  seznam_vrstic[[length(seznam_vrstic) + 1]] <- ena_vrstica # shrani vrstico
}
if (length(seznam_vrstic) == 0) { # preveri, ali seznam_vrstic ni prazen
  stop("Seznam vrstic je prazen. Noben veljaven vzorec ni bil pretvorjen v matriko.")
}

#preverjanje končne flux tabele
flux_matrika <- bind_rows(seznam_vrstic)
cat("\nDimenzije končne flux matrike:\n")
print(dim(flux_matrika))
cat("\nPreview:\n")
if (nrow(flux_matrika) > 0 && ncol(flux_matrika) > 0) {
  print(flux_matrika[1:min(5, nrow(flux_matrika)), 1:min(10, ncol(flux_matrika))])
} else {
  cat("flux_matrika je prazna.\n")
}
write.csv(flux_matrika, "flux_matrika_vsi_vzorci_NA.csv", row.names = FALSE)
cat("\nDatoteka shranjena kot: flux_matrika_vsi_vzorci_NA.csv\n")

#pretvorba NA v 0
flux_matrika[is.na(flux_matrika)] <- 0

library(dplyr)
library(tidyr)

#povprečenje fluxov tistih vzorcev z istim patient_id in type

# preštej vzorce na pacienta
patient_counts_flux <- flux_matrika %>%
  count(patient_id, type) %>%
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  )

# obdrži samo paciente, ki imajo vsaj en normal in vsaj en tumor
valid_patients_flux <- patient_counts_flux %>%
  filter(normal >= 1, tumor >= 1) %>%
  pull(patient_id)

# filtriraj samo veljavne paciente
flux_matrika_filtered <- flux_matrika %>%
  filter(patient_id %in% valid_patients_flux)

# imena reakcijskih stolpcev
reaction_cols <- colnames(flux_matrika_filtered)[4:ncol(flux_matrika_filtered)]

# združi po pacientu in tipu vzorca kot povprečje flux vrednosti
flux_clean <- flux_matrika_filtered %>%
  group_by(patient_id, type) %>%
  summarise(
    across(
      all_of(reaction_cols),
      ~ mean(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

# hiter izpis rezultatov
cat("Dimenzije očiščene flux matrike:\n")
print(dim(flux_clean))

cat("\nPreview:\n")
print(flux_clean[1:min(5, nrow(flux_clean)), 1:min(10, ncol(flux_clean))])

# shrani rezultat
write.csv(flux_clean, "flux_clean.csv", row.names = FALSE)
cat("\nDatoteka shranjena kot: flux_clean.csv\n")

#
human1 <- read_excel("Human-GEM.xlsx")

#_________________________________________________
#2. izračun evklidove distance fluxov za vsak podsistem (analogno log2TPM pipeline)
#_________________________________________________

library(dplyr)
library(ggplot2)

# reaction matrika (brez metadata stolpcev)
reaction_matrix <- flux_clean[, 3:ncol(flux_clean)]

# pripravi mapping reakcija -> podsistem
reaction_subsystem <- human1 %>%
  select(ID, SUBSYSTEM) %>%
  filter(!is.na(ID), !is.na(SUBSYSTEM)) %>%
  distinct()

# seznam reakcij po podsistemih
subsystem_list <- split(reaction_subsystem$ID, reaction_subsystem$SUBSYSTEM)

# unikatni pacienti
patients <- unique(flux_clean$patient_id)

# prazen dataframe za rezultate
results_flux <- data.frame()

for (p in patients) {
  
  # poišči vrstico normalnega in tumorskega vzorca za tega pacienta
  normal_row <- which(flux_clean$patient_id == p & flux_clean$type == "normal")
  tumor_row  <- which(flux_clean$patient_id == p & flux_clean$type == "tumor")
  
  # nadaljuj samo, če ima pacient točno en normal in en tumor
  if (length(normal_row) == 1 & length(tumor_row) == 1) {
    
    for (sub in names(subsystem_list)) {
      
      # reakcije v tem podsistemu
      reactions <- subsystem_list[[sub]]
      
      # obdrži samo reakcije, ki obstajajo v flux_clean
      reactions <- reactions[reactions %in% colnames(reaction_matrix)]
      
      # računaj samo, če sta v podsistemu vsaj 2 reakciji
      if (length(reactions) > 1) {
        
        normal_vec <- as.numeric(reaction_matrix[normal_row, reactions])
        tumor_vec  <- as.numeric(reaction_matrix[tumor_row, reactions])
        
        # evklidova distanca
        d <- sqrt(sum((tumor_vec - normal_vec)^2))
        
        results_flux <- rbind(
          results_flux,
          data.frame(
            patient = p,
            subsystem = sub,
            distance = d,
            n_reactions = length(reactions)
          )
        )
      }
    }
  }
}

# sortiraj po absolutni distanci
results_flux <- results_flux[order(-results_flux$distance), ]

# normalizacija glede na velikost podsistema
results_flux$distance_norm <- results_flux$distance / sqrt(results_flux$n_reactions)

# povzetek po podsistemih
summary_results_flux <- results_flux %>%
  group_by(subsystem) %>%
  summarise(
    mean_distance = mean(distance_norm, na.rm = TRUE),
    median_distance = median(distance_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_distance))

# graf top 15 podsistemov
ggplot(summary_results_flux[1:min(15, nrow(summary_results_flux)), ],
       aes(x = reorder(subsystem, median_distance),
           y = median_distance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    x = "Subsystem",
    y = "Median tumor-normal flux distance",
    title = "Top podsistemi po mediani normalizirane flux distance"
  )

# sortiraj še po normalizirani distanci
results_flux <- results_flux[order(-results_flux$distance_norm), ]
summary_results_flux <- summary_results_flux[order(-summary_results_flux$median_distance), ]

# shrani rezultate
write.csv2(results_flux, "distanca_pacientov_flux.csv", row.names = FALSE)
write.csv2(summary_results_flux, "mediana_distanc_podsistemov_flux.csv", row.names = FALSE)

# za nadaljnjo uporabo
mediana_podsistemov_flux <- summary_results_flux

#_____________________________________________________________________
#3. izračun AUC za top 10 podsistemov glede na normalizirano distanco fluxov
#______________________________________________________________________

library(dplyr)
library(pROC)

# 1) namesto n() vstavimo koliko top podsistemov želimo nadaljno analizirati, če se pusti n, se vzamejo vsi podsistemi
top10_subsystems_flux <- mediana_podsistemov_flux %>%
  arrange(desc(median_distance)) %>%
  slice(1:n()) %>%
  pull(subsystem)

# 2) reakcije iz teh podsistemov
top_reactions_df <- human1 %>%
  select(ID, SUBSYSTEM) %>%
  filter(!is.na(ID), !is.na(SUBSYSTEM)) %>%
  distinct() %>%
  filter(SUBSYSTEM %in% top10_subsystems_flux)

# 3) priprava labela
flux_auc <- flux_clean %>%
  filter(type %in% c("normal", "tumor")) %>%
  mutate(label = ifelse(type == "tumor", 1, 0))

# 4) izračun AUC za vsako reakcijo
reaction_auc_results <- data.frame()

for (i in seq_len(nrow(top_reactions_df))) {
  
  sub_i <- top_reactions_df$SUBSYSTEM[i]
  reaction_i <- top_reactions_df$ID[i]
  
  if (reaction_i %in% colnames(flux_auc)) {
    
    tmp <- flux_auc %>%
      select(patient_id, type, label, all_of(reaction_i)) %>%
      filter(!is.na(.data[[reaction_i]]))
    
    # AUC rabi vsaj 2 razreda
    if (length(unique(tmp$label)) == 2) {
      
      roc_obj <- roc(
        response = tmp$label,
        predictor = tmp[[reaction_i]],
        quiet = TRUE
      )
      
      auc_val <- as.numeric(auc(roc_obj))
      
      # obrni smer, da bo AUC vedno >= 0.5
      if (auc_val < 0.5) {
        auc_val <- 1 - auc_val
      }
      
      reaction_auc_results <- rbind(
        reaction_auc_results,
        data.frame(
          subsystem = sub_i,
          reaction_id = reaction_i,
          AUC = auc_val,
          n_samples = nrow(tmp)
        )
      )
    }
  }
}

# 5) uredi po AUC
reaction_auc_results <- reaction_auc_results %>%
  arrange(desc(AUC))

# 6) shrani
write.csv2(
  reaction_auc_results,
  "all_podsistemi_reaction_auc_flux.csv",
  row.names = FALSE
)
#_____________________________________________________________________
#4. Analiza rezultatov
#______________________________________________________________________
#izpis aktivnih reakcij v top n podsistemih
reaction_counts_data <- human1 %>%
  select(ID, SUBSYSTEM) %>%
  filter(!is.na(ID), !is.na(SUBSYSTEM)) %>%
  filter(SUBSYSTEM %in% top10_subsystems_flux) %>%
  filter(ID %in% colnames(flux_clean)) %>%
  distinct() %>%
  group_by(SUBSYSTEM) %>%
  summarise(n_reactions = n(), .groups = "drop") %>%
  arrange(desc(n_reactions))

reaction_counts_data

library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
#izbira reakcije in priprava podatkov
reaction_id <- "MAR09131"

# poišči ime reakcije in podsistem
reaction_info <- human1 %>%
  select(ID, NAME, SUBSYSTEM, EQUATION) %>%
  filter(ID == reaction_id)

print(reaction_info)

# pripravi podatke za grafiranje
reaction_data <- flux_clean %>%
  select(patient_id, type, all_of(reaction_id)) %>%
  rename(flux = all_of(reaction_id)) %>%
  mutate(type = factor(type, levels = c("normal", "tumor")))

summary(reaction_data$flux)
table(reaction_data$type)
#jitter boxplot

ggplot(reaction_data, aes(x = type, y = flux, color = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  labs(
    title = paste("Flux reaction", reaction_id),
    subtitle = paste(
      unique(reaction_info$NAME),
      "|",
      unique(reaction_info$SUBSYSTEM)
    ),
    x = "Tip vzorca",
    y = "Flux"
  ) +
  theme_minimal()
#density plot

ggplot(reaction_data, aes(x = flux, fill = type)) +
  geom_density(alpha = 0.4) +
  labs(
    title = paste("Gostota porazdelitve fluxa za", reaction_id),
    subtitle = paste(
      unique(reaction_info$NAME),
      "|",
      unique(reaction_info$SUBSYSTEM)
    ),
    x = "Flux",
    y = "Gostota"
  ) +
  theme_minimal()
#histogram
ggplot(reaction_data, aes(x = flux, fill = type)) +
  geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  labs(
    title = paste("Histogram fluxa za", reaction_id),
    x = "Flux",
    y = "Število vzorcev"
  ) +
  theme_minimal()
#ROC krivulja
roc_data <- reaction_data %>%
  mutate(label = ifelse(type == "tumor", 1, 0))

roc_obj <- roc(
  response = roc_data$label,
  predictor = roc_data$flux,
  quiet = TRUE
)

auc_val <- as.numeric(auc(roc_obj))
auc_val_display <- ifelse(auc_val < 0.5, 1 - auc_val, auc_val)

plot(
  roc_obj,
  main = paste("ROC krivulja za", reaction_id,
               "- AUC =", round(auc_val_display, 4))
)
abline(a = 0, b = 1, lty = 2, col = "grey")

#paired plot
ggplot(reaction_data, aes(x = type, y = flux, group = patient_id)) +
  geom_line(alpha = 0.35) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = paste("Paired flux for", reaction_id),
    subtitle = paste(
      unique(reaction_info$NAME),
      "|",
      unique(reaction_info$SUBSYSTEM)
    ),
    x = "Tip vzorca",
    y = "Flux"
  ) +
  theme_minimal()

#izračun wilcoxon za fluxe
#__________________________________________
# WILCOXON TEST ZA FLUX PODATKE
# VERZIJA 2:
# globalna BH-korekcija na unikatnih reakcijah
#__________________________________________

library(dplyr)


# ------------------------------------------------------------
# 1) VHODNI PODATKI
# ------------------------------------------------------------

# tabela z AUC-vrednostmi reakcij in pripadajočimi podsistemi
reaction_flux <- read.csv2(
  "all_podsistemi_reaction_auc_flux.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# očiščeni podatki metabolnih pretokov
flux_clean <- read.csv(
  "flux_clean.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 2) SEZNAM UNIKATNIH REAKCIJ
# ------------------------------------------------------------

# ista reakcija se lahko v reaction_flux pojavi večkrat,
# ker je lahko povezana z več podsistemi
unique_reactions <- reaction_flux %>%
  distinct(reaction_id)

cat(
  "Število unikatnih reakcij v tabeli:",
  nrow(unique_reactions),
  "\n"
)

cat(
  "Število unikatnih reakcij, prisotnih v flux_clean:",
  sum(unique_reactions$reaction_id %in% colnames(flux_clean)),
  "\n"
)


# ------------------------------------------------------------
# 3) PRIPRAVA PARNIH NORMALNIH IN TUMORSKIH VZORCEV
# ------------------------------------------------------------

normal_flux <- flux_clean %>%
  filter(type == "normal") %>%
  arrange(patient_id)

tumor_flux <- flux_clean %>%
  filter(type == "tumor") %>%
  arrange(patient_id)

# preveri, da je število normalnih in tumorskih vzorcev enako
stopifnot(
  nrow(normal_flux) == nrow(tumor_flux)
)

# preveri, da so pacienti v popolnoma enakem vrstnem redu
stopifnot(
  identical(
    as.character(normal_flux$patient_id),
    as.character(tumor_flux$patient_id)
  )
)

cat(
  "Število parnih pacientov:",
  nrow(normal_flux),
  "\n"
)


# ------------------------------------------------------------
# 4) PRAZNA TABELA REZULTATOV
# ena vrstica = ena unikatna reakcija
# ------------------------------------------------------------

reaction_flux_unique <- unique_reactions %>%
  mutate(
    n_pairs = NA_integer_,
    n_samples = NA_integer_,
    wilcox_p = NA_real_
  )


# ------------------------------------------------------------
# 5) PARNI WILCOXONOV TEST
# ------------------------------------------------------------

for (i in seq_len(nrow(reaction_flux_unique))) {
  
  reaction <- reaction_flux_unique$reaction_id[i]
  
  # test lahko izvedemo samo, če je reakcija stolpec v flux_clean
  if (reaction %in% colnames(flux_clean)) {
    
    normal_vals <- normal_flux[[reaction]]
    tumor_vals  <- tumor_flux[[reaction]]
    
    # obdržimo samo paciente, pri katerih sta prisotni obe vrednosti
    valid_pairs <- !is.na(normal_vals) & !is.na(tumor_vals)
    
    n_pairs <- sum(valid_pairs)
    
    reaction_flux_unique$n_pairs[i] <- n_pairs
    reaction_flux_unique$n_samples[i] <- 2 * n_pairs
    
    # test izvedemo, če imamo vsaj 5 veljavnih parov
    if (n_pairs >= 5) {
      
      normal_valid <- normal_vals[valid_pairs]
      tumor_valid  <- tumor_vals[valid_pairs]
      
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
        reaction_flux_unique$wilcox_p[i] <-
          wilcox_result$p.value
      }
    }
  }
}


# ------------------------------------------------------------
# 6) GLOBALNA BH-KOREKCIJA
# korekcija čez vse unikatne analizirane reakcije
# ------------------------------------------------------------

reaction_flux_unique$wilcox_padj <- p.adjust(
  reaction_flux_unique$wilcox_p,
  method = "BH"
)


# ------------------------------------------------------------
# 7) SHRANJEVANJE REZULTATOV NA RAVNI UNIKATNIH REAKCIJ
# ------------------------------------------------------------

write.csv2(
  reaction_flux_unique,
  "globalne_unikatne_reakcije_wilcox_flux.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 8) PRIKLJUČITEV REZULTATOV NA TABELO REAKCIJA–PODSISTEM
# ------------------------------------------------------------

# odstranimo morebitne stare statistične stolpce,
# da po joinu ne nastanejo wilcox_p.x in wilcox_p.y
reaction_flux <- reaction_flux %>%
  select(
    -any_of(
      c(
        "n_pairs",
        "n_samples",
        "wilcox_p",
        "wilcox_padj"
      )
    )
  ) %>%
  left_join(
    reaction_flux_unique %>%
      select(
        reaction_id,
        n_pairs,
        n_samples,
        wilcox_p,
        wilcox_padj
      ),
    by = "reaction_id"
  )


# ------------------------------------------------------------
# 9) IZRAČUN RANGA NORMALIZIRANIH DISTANC PODSISTEMOV
# ------------------------------------------------------------

summary_results_flux <- read.csv2(
  "mediana_distanc_podsistemov_flux.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# večja mediana normaliziranih distanc =
# višje uvrščen podsistem
summary_results_flux <- summary_results_flux %>%
  arrange(desc(median_distance)) %>%
  mutate(
    Rank_flux_norm_distance = row_number()
  )

write.csv2(
  summary_results_flux,
  "mediana_distanc_podsistemov_flux.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 10) PRIKLJUČITEV RANGA PODSISTEMA
# ------------------------------------------------------------

rank_flux <- summary_results_flux %>%
  select(
    subsystem,
    Rank_flux_norm_distance
  )

# prepreči nenamerno podvajanje vrstic pri left_join()
stopifnot(
  !anyDuplicated(rank_flux$subsystem)
)

# odstranimo morebitni stari stolpec z rangom
reaction_flux <- reaction_flux %>%
  select(
    -any_of("Rank_flux_norm_distance")
  ) %>%
  left_join(
    rank_flux,
    by = "subsystem"
  )


# ------------------------------------------------------------
# 11) SHRANJEVANJE KONČNE TABELE
# ------------------------------------------------------------

write.csv2(
  reaction_flux,
  "all_podsistemi_reaction_auc_flux.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 12) OSNOVNI PREGLED REZULTATOV
# ------------------------------------------------------------

cat(
  "Število unikatnih reakcij z izračunano p-vrednostjo:",
  sum(!is.na(reaction_flux_unique$wilcox_p)),
  "\n"
)

cat(
  "Število unikatnih reakcij s p < 0,01:",
  sum(
    reaction_flux_unique$wilcox_p < 0.01,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Število unikatnih reakcij s padj < 0,01:",
  sum(
    reaction_flux_unique$wilcox_padj < 0.01,
    na.rm = TRUE
  ),
  "\n"
)