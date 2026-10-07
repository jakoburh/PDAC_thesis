# ==============================================================================
# LASSO PIPELINE V2: GENI + METABOLNI PRETOKI + KONČNI LOGISTIČNI MODELI
# ==============================================================================
#
# Namen:
#   1) Iz integrirane preglednice izbere kandidatne gene in reakcije
#      z wilcox_padj < 0,01.
#   2) Za gene in reakcije ločeno izvede logistični LASSO:
#        - en podroben izračun pri seedu 123,
#        - stabilnostno analizo pri 100 različnih seedih.
#   3) Izračuna frekvenco izbire posameznih genov in reakcij pri lambda.1se.
#   4) Izbere vse gene in reakcije, ki so bili pri lambda.1se izbrani
#      v več kot 70 od 100 ponovitev.
#   5) Pri 100 delitvah pacientov preveri:
#        - model samo s stabilnimi geni,
#        - model samo s stabilnimi reakcijami,
#        - kombinirani model s stabilnimi geni in reakcijami.
#
# Pomembne lastnosti različice 2:
#   - vsi modeli uporabljajo isti nabor 52 skupnih parnih pacientov;
#   - normalni in tumorski vzorec istega pacienta vedno ostaneta skupaj;
#   - tudi notranji pregibi cv.glmnet so določeni na ravni pacientov;
#   - smer ROC je določena eksplicitno;
#   - preverjajo se podvojeni patient_id/type pari;
#   - vsi rezultati se shranijo v ločeno mapo LASSO_pipeline_v2.
#
# Zahtevane vhodne datoteke v delovni mapi:
#   - DESeq2_z_log2TPM_z_flux.csv
#   - CPTAC_clean_log2.csv
#   - flux_clean.csv
#
# ==============================================================================


# ==============================================================================
# 0) NASTAVITVE
# ==============================================================================

# Imena vhodnih datotek.
FILE_INTEGRATED <- "DESeq2_z_log2TPM_z_flux.csv"
FILE_GENES      <- "CPTAC_clean_log2.csv"
FILE_FLUX       <- "flux_clean.csv"

# Glavna izhodna mapa različice 2.
OUTPUT_DIR <- "LASSO_pipeline_v2_frequency_gt70"

# Prag za vključitev kandidatnih značilk.
PADJ_CUTOFF <- 0.01

# Število pacientov in velikost zunanje train/test delitve.
N_PATIENTS_TOTAL <- 52
N_TRAIN_PATIENTS <- 35

# Podroben enkratni izračun.
SINGLE_SEED <- 123

# Stabilnostna analiza in končna validacija.
STABILITY_SEEDS <- 1:100
VALIDATION_SEEDS <- 1:100

# Število notranjih pregibov pri cv.glmnet.
N_INNER_FOLDS <- 10

# Prag stabilnosti izbora pri lambda.1se.
# Ker uporabljamo operator >, mora biti značilka izbrana najmanj 71-krat.
SELECTION_FREQUENCY_CUTOFF <- 70

# Prag za pretvorbo napovedane verjetnosti v razred.
CLASS_THRESHOLD <- 0.5

# TRUE pomeni, da se značilke pred končnim navadnim glm modelom
# standardizirajo z uporabo povprečja in SD učnega nabora.
# FALSE ohrani logiko prvotnega pipeline-a.
STANDARDIZE_FINAL_GLM <- TRUE


# ==============================================================================
# 1) KNJIŽNICE
# ==============================================================================

# Funkcija preveri, ali so vsi zahtevani paketi nameščeni.
# Če kateri manjka, ustavi izvajanje in izpiše ukaz za namestitev.
check_packages <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  
  if (length(missing_packages) > 0) {
    stop(
      paste0(
        "Manjkajo paketi: ",
        paste(missing_packages, collapse = ", "),
        ". Namesti jih z: install.packages(c(",
        paste(sprintf('"%s"', missing_packages), collapse = ", "),
        "))"
      )
    )
  }
}

required_packages <- c("dplyr", "readr", "pROC", "glmnet", "tibble")
check_packages(required_packages)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pROC)
  library(glmnet)
  library(tibble)
})


# ==============================================================================
# 2) USTVARJANJE MAP
# ==============================================================================

# Funkcija ustvari mapo, če ta še ne obstaja.
# recursive = TRUE omogoči ustvarjanje več nivojev map naenkrat.
make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

make_dir(OUTPUT_DIR)
make_dir(file.path(OUTPUT_DIR, "00_input_checks"))
make_dir(file.path(OUTPUT_DIR, "01_single_seed"))
make_dir(file.path(OUTPUT_DIR, "02_gene_lasso", "seed_results"))
make_dir(file.path(OUTPUT_DIR, "03_reaction_lasso", "seed_results"))
make_dir(file.path(OUTPUT_DIR, "04_final_models", "stable_genes_gt70", "predictions"))
make_dir(file.path(OUTPUT_DIR, "04_final_models", "stable_reactions_gt70", "predictions"))
make_dir(file.path(OUTPUT_DIR, "04_final_models", "combined_stable_gt70", "predictions"))
make_dir(file.path(OUTPUT_DIR, "05_final_summary"))


# ==============================================================================
# 3) SPLOŠNE POMOŽNE FUNKCIJE
# ==============================================================================

# Funkcija pretvori števila, zapisana z decimalno vejico ali piko,
# v numerični vektor. Besedilo, ki ni število, postane NA.
as_numeric_comma <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}


# Funkcija preveri, ali podatkovni okvir vsebuje zahtevane stolpce.
# Če kateri stolpec manjka, izvajanje ustavi z jasnim sporočilom.
assert_columns <- function(data, required_columns, data_name) {
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        data_name,
        " nima zahtevanih stolpcev: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}


# Funkcija standardizira imena metapodatkovnih stolpcev na:
#   - patient_id
#   - type
# in obdrži samo normalne in tumorske vzorce.
normalise_sample_table <- function(data, patient_column, type_column, data_name) {
  assert_columns(data, c(patient_column, type_column), data_name)
  
  names(data)[names(data) == patient_column] <- "patient_id"
  names(data)[names(data) == type_column] <- "type"
  
  data <- data %>%
    mutate(
      patient_id = as.character(patient_id),
      type = tolower(trimws(as.character(type)))
    ) %>%
    filter(type %in% c("normal", "tumor"))
  
  data
}


# Funkcija vse stolpce razen metapodatkov pretvori v numerične.
# To je posebej pomembno pri flux_clean.csv, kjer so lahko decimalke zapisane z vejico.
convert_feature_columns_to_numeric <- function(
    data,
    metadata_columns = c("patient_id", "type")
) {
  feature_columns <- setdiff(colnames(data), metadata_columns)
  
  data[feature_columns] <- lapply(
    data[feature_columns],
    as_numeric_comma
  )
  
  data
}


# Funkcija preveri, ali ima vsak patient_id/type natanko eno vrstico.
# Če obstajajo replikati ali podvojene vrstice, jih izpiše in ustavi izvajanje.
check_unique_patient_type <- function(data, data_name) {
  duplicates <- data %>%
    count(patient_id, type, name = "n") %>%
    filter(n != 1)
  
  if (nrow(duplicates) > 0) {
    print(duplicates)
    
    stop(
      paste0(
        data_name,
        " nima natanko ene vrstice za vsak patient_id/type. ",
        "Pred nadaljevanjem je treba razrešiti replikate ali duplikate."
      )
    )
  }
  
  invisible(TRUE)
}


# Funkcija vrne paciente, ki imajo hkrati:
#   - en normalen vzorec,
#   - en tumorski vzorec.
get_paired_patients <- function(data) {
  data %>%
    distinct(patient_id, type) %>%
    count(patient_id, name = "n_types") %>%
    filter(n_types == 2) %>%
    pull(patient_id)
}


# Funkcija izbere skupne parne paciente, ki obstajajo tako v genski
# kot v reakcijski podatkovni tabeli.
# Če jih je več kot 52, reproducibilno izbere 52 pacientov.
select_common_patients <- function(
    gene_data,
    flux_data,
    n_patients = 52,
    seed = 123
) {
  paired_genes <- get_paired_patients(gene_data)
  paired_flux  <- get_paired_patients(flux_data)
  
  common_patients <- intersect(paired_genes, paired_flux)
  
  if (length(common_patients) < n_patients) {
    stop(
      paste0(
        "Na voljo je samo ",
        length(common_patients),
        " skupnih parnih pacientov, potrebnih pa je ",
        n_patients,
        "."
      )
    )
  }
  
  common_patients <- sort(unique(common_patients))
  
  if (length(common_patients) == n_patients) {
    return(common_patients)
  }
  
  set.seed(seed)
  sort(sample(common_patients, size = n_patients, replace = FALSE))
}


# Funkcija razdeli paciente na učni in testni nabor.
# Delitev poteka na ravni pacientov, zato oba vzorca istega pacienta
# vedno ostaneta v isti skupini.
make_patient_split <- function(
    patient_ids,
    seed,
    n_train_patients = 35
) {
  patient_ids <- sort(unique(as.character(patient_ids)))
  
  if (length(patient_ids) <= n_train_patients) {
    stop("Število pacientov mora biti večje od velikosti učnega nabora.")
  }
  
  set.seed(seed)
  
  train_patients <- sample(
    patient_ids,
    size = n_train_patients,
    replace = FALSE
  )
  
  test_patients <- setdiff(patient_ids, train_patients)
  
  list(
    train = sort(train_patients),
    test = sort(test_patients)
  )
}


# Funkcija izdela notranje CV-pregibe na ravni pacientov.
# Normalni in tumorski vzorec istega pacienta dobita isti foldid,
# s čimer preprečimo uhajanje informacij med notranjimi pregibi.
make_patient_foldid <- function(
    train_data,
    nfolds = 10,
    seed = 123
) {
  train_patients <- sort(unique(train_data$patient_id))
  nfolds_actual <- min(nfolds, length(train_patients))
  
  if (nfolds_actual < 3) {
    stop("Za notranjo navzkrižno validacijo so potrebni vsaj 3 pregibi.")
  }
  
  set.seed(seed)
  
  patient_folds <- sample(
    rep(seq_len(nfolds_actual), length.out = length(train_patients))
  )
  
  fold_lookup <- setNames(patient_folds, train_patients)
  
  unname(fold_lookup[train_data$patient_id])
}


# Funkcija varno izračuna AUC.
# Če ROC zaradi neustreznih podatkov ni mogoče izračunati, vrne NA.
safe_auc <- function(y_true, probabilities) {
  if (
    length(unique(y_true)) < 2 ||
    all(is.na(probabilities)) ||
    length(unique(probabilities[!is.na(probabilities)])) < 2
  ) {
    return(NA_real_)
  }
  
  tryCatch(
    {
      roc_object <- pROC::roc(
        response = y_true,
        predictor = probabilities,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      )
      
      as.numeric(pROC::auc(roc_object))
    },
    error = function(e) NA_real_
  )
}


# Funkcija iz napovedanih verjetnosti izračuna:
#   - točnost,
#   - občutljivost,
#   - specifičnost,
#   - TP, TN, FP in FN.
classification_metrics <- function(
    y_true,
    probabilities,
    threshold = 0.5
) {
  predicted_class <- ifelse(probabilities > threshold, 1, 0)
  
  confusion_matrix <- table(
    Predicted = factor(predicted_class, levels = c(0, 1)),
    Actual = factor(y_true, levels = c(0, 1))
  )
  
  TN <- as.numeric(confusion_matrix["0", "0"])
  FN <- as.numeric(confusion_matrix["0", "1"])
  FP <- as.numeric(confusion_matrix["1", "0"])
  TP <- as.numeric(confusion_matrix["1", "1"])
  
  tibble(
    accuracy = mean(predicted_class == y_true),
    sensitivity = ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN)),
    specificity = ifelse((TN + FP) == 0, NA_real_, TN / (TN + FP)),
    TN = TN,
    FN = FN,
    FP = FP,
    TP = TP
  )
}


# Funkcija varno izračuna minimum in maksimum ob prisotnih NA.
safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


# ==============================================================================
# 4) FUNKCIJE ZA PRIPRAVO X/Y MATRIK
# ==============================================================================

# Funkcija iz ene modalnosti (geni ali reakcije) pripravi:
#   - X_train,
#   - X_test,
#   - y_train,
#   - y_test,
#   - metapodatke testnega nabora.
#
# Odstrani značilke, ki:
#   - ne obstajajo v tabeli,
#   - vsebujejo NA v učnem ali testnem naboru,
#   - imajo ničelno varianco v učnem naboru.
prepare_xy <- function(
    data,
    candidate_features,
    split
) {
  train_data <- data %>%
    filter(patient_id %in% split$train) %>%
    arrange(patient_id, type)
  
  test_data <- data %>%
    filter(patient_id %in% split$test) %>%
    arrange(patient_id, type)
  
  expected_train_rows <- 2 * length(split$train)
  expected_test_rows  <- 2 * length(split$test)
  
  if (nrow(train_data) != expected_train_rows) {
    stop("Učni nabor nima pričakovanih dveh vzorcev na pacienta.")
  }
  
  if (nrow(test_data) != expected_test_rows) {
    stop("Testni nabor nima pričakovanih dveh vzorcev na pacienta.")
  }
  
  available_features <- intersect(
    candidate_features,
    intersect(colnames(train_data), colnames(test_data))
  )
  
  if (length(available_features) == 0) {
    stop("Nobena kandidatna značilka ni prisotna v učnem in testnem naboru.")
  }
  
  X_train <- as.matrix(
    train_data[, available_features, drop = FALSE]
  )
  
  X_test <- as.matrix(
    test_data[, available_features, drop = FALSE]
  )
  
  storage.mode(X_train) <- "double"
  storage.mode(X_test) <- "double"
  
  no_missing <- (
    colSums(is.na(X_train)) == 0 &
      colSums(is.na(X_test)) == 0
  )
  
  train_sd <- apply(X_train, 2, sd)
  non_constant <- !is.na(train_sd) & train_sd > 0
  
  good_features <- available_features[no_missing & non_constant]
  
  if (length(good_features) == 0) {
    stop(
      "Po odstranitvi značilk z NA ali ničelno varianco ni ostala nobena značilka."
    )
  }
  
  X_train <- X_train[, good_features, drop = FALSE]
  X_test  <- X_test[, good_features, drop = FALSE]
  
  y_train <- ifelse(train_data$type == "tumor", 1, 0)
  y_test  <- ifelse(test_data$type == "tumor", 1, 0)
  
  list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test,
    train_metadata = train_data %>% select(patient_id, type),
    test_metadata = test_data %>% select(patient_id, type),
    features = good_features,
    train_data = train_data,
    test_data = test_data
  )
}


# Funkcija iz objekta glmnet izlušči neničelne koeficiente
# pri izbrani vrednosti lambda.
extract_nonzero_coefficients <- function(
    cv_fit,
    lambda_name,
    feature_id_column = "feature_id"
) {
  coefficient_matrix <- as.matrix(
    coef(cv_fit, s = lambda_name)
  )
  
  result <- tibble(
    feature_id = rownames(coefficient_matrix),
    coefficient = as.numeric(coefficient_matrix[, 1])
  ) %>%
    filter(
      feature_id != "(Intercept)",
      coefficient != 0
    ) %>%
    arrange(desc(abs(coefficient)))
  
  names(result)[names(result) == "feature_id"] <- feature_id_column
  result
}


# ==============================================================================
# 5) LOGISTIČNI LASSO ZA ENO TRAIN/TEST DELITEV
# ==============================================================================

# Funkcija izvede celoten logistični LASSO pri eni zunanji delitvi:
#   1) pripravi učne in testne matrike;
#   2) izdela patient-level notranje pregibe;
#   3) s cv.glmnet izbere lambda;
#   4) izračuna napovedi pri lambda.1se in lambda.min;
#   5) vrne izbrane značilke, AUC in druge rezultate.
fit_lasso_one_split <- function(
    data,
    candidate_features,
    split,
    seed,
    n_inner_folds = 10
) {
  prepared <- prepare_xy(
    data = data,
    candidate_features = candidate_features,
    split = split
  )
  
  foldid <- make_patient_foldid(
    train_data = prepared$train_data,
    nfolds = n_inner_folds,
    seed = seed
  )
  
  set.seed(seed)
  
  cv_fit <- cv.glmnet(
    x = prepared$X_train,
    y = prepared$y_train,
    family = "binomial",
    alpha = 1,
    foldid = foldid,
    standardize = TRUE,
    type.measure = "auc"
  )
  
  probability_1se <- as.vector(
    predict(
      cv_fit,
      newx = prepared$X_test,
      s = "lambda.1se",
      type = "response"
    )
  )
  
  probability_min <- as.vector(
    predict(
      cv_fit,
      newx = prepared$X_test,
      s = "lambda.min",
      type = "response"
    )
  )
  
  coefficients_1se <- extract_nonzero_coefficients(
    cv_fit,
    lambda_name = "lambda.1se"
  )
  
  coefficients_min <- extract_nonzero_coefficients(
    cv_fit,
    lambda_name = "lambda.min"
  )
  
  list(
    cv_fit = cv_fit,
    prepared = prepared,
    probability_1se = probability_1se,
    probability_min = probability_min,
    auc_1se = safe_auc(prepared$y_test, probability_1se),
    auc_min = safe_auc(prepared$y_test, probability_min),
    coefficients_1se = coefficients_1se,
    coefficients_min = coefficients_min
  )
}


# ==============================================================================
# 6) SHRANJEVANJE PODROBNEGA IZRAČUNA PRI ENEM SEEDU
# ==============================================================================

# Funkcija shrani podrobne rezultate ene LASSO delitve:
#   - CV-graf,
#   - ROC-grafa,
#   - koeficiente pri lambda.1se in lambda.min,
#   - napovedi na testnem naboru,
#   - povzetek modela.
save_single_seed_lasso <- function(
    fit,
    seed,
    modality_name,
    output_directory
) {
  make_dir(output_directory)
  
  coefficients_1se <- fit$coefficients_1se
  coefficients_min <- fit$coefficients_min
  
  write.csv2(
    coefficients_1se,
    file.path(output_directory, paste0(modality_name, "_lambda1se_coefficients.csv")),
    row.names = FALSE
  )
  
  write.csv2(
    coefficients_min,
    file.path(output_directory, paste0(modality_name, "_lambdamin_coefficients.csv")),
    row.names = FALSE
  )
  
  predictions <- fit$prepared$test_metadata %>%
    mutate(
      y_true = fit$prepared$y_test,
      probability_lambda1se = fit$probability_1se,
      probability_lambdamin = fit$probability_min,
      predicted_lambda1se = ifelse(
        probability_lambda1se > CLASS_THRESHOLD, 1, 0
      ),
      predicted_lambdamin = ifelse(
        probability_lambdamin > CLASS_THRESHOLD, 1, 0
      )
    )
  
  write.csv2(
    predictions,
    file.path(output_directory, paste0(modality_name, "_test_predictions.csv")),
    row.names = FALSE
  )
  
  summary_table <- tibble(
    seed = seed,
    n_train_patients = length(unique(fit$prepared$train_metadata$patient_id)),
    n_test_patients = length(unique(fit$prepared$test_metadata$patient_id)),
    n_features_input = ncol(fit$prepared$X_train),
    lambda_min = fit$cv_fit$lambda.min,
    lambda_1se = fit$cv_fit$lambda.1se,
    max_internal_cv_auc = max(fit$cv_fit$cvm, na.rm = TRUE),
    test_auc_1se = fit$auc_1se,
    test_auc_min = fit$auc_min,
    n_selected_1se = nrow(coefficients_1se),
    n_selected_min = nrow(coefficients_min)
  )
  
  write.csv2(
    summary_table,
    file.path(output_directory, paste0(modality_name, "_single_seed_summary.csv")),
    row.names = FALSE
  )
  
  png(
    file.path(output_directory, paste0(modality_name, "_cv_curve.png")),
    width = 1200,
    height = 900,
    res = 140
  )
  plot(fit$cv_fit)
  title(main = paste0(modality_name, ": notranja CV, seed ", seed))
  dev.off()
  
  roc_1se <- tryCatch(
    pROC::roc(
      fit$prepared$y_test,
      fit$probability_1se,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
  
  roc_min <- tryCatch(
    pROC::roc(
      fit$prepared$y_test,
      fit$probability_min,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(roc_1se) && !is.null(roc_min)) {
    png(
      file.path(output_directory, paste0(modality_name, "_roc_comparison.png")),
      width = 1200,
      height = 900,
      res = 140
    )
    plot(
      roc_1se,
      col = "blue",
      main = paste0(modality_name, ": ROC primerjava, seed ", seed)
    )
    plot(roc_min, col = "red", add = TRUE)
    legend(
      "bottomright",
      legend = c(
        paste0("lambda.1se, AUC = ", round(fit$auc_1se, 3)),
        paste0("lambda.min, AUC = ", round(fit$auc_min, 3))
      ),
      col = c("blue", "red"),
      lwd = 2
    )
    dev.off()
  }
  
  invisible(summary_table)
}


# ==============================================================================
# 7) STABILNOSTNA LASSO ANALIZA ČEZ VEČ SEEDOV
# ==============================================================================

# Funkcija pri vseh podanih seedih:
#   - naredi novo patient-level train/test delitev;
#   - izvede logistični LASSO;
#   - shrani izbrane značilke za posamezen seed;
#   - zbere AUC in število izbranih značilk;
#   - izračuna frekvenco izbire vsake značilke pri lambda.1se in lambda.min.
run_lasso_stability <- function(
    data,
    candidate_features,
    patient_ids,
    seeds,
    modality_name,
    feature_id_column,
    output_directory
) {
  make_dir(output_directory)
  seed_directory <- file.path(output_directory, "seed_results")
  make_dir(seed_directory)
  
  candidate_features <- sort(unique(candidate_features))
  
  seed_summaries <- list()
  all_selected_coefficients <- list()
  
  for (seed in seeds) {
    cat("LASSO stabilnost -", modality_name, "- seed:", seed, "\n")
    
    split <- make_patient_split(
      patient_ids = patient_ids,
      seed = seed,
      n_train_patients = N_TRAIN_PATIENTS
    )
    
    fit <- tryCatch(
      fit_lasso_one_split(
        data = data,
        candidate_features = candidate_features,
        split = split,
        seed = seed,
        n_inner_folds = N_INNER_FOLDS
      ),
      error = function(e) {
        warning(
          paste0(
            modality_name,
            ", seed ",
            seed,
            ": ",
            conditionMessage(e)
          )
        )
        NULL
      }
    )
    
    if (is.null(fit)) next
    
    coefficients_1se <- fit$coefficients_1se %>%
      rename(
        coefficient_1se = coefficient
      )
    
    coefficients_min <- fit$coefficients_min %>%
      rename(
        coefficient_min = coefficient
      )
    
    seed_coefficients <- full_join(
      coefficients_1se,
      coefficients_min,
      by = "feature_id"
    ) %>%
      mutate(seed = seed) %>%
      select(seed, feature_id, coefficient_1se, coefficient_min) %>%
      arrange(feature_id)
    
    output_seed_table <- seed_coefficients
    names(output_seed_table)[
      names(output_seed_table) == "feature_id"
    ] <- feature_id_column
    
    write.csv2(
      output_seed_table,
      file.path(seed_directory, paste0(seed, ".csv")),
      row.names = FALSE
    )
    
    metrics_1se <- classification_metrics(
      fit$prepared$y_test,
      fit$probability_1se,
      CLASS_THRESHOLD
    )
    
    metrics_min <- classification_metrics(
      fit$prepared$y_test,
      fit$probability_min,
      CLASS_THRESHOLD
    )
    
    seed_summaries[[as.character(seed)]] <- tibble(
      seed = seed,
      n_features_available = ncol(fit$prepared$X_train),
      auc_1se = fit$auc_1se,
      auc_min = fit$auc_min,
      accuracy_1se = metrics_1se$accuracy,
      accuracy_min = metrics_min$accuracy,
      n_selected_1se = nrow(fit$coefficients_1se),
      n_selected_min = nrow(fit$coefficients_min),
      lambda_1se = fit$cv_fit$lambda.1se,
      lambda_min = fit$cv_fit$lambda.min
    )
    
    all_selected_coefficients[[as.character(seed)]] <- seed_coefficients
  }
  
  stability_summary <- bind_rows(seed_summaries)
  coefficient_long <- bind_rows(all_selected_coefficients)
  
  write.csv2(
    stability_summary,
    file.path(output_directory, paste0(modality_name, "_stability_summary.csv")),
    row.names = FALSE
  )
  
  write.csv2(
    coefficient_long,
    file.path(output_directory, paste0(modality_name, "_all_selected_coefficients.csv")),
    row.names = FALSE
  )
  
  frequency_summary <- tibble(feature_id = candidate_features) %>%
    left_join(
      coefficient_long %>%
        group_by(feature_id) %>%
        summarise(
          frequency = sum(!is.na(coefficient_1se)),
          frequency_min = sum(!is.na(coefficient_min)),
          minimal_coefficient_1se = safe_min(coefficient_1se),
          maximal_coefficient_1se = safe_max(coefficient_1se),
          median_coefficient_1se = ifelse(
            all(is.na(coefficient_1se)),
            NA_real_,
            median(coefficient_1se, na.rm = TRUE)
          ),
          .groups = "drop"
        ),
      by = "feature_id"
    ) %>%
    mutate(
      frequency = ifelse(is.na(frequency), 0, frequency),
      frequency_min = ifelse(is.na(frequency_min), 0, frequency_min)
    ) %>%
    arrange(desc(frequency), feature_id)
  
  names(frequency_summary)[
    names(frequency_summary) == "feature_id"
  ] <- feature_id_column
  
  write.csv2(
    frequency_summary,
    file.path(output_directory, paste0(modality_name, "_frequency_summary.csv")),
    row.names = FALSE
  )
  
  list(
    stability_summary = stability_summary,
    coefficient_long = coefficient_long,
    frequency_summary = frequency_summary
  )
}


# ==============================================================================
# 8) IZBOR STABILNIH ZNAČILK
# ==============================================================================

# Funkcija iz frekvenčne preglednice obdrži vse značilke,
# ki so bile pri lambda.1se izbrane pogosteje od podanega praga.
#
# Pri frequency_cutoff = 70 in pogoju frequency > 70
# se obdržijo značilke, izbrane najmanj 71-krat od 100 ponovitev.
#
# Funkcija vrne celotno filtrirano preglednico, ne le identifikatorjev,
# zato ostanejo na voljo tudi frekvence ter najmanjši, največji
# in mediani koeficienti.
select_stable_features <- function(
    frequency_summary,
    feature_id_column,
    frequency_cutoff
) {
  stable_table <- frequency_summary %>%
    filter(frequency > frequency_cutoff) %>%
    arrange(
      desc(frequency),
      .data[[feature_id_column]]
    )
  
  if (nrow(stable_table) == 0) {
    warning(
      paste0(
        "Nobena značilka ni imela frekvence > ",
        frequency_cutoff,
        "."
      )
    )
  }
  
  stable_table
}


# ==============================================================================
# 9) PRIPRAVA PODATKOV ZA KONČNE LOGISTIČNE MODELE
# ==============================================================================

# Funkcija združi top genske in reakcijske značilke v eno tabelo.
# Pri modelu samo z geni ali samo z reakcijami sprejme tudi prazen vektor
# druge vrste značilk.
build_final_model_table <- function(
    gene_data,
    flux_data,
    patient_ids,
    gene_features = character(0),
    reaction_features = character(0)
) {
  patient_ids <- sort(unique(as.character(patient_ids)))
  
  use_genes <- length(gene_features) > 0
  use_reactions <- length(reaction_features) > 0
  
  if (!use_genes && !use_reactions) {
    stop("Za končni model ni podana nobena značilka.")
  }
  
  if (length(intersect(gene_features, reaction_features)) > 0) {
    stop("Genski in reakcijski identifikatorji se prekrivajo.")
  }
  
  gene_part <- NULL
  reaction_part <- NULL
  
  if (use_genes) {
    missing_genes <- setdiff(gene_features, colnames(gene_data))
    
    if (length(missing_genes) > 0) {
      stop(
        paste0(
          "V genski tabeli manjkajo značilke: ",
          paste(missing_genes, collapse = ", ")
        )
      )
    }
    
    gene_part <- gene_data %>%
      filter(patient_id %in% patient_ids) %>%
      select(patient_id, type, all_of(gene_features))
  }
  
  if (use_reactions) {
    missing_reactions <- setdiff(reaction_features, colnames(flux_data))
    
    if (length(missing_reactions) > 0) {
      stop(
        paste0(
          "V reakcijski tabeli manjkajo značilke: ",
          paste(missing_reactions, collapse = ", ")
        )
      )
    }
    
    reaction_part <- flux_data %>%
      filter(patient_id %in% patient_ids) %>%
      select(patient_id, type, all_of(reaction_features))
  }
  
  if (use_genes && use_reactions) {
    model_data <- inner_join(
      gene_part,
      reaction_part,
      by = c("patient_id", "type")
    )
  } else if (use_genes) {
    model_data <- gene_part
  } else {
    model_data <- reaction_part
  }
  
  model_data <- model_data %>%
    arrange(patient_id, type)
  
  expected_rows <- 2 * length(patient_ids)
  
  if (nrow(model_data) != expected_rows) {
    stop(
      paste0(
        "Končna modelna tabela ima ",
        nrow(model_data),
        " vrstic, pričakovanih pa je ",
        expected_rows,
        "."
      )
    )
  }
  
  model_data
}


# Funkcija standardizira učne in testne značilke.
# Povprečje in standardni odklon se izračunata samo na učnem naboru,
# nato pa se enaka transformacija uporabi tudi na testnem naboru.
standardize_train_test <- function(X_train, X_test) {
  train_means <- colMeans(X_train)
  train_sds <- apply(X_train, 2, sd)
  
  if (any(is.na(train_sds) | train_sds == 0)) {
    stop("Standardizacija ni mogoča zaradi ničelne variance.")
  }
  
  X_train_scaled <- sweep(X_train, 2, train_means, "-")
  X_train_scaled <- sweep(X_train_scaled, 2, train_sds, "/")
  
  X_test_scaled <- sweep(X_test, 2, train_means, "-")
  X_test_scaled <- sweep(X_test_scaled, 2, train_sds, "/")
  
  list(
    X_train = X_train_scaled,
    X_test = X_test_scaled
  )
}


# Funkcija pripravi podatke za navadno končno logistično regresijo.
# Odstrani značilke z NA ali ničelno varianco ter po želji
# standardizira značilke na osnovi učnega nabora.
prepare_final_glm_data <- function(
    model_data,
    feature_names,
    split,
    standardize = FALSE
) {
  train_data <- model_data %>%
    filter(patient_id %in% split$train) %>%
    arrange(patient_id, type)
  
  test_data <- model_data %>%
    filter(patient_id %in% split$test) %>%
    arrange(patient_id, type)
  
  X_train <- as.matrix(
    train_data[, feature_names, drop = FALSE]
  )
  
  X_test <- as.matrix(
    test_data[, feature_names, drop = FALSE]
  )
  
  storage.mode(X_train) <- "double"
  storage.mode(X_test) <- "double"
  
  no_missing <- (
    colSums(is.na(X_train)) == 0 &
      colSums(is.na(X_test)) == 0
  )
  
  train_sd <- apply(X_train, 2, sd)
  non_constant <- !is.na(train_sd) & train_sd > 0
  
  good_features <- feature_names[no_missing & non_constant]
  
  if (length(good_features) == 0) {
    stop("Za končni glm model ni ostala nobena uporabna značilka.")
  }
  
  X_train <- X_train[, good_features, drop = FALSE]
  X_test <- X_test[, good_features, drop = FALSE]
  
  if (standardize) {
    scaled <- standardize_train_test(X_train, X_test)
    X_train <- scaled$X_train
    X_test <- scaled$X_test
  }
  
  y_train <- ifelse(train_data$type == "tumor", 1, 0)
  y_test <- ifelse(test_data$type == "tumor", 1, 0)
  
  train_frame <- as.data.frame(X_train, check.names = FALSE)
  test_frame <- as.data.frame(X_test, check.names = FALSE)
  
  train_frame$y <- y_train
  
  list(
    train_frame = train_frame,
    test_frame = test_frame,
    y_train = y_train,
    y_test = y_test,
    train_metadata = train_data %>% select(patient_id, type),
    test_metadata = test_data %>% select(patient_id, type),
    features = good_features
  )
}


# ==============================================================================
# 10) NAVADNA LOGISTIČNA REGRESIJA ZA ENO DELITEV
# ==============================================================================

# Funkcija postavi večspremenljivčno logistično regresijo glm,
# izračuna testne napovedi in vrne opozorila ter koeficiente.
fit_final_glm_one_split <- function(
    model_data,
    feature_names,
    split,
    standardize = FALSE
) {
  prepared <- prepare_final_glm_data(
    model_data = model_data,
    feature_names = feature_names,
    split = split,
    standardize = standardize
  )
  
  model_warnings <- character(0)
  
  model <- withCallingHandlers(
    glm(
      y ~ .,
      data = prepared$train_frame,
      family = binomial()
    ),
    warning = function(w) {
      model_warnings <<- c(model_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  
  probabilities <- as.vector(
    predict(
      model,
      newdata = prepared$test_frame,
      type = "response"
    )
  )
  
  model_coefficients <- tibble(
    term = names(coef(model)),
    coefficient = as.numeric(coef(model))
  )
  
  list(
    model = model,
    prepared = prepared,
    probabilities = probabilities,
    auc = safe_auc(prepared$y_test, probabilities),
    metrics = classification_metrics(
      prepared$y_test,
      probabilities,
      CLASS_THRESHOLD
    ),
    coefficients = model_coefficients,
    warnings = unique(model_warnings)
  )
}


# ==============================================================================
# 11) KONČNA VALIDACIJA MODELA ČEZ 100 SEEDOV
# ==============================================================================

# Funkcija pri več seedih preveri napovedno uspešnost končnega
# navadnega logističnega modela.
#
# Za vsak seed:
#   - naredi isto patient-level delitev 35/17;
#   - model nauči samo na učnih pacientih;
#   - napovedi izračuna na testnih pacientih;
#   - shrani AUC, točnost, občutljivost, specifičnost in koeficiente.
run_final_glm_validation <- function(
    model_data,
    feature_names,
    patient_ids,
    seeds,
    model_name,
    output_directory,
    standardize = FALSE
) {
  make_dir(output_directory)
  prediction_directory <- file.path(output_directory, "predictions")
  make_dir(prediction_directory)
  
  summary_list <- list()
  coefficient_list <- list()
  
  for (seed in seeds) {
    cat("Končna validacija -", model_name, "- seed:", seed, "\n")
    
    split <- make_patient_split(
      patient_ids = patient_ids,
      seed = seed,
      n_train_patients = N_TRAIN_PATIENTS
    )
    
    fit <- tryCatch(
      fit_final_glm_one_split(
        model_data = model_data,
        feature_names = feature_names,
        split = split,
        standardize = standardize
      ),
      error = function(e) {
        warning(
          paste0(
            model_name,
            ", seed ",
            seed,
            ": ",
            conditionMessage(e)
          )
        )
        NULL
      }
    )
    
    if (is.null(fit)) next
    
    prediction_table <- fit$prepared$test_metadata %>%
      mutate(
        seed = seed,
        y_true = fit$prepared$y_test,
        prob_tumor = fit$probabilities,
        pred_class = ifelse(
          prob_tumor > CLASS_THRESHOLD,
          1,
          0
        )
      ) %>%
      select(
        seed,
        patient_id,
        type,
        y_true,
        prob_tumor,
        pred_class
      )
    
    write.csv2(
      prediction_table,
      file.path(prediction_directory, paste0(seed, ".csv")),
      row.names = FALSE
    )
    
    summary_list[[as.character(seed)]] <- tibble(
      seed = seed,
      n_features_used = length(fit$prepared$features),
      features_used = paste(fit$prepared$features, collapse = ","),
      auc = fit$auc,
      accuracy = fit$metrics$accuracy,
      sensitivity = fit$metrics$sensitivity,
      specificity = fit$metrics$specificity,
      TN = fit$metrics$TN,
      FN = fit$metrics$FN,
      FP = fit$metrics$FP,
      TP = fit$metrics$TP,
      converged = isTRUE(fit$model$converged),
      warning_messages = paste(fit$warnings, collapse = " | ")
    )
    
    coefficient_list[[as.character(seed)]] <- fit$coefficients %>%
      mutate(seed = seed) %>%
      select(seed, term, coefficient)
  }
  
  summary_table <- bind_rows(summary_list)
  coefficient_table <- bind_rows(coefficient_list)
  
  write.csv2(
    summary_table,
    file.path(output_directory, paste0(model_name, "_summary_100seeds.csv")),
    row.names = FALSE
  )
  
  write.csv2(
    coefficient_table,
    file.path(output_directory, paste0(model_name, "_coefficients_100seeds.csv")),
    row.names = FALSE
  )
  
  aggregate_summary <- summary_table %>%
    summarise(
      n_successful_seeds = n(),
      median_auc = median(auc, na.rm = TRUE),
      iqr_auc = IQR(auc, na.rm = TRUE),
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      median_accuracy = median(accuracy, na.rm = TRUE),
      iqr_accuracy = IQR(accuracy, na.rm = TRUE),
      median_sensitivity = median(sensitivity, na.rm = TRUE),
      iqr_sensitivity = IQR(sensitivity, na.rm = TRUE),
      median_specificity = median(specificity, na.rm = TRUE),
      iqr_specificity = IQR(specificity, na.rm = TRUE),
      n_nonconverged = sum(!converged, na.rm = TRUE),
      n_with_warnings = sum(warning_messages != "", na.rm = TRUE)
    )
  
  write.csv2(
    aggregate_summary,
    file.path(output_directory, paste0(model_name, "_aggregate_summary.csv")),
    row.names = FALSE
  )
  
  list(
    summary = summary_table,
    coefficients = coefficient_table,
    aggregate = aggregate_summary
  )
}


# ==============================================================================
# 12) BRANJE IN PRIPRAVA VHODNIH PODATKOV
# ==============================================================================

cat("\n=== Branje vhodnih podatkov ===\n")

integrated_results <- read.csv2(
  FILE_INTEGRATED,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_data <- read.csv(
  FILE_GENES,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

flux_data <- read.csv(
  FILE_FLUX,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assert_columns(
  integrated_results,
  c(
    "gene_id",
    "reaction_id",
    "wilcox_padj_log2TPM",
    "wilcox_padj_flux"
  ),
  "Integrirana preglednica"
)

gene_data <- normalise_sample_table(
  gene_data,
  patient_column = "id_pacienta",
  type_column = "tip_vzorca",
  data_name = "CPTAC_clean_log2.csv"
)

flux_data <- normalise_sample_table(
  flux_data,
  patient_column = "patient_id",
  type_column = "type",
  data_name = "flux_clean.csv"
)

gene_data <- convert_feature_columns_to_numeric(gene_data)
flux_data <- convert_feature_columns_to_numeric(flux_data)

check_unique_patient_type(gene_data, "Genski podatki")
check_unique_patient_type(flux_data, "Podatki metabolnih pretokov")

integrated_results$wilcox_padj_log2TPM <- as_numeric_comma(
  integrated_results$wilcox_padj_log2TPM
)

integrated_results$wilcox_padj_flux <- as_numeric_comma(
  integrated_results$wilcox_padj_flux
)


# ==============================================================================
# 13) IZBOR KANDIDATNIH GENOV IN REAKCIJ
# ==============================================================================

# Kandidatni geni: unikatni gene_id z wilcox_padj_log2TPM < 0,01.
candidate_genes <- integrated_results %>%
  filter(
    !is.na(gene_id),
    gene_id != "",
    !is.na(wilcox_padj_log2TPM),
    wilcox_padj_log2TPM < PADJ_CUTOFF
  ) %>%
  distinct(gene_id) %>%
  pull(gene_id) %>%
  intersect(colnames(gene_data)) %>%
  sort()

# Kandidatne reakcije: unikatni reaction_id z wilcox_padj_flux < 0,01.
candidate_reactions <- integrated_results %>%
  filter(
    !is.na(reaction_id),
    reaction_id != "",
    !is.na(wilcox_padj_flux),
    wilcox_padj_flux < PADJ_CUTOFF
  ) %>%
  distinct(reaction_id) %>%
  pull(reaction_id) %>%
  intersect(colnames(flux_data)) %>%
  sort()

if (length(candidate_genes) == 0) {
  stop("Po filtriranju ni ostal noben kandidatni gen.")
}

if (length(candidate_reactions) == 0) {
  stop("Po filtriranju ni ostala nobena kandidatna reakcija.")
}

candidate_counts <- tibble(
  feature_type = c("gene", "reaction"),
  n_candidates = c(
    length(candidate_genes),
    length(candidate_reactions)
  ),
  padj_cutoff = PADJ_CUTOFF
)

write.csv2(
  candidate_counts,
  file.path(OUTPUT_DIR, "00_input_checks", "candidate_counts.csv"),
  row.names = FALSE
)

write.csv2(
  tibble(gene_id = candidate_genes),
  file.path(OUTPUT_DIR, "00_input_checks", "candidate_genes.csv"),
  row.names = FALSE
)

write.csv2(
  tibble(reaction_id = candidate_reactions),
  file.path(OUTPUT_DIR, "00_input_checks", "candidate_reactions.csv"),
  row.names = FALSE
)

cat("Število kandidatnih genov:", length(candidate_genes), "\n")
cat("Število kandidatnih reakcij:", length(candidate_reactions), "\n")


# ==============================================================================
# 14) IZBOR ISTIH 52 SKUPNIH PARNIH PACIENTOV
# ==============================================================================

selected_patients <- select_common_patients(
  gene_data = gene_data,
  flux_data = flux_data,
  n_patients = N_PATIENTS_TOTAL,
  seed = SINGLE_SEED
)

gene_data_52 <- gene_data %>%
  filter(patient_id %in% selected_patients) %>%
  arrange(patient_id, type)

flux_data_52 <- flux_data %>%
  filter(patient_id %in% selected_patients) %>%
  arrange(patient_id, type)

write.csv2(
  tibble(patient_id = selected_patients),
  file.path(OUTPUT_DIR, "00_input_checks", "selected_52_patients.csv"),
  row.names = FALSE
)

cat("Izbranih skupnih parnih pacientov:", length(selected_patients), "\n")


# ==============================================================================
# 15) PODROBEN LASSO PRI SEEDU 123
# ==============================================================================

cat("\n=== Podroben logistični LASSO pri seedu", SINGLE_SEED, "===\n")

single_split <- make_patient_split(
  patient_ids = selected_patients,
  seed = SINGLE_SEED,
  n_train_patients = N_TRAIN_PATIENTS
)

single_gene_fit <- fit_lasso_one_split(
  data = gene_data_52,
  candidate_features = candidate_genes,
  split = single_split,
  seed = SINGLE_SEED,
  n_inner_folds = N_INNER_FOLDS
)

save_single_seed_lasso(
  fit = single_gene_fit,
  seed = SINGLE_SEED,
  modality_name = "genes",
  output_directory = file.path(
    OUTPUT_DIR,
    "01_single_seed",
    "genes"
  )
)

single_reaction_fit <- fit_lasso_one_split(
  data = flux_data_52,
  candidate_features = candidate_reactions,
  split = single_split,
  seed = SINGLE_SEED,
  n_inner_folds = N_INNER_FOLDS
)

save_single_seed_lasso(
  fit = single_reaction_fit,
  seed = SINGLE_SEED,
  modality_name = "reactions",
  output_directory = file.path(
    OUTPUT_DIR,
    "01_single_seed",
    "reactions"
  )
)


# ==============================================================================
# 16) STABILNOSTNA ANALIZA GENOV PRI 100 SEEDIH
# ==============================================================================

cat("\n=== Stabilnostna LASSO analiza genov ===\n")

gene_stability <- run_lasso_stability(
  data = gene_data_52,
  candidate_features = candidate_genes,
  patient_ids = selected_patients,
  seeds = STABILITY_SEEDS,
  modality_name = "genes",
  feature_id_column = "gene_id",
  output_directory = file.path(
    OUTPUT_DIR,
    "02_gene_lasso"
  )
)


# ==============================================================================
# 17) STABILNOSTNA ANALIZA REAKCIJ PRI 100 SEEDIH
# ==============================================================================

cat("\n=== Stabilnostna LASSO analiza reakcij ===\n")

reaction_stability <- run_lasso_stability(
  data = flux_data_52,
  candidate_features = candidate_reactions,
  patient_ids = selected_patients,
  seeds = STABILITY_SEEDS,
  modality_name = "reactions",
  feature_id_column = "reaction_id",
  output_directory = file.path(
    OUTPUT_DIR,
    "03_reaction_lasso"
  )
)


# ==============================================================================
# 18) IZBOR STABILNIH GENOV IN REAKCIJ S FREKVENCO > 70
# ==============================================================================

# Iz frekvenčne preglednice obdržimo vse gene,
# ki so bili pri lambda.1se izbrani v več kot 70 od 100 ponovitev.
stable_gene_table <- select_stable_features(
  frequency_summary = gene_stability$frequency_summary,
  feature_id_column = "gene_id",
  frequency_cutoff = SELECTION_FREQUENCY_CUTOFF
)

# Iz frekvenčne preglednice obdržimo vse reakcije,
# ki so bile pri lambda.1se izbrane v več kot 70 od 100 ponovitev.
stable_reaction_table <- select_stable_features(
  frequency_summary = reaction_stability$frequency_summary,
  feature_id_column = "reaction_id",
  frequency_cutoff = SELECTION_FREQUENCY_CUTOFF
)

# Iz celotnih preglednic izluščimo vektorja identifikatorjev,
# ki ju bodo uporabili končni logistični modeli.
stable_genes <- stable_gene_table$gene_id
stable_reactions <- stable_reaction_table$reaction_id

# Če pri eni vrsti značilk nobena značilka ne preseže praga,
# modela te vrste in kombiniranega modela ni mogoče postaviti.
if (length(stable_genes) == 0) {
  stop(
    paste0(
      "Noben gen nima frekvence > ",
      SELECTION_FREQUENCY_CUTOFF,
      "."
    )
  )
}

if (length(stable_reactions) == 0) {
  stop(
    paste0(
      "Nobena reakcija nima frekvence > ",
      SELECTION_FREQUENCY_CUTOFF,
      "."
    )
  )
}

# Gene in reakcije združimo v eno preglednico.
# Ohranimo frekvenco ter najmanjši, največji in mediani koeficient
# pri lambda.1se.
selected_stable_feature_table <- bind_rows(
  stable_gene_table %>%
    transmute(
      feature_type = "gene",
      feature_id = gene_id,
      frequency = frequency,
      minimal_coefficient_1se = minimal_coefficient_1se,
      maximal_coefficient_1se = maximal_coefficient_1se,
      median_coefficient_1se = median_coefficient_1se
    ),
  
  stable_reaction_table %>%
    transmute(
      feature_type = "reaction",
      feature_id = reaction_id,
      frequency = frequency,
      minimal_coefficient_1se = minimal_coefficient_1se,
      maximal_coefficient_1se = maximal_coefficient_1se,
      median_coefficient_1se = median_coefficient_1se
    )
) %>%
  arrange(
    desc(frequency),
    feature_type,
    feature_id
  )

# Shranimo skupno preglednico stabilnih genov in reakcij.
write.csv2(
  selected_stable_feature_table,
  file.path(
    OUTPUT_DIR,
    "05_final_summary",
    "stable_features_frequency_gt70.csv"
  ),
  row.names = FALSE
)

cat(
  "Število genov s frekvenco >",
  SELECTION_FREQUENCY_CUTOFF,
  ":",
  length(stable_genes),
  "\n"
)

cat(
  "Število reakcij s frekvenco >",
  SELECTION_FREQUENCY_CUTOFF,
  ":",
  length(stable_reactions),
  "\n"
)

cat("Stabilni geni:\n")
print(stable_genes)

cat("Stabilne reakcije:\n")
print(stable_reactions)


# ==============================================================================
# 19) PRIPRAVA TREH KONČNIH MODELOV
# ==============================================================================

# Pripravimo podatkovno tabelo samo s stabilnimi geni.
stable_gene_model_data <- build_final_model_table(
  gene_data = gene_data_52,
  flux_data = flux_data_52,
  patient_ids = selected_patients,
  gene_features = stable_genes,
  reaction_features = character(0)
)

# Pripravimo podatkovno tabelo samo s stabilnimi reakcijami.
stable_reaction_model_data <- build_final_model_table(
  gene_data = gene_data_52,
  flux_data = flux_data_52,
  patient_ids = selected_patients,
  gene_features = character(0),
  reaction_features = stable_reactions
)

# Pripravimo kombinirano podatkovno tabelo s stabilnimi geni in reakcijami.
combined_model_data <- build_final_model_table(
  gene_data = gene_data_52,
  flux_data = flux_data_52,
  patient_ids = selected_patients,
  gene_features = stable_genes,
  reaction_features = stable_reactions
)


# ==============================================================================
# 20) VALIDACIJA MODELA S STABILNIMI GENI
# ==============================================================================

cat("\n=== Validacija modela s stabilnimi geni (frequency > 70) ===\n")

stable_gene_validation <- run_final_glm_validation(
  model_data = stable_gene_model_data,
  feature_names = stable_genes,
  patient_ids = selected_patients,
  seeds = VALIDATION_SEEDS,
  model_name = "stable_genes_gt70",
  output_directory = file.path(
    OUTPUT_DIR,
    "04_final_models",
    "stable_genes_gt70"
  ),
  standardize = STANDARDIZE_FINAL_GLM
)


# ==============================================================================
# 21) VALIDACIJA MODELA S STABILNIMI REAKCIJAMI
# ==============================================================================

cat("\n=== Validacija modela s stabilnimi reakcijami (frequency > 70) ===\n")

stable_reaction_validation <- run_final_glm_validation(
  model_data = stable_reaction_model_data,
  feature_names = stable_reactions,
  patient_ids = selected_patients,
  seeds = VALIDATION_SEEDS,
  model_name = "stable_reactions_gt70",
  output_directory = file.path(
    OUTPUT_DIR,
    "04_final_models",
    "stable_reactions_gt70"
  ),
  standardize = STANDARDIZE_FINAL_GLM
)


# ==============================================================================
# 22) VALIDACIJA KOMBINIRANEGA MODELA
# ==============================================================================

cat("\n=== Validacija kombiniranega modela s stabilnimi značilkami ===\n")

# Kombinirani model vsebuje vse stabilne gene in vse stabilne reakcije.
combined_features <- c(
  stable_genes,
  stable_reactions
)

combined_validation <- run_final_glm_validation(
  model_data = combined_model_data,
  feature_names = combined_features,
  patient_ids = selected_patients,
  seeds = VALIDATION_SEEDS,
  model_name = "combined_stable_gt70",
  output_directory = file.path(
    OUTPUT_DIR,
    "04_final_models",
    "combined_stable_gt70"
  ),
  standardize = STANDARDIZE_FINAL_GLM
)


# ==============================================================================
# 23) SKUPNI POVZETEK KONČNIH MODELOV
# ==============================================================================

# Združimo agregirane rezultate vseh treh končnih modelov
# v eno preglednico za neposredno primerjavo.
final_model_comparison <- bind_rows(
  stable_gene_validation$aggregate %>%
    mutate(
      model = "stable_genes_gt70",
      n_selected_features = length(stable_genes)
    ),
  
  stable_reaction_validation$aggregate %>%
    mutate(
      model = "stable_reactions_gt70",
      n_selected_features = length(stable_reactions)
    ),
  
  combined_validation$aggregate %>%
    mutate(
      model = "combined_stable_gt70",
      n_selected_features = length(combined_features)
    )
) %>%
  select(
    model,
    n_selected_features,
    everything()
  )

write.csv2(
  final_model_comparison,
  file.path(
    OUTPUT_DIR,
    "05_final_summary",
    "final_model_comparison.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 24) KONEC
# ==============================================================================

cat("\n============================================================\n")
cat("LASSO pipeline V2 z izborom frequency > 70 je končan.\n")
cat("Vsi rezultati so shranjeni v mapi:\n")
cat(normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n")
cat("============================================================\n")