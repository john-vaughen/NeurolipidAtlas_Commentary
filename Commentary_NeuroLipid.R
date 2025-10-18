####Libraries and function####
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)

#### Merging datasets ####
setwd("/Users/jvaughen/Dropbox/Lipidomics/Processed")
ad_path    <- "CSV_lipidomics/AD_neurolipid.csv"
ftd_path   <- "CSV_lipidomics/ftd_neurolipid.csv"
ineu1_path <- "CSV_lipidomics/Isaacs_ineuron.csv"
ineu2_path <- "CSV_lipidomics/VDK_ineuron.csv"

AD    <- read_csv(ad_path,    show_col_types = FALSE)
FTD   <- read_csv(ftd_path,   show_col_types = FALSE)
INEU1 <- read_csv(ineu1_path, show_col_types = FALSE)
INEU2 <- read_csv(ineu2_path, show_col_types = FALSE)

normalize_meta_names <- function(nm) {
  nm0 <- nm
  nm0 <- gsub("\\s+", "_", nm0)
  nm0 <- gsub("[.-]", "_", nm0)
  nm0 <- gsub("__+", "_", nm0)
  nm0 <- trimws(nm0)
  # Canonicalize common ones
  nm0 <- dplyr::case_when(
    tolower(nm0) %in% c("genotype_ftd","genotypeftd","genotype_ftld","genotypeftld") ~ "Genotype_FTD",
    tolower(nm0) %in% c("genotype_ad","genotypead")                                   ~ "Genotype_AD",
    tolower(nm0) == "diagnosis"                                                       ~ "Diagnosis",
    tolower(nm0) == "tissue"                                                          ~ "Tissue",
    tolower(nm0) == "genotype"                                                        ~ "Genotype",
    tolower(nm0) == "lab"                                                             ~ "Lab",
    tolower(nm0) == "cell"                                                            ~ "Cell",
    tolower(nm0) == "source"                                                          ~ "Source",
    TRUE ~ nm0
  )
  nm0
}
normalize_header_dot <- function(h) {
  if (is.na(h) || !nzchar(h)) return(h)
  x <- h
  
  # trim + drop obvious tails
  x <- stringr::str_squish(x)
  x <- sub("\\s*;.*$", "", x, perl = TRUE)
  x <- sub("\\s*\\[[^\\]]*\\]$", "", x, perl = TRUE)
  x <- sub("\\s*\\([^)]*adduct[^)]*\\)$", "", x, ignore.case = TRUE)
  
  # ensure class prefix separated from next token with a dot: "CE 14:0"/"CE.14.0" both fine
  x <- sub("^([A-Za-z]{1,6})[\\s._/:-]+", "\\1.", x, perl = TRUE)
  
  # sphingoid bases d/t/m: "d18.1" or "d18:1" or "d18_1" -> "d18.1"
  x <- gsub("([dtmDTM])(\\d{1,2})[.:_/](\\d{1,2})", "\\1\\2.\\3", x, perl = TRUE)
  
  # generic FA tokens: "18:1" or "18_1" -> "18.1"
  x <- gsub("(\\b\\d{1,2})[:_/](\\d{1,2}\\b)", "\\1.\\2", x, perl = TRUE)
  
  # all remaining separators (space/underscore/slash/colon) -> "."
  x <- gsub("[\\s_/:]+", ".", x, perl = TRUE)
  
  # collapse repeated dots, trim edge dots
  x <- gsub("\\.+", ".", x)
  x <- sub("^\\.", "", x); x <- sub("\\.$", "", x)
  
  # uppercase lipid class prefix (first token)
  x <- sub("^([a-z]{1,6})(\\.)", "\\U\\1\\E\\2", x, perl = TRUE)
  
  x
}
meta_candidates <- c(
  "Tissue","Diagnosis","Genotype_AD","Genotype_FTD","Genotype","Lab","Cell","Source",
  "Subject","Subject_ID","Sample","Sample_ID","ID","Sex","Age","Batch","Plate","Well",
  "Region","Area","Notes","Comment","Group","Class","Filename","File","Run",
  "Acquisition","RetentionTime","RT","mz","m_z","m.z","Unnamed_0","Unnamed: 0"
)

prepare_one <- function(df, source_label) {
  if (!("Source" %in% names(df))) df <- dplyr::mutate(df, Source = source_label, .before = 1)
  names(df) <- normalize_meta_names(names(df))
  
  nm_raw    <- names(df)
  meta_keep <- intersect(meta_candidates, nm_raw)
  lipid_raw <- setdiff(nm_raw, meta_keep)
  
  # raw -> dot-canonical
  lipid_clean <- vapply(lipid_raw, normalize_header_dot, character(1))
  
  tibble(
    raw   = lipid_raw,
    clean = lipid_clean
  ) -> map
  
  list(df = df, meta_keep = meta_keep, map = map)
}

P_AD  <- prepare_one(AD,   "AD")
P_FTD <- prepare_one(FTD,  "FTD")
P_I1  <- prepare_one(INEU1,"iNeuron-I")
P_I2  <- prepare_one(INEU2,"iNeuron-II")

all_lipids <- sort(unique(c(P_AD$map$clean, P_FTD$map$clean, P_I1$map$clean, P_I2$map$clean)))
meta_force <- c("Source","Tissue","Diagnosis","Genotype_AD","Genotype_FTD")
meta_union <- sort(unique(c(P_AD$meta_keep, P_FTD$meta_keep, P_I1$meta_keep, P_I2$meta_keep, meta_force)))

num_block_from_map <- function(df, map, clean_names) {
  out <- matrix(NA_real_, nrow = nrow(df), ncol = length(clean_names),
                dimnames = list(NULL, clean_names))
  # group raw columns by cleaned name
  groups <- split(map$raw, map$clean)  # list: clean -> c(raw1, raw2, ...)
  for (cn in clean_names) {
    raws <- groups[[cn]]
    if (is.null(raws)) next
    raws <- intersect(raws, names(df))
    if (!length(raws)) next
    
    # coerce and sum across duplicates
    mat <- sapply(raws, function(col) {
      v <- df[[col]]
      if (is.character(v)) v <- gsub(",", "", trimws(v))
      suppressWarnings(as.numeric(v))
    })
    if (!is.matrix(mat)) {
      out[, cn] <- mat
    } else {
      out[, cn] <- rowSums(mat, na.rm = TRUE)
      # if a row was all NA across duplicates, turn 0 back into NA
      all_na <- apply(is.na(mat), 1, all)
      out[all_na, cn] <- NA_real_
    }
  }
  as_tibble(out)
}

ensure_meta <- function(df, meta_cols) {
  miss <- setdiff(meta_cols, names(df))
  if (length(miss)) df[miss] <- NA_character_
  df[, meta_cols, drop = FALSE]
}

W_AD  <- bind_cols(ensure_meta(P_AD$df,  meta_union), num_block_from_map(P_AD$df,  P_AD$map,  all_lipids))
W_FTD <- bind_cols(ensure_meta(P_FTD$df, meta_union), num_block_from_map(P_FTD$df, P_FTD$map, all_lipids))
W_I1  <- bind_cols(ensure_meta(P_I1$df,  meta_union), num_block_from_map(P_I1$df,  P_I1$map,  all_lipids))
W_I2  <- bind_cols(ensure_meta(P_I2$df,  meta_union), num_block_from_map(P_I2$df,  P_I2$map,  all_lipids))

neurolipid_wide <- bind_rows(W_AD, W_FTD, W_I1, W_I2)

neurolipid_wide <- neurolipid_wide %>%
  mutate(
    Tissue = if ("Tissue" %in% names(.)) str_squish(Tissue) else NA_character_,
    Tissue = case_when(
      !is.na(Tissue) & Tissue != "" ~ Tissue,
      is.na(Tissue) & Source == "iNeuron-I"  ~ "iNeuron-I",
      is.na(Tissue) & Source == "iNeuron-II" ~ "iNeuron-II",
      TRUE ~ Tissue
    )
  )

meta_order <- c("Source","Tissue","Diagnosis","Genotype_AD","Genotype_FTD")
neurolipid_wide <- neurolipid_wide %>% relocate(any_of(meta_order), .before = 1)

dir.create("CSV_lipidomics", showWarnings = FALSE, recursive = TRUE)
write_csv(neurolipid_wide, "CSV_lipidomics/Merge_neurolipid.csv")

#### confirming data is raw and not yet normalized ####
merged_path <- "CSV_lipidomics/Merge_neurolipid.csv"
neu <- read_csv(merged_path, show_col_types = FALSE) %>%
  mutate(
    Tissue       = str_squish(Tissue),
    Diagnosis    = str_squish(Diagnosis),
    Diagnosis_up = toupper(Diagnosis),
    Condition = case_when(
      Source == "AD"  & Diagnosis_up == "CONTROL" & Tissue == "PFC (grey)"  ~ "PFC (grey)",
      Source == "AD"  & Diagnosis_up == "CONTROL" & Tissue == "PFC (white)" ~ "PFC (white)",
      Source == "FTD" & Diagnosis_up %in% c("FTD","FTLD") & Tissue == "PFC (grey)" ~ "FTLD (PFC grey)",
      Source == "AD"  & Diagnosis_up == "AD"     & Tissue == "PFC (grey)"  ~ "AD (PFC grey)",
      Source == "AD"  & Diagnosis_up == "AD"     & Tissue == "PFC (white)" ~ "AD (PFC white)",
      Source == "iNeuron-I"  ~ "iNeuron-I",
      Source == "iNeuron-II" ~ "iNeuron-II",
      TRUE ~ NA_character_
    )
  ) %>% filter(!is.na(Condition))

meta_cols <- intersect(c(
  "Tissue","Diagnosis","Diagnosis_up","Genotype_AD","Genotype_FTD","Genotype",
  "Lab","Cell","Source","Condition","Subject","Subject_ID","Sample","Sample_ID","ID",
  "Sex","Age","Batch","Plate","Well","Region","Area","Notes","Comment","Group","Class",
  "Filename","File","Run","Acquisition","RetentionTime","RT","mz","m.z","Unnamed: 0"
), names(neu))
all_features <- setdiff(names(neu), meta_cols)

# Convenience coercion
as_num_mat <- function(df, cols) {
  M <- df %>% select(all_of(cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  M
}

gpl_regex      <- "(?i)^(PC|PE|PI|PS)\\."                     # all GPL
ether_only     <- "(?i)^(PC|PE|PI|PS)\\.(?:O|P)(?:\\.|$)"     # ether forms
sphingo_regex  <- "(?i)^(HexCer|GlcCer|GalCer|LacCer|SM|Cer(?:amide)?)\\."
neutral_regex  <- "(?i)^(DG|TAG|TG|Triacylglycer|CE|Cholesteryl\\.Ester)\\."

pick_features <- function(include_regex, exclude_regex = NULL) {
  cols <- grep(include_regex, all_features, value = TRUE, perl = TRUE)
  if (!is.null(exclude_regex)) cols <- cols[!grepl(exclude_regex, cols, perl = TRUE)]
  cols
}

subsets <- list(
  ALL           = list(cols = all_features,                               label = "All features"),
  GPL_nonether  = list(cols = pick_features(gpl_regex, ether_only),        label = "GPL (non-ether)"),
  Ether_GPL     = list(cols = pick_features(ether_only),                   label = "Ether GPL"),
  Sphingolipids = list(cols = pick_features(sphingo_regex),                label = "Sphingolipids"),
  Neutral       = list(cols = pick_features(neutral_regex),                label = "Neutral lipids")
)

check_closure <- function(M, name, tol = 1e-6) {
  rs <- rowSums(M, na.rm = TRUE)
  # Basic stats
  rng  <- range(rs, na.rm = TRUE)
  sdv  <- sd(rs, na.rm = TRUE)
  mn   <- mean(rs, na.rm = TRUE)
  cv   <- sdv / mn
  n    <- sum(is.finite(rs))
  
  # Are many rows ~1 or ~100 (common pre-normalized targets)?
  near1   <- mean(abs(rs - 1)   <= 1e-3, na.rm = TRUE)   # fraction near 1
  near100 <- mean(abs(rs - 100) <= 1e-1, na.rm = TRUE)   # fraction near 100
  
  cat(sprintf("\n[%s] rows=%d | sum range = [%.6f, %.6f] | mean=%.6f sd=%.6f CV=%.4f\n",
              name, n, rng[1], rng[2], mn, sdv, cv))
  cat(sprintf("   Fraction near 1.0:   %.2f\n", near1))
  cat(sprintf("   Fraction near 100.0: %.2f\n", near100))
  
  invisible(list(row_sums = rs, stats = c(mean = mn, sd = sdv, cv = cv)))
}


for (nm in names(subsets)) {
  cols <- subsets[[nm]]$cols
  if (!length(cols)) {
    cat(sprintf("\n[%s] No matching columns.\n", nm)); next
  }
  M <- as_num_mat(neu, cols)
  check_closure(M, nm)
}


plot_rowsums <- function(M, title) {
  df <- tibble(RowSum = rowSums(M, na.rm = TRUE))
  ggplot(df, aes(RowSum)) +
    geom_histogram(bins = 40, alpha = 0.85) +
    labs(x = "Row total", y = "Count", title = paste("Row-sum distribution —", title)) +
    theme_classic()
}


p_all <- plot_rowsums(as_num_mat(neu, subsets$ALL$cols), "All features")
p_gpl <- plot_rowsums(as_num_mat(neu, subsets$GPL_nonether$cols), "GPL (non-ether)")
print(p_all); print(p_gpl)

summarize_totals <- function(M, by = c("Source","Condition")) {
  rs <- rowSums(M, na.rm = TRUE)
  tibble(Total = rs) %>%
    bind_cols(neu[, by, drop = FALSE]) %>%
    group_by(across(all_of(by))) %>%
    summarize(n = n(), mean = mean(Total, na.rm = TRUE),
              sd = sd(Total, na.rm = TRUE),
              cv = sd/mean, .groups = "drop")
}

tot_all  <- summarize_totals(as_num_mat(neu, subsets$ALL$cols))
tot_gpl  <- summarize_totals(as_num_mat(neu, subsets$GPL_nonether$cols))
tot_neut <- summarize_totals(as_num_mat(neu, subsets$Neutral$cols))

print(tot_all, n = Inf)
print(tot_gpl, n = Inf)
print(tot_neut, n = Inf)





#### PCA — different lipid classes, testing 3 normalizations (PQN / SUM / CLR) ####

neu <- read.csv("CSV_lipidomics/Merge_neurolipid.csv", check.names = FALSE) %>%
  mutate(
    Tissue       = stringr::str_squish(Tissue),
    Diagnosis    = stringr::str_squish(Diagnosis),
    Diagnosis_up = toupper(Diagnosis)
  ) %>%
  # --- NEW: drop Cerebellum entirely ---
  dplyr::filter(!grepl("(?i)cerebell", Tissue)) %>%
  # --- UPDATED: label FTLD controls explicitly ---
  mutate(
    Condition = dplyr::case_when(
      Source == "AD"  & Diagnosis_up == "CONTROL" & Tissue == "PFC (grey)"  ~ "PFC grey Control",
      Source == "AD"  & Diagnosis_up == "CONTROL" & Tissue == "PFC (white)" ~ "PFC white Control",
      Source == "FTD" & Diagnosis_up %in% c("FTD","FTLD") & Tissue == "PFC (grey)" ~ "PFC grey FTLD",
      Source == "FTD" & Diagnosis_up == "CONTROL" & Tissue %in% c("PFC (grey)","PFC (white)") ~ "PFC grey Control1",
      Source == "AD"  & Diagnosis_up == "AD"     & Tissue == "PFC (grey)"  ~ "PFC grey AD",
      Source == "AD"  & Diagnosis_up == "AD"     & Tissue == "PFC (white)" ~ "PFC white AD",
      Source == "iNeuron-I"  ~ "iNeuron-I",
      Source == "iNeuron-II" ~ "iNeuron-II",
      TRUE ~ NA_character_
    )
  )

#if want to exclude any conditions in plot, do this:
#neu <- neu %>% dplyr::filter(Condition != "PFC grey Control1")


meta_cols <- intersect(c(
  "Tissue","Diagnosis","Diagnosis_up","Genotype_AD","Genotype_FTD","Genotype",
  "Lab","Cell","Source","Condition","Subject","Subject_ID","Sample","Sample_ID","ID",
  "Sex","Age","Batch","Plate","Well","Region","Area","Notes","Comment","Group","Class",
  "Filename","File","Run","Acquisition","RetentionTime","RT","mz","m.z","Unnamed: 0"
), names(neu))
all_features <- setdiff(names(neu), meta_cols)

sphingo_regex <- "(?i)^(HexCer|GlcCer|GalCer|LacCer|SM|Cer(?:amide)?)\\b"
neutral_regex <- "(?i)^(DG|Diacylglycer|TG|TAG|Triacylglycer|CE|Cholesteryl\\s*Ester)\\b"
gpl_regex     <- "(?i)^(PC|PE|PI|PS)\\b"
ether_regex   <- "(?i)^(PC|PE|PI|PS)\\.[OP](?:\\.|\\d)"  # PC.O... / PC.P...

# Column selections
sphingo_cols <- grep(sphingo_regex, all_features, value = TRUE, perl = TRUE)
neutral_cols <- grep(neutral_regex, all_features, value = TRUE, perl = TRUE)
gpl_all      <- grep(gpl_regex,   all_features, value = TRUE, perl = TRUE)
gpl_ether    <- grep(ether_regex, all_features, value = TRUE, perl = TRUE)
gpl_nonether <- setdiff(gpl_all, gpl_ether)

# Subset definitions (titles are correct here)
subsets <- list(
  Sphingo = list(cols = sphingo_cols, stub = "PCA_sphingolipids", title = "Sphingolipids"),
  Neutral = list(cols = neutral_cols, stub = "PCA_neutral",       title = "Neutral lipids"),
  GPL     = list(cols = gpl_nonether, stub = "PCA_GPL",           title = "Glycerophospholipids"),
  Ether   = list(cols = gpl_ether,    stub = "PCA_ether",         title = "Ether phospholipids")
)

# --- UPDATED: include FTLD control in preferred legend order ---
pref <- c("PFC (grey)", "AD (PFC grey)", "FTLD (PFC grey)", "FTLD control",
          "PFC (white)", "AD (PFC white)", "iNeuron-I", "iNeuron-II")
present <- intersect(pref, unique(neu$Condition))
extras  <- setdiff(unique(neu$Condition), present)
cond_order <- c(present, sort(extras))

# --- UPDATED: add a dedicated light grey for FTLD control ---
#base_cols <- c(
#  "PFC (grey)"      = "#7FB3D5",
#  "AD (PFC grey)"   = "#3F88C5",
#  "FTLD (PFC grey)" = "#1F5AA6",
#  "FTLD control"    = "#B0B0B0",   # light grey for FTLD controls
#  "PFC (white)"     = "#7FC8A9",
#  "AD (PFC white)"  = "#2CA25F",
#  "iNeuron-I"       = "#F3A530",
#  "iNeuron-II"      = "#D97A00"
#)
base_cols <- c(
  # Grey matter (controls lighter, disease darker)
  "PFC grey Control"      = "#56B4E9",  # light blue (control)
  "PFC grey AD"   = "#0072B2",  # dark blue (AD)
  "PFC grey FTLD" = "#B2ABD2",  # purple (FTLD case)
  "PFC grey Control1"    = "#CC79A7",  # light lavender (FTLD control)
  
  # White matter (greenish)
  "PFC white Control"     = "#7BC8B3",  # lighter teal/green (control)
  "PFC white AD"  = "#009E73",  # darker bluish-green (AD)
  
  # iPSC / iNeuron lines (oranges)
  "iNeuron-I"       = "#E69F00",  # orange
  "iNeuron-II"      = "#D55E00"   # vermilion
)


missing <- setdiff(cond_order, names(base_cols))
cond_colors <- c(base_cols, setNames(rep("#777777", length(missing)), missing))
cond_colors <- cond_colors[cond_order]

sanitize_for_pca <- function(Z, max_na_frac_col = 0.98, nzv_eps = 1e-12) {
  Z <- as.matrix(Z)
  Z[!is.finite(Z)] <- NA_real_
  all_nonfinite <- apply(Z, 2, function(v) all(is.na(v)))
  too_many_na   <- apply(Z, 2, function(v) mean(is.na(v)) > max_na_frac_col)
  nzv_col <- apply(Z, 2, function(v) {
    vv <- v[is.finite(v)]
    if (length(vv) < 2) return(TRUE)
    (sd(vv) <= nzv_eps) || (length(unique(vv)) < 2)
  })
  keep_cols <- !(all_nonfinite | too_many_na | nzv_col)
  Zc <- Z[, keep_cols, drop = FALSE]
  if (!ncol(Zc)) stop("All features invalid after column sanitization.")
  if (anyNA(Zc)) {
    for (j in seq_len(ncol(Zc))) {
      v <- Zc[, j]
      if (anyNA(v)) {
        m <- stats::median(v, na.rm = TRUE)
        if (!is.finite(m)) m <- 0
        v[is.na(v)] <- m
        Zc[, j] <- v
      }
    }
  }
  keep_rows <- apply(Zc, 1, function(r) all(is.finite(r)))
  Zcr <- Zc[keep_rows, , drop = FALSE]
  if (!nrow(Zcr)) stop("All samples invalid after row sanitization.")
  list(Z = Zcr, keep_rows = keep_rows, keep_cols = keep_cols)
}

normalize_block <- function(X, method = c("clr","pqn","median","sum","percent","quantile","vsn","none"),
                            ref = c("median","first","given"), ref_vector = NULL) {
  method <- match.arg(method)
  ref    <- match.arg(ref)
  X <- as.matrix(X); storage.mode(X) <- "double"
  all_na <- apply(X, 2, function(v) all(is.na(v)))
  if (any(all_na)) X <- X[, !all_na, drop = FALSE]
  X_na <- X
  X[is.na(X)] <- 0
  
  get_ref <- function(M) {
    if (ref == "median") {
      rv <- apply(M, 2, function(v) median(v, na.rm = TRUE))
    } else if (ref == "first") {
      rv <- M[1, ]
    } else {
      if (is.null(ref_vector)) stop("ref_vector must be provided when ref='given'")
      rv <- ref_vector
    }
    eps <- 1e-12
    rv[!is.finite(rv)] <- NA_real_
    if (anyNA(rv)) {
      med <- apply(X_na, 2, function(v) median(v, na.rm = TRUE))
      idx <- which(is.na(rv)); rv[idx] <- med[idx]
    }
    rv[!is.finite(rv) | rv == 0] <- eps
    rv
  }
  
  if (method == "none") return(X_na)
  
  if (method %in% c("sum","percent")) {
    rs <- rowSums(X)                # NAs treated as 0
    rs_safe <- ifelse(rs > 0, rs, NA_real_)
    Z <- sweep(X, 1, rs_safe, "/")
    if (method == "percent") Z <- 100 * Z
    return(Z)
  }
  
  if (method %in% c("median","pqn")) {
    ref_s <- get_ref(X_na)
    Z <- X_na
    for (i in seq_len(nrow(Z))) {
      v <- Z[i, ]
      r <- v / ref_s
      r <- r[is.finite(r)]
      if (!length(r)) { Z[i, ] <- v; next }
      f <- stats::median(r, na.rm = TRUE)
      if (!is.finite(f) || f == 0) f <- 1
      Z[i, ] <- v / f
    }
    return(Z)
  }
  
  if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE))
      stop("Install preprocessCore for quantile normalization: install.packages('preprocessCore')")
    keep_cols <- apply(X_na, 2, function(v) any(is.finite(v)))
    Zt <- preprocessCore::normalize.quantiles(t(X_na[, keep_cols, drop = FALSE]))
    Z  <- matrix(NA_real_, nrow = nrow(X_na), ncol = ncol(X_na), dimnames = dimnames(X_na))
    Z[, keep_cols] <- t(Zt)
    return(Z)
  }
  
  if (method == "vsn") {
    if (!requireNamespace("vsn", quietly = TRUE))
      stop("Install vsn: install.packages('vsn')")
    keep_cols <- apply(X_na, 2, function(v) any(is.finite(v)))
    fit <- vsn::vsn2(X_na[, keep_cols, drop = FALSE])
    Z   <- matrix(NA_real_, nrow = nrow(X_na), ncol = ncol(X_na), dimnames = dimnames(X_na))
    Z[, keep_cols] <- vsn::predict(fit, X_na[, keep_cols, drop = FALSE])
    return(Z)
  }
  
  if (method == "clr") {
    Xc <- X_na
    nonpos_all <- apply(Xc, 2, function(v) all(!is.finite(v) | v <= 0, na.rm = TRUE))
    if (any(nonpos_all)) Xc <- Xc[, !nonpos_all, drop = FALSE]
    Xc[is.na(Xc)] <- 0
    min_pos <- suppressWarnings(min(Xc[Xc > 0], na.rm = TRUE))
    pseudo  <- if (is.finite(min_pos)) min_pos * 0.5 else 1e-9
    Xc[Xc <= 0] <- pseudo
    rs <- rowSums(Xc)
    rs[!is.finite(rs) | rs <= 0] <- NA_real_
    P  <- sweep(Xc, 1, rs, "/")
    LP <- log(P)
    Z  <- sweep(LP, 1, rowMeans(LP, na.rm = TRUE), "-")
    return(Z)
  }
  
  stop("Unknown method.")
}

run_pca_norm <- function(methods = c("pqn","sum","clr"), ref = "median", ref_vector = NULL) {
  
  outdir <- "Graphs/2025/PUFA_Paper/Commentary/PCA_norms"
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  for (sname in names(subsets)) {
    feat_cols <- subsets[[sname]]$cols
    if (!length(feat_cols)) { message(sname, ": no matching columns."); next }
    
    X <- neu %>%
      dplyr::select(dplyr::all_of(feat_cols)) %>%
      dplyr::mutate(dplyr::across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
      as.matrix()
    
    # Pre-filter: drop columns that are ALL NA or ALL zero
    all_na   <- apply(X, 2, function(v) all(is.na(v)))
    all_zero <- apply(X, 2, function(v) all(is.na(v) | v == 0))
    keep0    <- !(all_na | all_zero)
    X        <- X[, keep0, drop = FALSE]
    if (!ncol(X)) { message(sname, ": no usable features after prefilter."); next }
    
    readr::write_lines(colnames(X), file.path(outdir, paste0(subsets[[sname]]$stub, "_features.txt")))
    
    for (m in methods) {
      Z0 <- normalize_block(X, method = m, ref = ref, ref_vector = ref_vector)
      
      cat(sprintf("\n[%s / %s] pre-sanitize: rows=%d cols=%d | any nonfinite? %s\n",
                  sname, toupper(m), nrow(Z0), ncol(Z0), any(!is.finite(Z0))))
      
      S <- sanitize_for_pca(Z0)
      Z <- S$Z
      
      cat(sprintf("[%s / %s] post-sanitize: rows=%d cols=%d | any nonfinite? %s\n",
                  sname, toupper(m), nrow(Z), ncol(Z), any(!is.finite(Z))))
      
      # Guarded Condition vector (handles missing or filtered rows safely)
      cond_all <- if ("Condition" %in% names(neu)) neu$Condition else rep(NA_character_, nrow(neu))
      cond_vec <- factor(cond_all[S$keep_rows], levels = cond_order)
      
      pc <- prcomp(Z, center = TRUE, scale. = TRUE)
      var_exp  <- (pc$sdev^2) / sum(pc$sdev^2)
      pc12_lab <- paste0("PC1 (", scales::percent(var_exp[1], accuracy = 0.1), ")")
      pc22_lab <- paste0("PC2 (", scales::percent(var_exp[2], accuracy = 0.1), ")")
      
      scores <- as.data.frame(pc$x[, 1:2])
      names(scores) <- c("PC1","PC2")
      scores$Condition <- cond_vec
      
      # ---- (optional) quick diagnostics for NA color points ----
      diag_cols <- intersect(c(
        "Source","Tissue","Diagnosis","Diagnosis_up",
        "APOE","Genotype","Genotype_AD","Genotype_FTD",
        "Subject","Subject_ID","Sample","Sample_ID","ID",
        "Region","Area","Batch","Plate","Well","Lab","Notes","Comment"
      ), names(neu))
      kept_meta <- neu[S$keep_rows, diag_cols, drop = FALSE]
      diag_df <- kept_meta %>%
        dplyr::mutate(Condition = as.character(cond_vec)) %>%
        dplyr::bind_cols(scores[, c("PC1","PC2")])
      grey_df <- diag_df %>% dplyr::filter(is.na(Condition))
      message(sprintf("[Diag] Grey points in %s / %s: %d", sname, toupper(m), nrow(grey_df)))
      if (nrow(grey_df) > 0) {
        readr::write_csv(
          grey_df,
          file.path(outdir, paste0(subsets[[sname]]$stub, "_", toupper(m), "_GREY_points_metadata.csv"))
        )
      }
      # ---- end diagnostics ----
      
      # Title override (typo-proof)
      pretty_titles <- c(GPL = "Phospholipids", Ether = "Ether lipids",
                         Neutral = "Neutral lipids", Sphingo = "Sphingolipids")
      title_raw  <- subsets[[sname]]$title
      title_safe <- if (sname %in% names(pretty_titles)) pretty_titles[[sname]] else
        if (is.character(title_raw) && nzchar(title_raw)) title_raw else sname
      
      p <- ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2, color = Condition)) +
        ggplot2::geom_point(size = 2.8, alpha = 0.9) +
        ggplot2::scale_color_manual(values = cond_colors, limits = cond_order, breaks = cond_order, drop = FALSE, name = NULL) +
        ggplot2::labs(x = pc12_lab, y = pc22_lab, title = paste0(title_safe, " — ", toupper(m))) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill = "white", colour = NA),
          plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
          panel.grid       = ggplot2::element_blank(),
          axis.ticks       = ggplot2::element_blank(),
          axis.text.x      = ggplot2::element_blank(),
          axis.text.y      = ggplot2::element_blank(),
          plot.title       = ggplot2::element_text(size = 26, color = "black", face = "plain"),
          axis.title.x     = ggplot2::element_text(size = 28, color = "black", face = "plain", margin = ggplot2::margin(t = 8)),
          axis.title.y     = ggplot2::element_text(size = 28, color = "black", face = "plain", margin = ggplot2::margin(r = 8)),
          legend.background= ggplot2::element_rect(fill = "white", colour = NA),
          legend.key       = ggplot2::element_rect(fill = "white", colour = NA),
          legend.position  = "none"
        )
      
      pngfile <- file.path(outdir, paste0(subsets[[sname]]$stub, "_", toupper(m), "_pc1_pc2.png"))
      ggplot2::ggsave(pngfile, p, width = 5, height = 5, dpi = 300)
      print(p)
      
      load_pc1 <- tibble::tibble(Feature = colnames(Z), PC1_Loading = pc$rotation[, 1]) %>%
        dplyr::mutate(Abs_PC1 = abs(PC1_Loading)) %>%
        dplyr::arrange(dplyr::desc(Abs_PC1))
      readr::write_csv(head(load_pc1, 10),
                       file.path(outdir, paste0(subsets[[sname]]$stub, "_", toupper(m), "_top10_PC1.csv")))
      readr::write_csv(load_pc1,
                       file.path(outdir, paste0(subsets[[sname]]$stub, "_", toupper(m), "_PC1_loadings_sorted.csv")))
    }
  }
}

theme_set(theme_get() + theme(
  legend.text  = element_text(size = 12),  # bigger legend labels
  legend.key.size = unit(6, "pt")
))

# Run
run_pca_norm(methods = c("pqn","sum","clr"), ref = "median")

#### DHA % plotting  ####


plot_metric_points_sem_by_region <- function(
    long_df,               # columns: Region, METRIC, LIPID_CLASS, Value
    metric,
    custom_colors,
    regions_order,
    x_axis_labels = NULL,  # named vector: names = raw region keys, values = pretty labels
    y_mode = c("smart","zero","free"),
    sem_mult = 1.5,
    y_pad_frac = 0.08,
    point_size = 3,
    title_map = c(
      Prop_Omega3="% ω3", Prop_Omega6="% ω6", Omega3_to_Omega6="ω3/ω6",
      Prop_PUFA="% PUFA", Prop_MUFA="% MUFA", Prop_SFA="% SFA",
      Weighted_DB_Sum="Double bonds", Weighted_Chain_Sum="Chain length"
    ),
    y_label_map = c(
      Prop_Omega3="%", Prop_Omega6="%", Omega3_to_Omega6="ω3/ω6",
      Prop_PUFA="%", Prop_MUFA="%", Prop_SFA="%",
      Weighted_DB_Sum="#", Weighted_Chain_Sum="Length"
    )
){
  library(dplyr); library(ggplot2)
  
  y_mode <- match.arg(y_mode)
  df <- long_df %>% filter(METRIC == metric)
  
  # enforce region order
  df$Region <- factor(df$Region, levels = regions_order, ordered = TRUE)
  
  # summarise mean & SEM per class per region
  grand <- df %>%
    group_by(LIPID_CLASS, Region) %>%
    summarise(mean = mean(Value, na.rm = TRUE),
              sem  = sd(Value, na.rm = TRUE)/sqrt(n()),
              .groups = "drop") %>%
    filter(!is.na(mean))
  
  if (!nrow(grand)) {
    the_title  <- if (metric %in% names(title_map)) title_map[[metric]] else gsub("_"," ", metric)
    the_ylabel <- if (metric %in% names(y_label_map)) y_label_map[[metric]] else ""
    return(ggplot() + theme_void() + labs(title = paste0(the_title, " — no data"), y = the_ylabel, x = ""))
  }
  
  # colors
  present_classes <- sort(unique(grand$LIPID_CLASS))
  cols_use <- custom_colors[intersect(names(custom_colors), present_classes)]
  if (length(setdiff(present_classes, names(cols_use))) > 0) {
    missing <- setdiff(present_classes, names(cols_use))
    cols_use <- c(cols_use, setNames(rep("#333333", length(missing)), missing))
  }
  
  # x labels
  lvl <- levels(grand$Region); if (is.null(lvl)) lvl <- unique(grand$Region)
  xlab_vec <- if (!is.null(x_axis_labels)) unname(x_axis_labels[lvl]) else as.character(lvl)
  xlab_vec[is.na(xlab_vec)] <- as.character(lvl[is.na(xlab_vec)])
  
  # y-limits
  get_smart_limits <- function(grand_df) {
    mu <- grand_df$mean; se <- grand_df$sem; se[is.na(se)] <- 0
    lo <- suppressWarnings(min(mu - sem_mult * se, na.rm = TRUE))
    hi <- suppressWarnings(max(mu + sem_mult * se, na.rm = TRUE))
    if (!is.finite(lo) || !is.finite(hi)) {
      lo <- suppressWarnings(min(mu, na.rm = TRUE)); hi <- suppressWarnings(max(mu, na.rm = TRUE))
    }
    rng <- hi - lo
    if (!is.finite(rng) || rng == 0) {
      pad <- max(0.05 * abs(hi), 0.5); return(c(lo - pad, hi + pad))
    }
    c(lo - y_pad_frac * rng, hi + y_pad_frac * rng)
  }
  is_percent <- grepl("^Prop_", metric)
  is_ratio   <- grepl("Omega3_to_Omega6", metric)
  is_weight  <- grepl("^Weighted_", metric)
  y_limits <- switch(y_mode,
                     zero  = c(0, NA),
                     free  = NULL,
                     smart = if (is_weight) get_smart_limits(grand) else c(0, NA)
  )
  
  the_title  <- if (metric %in% names(title_map)) title_map[[metric]] else gsub("_"," ", metric)
  the_ylabel <- if (metric %in% names(y_label_map)) y_label_map[[metric]] else ""
  
  ggplot() +
    geom_errorbar(
      data = grand,
      aes(x = Region, ymin = mean - sem, ymax = mean + sem, color = LIPID_CLASS, group = LIPID_CLASS),
      width = 0.18, linewidth = 0.6
    ) +
    geom_point(
      data = grand,
      aes(x = Region, y = mean, color = LIPID_CLASS, group = LIPID_CLASS),
      size = point_size
    ) +
    scale_x_discrete(drop = FALSE, labels = setNames(xlab_vec, lvl)) +
    scale_color_manual(values = cols_use, limits = present_classes, drop = FALSE) +
    labs(x = "", y = the_ylabel, title = the_title) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 26, face="plain"),
      axis.text.y = element_text(size = 26, face="plain"),
      plot.title  = element_text(hjust = 0.5, size = 30, face="plain"),
      axis.title.y = element_text(size = 30, face="plain")
    ) +
    { if (!is.null(y_limits)) coord_cartesian(ylim = y_limits) }
}

plot_gpl_metrics_points_many_by_region <- function(
    long_df, metrics = NULL, custom_colors, regions_order, x_axis_labels = NULL,
    y_mode = "smart", save_dir = NULL, width = 7, height = 5, dpi = 300
){
  if (is.null(metrics)) {
    metrics <- sort(unique(long_df$METRIC))
  }
  if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  plots <- list()
  for (m in metrics) {
    p <- plot_metric_points_sem_by_region(
      long_df, metric = m, custom_colors = custom_colors,
      regions_order = regions_order, x_axis_labels = x_axis_labels,
      y_mode = y_mode
    ) + theme(axis.text.x = element_text(size = 22, face = "plain"))
    plots[[m]] <- p
    if (!is.null(save_dir)) ggsave(file.path(save_dir, paste0(m, ".png")),
                                   p, width = width, height = height, dpi = dpi)
  }
  plots
}

# 1) Helpers that parse current wide colnames
is_chain_species_name <- function(nm) {
  nm <- trimws(nm)
  # Form A: "PC 16:0_22:6" or with dots "PC 16.0_22.6"
  a <- grepl("^([A-Za-z]{2,3})\\s*[0-9]{1,2}[:\\.]\\d{1,2}[_/][0-9]{1,2}[:\\.]\\d{1,2}$", nm, perl = TRUE)
  # Form B: "PC.16.0.22.6"
  b <- grepl("^([A-Za-z]{2,3})\\.\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}$", nm, perl = TRUE)
  a || b
}

lipid_class_of_col <- function(nm) toupper(sub("^([A-Za-z]{2,3}).*$", "\\1", nm, perl = TRUE))

extract_two_chains <- function(nm) {
  nm <- trimws(nm)
  # A) "PC 16:0_22:6" (allow dots for ':')
  m1 <- regexec("^([A-Za-z]{2,3})\\s*([0-9]{1,2}[:\\.]\\d{1,2})[_/]([0-9]{1,2}[:\\.]\\d{1,2})$", nm, perl = TRUE)
  mm1 <- regmatches(nm, m1)[[1]]
  if (length(mm1) == 4) {
    ch1 <- gsub("\\.", ":", mm1[2 + 0]) # chains positions 2 & 3 in mm1
    ch2 <- gsub("\\.", ":", mm1[2 + 1])
    return(c(ch1, ch2))
  }
  # B) "PC.16.0.22.6"
  m2 <- regexec("^([A-Za-z]{2,3})\\.(\\d{1,2})\\.(\\d{1,2})\\.(\\d{1,2})\\.(\\d{1,2})$", nm, perl = TRUE)
  mm2 <- regmatches(nm, m2)[[1]]
  if (length(mm2) == 6) {
    ch1 <- sprintf("%s:%s", mm2[3], mm2[4])
    ch2 <- sprintf("%s:%s", mm2[5], mm2[6])
    return(c(ch1, ch2))
  }
  c(NA_character_, NA_character_)
}

is_dha_colname <- function(nm) {
  if (!is_chain_species_name(nm)) return(FALSE)
  ch <- extract_two_chains(nm)
  any(ch %in% c("22:6", "22.6", "22:06", "22.06")) || any(gsub("\\.", ":", ch) == "22:6")
}

# 2) Work directly from neurolipid_wide
stopifnot(exists("neurolipid_wide"))

all_cols <- names(neurolipid_wide)
gpl_cols <- grep("^(?i)(PC|PE|PI|PS)", all_cols, value = TRUE, perl = TRUE)

# keep only chain-resolved species (skip sums like PC(38:6))
chain_cols <- gpl_cols[sapply(gpl_cols, is_chain_species_name)]
if (!length(chain_cols)) stop("No chain-resolved GPL columns found in neurolipid_wide.")

classes_vec <- setNames(vapply(chain_cols, lipid_class_of_col, character(1)), chain_cols)
dha_mask    <- setNames(vapply(chain_cols, is_dha_colname, logical(1)), chain_cols)

class_levels <- c("PC","PE","PI","PS")

# 3) Compute %DHA per class for each row (safe numeric conversion + NA handling)
vals_mat <- neurolipid_wide[, chain_cols, drop = FALSE]
for (j in seq_along(chain_cols)) {
  v <- vals_mat[[j]]
  if (is.character(v)) v <- gsub(",", "", v)
  vals_mat[[j]] <- suppressWarnings(as.numeric(v))
}

pctDHA_df <- lapply(class_levels, function(cls) {
  in_cls <- names(classes_vec)[classes_vec == cls]
  if (!length(in_cls)) return(rep(NA_real_, nrow(vals_mat)))
  den <- rowSums(as.matrix(vals_mat[, in_cls, drop = FALSE]), na.rm = TRUE)
  num <- rowSums(as.matrix(vals_mat[, in_cls[dha_mask[in_cls]], drop = FALSE]), na.rm = TRUE)
  out <- ifelse(den > 0, 100 * num / den, NA_real_)
  out
})
names(pctDHA_df) <- paste0(class_levels, "_pctDHA")
pct_tbl <- as_tibble(pctDHA_df)

# 4) Build metrics table and flags (using only columns that exist)
safe_col <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NA
neurolipid_metrics <- bind_cols(
  tibble(
    Tissue       = safe_col(neurolipid_wide,"Tissue"),
    Diagnosis    = safe_col(neurolipid_wide,"Diagnosis"),
    Genotype_AD  = safe_col(neurolipid_wide,"Genotype_AD"),
    Genotype_FTD = safe_col(neurolipid_wide,"Genotype_FTD"),
    Genotype     = safe_col(neurolipid_wide,"Genotype"),
    Lab          = safe_col(neurolipid_wide,"Lab"),
    Cell         = safe_col(neurolipid_wide,"Cell"),
    Source       = safe_col(neurolipid_wide,"Source")
  ),
  pct_tbl
)

# 5) Region flags as you had (lightly robust)
neu_flags <- neurolipid_metrics %>%
  mutate(
    Genotype_AD = str_squish(Genotype_AD),
    Tissue      = str_squish(Tissue),
    Diagnosis_norm = case_when(
      !is.na(Diagnosis) & toupper(str_squish(Diagnosis)) %in% c("CONTROL","CTRL") ~ "Control",
      !is.na(Diagnosis) & toupper(str_squish(Diagnosis)) %in% c("FTD","FTLD")    ~ "FTD",
      TRUE ~ NA_character_
    ),
    grp_ctrl_grey  = (Source == "AD"  & Diagnosis_norm == "Control" & Tissue == "PFC (grey)"),
    grp_ctrl_white = (Source == "AD"  & Diagnosis_norm == "Control" & Tissue == "PFC (white)"),
    grp_ftld_grey  = (Source == "FTD" & Diagnosis_norm == "FTD"     & Tissue == "PFC (grey)"),
    grp_apoe2_grey = (Source == "AD" & !is.na(Genotype_AD) &
                        str_detect(Genotype_AD, regex("\\bApoE2\\b", ignore_case = TRUE)) &
                        Tissue == "PFC (grey)"),
    grp_apoe4_grey = (Source == "AD" & !is.na(Genotype_AD) &
                        !str_detect(Genotype_AD, regex("\\bApoE2\\b", ignore_case = TRUE)) &
                        str_detect(Genotype_AD, regex("^ApoE\\s*(3/4|4/4)\\s*$", ignore_case = TRUE)) &
                        Tissue == "PFC (grey)"),
    grp_i1 = (Source == "iNeuron-I"),
    grp_i2 = (Source == "iNeuron-II")
  )

neu_expand <- neu_flags %>%
  tidyr::pivot_longer(starts_with("grp_"), names_to = "RegionFlag", values_to = "in_group") %>%
  filter(in_group) %>%
  mutate(
    Region = dplyr::recode(RegionFlag,
                           grp_ctrl_grey  = "PFC (grey)",
                           grp_ctrl_white = "PFC (white)",
                           grp_ftld_grey  = "FTLD (PFC grey)",
                           grp_apoe2_grey = "APOE2 (PFC grey)",
                           grp_apoe4_grey = "APOE4 (PFC grey)",
                           grp_i1         = "iNeuron-I",
                           grp_i2         = "iNeuron-II"
    )
  )

# 6) Long plotting DF for %DHA
long_df <- neu_expand %>%
  transmute(
    Region,
    PC = PC_pctDHA, PE = PE_pctDHA, PI = PI_pctDHA, PS = PS_pctDHA
  ) %>%
  tidyr::pivot_longer(c(PC,PE,PI,PS), names_to = "LIPID_CLASS", values_to = "Value") %>%
  mutate(
    LIPID_CLASS = factor(LIPID_CLASS, levels = class_levels),
    METRIC = "% DHA"
  ) %>%
  filter(!is.na(Value))

# 7) Your plotting call (unchanged)
gpl_colors   <- c(PC="#88CCEE", PE="#44AA99", PI="#DDCC77", PS="#AA4499")
regions_pref <- c("PFC (grey)","APOE2 (PFC grey)","APOE4 (PFC grey)",
                  "FTLD (PFC grey)","PFC (white)","iNeuron-I","iNeuron-II")
regions_order <- intersect(regions_pref, unique(long_df$Region))
x_labels <- setNames(regions_order, regions_order)


plots <- plot_gpl_metrics_points_many_by_region(
  long_df,
  metrics       = c("% DHA"),
  custom_colors = gpl_colors,
  regions_order = regions_order,
  x_axis_labels = x_labels,
  y_mode        = "smart",
  save_dir      = "Graphs/2025/PUFA_Paper/Commentary/",
  width = 6.1, height = 6, dpi = 300
)


#### Svennerholm human grey and white cortex####
hum <- read.csv("CSV_lipidomics/Svennerholm.csv", check.names = FALSE)
names(hum) <- trimws(names(hum))

meta_cols <- c("AGE","Region","Age (Days from conception)")
classes   <- c("PC","PE","PS","PI")
dha_cols  <- paste0(classes, ".22.6(n3)")

# 0) sanity
need <- c(meta_cols, dha_cols)
miss <- setdiff(need, names(hum))
if (length(miss)) stop("Missing expected columns: ", paste(miss, collapse=", "))

# 1) restrict to PFC grey / white
hum_grey  <- hum %>% filter(Region == "Cerebral_Cortex_Grey")
hum_white <- hum %>% filter(Region == "Cerebral_White")

# 2) build %DHA tables (already 0–100 in Svennerholm sheet)
make_dha <- function(df) {
  out <- df %>% select(all_of(meta_cols), any_of(dha_cols))
  # coerce to numeric safely
  for (c in dha_cols) out[[c]] <- suppressWarnings(as.numeric(out[[c]]))
  out
}

comp_grey  <- make_dha(hum_grey)
comp_white <- make_dha(hum_white)

# 3) pivot to long (LIPID_CLASS, METRIC="% DHA", Value)
to_long_dha <- function(comp_df) {
  comp_df %>%
    pivot_longer(cols = all_of(dha_cols),
                 names_to = "LIPID_CLASS", values_to = "Value") %>%
    mutate(
      LIPID_CLASS = sub("\\.22\\.6\\(n3\\)$", "", LIPID_CLASS),
      LIPID_CLASS = factor(LIPID_CLASS, levels = classes, ordered = TRUE),
      METRIC = "% DHA"
    ) %>%
    filter(!is.na(Value))
}

long_grey  <- to_long_dha(comp_grey)  %>% select(all_of(meta_cols), LIPID_CLASS, METRIC, Value)
long_white <- to_long_dha(comp_white) %>% select(all_of(meta_cols), LIPID_CLASS, METRIC, Value)

# 4) small plotting helper (mean±SEM across AGE order)
if (!exists("plot_gpl_metrics_many")) {
  plot_gpl_metrics_many <- function(long_df, metrics, custom_colors, myages,
                                    y_mode="smart", save_dir=NULL, width=8, height=6, dpi=300) {
    plots <- list()
    for (m in metrics) {
      df <- long_df %>% filter(METRIC == m)
      df$AGE <- factor(df$AGE, levels = myages, ordered = TRUE)
      
      sumdf <- df %>%
        group_by(LIPID_CLASS, AGE) %>%
        summarise(mean = mean(Value, na.rm = TRUE),
                  sem  = sd(Value,  na.rm = TRUE)/sqrt(n()),
                  .groups = "drop") %>%
        filter(!is.na(mean))
      
      p <- ggplot(sumdf, aes(AGE, mean, color = LIPID_CLASS, group = LIPID_CLASS)) +
        geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.18, linewidth = 0.6) +
        geom_point(size = 3) +
        scale_color_manual(values = custom_colors, limits = levels(df$LIPID_CLASS), drop = FALSE) +
        labs(x = "", y = "%", title = m) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 26, face="plain"),
          axis.text.y = element_text(size = 26),
          plot.title  = element_text(hjust = 0.5, size = 30),
          axis.title.y = element_text(size = 30)
        )
      if (!is.null(save_dir) && !dir.exists(save_dir)) dir.create(save_dir, TRUE)
      if (!is.null(save_dir)) ggsave(file.path(save_dir, paste0(m, ".png")), p, width=width, height=height, dpi=dpi)
      plots[[m]] <- p
    }
    plots
  }
}

# 5) colors + age order per region
gpl_colors   <- c(PC="#88CCEE", PE="#44AA99", PI="#DDCC77", PS="#AA4499")

myages_grey  <- hum_grey  %>% arrange(`Age (Days from conception)`) %>% pull(AGE) %>% unique()
myages_white <- hum_white %>% arrange(`Age (Days from conception)`) %>% pull(AGE) %>% unique()

# 6) plot and save
theme_set(theme_get() + theme(axis.text.x = element_text(face = "plain")))
theme_set(theme_get() + theme(axis.text.y = element_text(face = "plain")))
theme_set(theme_get() + theme(axis.title.y = element_text(face = "plain")))
theme_set(theme_get() +theme(plot.title = element_text(face = "plain")))

plots_grey_dha <- plot_gpl_metrics_many(
  long_grey,
  metrics = c("% DHA"),
  custom_colors = gpl_colors,
  myages = myages_grey,
  y_mode = "smart",
  save_dir = "Graphs/2025/PUFA_Paper/Commentary/grey/",
  width = 6, height = 6, dpi = 300
)

plots_white_dha <- plot_gpl_metrics_many(
  long_white,
  metrics = c("% DHA"),
  custom_colors = gpl_colors,
  myages = myages_white,
  y_mode = "smart",
  save_dir = "Graphs/2025/PUFA_Paper/Commentary/white/",
  width = 6, height = 6, dpi = 300
)
