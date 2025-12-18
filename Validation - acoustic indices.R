# ============================================================
# HawkEars: species-specific threshold calibration 
## with inspiration from: Tseng, S., Hodder, D.P. & Otter, K.A. 
#Setting BirdNET confidence thresholds: 
#species-specific vs. universal approaches. J Ornithol 166, 1123–1135 (2025). 
#https://doi.org/10.1007/s10336-025-02260-w
# ============================================================

library(dplyr)
library(purrr)
library(tidyr)

# ---- CSVs ----
validation<- read.csv("E:/BAR-LT_LocalizationProject/localization_05312025/hawkears_lowthresh/HawkEars_validation_results.csv")
labels <- read.csv("E:/BAR-LT_LocalizationProject/localization_05312025/hawkears_lowthresh/HawkEars_labels.csv")



# -------------------------
# Settings (match Tseng)
# -------------------------
min_conf <- 0.10
max_conf <- 1.00
step     <- 0.05
target_precision <- 0.90
fallback_threshold <- 0.95  # Tseng fallback when 0.9 can't be achieved 
threshold_grid <- seq(min_conf, fallback_threshold, by = step)

# -------------------------
# 1) Clean validation data: TP (1) vs FP (0)
# -------------------------
val <- as_tibble(validation) %>%
  mutate(
    species = class_code,  # use class_code to avoid name mismatches; swap to class_name if you prefer
    score   = as.numeric(score),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf)

# Fit logistic regression per species (Tseng: GLM binomial, score predictor) :contentReference[oaicite:6]{index=6}
fits <- val %>%
  group_by(species) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(tp ~ score, family = binomial(), data = .x)),
    n_val = map_int(data, nrow)
  ) %>%
  select(species, model, n_val)

# -------------------------
# 2) Bin full labels into 0.05 confidence classes and count N_i per species
# -------------------------
lab <- as_tibble(labels) %>%
  mutate(
    species = class_code,       # must match val species choice
    score   = as.numeric(score)
  ) %>%
  filter(!is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf) %>%
  semi_join(fits %>% select(species), by = "species")  # only species you validated

# Fast numeric binning to the LOWER edge of the 0.05 bin:
lab_bins <- lab %>%
  mutate(
    bin_lower = min_conf + floor((score - min_conf) / step) * step,
    bin_lower = pmin(fallback_threshold, pmax(min_conf, bin_lower)),
    bin_mid   = bin_lower + step/2
  ) %>%
  count(species, bin_lower, bin_mid, name = "N_i")

# -------------------------
# 3) For each species/bin, predict TPR_i from GLM and compute NTP_i / NFP_i
#    Then compute precision(T) over thresholds and pick minimum T with precision>=0.9
# -------------------------
compute_tseng_threshold <- function(df_bins, model) {
  # df_bins has: bin_lower, bin_mid, N_i
  df_bins <- df_bins %>%
    mutate(
      TPR_i = predict(model, newdata = data.frame(score = bin_mid), type = "response"),
      NTP_i = TPR_i * N_i,                           # :contentReference[oaicite:7]{index=7}
      NFP_i = (1 - TPR_i) * N_i
    )
  
  curve <- map_dfr(threshold_grid, function(T) {
    kept <- df_bins %>% filter(bin_lower >= T)
    denom <- sum(kept$NTP_i + kept$NFP_i)
    prec  <- if (denom == 0) NA_real_ else sum(kept$NTP_i) / denom   # :contentReference[oaicite:8]{index=8}
    retained <- if (sum(df_bins$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(df_bins$N_i)
    tibble(threshold = T, precision = prec, prop_retained = retained, n_retained = sum(kept$N_i))
  })
  
  hit <- curve %>% filter(!is.na(precision), precision >= target_precision) %>% arrange(threshold) %>% slice(1)
  
  if (nrow(hit) == 0) {
    # Tseng fallback when 0.9 not achievable :contentReference[oaicite:9]{index=9}
    hit <- curve %>% filter(threshold == fallback_threshold)
  }
  
  list(summary = hit, curve = curve, bins = df_bins)
}

results <- fits %>%
  left_join(lab_bins, by = "species") %>%
  group_by(species, model, n_val) %>%
  summarise(
    out = list(compute_tseng_threshold(cur_data_all() %>% select(bin_lower, bin_mid, N_i), model[[1]])),
    total_detections = sum(N_i, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    threshold = map_dbl(out, ~ .x$summary$threshold),
    precision = map_dbl(out, ~ .x$summary$precision),
    prop_retained = map_dbl(out, ~ .x$summary$prop_retained)
  ) %>%
  select(species, n_val, total_detections, threshold, precision, prop_retained)

print(results)

# -------------------------
# Optional: plot the precision curve for one species (Tseng Fig.5-style)
# -------------------------
plot_species <- function(species_code) {
  sp_fit <- fits %>% filter(species == species_code)
  sp_bins <- lab_bins %>% filter(species == species_code)
  if (nrow(sp_fit) == 0 || nrow(sp_bins) == 0) stop("No fit/bins for that species.")
  
  out <- compute_tseng_threshold(sp_bins %>% select(bin_lower, bin_mid, N_i), sp_fit$model[[1]])
  
  plot(out$curve$threshold, out$curve$precision, type = "l",
       xlab = "Confidence threshold (T)", ylab = "Precision of retained detections")
  abline(h = target_precision, lty = 2)
  abline(v = out$summary$threshold, lty = 3)
  title(main = paste0(species_code, " — Tseng-style threshold"))
}

# Example:
plot_species("VEER")
