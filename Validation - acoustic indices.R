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
validation<- read.csv("D:/BARLT Localization Project/localization_05312025/hawkears_lowthresh/HawkEars_validation_results.csv")
labels <- read.csv("D:/BARLT Localization Project/localization_05312025/hawkears_lowthresh/HawkEars_labels.csv")



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
    species = class_code,
    score   = as.numeric(score),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf)

val_models <- val %>%
  group_by(species) %>%
  nest() %>%
  mutate(
    model = map(data, ~ glm(tp ~ score, family = binomial(), data = .x)),
    n_val = map_int(data, nrow)
  ) %>%
  select(species, model, n_val)
fits <- val_models



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
  semi_join(fits %>% select(species), by = "species")  # only species  validated

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
      NTP_i = TPR_i * N_i,                          
      NFP_i = (1 - TPR_i) * N_i
    )
  
  curve <- map_dfr(threshold_grid, function(T) {
    kept <- df_bins %>% filter(bin_lower >= T)
    denom <- sum(kept$NTP_i + kept$NFP_i)
    prec  <- if (denom == 0) NA_real_ else sum(kept$NTP_i) / denom  
    retained <- if (sum(df_bins$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(df_bins$N_i)
    tibble(threshold = T, precision = prec, prop_retained = retained, n_retained = sum(kept$N_i))
  })
  
  hit <- curve %>% filter(!is.na(precision), precision >= target_precision) %>% arrange(threshold) %>% slice(1)
  
  if (nrow(hit) == 0) {
    # Tseng fallback when 0.9 not achievable 
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

curves <- fits %>%
  select(species, model) %>%
  left_join(lab_bins, by="species") %>%
  group_by(species) %>%
  summarise(
    model = list(first(model)),
    bins = list(cur_data() %>% select(bin_lower, bin_mid, N_i)),
    .groups="drop"
  ) %>%
  mutate(curve = map2(bins, model, ~ compute_tseng_threshold(.x, .y)$curve)) %>%
  select(species, curve) %>%
  unnest(curve)

plot(curves)




# -------------------------
#  plot the precision curve for one species (Tseng Fig.5-style) - UGLY RIGHT NOW!!!! 
# -------------------------
plot_species <- function(species_code, target = 0.90) {
  sp_fit  <- fits %>% filter(species == species_code)
  sp_bins <- lab_bins %>% filter(species == species_code)
  if (nrow(sp_fit) == 0 || nrow(sp_bins) == 0) stop("No fit/bins for that species.")
  
  out <- compute_tseng_threshold(sp_bins %>% select(bin_lower, bin_mid, N_i), sp_fit$model[[1]])
  
  d <- out$curve %>%
    arrange(threshold)
  

  plot(d$threshold, d$precision, type = "l",
       xlab = "Confidence threshold (T)",
       ylab = "Precision of retained detections",
       ylim = c(0, 1),
       main = paste0(species_code, " — Tseng-style threshold (precision + retained)"))
  
  abline(h = target, lty = 2)
  abline(v = out$summary$threshold, lty = 3)
  
  # overlay proportion retained
  par(new = TRUE)
  plot(d$threshold, d$prop_retained, type = "l",
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(0, 1))
  
  axis(side = 4, at = seq(0, 1, by = 0.2))
  mtext("Proportion of detections retained", side = 4, line = 3)
  
  legend("topright",
         legend = c("Precision", "Proportion retained", "Target precision", "Chosen threshold"),
         lty = c(1, 1, 2, 3),
         bty = "n")
}

# Example:
plot_species("VEER", target=0.9)
plot_species("COYE", target=0.9)
plot_species("MAWA", target=0.9)



## ggplot version ## 
plot_species_gg_onepanel <- function(species_code, target = 0.90) {
  sp_fit  <- fits %>% dplyr::filter(species == species_code)
  sp_bins <- lab_bins %>% dplyr::filter(species == species_code)
  if (nrow(sp_fit) == 0 || nrow(sp_bins) == 0) stop("No fit/bins for that species.")

  out <- compute_tseng_threshold(
    sp_bins %>% dplyr::select(bin_lower, bin_mid, N_i),
    sp_fit$model[[1]]
  )

  d <- out$curve %>% dplyr::arrange(threshold)
  thr_star <- out$summary$threshold

  # Build long data for clean ggplot mapping
  d_long <- d %>%
    dplyr::select(threshold, precision, prop_retained) %>%
    tidyr::pivot_longer(
      cols = c(precision, prop_retained),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(metric,
        precision = "Precision",
        prop_retained = "Proportion retained"
      )
    )

  ggplot2::ggplot(d_long, ggplot2::aes(x = threshold, y = value, linetype = metric)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 1.6, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = target, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = thr_star, linetype = "dotdash", linewidth = 0.8) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      name = "Precision of retained detections",
      sec.axis = ggplot2::sec_axis(~ ., name = "Proportion of detections retained")
    ) +
    ggplot2::scale_x_continuous(name = "Confidence threshold (T)") +
    ggplot2::scale_linetype_manual(values = c("Precision" = "solid", "Proportion retained" = "solid")) +
    ggplot2::labs(
      title = paste0(species_code, " — Tseng-style threshold (precision + retained)"),
      subtitle = paste0("Chosen threshold T* = ", format(thr_star, digits = 3),
                        " | Target precision = ", target),
      linetype = NULL
    ) +
    ggplot2::annotate(
      "text",
      x = thr_star,
      y = 0.98,
      label = paste0("T* = ", format(thr_star, digits = 3)),
      hjust = -0.05,
      vjust = 1,
      size = 3.6
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 10)),
      legend.position = "top"
    )
}

# Examples:
plot_species_gg_onepanel("VEER", target = 0.9)
plot_species_gg_onepanel("COYE", target = 0.9)
plot_species_gg_onepanel("MAWA", target = 0.9)














## table to show the results quantitatively 

compute_tseng_threshold <- function(df_bins, model) {
  df_bins <- df_bins %>%
    mutate(
      TPR_i = predict(model, newdata = data.frame(score = bin_mid), type = "response"),
      NTP_i = TPR_i * N_i,
      NFP_i = (1 - TPR_i) * N_i
    )

  curve <- purrr::map_dfr(threshold_grid, function(T) {
    kept <- df_bins %>% dplyr::filter(bin_lower >= T)
    denom <- sum(kept$NTP_i + kept$NFP_i)
    prec  <- if (denom == 0) NA_real_ else sum(kept$NTP_i) / denom
    retained <- if (sum(df_bins$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(df_bins$N_i)

    tibble::tibble(
      threshold = T,
      precision = prec,
      prop_retained = retained,
      n_retained = sum(kept$N_i)
    )
  })

  hit <- curve %>%
    dplyr::filter(!is.na(precision), precision >= target_precision) %>%
    dplyr::arrange(threshold) %>%
    dplyr::slice(1)

  status <- "hit_target"
  if (nrow(hit) == 0) {
    hit <- curve %>% dplyr::filter(threshold == fallback_threshold)
    status <- "fallback"
  }

  hit <- hit %>% dplyr::mutate(status = status)

  list(summary = hit, curve = curve, bins = df_bins)
}

threshold_table <- fits %>%
  dplyr::select(species, model, n_val) %>%
  dplyr::left_join(lab_bins %>% dplyr::group_by(species) %>%
                     dplyr::summarise(total_detections = sum(N_i), .groups = "drop"),
                   by = "species") %>%
  dplyr::left_join(lab_bins, by = "species") %>%
  dplyr::group_by(species, model, n_val, total_detections) %>%
  dplyr::summarise(
    out = list(
      compute_tseng_threshold(
        dplyr::cur_data_all() %>% dplyr::select(bin_lower, bin_mid, N_i),
        model[[1]]
      )
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    chosen_threshold = purrr::map_dbl(out, ~ .x$summary$threshold),
    precision_at_threshold = purrr::map_dbl(out, ~ .x$summary$precision),
    prop_retained = purrr::map_dbl(out, ~ .x$summary$prop_retained),
    n_retained = purrr::map_int(out, ~ .x$summary$n_retained),
    status = purrr::map_chr(out, ~ .x$summary$status),
    n_dropped = total_detections - n_retained
  ) %>%
  dplyr::select(
    species,
    n_val,
    total_detections,
    chosen_threshold,
    precision_at_threshold,
    prop_retained,
    n_retained,
    n_dropped,
    status
  ) %>%
  dplyr::arrange(dplyr::desc(chosen_threshold), species)

threshold_table



threshold_table_pretty <- threshold_table %>%
  dplyr::mutate(
    chosen_threshold = round(chosen_threshold, 2),
    precision_at_threshold = round(precision_at_threshold, 3),
    prop_retained = round(prop_retained, 3)
  )

readr::write_csv(threshold_table_pretty, "HawkEars_thresholds_precision_0.90.csv")





# ============================================================
# Bioacoustic Index (BI) +  BI-aware thresholding
# ============================================================

# ---- new packages ----
#install.packages(c("tuneR","soundecology","readr","ggplot2"))  # if needed
library(tuneR)
library(soundecology)
library(readr)
library(ggplot2)

# -------------------------
# BI settings (simple + "canned")
# -------------------------
audio_root <- "D:/BARLT Localization Project/localization_05312025/localizationtrim_new"  # CHANGE if needed
bi_win_sec <- 6                         # window length around clip midpoint
bi_min_freq <- 2000                     # typical bird band for BI
bi_max_freq <- 8000
bi_fft_w <- 512

# -------------------------
# Helper: safe z-score
# -------------------------
z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# -------------------------
# Helper: read window centered on clip, shift at file edges
# -------------------------
safe_read_window <- function(path, mid_sec, win_sec) {
  start1 <- mid_sec - win_sec/2
  end1   <- mid_sec + win_sec/2

  try_read <- function(from_s, to_s) {
    tryCatch({
      w <- tuneR::readWave(path, from = max(0, from_s), to = max(0, to_s), units = "seconds")
      if (w@stereo) w <- tuneR::mono(w, which = "left")
      if (length(w@left) < 64) stop("too_short")
      w
    }, error = function(e) NULL)
  }

  # 1) centered
  w <- try_read(start1, end1)
  if (!is.null(w)) return(list(wave = w, win_type = "centered"))

  # 2) shift left (window ends at mid)
  w <- try_read(mid_sec - win_sec, mid_sec)
  if (!is.null(w)) return(list(wave = w, win_type = "left_shift"))

  # 3) shift right (window starts at mid)
  w <- try_read(mid_sec, mid_sec + win_sec)
  if (!is.null(w)) return(list(wave = w, win_type = "right_shift"))

  # 4) fallback (file start)
  w <- try_read(0, win_sec)
  if (!is.null(w)) return(list(wave = w, win_type = "file_start"))

  list(wave = NULL, win_type = "failed")
}

# -------------------------
# Helper: compute BI from soundecology
# -------------------------
calc_bi <- function(w) {
  bi <- tryCatch(
    soundecology::bioacoustic_index(w, min_freq = bi_min_freq, max_freq = bi_max_freq, fft_w = bi_fft_w),
    error = function(e) NULL
  )
  if (is.null(bi)) return(NA_real_)
  as.numeric(bi$left_area)
}

# -------------------------
# 1) Re-read validation with timing + join to wav paths
#    (assumes validation has filename + start_time + end_time)
# -------------------------
wav_index <- tibble::tibble(
  wav_path = list.files(audio_root, pattern = "\\.wav$", recursive = TRUE, full.names = TRUE)
) %>%
  dplyr::mutate(filename = basename(wav_path)) %>%
  dplyr::distinct(filename, .keep_all = TRUE)

val_bi <- as_tibble(validation) %>%
  mutate(
    species = class_code,
    score   = as.numeric(score),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    ),
    start_time = as.numeric(start_time),
    end_time   = as.numeric(end_time),
    clip_mid   = (start_time + end_time) / 2
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species), !is.na(clip_mid)) %>%
  filter(score >= min_conf, score <= max_conf) %>%
  left_join(wav_index, by = "filename")

if (any(is.na(val_bi$wav_path))) {
  print(val_bi %>% filter(is.na(wav_path)) %>% distinct(filename) %>% head(50))
  stop("Some validation filenames were not found under audio_root.")
}

# -------------------------
# 2) Compute BI per validation row
# -------------------------
val_bi <- val_bi %>%
  mutate(win = purrr::map2(wav_path, clip_mid, ~ safe_read_window(.x, .y, bi_win_sec))) %>%
  mutate(
    win_type = purrr::map_chr(win, "win_type"),
    wave_obj = purrr::map(win, "wave"),
    BI_raw   = purrr::map_dbl(wave_obj, ~ if (is.null(.x)) NA_real_ else calc_bi(.x))
  ) %>%
  select(-win, -wave_obj) %>%
  filter(!is.na(BI_raw))

# Optional QC: how often do we have to shift the window?
val_bi %>%
  count(win_type) %>%
  arrange(desc(n)) %>%
  print()

# -------------------------
# 3) Fit species-specific BI models: tp ~ score + BI_z (+ interaction)
# -------------------------
fits_bi <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  nest() %>%
  mutate(
    model_bi = map(data, ~ glm(tp ~ score + BI_z + score:BI_z, family = binomial(), data = .x)),
    n_val_bi = map_int(data, nrow)
  ) %>%
  select(species, model_bi, n_val_bi)

# -------------------------
# 4) Make BI-conditional threshold table
#    Rule: pick thresholds at BI_z = median (and optionally q25/q75)
# -------------------------
bi_levels <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  summarise(
    BI_z_q25 = quantile(BI_z, 0.25, na.rm = TRUE),
    BI_z_med = median(BI_z, na.rm = TRUE),
    BI_z_q75 = quantile(BI_z, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(BI_z_q25, BI_z_med, BI_z_q75),
                      names_to = "BI_level",
                      values_to = "BI_z") %>%
  mutate(BI_level = recode(BI_level,
                           BI_z_q25 = "Low BI (25%)",
                           BI_z_med = "Median BI",
                           BI_z_q75 = "High BI (75%)"))

compute_tseng_threshold_bi <- function(df_bins, model, BI_z_fixed) {
  # df_bins has: bin_lower, bin_mid, N_i
  df_bins <- df_bins %>%
    mutate(
      TPR_i = predict(model,
                      newdata = data.frame(score = bin_mid, BI_z = BI_z_fixed),
                      type = "response"),
      NTP_i = TPR_i * N_i,
      NFP_i = (1 - TPR_i) * N_i
    )

  curve <- purrr::map_dfr(threshold_grid, function(T) {
    kept <- df_bins %>% dplyr::filter(bin_lower >= T)
    denom <- sum(kept$NTP_i + kept$NFP_i)
    prec  <- if (denom == 0) NA_real_ else sum(kept$NTP_i) / denom
    retained <- if (sum(df_bins$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(df_bins$N_i)

    tibble::tibble(
      threshold = T,
      precision = prec,
      prop_retained = retained,
      n_retained = sum(kept$N_i)
    )
  })

  hit <- curve %>%
    dplyr::filter(!is.na(precision), precision >= target_precision) %>%
    dplyr::arrange(threshold) %>%
    dplyr::slice(1)

  status <- "hit_target"
  if (nrow(hit) == 0) {
    hit <- curve %>% dplyr::filter(threshold == fallback_threshold)
    status <- "fallback"
  }

  hit <- hit %>% dplyr::mutate(status = status)

  list(summary = hit, curve = curve, bins = df_bins)
}

threshold_table_bi <- fits_bi %>%
  left_join(bi_levels, by = "species") %>%
  left_join(lab_bins %>%
              group_by(species) %>%
              summarise(total_detections = sum(N_i), .groups = "drop"),
            by = "species") %>%
  left_join(lab_bins, by = "species") %>%
  group_by(species, BI_level, BI_z, model_bi, n_val_bi, total_detections) %>%
  summarise(
    out = list(
      compute_tseng_threshold_bi(
        dplyr::cur_data_all() %>% dplyr::select(bin_lower, bin_mid, N_i),
        model_bi[[1]],
        BI_z_fixed = BI_z
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    chosen_threshold = purrr::map_dbl(out, ~ .x$summary$threshold),
    precision_at_threshold = purrr::map_dbl(out, ~ .x$summary$precision),
    prop_retained = purrr::map_dbl(out, ~ .x$summary$prop_retained),
    n_retained = purrr::map_int(out, ~ .x$summary$n_retained),
    status = purrr::map_chr(out, ~ .x$summary$status),
    n_dropped = total_detections - n_retained
  ) %>%
  select(
    species,
    BI_level,
    n_val_bi,
    total_detections,
    chosen_threshold,
    precision_at_threshold,
    prop_retained,
    n_retained,
    n_dropped,
    status
  ) %>%
  arrange(species, BI_level)

threshold_table_bi

# Save a "pretty" version
threshold_table_bi_pretty <- threshold_table_bi %>%
  mutate(
    chosen_threshold = round(chosen_threshold, 2),
    precision_at_threshold = round(precision_at_threshold, 3),
    prop_retained = round(prop_retained, 3)
  )
threshold_table_bi_pretty

readr::write_csv(threshold_table_bi_pretty, "HawkEars_thresholds_precision_0.90_BIconditional.csv")

# -------------------------
# 5) Optional quick plot: one species, BI-conditional curves (precision only)
# -------------------------
curves_bi <- fits_bi %>%
  left_join(bi_levels, by = "species") %>%
  left_join(lab_bins, by = "species") %>%
  group_by(species, BI_level, BI_z) %>%
  summarise(
    model_bi = list(fits_bi$model_bi[match(first(species), fits_bi$species)][[1]]),
    bins = list(cur_data() %>% select(bin_lower, bin_mid, N_i)),
    .groups = "drop"
  ) %>%
  mutate(curve = pmap(list(bins, model_bi, BI_z),
                      ~ compute_tseng_threshold_bi(..1, ..2, ..3)$curve)) %>%
  select(species, BI_level, curve) %>%
  unnest(curve)

# Example plot for one species (change code)
one_sp <- "MAWA"
ggplot(curves_bi %>% filter(species == one_sp),
       aes(x = threshold, y = precision, linetype = BI_level)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = target_precision, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = "Confidence threshold (T)",
    y = "Precision of retained detections",
    linetype = "BI condition",
    title = paste0(one_sp, " — Tseng-weighted precision vs threshold (BI-conditional)")
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")



# ============================================================
# END-OF-SCRIPT: Model diagnostics + BI vs non-BI comparison tables
# ============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(broom)



# ---- helper: safe z-score (you already have this; keep ONE copy in your script) ----
z <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---- helper: metrics for binomial GLM ----
calc_glm_metrics <- function(mod, dat) {
  dat <- dat %>% dplyr::filter(!is.na(tp))
  p <- stats::predict(mod, type = "response")
  y <- dat$tp

  # AIC
  aic <- stats::AIC(mod)

  # McFadden pseudo-R2
  null_mod <- stats::glm(tp ~ 1, family = binomial(), data = dat)
  pseudoR2 <- 1 - as.numeric(stats::logLik(mod) / stats::logLik(null_mod))

  # AUC (rank-based; no extra packages)
  pos <- p[y == 1]
  neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) {
    auc <- NA_real_
  } else {
    r <- rank(c(pos, neg), ties.method = "average")
    r_pos_sum <- sum(r[seq_along(pos)])
    n_pos <- length(pos); n_neg <- length(neg)
    auc <- (r_pos_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
  }

  # Brier score
  brier <- mean((p - y)^2, na.rm = TRUE)

  tibble::tibble(AIC = aic, pseudoR2 = pseudoR2, AUC = auc, Brier = brier)
}

# ------------------------------------------------------------
# OPTIONAL BUT STRONGLY RECOMMENDED:
# Compare on the same rows (only rows with BI computed).
# This avoids "BI model looks better" just because it's evaluated on an easier subset.
# ------------------------------------------------------------
val_bi_eval <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  ungroup()

# define the matched set of validation rows (same filenames/times as BI-evaluable rows)
# if you have an ID column, use it; otherwise fall back to these columns:
join_keys <- intersect(c("filename","start_time","end_time","species","score","tp"), names(val_bi_eval))

val_nonbi_matched <- val %>%
  # bring in matching keys from the BI-evaluable dataset
  inner_join(val_bi_eval %>% dplyr::select(all_of(join_keys)) %>% distinct(), by = join_keys)

# If the join keys end up empty (because columns differ), just evaluate non-BI on full val:
if (length(join_keys) == 0 || nrow(val_nonbi_matched) == 0) {
  val_nonbi_matched <- val
}

# -------------------------
# 1) NON-BI model results: tp ~ score
# -------------------------
nonbi_model_results <- val_models %>%         # species, model, n_val
  mutate(
    metrics = purrr::map2(model, species, ~ calc_glm_metrics(.x, val_nonbi_matched %>% filter(species == .y)))
  ) %>%
  tidyr::unnest(metrics) %>%
  select(species, n_val, AIC, pseudoR2, AUC, Brier) %>%
  arrange(species)

# -------------------------
# 2) BI model results: tp ~ score + BI_z + score:BI_z
# -------------------------
bi_model_results <- fits_bi %>%               # species, model_bi, n_val_bi
  mutate(
    metrics = purrr::map2(model_bi, species, ~ calc_glm_metrics(.x, val_bi_eval %>% filter(species == .y)))
  ) %>%
  tidyr::unnest(metrics) %>%
  select(species, n_val_bi, AIC, pseudoR2, AUC, Brier) %>%
  arrange(species)

# -------------------------
# 3) Comparison table
# -------------------------
model_comparison <- nonbi_model_results %>%
  transmute(
    species,
    AIC_nonBI = AIC, pseudoR2_nonBI = pseudoR2, AUC_nonBI = AUC, Brier_nonBI = Brier
  ) %>%
  left_join(
    bi_model_results %>%
      transmute(
        species,
        AIC_BI = AIC, pseudoR2_BI = pseudoR2, AUC_BI = AUC, Brier_BI = Brier
      ),
    by = "species"
  ) %>%
  mutate(
    dAIC = AIC_BI - AIC_nonBI,                 # negative = BI better
    dPseudoR2 = pseudoR2_BI - pseudoR2_nonBI,  # positive = BI better
    dAUC = AUC_BI - AUC_nonBI,                 # positive = BI better
    dBrier = Brier_BI - Brier_nonBI            # negative = BI better
  ) %>%
  arrange(species)

# -------------------------
# 4) Print + save
# -------------------------
print(nonbi_model_results)
print(bi_model_results)
print(model_comparison)


readr::write_csv(nonbi_model_results, "HawkEars_model_metrics_nonBI.csv")
readr::write_csv(bi_model_results,   "HawkEars_model_metrics_BI.csv")
readr::write_csv(model_comparison,   "HawkEars_model_comparison_BI_vs_nonBI.csv")


#
#dAIC = AIC_BI − AIC_nonBI → negative means BI model is better
#dAUC = AUC_BI − AUC_nonBI → positive means BI model is better
#dBrier = Brier_BI − Brier_nonBI → negative means BI model is better
#dPseudoR2 → positive means BI model is better

final_summary_table <- threshold_table_pretty %>%
  rename(
    chosen_threshold_nonBI = chosen_threshold,
    precision_at_threshold_nonBI = precision_at_threshold,
    prop_retained_nonBI = prop_retained,
    n_retained_nonBI = n_retained,
    n_dropped_nonBI = n_dropped,
    status_nonBI = status
  ) %>%
  left_join(
    threshold_table_bi_pretty %>%
      filter(BI_level == "Median BI") %>%   # pick one BI condition for the “headline” table
      select(species, BI_level,
             chosen_threshold, precision_at_threshold, prop_retained,
             n_retained, n_dropped, status) %>%
      rename(
        chosen_threshold_BI = chosen_threshold,
        precision_at_threshold_BI = precision_at_threshold,
        prop_retained_BI = prop_retained,
        n_retained_BI = n_retained,
        n_dropped_BI = n_dropped,
        status_BI = status
      ),
    by = "species"
  ) %>%
  left_join(model_comparison, by = "species") %>%
  arrange(species)

print(final_summary_table)
readr::write_csv(final_summary_table, "HawkEars_FINAL_thresholds_plus_model_comparison.csv")
