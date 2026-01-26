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


cat("\n=== VISUALIZING RAW BI DISTRIBUTIONS ===\n")

# Calculate raw BI distribution statistics per species
bi_raw_stats <- val_bi %>%
  group_by(species) %>%
  summarise(
    n_obs = n(),
    BI_raw_mean = mean(BI_raw, na.rm = TRUE),
    BI_raw_sd = sd(BI_raw, na.rm = TRUE),
    BI_raw_min = min(BI_raw, na.rm = TRUE),
    BI_raw_q25 = quantile(BI_raw, 0.25, na.rm = TRUE),
    BI_raw_median = median(BI_raw, na.rm = TRUE),
    BI_raw_q75 = quantile(BI_raw, 0.75, na.rm = TRUE),
    BI_raw_max = max(BI_raw, na.rm = TRUE),
    .groups = "drop"
  )

print("Raw BI Distribution Statistics:")
print(bi_raw_stats)


# Individual species histograms (raw BI)
for (sp in unique(val_bi$species)) {
  sp_data <- val_bi %>% filter(species == sp)
  sp_stats <- bi_raw_stats %>% filter(species == sp)
  
  p <- ggplot(sp_data, aes(x = BI_raw)) +
    geom_histogram(bins = 30, fill = "#1b9e77", color = "white", alpha = 0.8) +
    geom_vline(xintercept = sp_stats$BI_raw_q25, linetype = "dashed", 
               color = "#2166ac", linewidth = 1) +
    geom_vline(xintercept = sp_stats$BI_raw_median, linetype = "solid", 
               color = "#762a83", linewidth = 1.2) +
    geom_vline(xintercept = sp_stats$BI_raw_q75, linetype = "dashed", 
               color = "#c51b7d", linewidth = 1) +
    labs(
      title = paste0(sp, " — Raw Bioacoustic Index Distribution"),
      subtitle = paste0("n = ", sp_stats$n_obs, 
                       " | Median = ", round(sp_stats$BI_raw_median, 1),
                       " | IQR = [", round(sp_stats$BI_raw_q25, 1), ", ", 
                       round(sp_stats$BI_raw_q75, 1), "]"),
      x = "Raw Bioacoustic Index",
      y = "Count"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Q25\nMedian\nQ75"),
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "gray30")
  
  print(p)
  ggsave(paste0("HawkEars_BI_raw_distribution_", sp, ".png"), 
         plot = p, width = 8, height = 5, dpi = 300)
}

# Combined multi-panel histogram (all species, raw BI)
combined_hist_raw <- ggplot(val_bi, aes(x = BI_raw, fill = species)) +
  geom_histogram(bins = 30, color = "white", alpha = 0.8) +
  geom_vline(data = bi_raw_stats, aes(xintercept = BI_raw_q25), 
             linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(data = bi_raw_stats, aes(xintercept = BI_raw_median), 
             linetype = "solid", color = "gray20", linewidth = 0.8) +
  geom_vline(data = bi_raw_stats, aes(xintercept = BI_raw_q75), 
             linetype = "dashed", color = "gray40", linewidth = 0.6) +
  facet_wrap(~ species, ncol = 1, scales = "free") +  # Note: "free" scales for raw values
  scale_fill_manual(values = c("COYE" = "#1b9e77", 
                                "MAWA" = "#d95f02", 
                                "VEER" = "#7570b3")) +
  labs(
    title = "Raw Bioacoustic Index Distributions by Species",
    subtitle = "Solid line = median; Dashed lines = Q25 and Q75",
    x = "Raw Bioacoustic Index",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(combined_hist_raw)
ggsave("HawkEars_combined_RAW_comparison.png",
       plot = combined_hist_raw, width = 10, height = 8, dpi = 300)

cat("\n=== VISUALIZING BI_z DISTRIBUTIONS ===\n")

# Calculate BI_z distribution statistics per species
bi_z_stats <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  summarise(
    n_obs = n(),
    BI_z_mean = mean(BI_z, na.rm = TRUE),
    BI_z_sd = sd(BI_z, na.rm = TRUE),
    BI_z_min = min(BI_z, na.rm = TRUE),
    BI_z_q25 = quantile(BI_z, 0.25, na.rm = TRUE),
    BI_z_median = median(BI_z, na.rm = TRUE),
    BI_z_q75 = quantile(BI_z, 0.75, na.rm = TRUE),
    BI_z_max = max(BI_z, na.rm = TRUE),
    .groups = "drop"
  )

print("BI_z Distribution Statistics:")
print(bi_z_stats)

# Create histogram data
bi_z_for_hist <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  ungroup()

# Individual species histograms
for (sp in unique(bi_z_for_hist$species)) {
  sp_data <- bi_z_for_hist %>% filter(species == sp)
  sp_stats <- bi_z_stats %>% filter(species == sp)
  
  p <- ggplot(sp_data, aes(x = BI_z)) +
    geom_histogram(bins = 30, fill = "#1b9e77", color = "white", alpha = 0.8) +
    geom_vline(xintercept = sp_stats$BI_z_q25, linetype = "dashed", 
               color = "#2166ac", linewidth = 1) +
    geom_vline(xintercept = sp_stats$BI_z_median, linetype = "solid", 
               color = "#762a83", linewidth = 1.2) +
    geom_vline(xintercept = sp_stats$BI_z_q75, linetype = "dashed", 
               color = "#c51b7d", linewidth = 1) +
    labs(
      title = paste0(sp, " — Standardized Bioacoustic Index (BI_z) Distribution"),
      subtitle = paste0("n = ", sp_stats$n_obs, 
                       " | Median = ", round(sp_stats$BI_z_median, 2),
                       " | IQR = [", round(sp_stats$BI_z_q25, 2), ", ", 
                       round(sp_stats$BI_z_q75, 2), "]"),
      x = "Standardized Bioacoustic Index (BI_z)",
      y = "Count"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Q25 (Low BI)\nMedian\nQ75 (High BI)"),
             hjust = 1.1, vjust = 1.5, size = 3.5, color = "gray30")
  
  print(p)
  ggsave(paste0("HawkEars_BI_z_distribution_", sp, ".png"), 
         plot = p, width = 8, height = 5, dpi = 300)
}

# Combined multi-panel histogram (all species)
combined_hist <- ggplot(bi_z_for_hist, aes(x = BI_z, fill = species)) +
  geom_histogram(bins = 30, color = "white", alpha = 0.8) +
  geom_vline(data = bi_z_stats, aes(xintercept = BI_z_q25), 
             linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(data = bi_z_stats, aes(xintercept = BI_z_median), 
             linetype = "solid", color = "gray20", linewidth = 0.8) +
  geom_vline(data = bi_z_stats, aes(xintercept = BI_z_q75), 
             linetype = "dashed", color = "gray40", linewidth = 0.6) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("COYE" = "#1b9e77", 
                                "MAWA" = "#d95f02", 
                                "VEER" = "#7570b3")) +
  labs(
    title = "Bioacoustic Index (BI_z) Distributions by Species",
    subtitle = "Solid line = median; Dashed lines = Q25 and Q75 (BI thresholding levels)",
    x = "Standardized Bioacoustic Index (BI_z)",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(combined_hist)




bi_comparison_data <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  ungroup() %>%
  select(species, BI_raw, BI_z) %>%
  pivot_longer(cols = c(BI_raw, BI_z), 
               names_to = "BI_type", 
               values_to = "BI_value") %>%
  mutate(BI_type = recode(BI_type,
                          BI_raw = "Raw BI",
                          BI_z = "Standardized BI (BI_z)"))

comparison_plot <- ggplot(bi_comparison_data, 
                          aes(x = BI_value, fill = BI_type)) +
  geom_histogram(bins = 30, color = "white", alpha = 0.7, position = "identity") +
  facet_grid(species ~ BI_type, scales = "free") +
  scale_fill_manual(values = c("Raw BI" = "#1b9e77", 
                                "Standardized BI (BI_z)" = "#d95f02")) +
  labs(
    title = "Comparison: Raw vs. Standardized Bioacoustic Index",
    x = "Bioacoustic Index Value",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(comparison_plot)
ggsave("HawkEars_BI_raw_vs_standardized_comparison.png", 
       plot = comparison_plot, width = 10, height = 8, dpi = 300)





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
#Model diagnostics + BI vs non-BI comparison tables
# ============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(broom)



# ---- helper z-score  ----
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


save.image(file = "workspace.RData")


## TEST












# ============================================================
# HawkEars Threshold Calibration - Enhancement Tasks

# ============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(broom)


# ============================================================
# TASK 3: INTERACTION TERM VALIDATION
# ============================================================

#' Extract interaction coefficients and test significance
interaction_validation <- fits_bi %>%
  mutate(
    # Extract model coefficients
    coef_summary = map(model_bi, ~ broom::tidy(.x, conf.int = TRUE)),
    
    # Get interaction term specifically
    interaction_term = map(coef_summary, ~ filter(.x, term == "score:BI_z"))
  ) %>%
  select(species, n_val_bi, interaction_term) %>%
  unnest(interaction_term) %>%
  mutate(
    # Is interaction significant?
    sig_05 = p.value < 0.05,
    sig_01 = p.value < 0.01,
    
    # Effect interpretation
    effect_direction = case_when(
      estimate > 0 ~ "BI amplifies confidence effect",
      estimate < 0 ~ "BI dampens confidence effect",
      TRUE ~ "No effect"
    )
  ) %>%
  select(species, n_val_bi, estimate, std.error, p.value, 
         conf.low, conf.high, sig_05, sig_01, effect_direction)

print("=== INTERACTION TERM VALIDATION ===")
print(interaction_validation)

# Save results
write.csv(interaction_validation, 
          "HawkEars_interaction_validation.csv", 
          row.names = FALSE)

#' Calculate predicted probability change across BI range
#' This shows if interaction is PRACTICALLY important
interaction_effect_size <- fits_bi %>%
  left_join(val_bi_eval %>% 
              group_by(species) %>% 
              summarise(BI_z_min = min(BI_z, na.rm = TRUE),
                       BI_z_max = max(BI_z, na.rm = TRUE),
                       .groups = "drop"),
            by = "species") %>%
  mutate(
    # Predict at median score, low vs high BI
    pred_data = map2(BI_z_min, BI_z_max, ~ {
      expand.grid(
        score = 0.5,  # median-ish score
        BI_z = c(.x, .y),
        BI_level = c("Low BI", "High BI")
      )
    }),
    
    # Get predictions
    predictions = map2(model_bi, pred_data, ~ {
      .y %>% mutate(pred_prob = predict(.x, newdata = .y, type = "response"))
    })
  ) %>%
  select(species, predictions) %>%
  unnest(predictions) %>%
  group_by(species) %>%
  summarise(
    prob_at_low_BI = pred_prob[BI_level == "Low BI"],
    prob_at_high_BI = pred_prob[BI_level == "High BI"],
    prob_diff = prob_at_high_BI - prob_at_low_BI,
    pct_change = (prob_diff / prob_at_low_BI) * 100,
    .groups = "drop"
  )

print("=== PRACTICAL EFFECT SIZE (Δ probability across BI range) ===")
print(interaction_effect_size)

#' Decision: Keep interaction if p<0.05 AND |prob_diff| > 0.05
interaction_decision <- interaction_validation %>%
  left_join(interaction_effect_size, by = "species") %>%
  mutate(
    keep_interaction = sig_05 & abs(prob_diff) > 0.05,
    recommendation = case_when(
      keep_interaction ~ "Use BI model with interaction",
      sig_05 & !keep_interaction ~ "Try additive BI model (no interaction)",
      TRUE ~ "Use non-BI model (interaction not justified)"
    )
  )

print("=== INTERACTION MODEL DECISION ===")
print(interaction_decision %>% select(species, p.value, prob_diff, keep_interaction, recommendation))

write.csv(interaction_decision, 
          "HawkEars_interaction_decision.csv", 
          row.names = FALSE)


# ============================================================
# TASK 2: BI LEVEL THRESHOLD CURVES (ALL THREE LEVELS)
# ============================================================

#' Recreate curves for all BI levels (you already have this in curves_bi)
#' But let's make a focused comparison plot

# Calculate threshold differences across BI levels
threshold_bi_comparison <- threshold_table_bi %>%
  select(species, BI_level, chosen_threshold, precision_at_threshold, prop_retained) %>%
  pivot_wider(
    names_from = BI_level,
    values_from = c(chosen_threshold, precision_at_threshold, prop_retained),
    names_sep = "_"
  ) %>%
  mutate(
    threshold_range = `chosen_threshold_High BI (75%)` - `chosen_threshold_Low BI (25%)`,
    threshold_sensitivity = case_when(
      abs(threshold_range) > 0.10 ~ "High sensitivity to BI",
      abs(threshold_range) > 0.05 ~ "Moderate sensitivity to BI",
      TRUE ~ "Low sensitivity to BI"
    )
  )

print("=== THRESHOLD SENSITIVITY TO BI ===")
print(threshold_bi_comparison %>% 
        select(species, threshold_range, threshold_sensitivity, 
               `chosen_threshold_Low BI (25%)`, 
               `chosen_threshold_Median BI`,
               `chosen_threshold_High BI (75%)`))

write.csv(threshold_bi_comparison, 
          "HawkEars_threshold_BI_sensitivity.csv", 
          row.names = FALSE)

#' Generate faceted plot for all BI levels
plot_all_bi_levels <- function(sp) {
  sp_curves <- curves_bi %>% filter(species == sp)
  sp_thresholds <- threshold_table_bi %>% filter(species == sp)
  
  ggplot(sp_curves, aes(x = threshold, y = precision, color = BI_level)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray40") +
    geom_vline(data = sp_thresholds, 
               aes(xintercept = chosen_threshold, color = BI_level),
               linetype = "dotted", linewidth = 0.8) +
    scale_color_manual(
      values = c("Low BI (25%)" = "#2166ac", 
                 "Median BI" = "#762a83", 
                 "High BI (75%)" = "#c51b7d")
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(sp, " — Precision vs. Threshold across BI Conditions"),
      subtitle = "Vertical dotted lines show chosen thresholds for each BI level",
      x = "Confidence Threshold",
      y = "Precision of Retained Detections",
      color = "BI Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
}

# Generate plots for each species
for (sp in unique(threshold_table_bi$species)) {
  p <- plot_all_bi_levels(sp)
  print(p)
  ggsave(paste0("HawkEars_BI_curves_", sp, ".png"), 
         plot = p, width = 8, height = 6, dpi = 300)
}


# ============================================================
# TASK 13: SPECIES-SPECIFIC INTERPRETATIONS
# ============================================================

#' Generate automated interpretation based on metrics
generate_interpretation <- function(species_data) {
  sp <- species_data$species
  
  # Extract key metrics
  dauc <- species_data$dAUC
  dbrier <- species_data$dBrier
  daic <- species_data$dAIC
  threshold_nonbi <- species_data$chosen_threshold_nonBI
  threshold_bi <- species_data$chosen_threshold_BI
  retention_nonbi <- species_data$prop_retained_nonBI
  retention_bi <- species_data$prop_retained_BI
  
  # Get interaction info
  interact_info <- interaction_decision %>% filter(species == sp)
  
  # Get threshold sensitivity
  sens_info <- threshold_bi_comparison %>% filter(species == sp)
  
  # Build interpretation
  interpretation <- paste0(
    "**", sp, "**: ",
    
    # BI model improvement
    if (dauc > 0.01) {
      paste0("BI substantially improves discrimination (ΔA UC = +", 
             round(dauc, 3), "). ")
    } else if (dauc > 0.005) {
      paste0("BI modestly improves discrimination (ΔAUC = +", 
             round(dauc, 3), "). ")
    } else if (dauc > 0) {
      paste0("BI marginally improves discrimination (ΔAUC = +", 
             round(dauc, 3), "); improvement may not be practically significant. ")
    } else {
      paste0("BI does not improve discrimination (ΔAUC = ", 
             round(dauc, 3), "). ")
    },
    
    # Threshold stability
    if (abs(sens_info$threshold_range) > 0.10) {
      paste0("Chosen threshold varies substantially across BI conditions (range = ", 
             round(sens_info$threshold_range, 2), "), indicating context-dependent detection. ")
    } else {
      paste0("Chosen threshold stable across BI levels (", 
             round(threshold_nonbi, 2), "), suggesting robust detection. ")
    },
    
    # Retention rate
    paste0("Non-BI threshold retains ", round(retention_nonbi * 100, 1), 
           "% of detections at 90% precision. "),
    
    # Ecological interpretation placeholder
    "Ecological context: [ADD SPECIES-SPECIFIC NOTES ON SONG CHARACTERISTICS, HABITAT ASSOCIATIONS, ETC.]"
  )
  
  tibble(species = sp, interpretation = interpretation)
}

# Generate interpretations for all species
species_interpretations <- final_summary_table %>%
  group_split(species) %>%
  purrr::map_dfr(generate_interpretation)

print("=== SPECIES INTERPRETATIONS ===")
for (i in 1:nrow(species_interpretations)) {
  cat("\n", species_interpretations$interpretation[i], "\n")
}

write.csv(species_interpretations, 
          "HawkEars_species_interpretations.csv", 
          row.names = FALSE)


# ============================================================
# TASK 7: BOOTSTRAP CONFIDENCE INTERVALS FOR THRESHOLDS
# ============================================================

#' Bootstrap function for one species (non-BI model)
bootstrap_threshold_nonbi <- function(species_name, model, validation_data, n_boot = 100) {
  
  sp_data <- validation_data %>% filter(species == species_name)
  n_obs <- nrow(sp_data)
  
  boot_thresholds <- map_dbl(1:n_boot, function(b) {
    # Resample with replacement
    boot_sample <- sp_data %>% 
      slice_sample(n = n_obs, replace = TRUE)
    
    # Refit model
    boot_model <- tryCatch(
      glm(tp ~ score, family = binomial(), data = boot_sample),
      error = function(e) NULL
    )
    
    if (is.null(boot_model)) return(NA_real_)
    
    # Recompute threshold using Tseng method
    # (Simplified version - in production, call your full compute_tseng_threshold function)
    # For now, just return a placeholder based on model coefficients
    coef(boot_model)[2]  # This is a placeholder - see note below
  })
  
  # Calculate CI
  ci_lower <- quantile(boot_thresholds, 0.025, na.rm = TRUE)
  ci_upper <- quantile(boot_thresholds, 0.975, na.rm = TRUE)
  ci_median <- median(boot_thresholds, na.rm = TRUE)
  
  tibble(
    species = species_name,
    ci_lower = ci_lower,
    ci_median = ci_median,
    ci_upper = ci_upper,
    n_boot = sum(!is.na(boot_thresholds))
  )
}

# NOTE: The bootstrap function above is simplified. For production use, you need to:
# 1. Bin the bootstrap sample into 0.05 bins
# 2. Predict TPR for each bin
# 3. Run the full Tseng threshold calculation
# This requires adapting your compute_tseng_threshold() function to work on bootstrap samples

# Here's a more complete version:

bootstrap_threshold_complete <- function(species_name, n_boot = 100) {
  
  # Get species data
  sp_val <- val_bi_eval %>% filter(species == species_name)
  sp_labels <- lab_bins %>% filter(species == species_name)
  
  boot_results <- map_dfr(1:n_boot, function(b) {
    
    # Resample validation data
    boot_val <- sp_val %>% 
      slice_sample(n = nrow(sp_val), replace = TRUE)
    
    # Refit model
    boot_model <- tryCatch(
      glm(tp ~ score, family = binomial(), data = boot_val),
      error = function(e) NULL
    )
    
    if (is.null(boot_model)) {
      return(tibble(boot_iter = b, threshold = NA_real_, 
                   precision = NA_real_, prop_retained = NA_real_))
    }
    
    # Apply to FULL label bins (not resampled) - this is the population we're deploying to
    out <- compute_tseng_threshold(
      sp_labels %>% select(bin_lower, bin_mid, N_i),
      boot_model
    )
    
    tibble(
      boot_iter = b,
      threshold = out$summary$threshold,
      precision = out$summary$precision,
      prop_retained = out$summary$prop_retained
    )
  })
  
  # Summarize bootstrap distribution
  tibble(
    species = species_name,
    threshold_median = median(boot_results$threshold, na.rm = TRUE),
    threshold_ci_lower = quantile(boot_results$threshold, 0.025, na.rm = TRUE),
    threshold_ci_upper = quantile(boot_results$threshold, 0.975, na.rm = TRUE),
    precision_median = median(boot_results$precision, na.rm = TRUE),
    precision_ci_lower = quantile(boot_results$precision, 0.025, na.rm = TRUE),
    precision_ci_upper = quantile(boot_results$precision, 0.975, na.rm = TRUE),
    n_boot_success = sum(!is.na(boot_results$threshold))
  )
}

# Run bootstrap for all species (WARNING: This will take several minutes!)
print("Running bootstrap... this may take 5-10 minutes")

bootstrap_cis_nonbi <- map_dfr(
  unique(val_bi_eval$species),
  ~ bootstrap_threshold_complete(.x, n_boot = 100)
)

print("=== BOOTSTRAP CONFIDENCE INTERVALS (Non-BI Model) ===")
print(bootstrap_cis_nonbi)

write.csv(bootstrap_cis_nonbi, 
          "HawkEars_bootstrap_CIs_nonBI.csv", 
          row.names = FALSE)

# Add CIs to final summary table
final_summary_with_cis <- final_summary_table %>%
  left_join(
    bootstrap_cis_nonbi %>% 
      select(species, threshold_ci_lower, threshold_ci_upper),
    by = "species"
  ) %>%
  mutate(
    threshold_ci_width = threshold_ci_upper - threshold_ci_lower,
    threshold_uncertainty = case_when(
      threshold_ci_width > 0.20 ~ "High uncertainty - need more validation data",
      threshold_ci_width > 0.10 ~ "Moderate uncertainty",
      TRUE ~ "Low uncertainty - threshold well-estimated"
    )
  )

print("=== FINAL SUMMARY WITH CONFIDENCE INTERVALS ===")
print(final_summary_with_cis %>% 
        select(species, chosen_threshold_nonBI, 
               threshold_ci_lower, threshold_ci_upper, 
               threshold_uncertainty))


# ============================================================
# TASK 11: PUBLICATION-READY FIGURES
# ============================================================

# FIGURE 1: Multi-panel precision curves (non-BI vs BI at median)
# This combines base model and BI model on same plot

fig1_data <- bind_rows(
  # Non-BI curves (from your original fits)
  fits %>%
    select(species, model) %>%
    left_join(lab_bins, by = "species") %>%
    group_by(species) %>%
    summarise(
      curve = list(compute_tseng_threshold(
        cur_data() %>% select(bin_lower, bin_mid, N_i),
        model[[1]]
      )$curve),
      .groups = "drop"
    ) %>%
    unnest(curve) %>%
    mutate(model_type = "Non-BI Model"),
  
  # BI curves at median
  curves_bi %>%
    filter(BI_level == "Median BI") %>%
    mutate(model_type = "BI Model (Median)")
) %>%
  left_join(
    # Add chosen thresholds for vertical lines
    bind_rows(
      threshold_table %>% 
        select(species, chosen_threshold) %>% 
        mutate(model_type = "Non-BI Model"),
      threshold_table_bi %>% 
        filter(BI_level == "Median BI") %>%
        select(species, chosen_threshold) %>%
        mutate(model_type = "BI Model (Median)")
    ),
    by = c("species", "model_type")
  )

figure1 <- ggplot(fig1_data, 
                  aes(x = threshold, y = precision, 
                      color = model_type, linetype = model_type)) +
  geom_line(linewidth = 1) +
  geom_vline(aes(xintercept = chosen_threshold, color = model_type),
             linetype = "dotted", linewidth = 0.6, alpha = 0.8) +
  geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  facet_wrap(~ species, ncol = 1, scales = "free_x") +
  scale_color_manual(
    values = c("Non-BI Model" = "#1b9e77", "BI Model (Median)" = "#d95f02"),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c("Non-BI Model" = "solid", "BI Model (Median)" = "dashed"),
    name = NULL
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title = "Species-Specific Threshold Calibration: BI vs. Non-BI Models",
    subtitle = "Dashed horizontal line = target precision (0.90); Dotted vertical lines = chosen thresholds",
    x = "Confidence Threshold",
    y = "Precision of Retained Detections"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 10)
  )

print(figure1)
ggsave("Figure1_Precision_Curves_BI_vs_nonBI.png", 
       plot = figure1, width = 8, height = 10, dpi = 300)


# FIGURE 2: Effect size forest plot (dAUC with CIs)
# For this we'd need bootstrap CIs for dAUC - simplified version here

figure2_data <- final_summary_table %>%
  mutate(
    dAUC_direction = ifelse(dAUC > 0, "BI Improves", "BI Worsens"),
    species_ordered = forcats::fct_reorder(species, dAUC)
  )

figure2 <- ggplot(figure2_data, 
                  aes(x = dAUC, y = species_ordered, color = dAUC_direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dotted", 
             color = "gray60", alpha = 0.5) +
  geom_point(size = 4) +
  # Add error bars if you have bootstrap CIs for AUC
  # geom_errorbarh(aes(xmin = dAUC_ci_lower, xmax = dAUC_ci_upper), height = 0.2) +
  scale_color_manual(
    values = c("BI Improves" = "#1b9e77", "BI Worsens" = "#d95f02"),
    name = NULL
  ) +
  labs(
    title = "Effect of BI Model on Discrimination Performance",
    subtitle = "Dotted lines mark ±0.01 AUC (practical significance threshold)",
    x = "Change in AUC (BI Model − Non-BI Model)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(figure2)
ggsave("Figure2_Effect_Sizes_dAUC.png", 
       plot = figure2, width = 7, height = 5, dpi = 300)


# FIGURE 3: Threshold heatmap (species × BI level)

figure3_data <- threshold_table_bi %>%
  mutate(
    BI_level = factor(BI_level, 
                     levels = c("Low BI (25%)", "Median BI", "High BI (75%)"))
  )

figure3 <- ggplot(figure3_data, 
                  aes(x = BI_level, y = species, fill = chosen_threshold)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = round(chosen_threshold, 2)), 
            color = "white", fontface = "bold", size = 5) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
    midpoint = 0.5,
    name = "Chosen\nThreshold",
    limits = c(0, 1)
  ) +
  labs(
    title = "Species-Specific Thresholds across BI Conditions",
    subtitle = "Lower thresholds (blue) = more permissive; Higher thresholds (red) = more conservative",
    x = "Bioacoustic Index Condition",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid = element_blank()
  )

print(figure3)
ggsave("Figure3_Threshold_Heatmap_BI_levels.png", 
       plot = figure3, width = 7, height = 4, dpi = 300)


# ============================================================
# SUMMARY: What to check in your outputs
# ============================================================

cat("\n=== TASK COMPLETION SUMMARY ===\n")
cat("✓ Task 3: Check 'HawkEars_interaction_decision.csv' - use recommendations\n")
cat("✓ Task 2: Check threshold sensitivity in 'HawkEars_threshold_BI_sensitivity.csv'\n")
cat("✓ Task 13: Review 'HawkEars_species_interpretations.csv' and add ecological context\n")
cat("✓ Task 7: Bootstrap CIs in 'HawkEars_bootstrap_CIs_nonBI.csv' - flag wide CIs\n")
cat("✓ Task 11: Three publication figures saved as PNG files\n")
cat("\nNEXT STEPS:\n")
cat("1. Review interaction decisions - consider refitting additive models where recommended\n")
cat("2. Add species-specific ecological interpretations to the auto-generated text\n")
cat("3. If bootstrap CIs are wide (>0.20), consider collecting more validation data\n")
cat("4. Refine figures for publication (adjust colors, add annotations, etc.)\n")



# ============================================================
# REFIT BI MODELS WITHOUT INTERACTION TERM
# Compare: Non-BI vs Additive BI vs Interaction BI
# ============================================================

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(broom)

cat("\n=== REFITTING BI MODELS (ADDITIVE ONLY) ===\n")

# -------------------------
# 1) Fit additive BI models: tp ~ score + BI_z (NO interaction)
# -------------------------

fits_bi_additive <- val_bi %>%
  group_by(species) %>%
  mutate(BI_z = z(BI_raw)) %>%
  nest() %>%
  mutate(
    model_bi_additive = map(data, ~ glm(tp ~ score + BI_z, family = binomial(), data = .x)),
    n_val_bi = map_int(data, nrow)
  ) %>%
  select(species, model_bi_additive, n_val_bi)

cat("✓ Additive BI models fitted for", nrow(fits_bi_additive), "species\n")

# -------------------------
# 2) Calculate metrics for additive models
# -------------------------

additive_model_results <- fits_bi_additive %>%
  mutate(
    metrics = map2(model_bi_additive, species, ~ calc_glm_metrics(.x, val_bi_eval %>% filter(species == .y)))
  ) %>%
  unnest(metrics) %>%
  select(species, n_val_bi, AIC, pseudoR2, AUC, Brier) %>%
  arrange(species)

cat("✓ Additive model metrics calculated\n")

# -------------------------
# 3) Three-way model comparison: Non-BI vs Additive BI vs Interaction BI
# -------------------------

model_comparison_full <- nonbi_model_results %>%
  transmute(
    species,
    AIC_nonBI = AIC, 
    pseudoR2_nonBI = pseudoR2, 
    AUC_nonBI = AUC, 
    Brier_nonBI = Brier
  ) %>%
  left_join(
    additive_model_results %>%
      transmute(
        species,
        AIC_additive = AIC, 
        pseudoR2_additive = pseudoR2, 
        AUC_additive = AUC, 
        Brier_additive = Brier
      ),
    by = "species"
  ) %>%
  left_join(
    bi_model_results %>%
      transmute(
        species,
        AIC_interaction = AIC, 
        pseudoR2_interaction = pseudoR2, 
        AUC_interaction = AUC, 
        Brier_interaction = Brier
      ),
    by = "species"
  ) %>%
  mutate(
    # Compare additive to non-BI
    dAIC_additive = AIC_additive - AIC_nonBI,
    dAUC_additive = AUC_additive - AUC_nonBI,
    dBrier_additive = Brier_additive - Brier_nonBI,
    
    # Compare interaction to additive
    dAIC_interaction_vs_additive = AIC_interaction - AIC_additive,
    dAUC_interaction_vs_additive = AUC_interaction - AUC_additive,
    
    # Best model by AIC (lower is better)
    best_model_AIC = case_when(
      AIC_nonBI == pmin(AIC_nonBI, AIC_additive, AIC_interaction) ~ "Non-BI",
      AIC_additive == pmin(AIC_nonBI, AIC_additive, AIC_interaction) ~ "Additive BI",
      TRUE ~ "Interaction BI"
    ),
    
    # Best model by AUC (higher is better)
    best_model_AUC = case_when(
      AUC_nonBI == pmax(AUC_nonBI, AUC_additive, AUC_interaction) ~ "Non-BI",
      AUC_additive == pmax(AUC_nonBI, AUC_additive, AUC_interaction) ~ "Additive BI",
      TRUE ~ "Interaction BI"
    ),
    
    # Substantive improvement? (dAIC < -2 is standard threshold; dAUC > 0.01 is practical)
    additive_better_than_nonBI = dAIC_additive < -2 | dAUC_additive > 0.01,
    interaction_better_than_additive = dAIC_interaction_vs_additive < -2 | dAUC_interaction_vs_additive > 0.01,
    
    # Final recommendation
    recommended_model = case_when(
      !additive_better_than_nonBI ~ "Non-BI (BI not helpful)",
      additive_better_than_nonBI & interaction_better_than_additive ~ "Interaction BI",
      additive_better_than_nonBI & !interaction_better_than_additive ~ "Additive BI",
      TRUE ~ "Non-BI"
    )
  ) %>%
  arrange(species)

print("=== FULL THREE-WAY MODEL COMPARISON ===")
print(model_comparison_full %>% 
        select(species, AIC_nonBI, AIC_additive, AIC_interaction, 
               AUC_nonBI, AUC_additive, AUC_interaction,
               recommended_model))

write.csv(model_comparison_full, 
          "HawkEars_model_comparison_FULL_threeway.csv", 
          row.names = FALSE)

# -------------------------
# 4) Visual comparison of model performance
# -------------------------

# Reshape for plotting
model_performance_long <- model_comparison_full %>%
  select(species, AIC_nonBI, AIC_additive, AIC_interaction,
         AUC_nonBI, AUC_additive, AUC_interaction) %>%
  pivot_longer(
    cols = -species,
    names_to = c("metric", "model"),
    names_pattern = "(.*)_(.*)",
    values_to = "value"
  ) %>%
  mutate(
    model = factor(model, levels = c("nonBI", "additive", "interaction"),
                   labels = c("Non-BI", "Additive BI", "Interaction BI"))
  )

# AIC comparison plot (lower is better)
aic_plot <- model_performance_long %>%
  filter(metric == "AIC") %>%
  ggplot(aes(x = model, y = value, color = model, group = species)) +
  geom_line(color = "gray60") +
  geom_point(size = 3) +
  facet_wrap(~ species, scales = "free_y") +
  scale_color_manual(values = c("Non-BI" = "#1b9e77", 
                                 "Additive BI" = "#d95f02",
                                 "Interaction BI" = "#7570b3")) +
  labs(
    title = "Model Comparison: AIC (Lower is Better)",
    x = NULL,
    y = "AIC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(aic_plot)
ggsave("HawkEars_AIC_comparison_threeway.png", 
       plot = aic_plot, width = 8, height = 4, dpi = 300)

# AUC comparison plot (higher is better)
auc_plot <- model_performance_long %>%
  filter(metric == "AUC") %>%
  ggplot(aes(x = model, y = value, color = model, group = species)) +
  geom_line(color = "gray60") +
  geom_point(size = 3) +
  facet_wrap(~ species, scales = "free_y") +
  scale_color_manual(values = c("Non-BI" = "#1b9e77", 
                                 "Additive BI" = "#d95f02",
                                 "Interaction BI" = "#7570b3")) +
  labs(
    title = "Model Comparison: AUC (Higher is Better)",
    x = NULL,
    y = "AUC"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(auc_plot)
ggsave("HawkEars_AUC_comparison_threeway.png", 
       plot = auc_plot, width = 8, height = 4, dpi = 300)

# -------------------------
# 5) Compute thresholds using ADDITIVE BI models
# -------------------------

cat("\n=== COMPUTING THRESHOLDS WITH ADDITIVE BI MODELS ===\n")

# Recompute BI-conditional thresholds using additive models
threshold_table_bi_additive <- fits_bi_additive %>%
  left_join(bi_levels, by = "species") %>%
  left_join(lab_bins %>%
              group_by(species) %>%
              summarise(total_detections = sum(N_i), .groups = "drop"),
            by = "species") %>%
  left_join(lab_bins, by = "species") %>%
  group_by(species, BI_level, BI_z, model_bi_additive, n_val_bi, total_detections) %>%
  summarise(
    out = list(
      compute_tseng_threshold_bi(
        cur_data_all() %>% select(bin_lower, bin_mid, N_i),
        model_bi_additive[[1]],
        BI_z_fixed = BI_z
      )
    ),
    .groups = "drop"
  ) %>%
  mutate(
    chosen_threshold = map_dbl(out, ~ .x$summary$threshold),
    precision_at_threshold = map_dbl(out, ~ .x$summary$precision),
    prop_retained = map_dbl(out, ~ .x$summary$prop_retained),
    n_retained = map_int(out, ~ .x$summary$n_retained),
    status = map_chr(out, ~ .x$summary$status),
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

threshold_table_bi_additive_pretty <- threshold_table_bi_additive %>%
  mutate(
    chosen_threshold = round(chosen_threshold, 2),
    precision_at_threshold = round(precision_at_threshold, 3),
    prop_retained = round(prop_retained, 3)
  )

print("=== ADDITIVE BI MODEL THRESHOLDS ===")
print(threshold_table_bi_additive_pretty)

write.csv(threshold_table_bi_additive_pretty, 
          "HawkEars_thresholds_BI_ADDITIVE.csv", 
          row.names = FALSE)

# -------------------------
# 6) Compare thresholds: Non-BI vs Additive BI vs Interaction BI
# -------------------------

threshold_comparison <- threshold_table_pretty %>%
  select(species, chosen_threshold_nonBI = chosen_threshold, 
         prop_retained_nonBI = prop_retained) %>%
  left_join(
    threshold_table_bi_additive_pretty %>%
      filter(BI_level == "Median BI") %>%
      select(species, chosen_threshold_additive = chosen_threshold,
             prop_retained_additive = prop_retained),
    by = "species"
  ) %>%
  left_join(
    threshold_table_bi_pretty %>%
      filter(BI_level == "Median BI") %>%
      select(species, chosen_threshold_interaction = chosen_threshold,
             prop_retained_interaction = prop_retained),
    by = "species"
  ) %>%
  mutate(
    threshold_diff_additive = chosen_threshold_additive - chosen_threshold_nonBI,
    threshold_diff_interaction = chosen_threshold_interaction - chosen_threshold_nonBI,
    retention_diff_additive = prop_retained_additive - prop_retained_nonBI,
    retention_diff_interaction = prop_retained_interaction - prop_retained_nonBI
  )

print("=== THRESHOLD COMPARISON (at Median BI) ===")
print(threshold_comparison)

write.csv(threshold_comparison, 
          "HawkEars_threshold_comparison_threeway.csv", 
          row.names = FALSE)

# -------------------------
# 7) Generate curves for additive BI models
# -------------------------

curves_bi_additive <- fits_bi_additive %>%
  left_join(bi_levels, by = "species") %>%
  left_join(lab_bins, by = "species") %>%
  group_by(species, BI_level, BI_z) %>%
  summarise(
    model_bi_additive = list(fits_bi_additive$model_bi_additive[match(first(species), fits_bi_additive$species)][[1]]),
    bins = list(cur_data() %>% select(bin_lower, bin_mid, N_i)),
    .groups = "drop"
  ) %>%
  mutate(curve = pmap(list(bins, model_bi_additive, BI_z),
                      ~ compute_tseng_threshold_bi(..1, ..2, ..3)$curve)) %>%
  select(species, BI_level, curve) %>%
  unnest(curve)

# Plot additive BI curves for comparison
for (sp in unique(threshold_table_bi_additive$species)) {
  sp_curves <- curves_bi_additive %>% filter(species == sp)
  sp_thresholds <- threshold_table_bi_additive %>% filter(species == sp)
  
  p <- ggplot(sp_curves, aes(x = threshold, y = precision, color = BI_level)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray40") +
    geom_vline(data = sp_thresholds, 
               aes(xintercept = chosen_threshold, color = BI_level),
               linetype = "dotted", linewidth = 0.8) +
    scale_color_manual(
      values = c("Low BI (25%)" = "#2166ac", 
                 "Median BI" = "#762a83", 
                 "High BI (75%)" = "#c51b7d")
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(sp, " — ADDITIVE BI Model Thresholds"),
      subtitle = "Vertical dotted lines show chosen thresholds for each BI level",
      x = "Confidence Threshold",
      y = "Precision of Retained Detections",
      color = "BI Condition"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p)
  ggsave(paste0("HawkEars_BI_curves_ADDITIVE_", sp, ".png"), 
         plot = p, width = 8, height = 6, dpi = 300)
}

# -------------------------
# 8) Side-by-side comparison plot: Additive vs Interaction BI
# -------------------------

curves_comparison <- bind_rows(
  curves_bi_additive %>% mutate(model_type = "Additive BI"),
  curves_bi %>% mutate(model_type = "Interaction BI")
)

for (sp in unique(curves_comparison$species)) {
  sp_curves <- curves_comparison %>% filter(species == sp, BI_level == "Median BI")
  
  p <- ggplot(sp_curves, aes(x = threshold, y = precision, 
                              color = model_type, linetype = model_type)) +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = c("Additive BI" = "#d95f02", 
                                   "Interaction BI" = "#7570b3")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(sp, " — Additive vs Interaction BI Models (Median BI)"),
      x = "Confidence Threshold",
      y = "Precision of Retained Detections",
      color = "Model Type",
      linetype = "Model Type"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p)
  ggsave(paste0("HawkEars_additive_vs_interaction_", sp, ".png"), 
         plot = p, width = 8, height = 6, dpi = 300)
}

# -------------------------
# 9) FINAL SUMMARY TABLE with all models
# -------------------------

final_summary_all_models <- threshold_table_pretty %>%
  rename(
    chosen_threshold_nonBI = chosen_threshold,
    precision_at_threshold_nonBI = precision_at_threshold,
    prop_retained_nonBI = prop_retained,
    n_retained_nonBI = n_retained,
    n_dropped_nonBI = n_dropped,
    status_nonBI = status
  ) %>%
  left_join(
    threshold_table_bi_additive_pretty %>%
      filter(BI_level == "Median BI") %>%
      select(species, 
             chosen_threshold_additive = chosen_threshold, 
             precision_at_threshold_additive = precision_at_threshold, 
             prop_retained_additive = prop_retained,
             n_retained_additive = n_retained, 
             n_dropped_additive = n_dropped, 
             status_additive = status),
    by = "species"
  ) %>%
  left_join(
    threshold_table_bi_pretty %>%
      filter(BI_level == "Median BI") %>%
      select(species, 
             chosen_threshold_interaction = chosen_threshold, 
             precision_at_threshold_interaction = precision_at_threshold, 
             prop_retained_interaction = prop_retained,
             n_retained_interaction = n_retained, 
             n_dropped_interaction = n_dropped, 
             status_interaction = status),
    by = "species"
  ) %>%
  left_join(model_comparison_full, by = "species") %>%
  arrange(species)

print("=== FINAL SUMMARY: ALL THREE MODELS ===")
print(final_summary_all_models %>% 
        select(species, 
               chosen_threshold_nonBI, chosen_threshold_additive, chosen_threshold_interaction,
               AUC_nonBI, AUC_additive, AUC_interaction,
               recommended_model))

write.csv(final_summary_all_models, 
          "HawkEars_FINAL_summary_all_models.csv", 
          row.names = FALSE)

# -------------------------
# 10) DECISION SUMMARY with clear recommendations
# -------------------------

decision_summary <- model_comparison_full %>%
  select(species, recommended_model, 
         dAIC_additive, dAUC_additive, 
         dAIC_interaction_vs_additive, dAUC_interaction_vs_additive) %>%
  left_join(
    threshold_comparison %>%
      select(species, chosen_threshold_nonBI, 
             chosen_threshold_additive, chosen_threshold_interaction),
    by = "species"
  ) %>%
  mutate(
    rationale = case_when(
      recommended_model == "Non-BI (BI not helpful)" ~ 
        "BI provides no meaningful improvement in discrimination (dAUC < 0.01 and dAIC > -2). Use simple non-BI threshold.",
      
      recommended_model == "Additive BI" ~ 
        paste0("BI improves discrimination (dAUC = +", round(dAUC_additive, 3), 
               ") but interaction term not justified. Use additive BI model with context-dependent thresholds."),
      
      recommended_model == "Interaction BI" ~ 
        paste0("Both BI and BI×score interaction improve discrimination. ",
               "Use interaction BI model with context-dependent thresholds."),
      
      TRUE ~ "Review metrics manually"
    ),
    
    deployment_threshold = case_when(
      recommended_model == "Non-BI (BI not helpful)" ~ chosen_threshold_nonBI,
      recommended_model == "Additive BI" ~ chosen_threshold_additive,
      recommended_model == "Interaction BI" ~ chosen_threshold_interaction,
      TRUE ~ chosen_threshold_nonBI
    )
  )

print("=== DEPLOYMENT DECISION SUMMARY ===")
print(decision_summary %>% select(species, recommended_model, deployment_threshold, rationale))

write.csv(decision_summary, 
          "HawkEars_DEPLOYMENT_decisions.csv", 
          row.names = FALSE)

# -------------------------
# SUMMARY OUTPUT
# -------------------------

cat("\n")
cat("============================================================\n")
cat("  ADDITIVE BI MODEL ANALYSIS COMPLETE\n")
cat("============================================================\n")

for (i in 1:nrow(decision_summary)) {
  cat(sprintf("  %s: %s (threshold = %.2f)\n", 
              decision_summary$species[i], 
              decision_summary$recommended_model[i],
              decision_summary$deployment_threshold[i]))
}


## N per combo of score and BI_z 

r# ============================================================
# SAMPLE SIZE DISTRIBUTION: Score × BI_z Combinations
# ============================================================

cat("\n=== ANALYZING SAMPLE SIZE ACROSS SCORE × BI_z SPACE ===\n")

# Create bins for score and BI_z
sample_coverage <- val_bi_eval %>%
  mutate(
    # Bin score into 0.1 intervals for visualization
    score_bin = cut(score, 
                    breaks = seq(0, 1, by = 0.1),
                    include.lowest = TRUE,
                    labels = paste0(seq(0, 0.9, by = 0.1), "-", seq(0.1, 1.0, by = 0.1))),
    
    # Bin BI_z into quartile-based categories
    BI_z_category = cut(BI_z,
                        breaks = c(-Inf, -0.674, 0, 0.674, Inf),  # approximate quartiles for std normal
                        labels = c("Very Low BI", "Low BI", "High BI", "Very High BI"))
  ) %>%
  count(species, score_bin, BI_z_category, name = "n_obs") %>%
  filter(!is.na(score_bin), !is.na(BI_z_category))

print("Sample Size Summary:")
print(sample_coverage %>% 
        group_by(species) %>% 
        summarise(
          total_obs = sum(n_obs),
          mean_per_combo = mean(n_obs),
          min_per_combo = min(n_obs),
          max_per_combo = max(n_obs),
          n_combos = n(),
          .groups = "drop"
        ))

# PLOT 1: Heatmap of sample sizes (score × BI_z)
heatmap_plot <- ggplot(sample_coverage, 
                       aes(x = score_bin, y = BI_z_category, fill = n_obs)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_obs), color = "white", size = 3, fontface = "bold") +
  facet_wrap(~ species, ncol = 1) +
  scale_fill_viridis_c(option = "plasma", name = "Sample\nSize (n)") +
  labs(
    title = "Sample Size Distribution: Confidence Score × BI Condition",
    subtitle = "Each cell shows number of validation observations",
    x = "Confidence Score Bin",
    y = "Bioacoustic Index Category"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(heatmap_plot)
ggsave("HawkEars_sample_coverage_heatmap.png", 
       plot = heatmap_plot, width = 10, height = 8, dpi = 300)

# PLOT 2: Stacked bar chart by species
stacked_bar <- ggplot(sample_coverage, 
                      aes(x = score_bin, y = n_obs, fill = BI_z_category)) +
  geom_col(position = "stack", width = 0.8) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(
    values = c("Very Low BI" = "#2166ac", 
               "Low BI" = "#92c5de", 
               "High BI" = "#f4a582", 
               "Very High BI" = "#b2182b"),
    name = "BI Condition"
  ) +
  labs(
    title = "Sample Size Distribution Across Confidence Score Bins",
    subtitle = "Stacked by BI condition",
    x = "Confidence Score Bin",
    y = "Number of Observations"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

print(stacked_bar)
ggsave("HawkEars_sample_coverage_stacked_bar.png", 
       plot = stacked_bar, width = 10, height = 8, dpi = 300)

# PLOT 3: Fine-grained scatter plot (actual continuous values)
scatter_coverage <- val_bi_eval %>%
  mutate(score_round = round(score, 2)) %>%
  count(species, score_round, BI_z, name = "n_combo") %>%
  mutate(n_category = cut(n_combo, 
                          breaks = c(0, 1, 3, 5, 10, Inf),
                          labels = c("1", "2-3", "4-5", "6-10", ">10")))

scatter_plot <- ggplot(scatter_coverage, 
                       aes(x = score_round, y = BI_z, size = n_combo, color = n_category)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ species, ncol = 1) +
  scale_size_continuous(range = c(1, 8), name = "Sample Size") +
  scale_color_viridis_d(option = "rocket", name = "Sample\nSize\nCategory") +
  labs(
    title = "Validation Data Coverage: Score × BI_z Space",
    subtitle = "Point size and color indicate number of observations per combination",
    x = "Confidence Score",
    y = "Standardized BI (BI_z)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(scatter_plot)
ggsave("HawkEars_sample_coverage_scatter.png", 
       plot = scatter_plot, width = 10, height = 8, dpi = 300)

# PLOT 4: Marginal distributions to show where data is concentrated
marginal_score <- val_bi_eval %>%
  ggplot(aes(x = score, fill = species)) +
  geom_histogram(bins = 20, color = "white", alpha = 0.7) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("COYE" = "#1b9e77", 
                                "MAWA" = "#d95f02", 
                                "VEER" = "#7570b3")) +
  labs(
    title = "Marginal Distribution: Confidence Scores",
    x = "Confidence Score",
    y = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

print(marginal_score)
ggsave("HawkEars_marginal_score.png", 
       plot = marginal_score, width = 8, height = 6, dpi = 300)

# Check for sparse regions (combinations with <5 observations)
sparse_regions <- sample_coverage %>%
  filter(n_obs < 5) %>%
  arrange(species, n_obs)

cat("\n=== SPARSE REGIONS (n < 5) ===\n")
print(sparse_regions)

if (nrow(sparse_regions) > 0) {
  cat("\nWARNING: Found", nrow(sparse_regions), "score×BI combinations with <5 observations.\n")
  cat("These regions may have unstable model predictions.\n")
}

# Summary statistics per species
coverage_summary <- val_bi_eval %>%
  group_by(species) %>%
  summarise(
    n_total = n(),
    score_min = min(score),
    score_max = max(score),
    score_median = median(score),
    BI_z_min = min(BI_z),
    BI_z_max = max(BI_z),
    BI_z_median = median(BI_z),
    # Count unique score bins (0.05 resolution)
    n_unique_score_bins = n_distinct(round(score / 0.05) * 0.05),
    # Correlation between score and BI_z (to check if they're confounded)
    score_BI_cor = cor(score, BI_z, use = "complete.obs"),
    .groups = "drop"
  )

cat("\n=== COVERAGE SUMMARY ===\n")
print(coverage_summary)

write.csv(coverage_summary, 
          "HawkEars_sample_coverage_summary.csv", 
          row.names = FALSE)

write.csv(sample_coverage, 
          "HawkEars_sample_coverage_by_bin.csv", 
          row.names = FALSE)

# Check for correlation between score and BI_z (important diagnostic!)
cor_plot <- ggplot(val_bi_eval, aes(x = score, y = BI_z)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
  facet_wrap(~ species) +
  labs(
    title = "Correlation Between Confidence Score and BI_z",
    subtitle = "Red line = linear fit; Check for confounding",
    x = "Confidence Score",
    y = "Standardized BI (BI_z)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

print(cor_plot)
ggsave("HawkEars_score_BI_correlation.png", 
       plot = cor_plot, width = 10, height = 4, dpi = 300)

