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
    species = class_code,  # use class_code to avoid name mismatches
    score   = as.numeric(score),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species)) %>%
  filter(score >= min_conf, score <= max_conf)

# Fit logistic regression per species (Tseng: GLM binomial, score predictor) 
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




##### ok now building on this GLM with bioacoustic indices ####
## does soundscape interference matter?? #

library(tidyverse)
library(tuneR)
library(seewave)

# -------------------------
# Paths
# -------------------------
val_path   <- "E:/BAR-LT_LocalizationProject/localization_05312025/hawkears_lowthresh/HawkEars_validation_results.csv"
audio_root <- "E:/BAR-LT_LocalizationProject/localization_05312025/localizationtrim_new"

# -------------------------
# Noise metric band definitions (can edit if we think of other things)
# -------------------------
LF_band <- c(0, 1000)          # wind/handling/engines often live here
HF_band <- c(6000, 10000)      # insects/hiss often live here
entropy_band <- c(200, 10000)  # avoid DC/very low junk

# -------------------------
# OPTIONAL BUT HIGHLY RECOMMENDED (MEG THINKS): species call-frequency windows (Hz)

# -------------------------
call_bands <- tribble(
  ~species, ~call_low, ~call_high,
  "VEER",   2000,      6000
  # "ADD A NICE BIRD HERE",   300,       900,
  # "ADD ANOTHER COOL BIRD, EWPW PERHAPS?",   3000,      8000,
)

# -------------------------
# 1) Read validation + build wav index
# -------------------------
validation <- readr::read_csv(val_path, show_col_types = FALSE) %>%
  mutate(
    score = as.numeric(score),
    species = as.character(class_code),
    tp = case_when(
      tolower(label) %in% c("yes","y","true","1","present") ~ 1L,
      tolower(label) %in% c("no","n","false","0","absent")  ~ 0L,
      TRUE ~ NA_integer_
    ),
    dur = as.numeric(end_time - start_time)
  ) %>%
  filter(!is.na(tp), !is.na(score), !is.na(species)) %>%
  filter(dur > 0)

wav_index <- tibble(
  wav_path = list.files(audio_root, pattern = "\\.wav$", recursive = TRUE, full.names = TRUE)
) %>%
  mutate(filename = basename(wav_path)) %>%
  distinct(filename, .keep_all = TRUE)  # assumes basename is unique

val2 <- validation %>%
  left_join(wav_index, by = "filename") %>%
  left_join(call_bands, by = c("species" = "species"))

if (any(is.na(val2$wav_path))) {
  print(val2 %>% filter(is.na(wav_path)) %>% distinct(filename) %>% head(50))
  stop("Some validation filenames were not found under audio_root.")
}

# -------------------------
# 2) Audio reading: exact 3-sec window
# -------------------------
read_clip <- function(path, start_sec, end_sec) {
  w <- tuneR::readWave(path, from = start_sec, to = end_sec, units = "seconds")
  if (w@stereo) w <- tuneR::mono(w, which = "left")
  w
}

# -------------------------
# 3) Spectral metrics from a Wave object
#    - band powers from a power spectrum
#    - Shannon spectral entropy (0..1)
#    -  call band power + simple SNR-like ratio
# -------------------------
power_spectrum_fft <- function(w) {
  if (w@stereo) w <- tuneR::mono(w, which = "left")
  
  x <- w@left / (2^(w@bit - 1))      # scale to ~[-1,1]
  n <- length(x)
  if (n < 64) return(NULL)
  
  # Hann window
  win <- 0.5 - 0.5 * cos(2*pi*(0:(n-1))/(n-1))
  xw  <- x * win
  
  X <- fft(xw)
  
  # one-sided spectrum
  kmax <- floor(n/2) + 1
  P <- Mod(X[1:kmax])^2
  f <- (0:(kmax-1)) * (w@samp.rate / n)
  
  list(freq = f, power = P)
}

bandpower_from_spec <- function(freq, power, f_lo, f_hi) {
  if (!is.finite(f_lo) || !is.finite(f_hi) || f_hi <= f_lo) return(NA_real_)
  
  keep <- is.finite(freq) & is.finite(power) & (freq >= f_lo) & (freq < f_hi)
  if (sum(keep) == 0) return(0)        # 0 is often nicer than NA for "no energy in band"
  sum(power[keep], na.rm = TRUE)
}

spectral_entropy <- function(freq, power, f_lo, f_hi) {
  if (!is.finite(f_lo) || !is.finite(f_hi) || f_hi <= f_lo) return(NA_real_)
  
  keep <- is.finite(freq) & is.finite(power) & (freq >= f_lo) & (freq < f_hi)
  if (sum(keep) < 10) return(NA_real_)
  
  p <- power[keep]
  p <- p[p > 0]
  if (length(p) < 10) return(NA_real_)
  
  p <- p / sum(p)
  H <- -sum(p * log(p))
  H / log(length(p))  # normalized 0..1
}

calc_metrics_3sec <- function(path, start_sec, end_sec,
                              LF_band, HF_band, entropy_band,
                              call_low = NA_real_, call_high = NA_real_) {
  out <- tryCatch({
    w <- read_clip(path, start_sec, end_sec)
    
    nyq <- w@samp.rate / 2
    lf <- c(max(0, LF_band[1]), min(nyq, LF_band[2]))
    hf <- c(max(0, HF_band[1]), min(nyq, HF_band[2]))
    eb <- c(max(0, entropy_band[1]), min(nyq, entropy_band[2]))
    
    sp <- power_spectrum_fft(w)
    if (is.null(sp)) stop("spectrum_failed")
    
    freq  <- sp$freq
    power <- sp$power
    
    LF_power <- bandpower_from_spec(freq, power, lf[1], lf[2])
    HF_power <- bandpower_from_spec(freq, power, hf[1], hf[2])
    entropy  <- spectral_entropy(freq, power, eb[1], eb[2])
    
    CALL_power <- NA_real_
    OFF_power  <- NA_real_
    CALL_SNR   <- NA_real_
    
    if (is.finite(call_low) && is.finite(call_high)) {
      cl <- c(max(0, call_low), min(nyq, call_high))
      if (cl[2] > cl[1]) {
        CALL_power <- bandpower_from_spec(freq, power, cl[1], cl[2])
        
        # Off-band = entropy band excluding call band
        in_call <- is.finite(freq) & (freq >= cl[1]) & (freq < cl[2])
        in_ent  <- is.finite(freq) & (freq >= eb[1]) & (freq < eb[2])
        off_idx <- in_ent & !in_call
        
        OFF_power <- if (any(off_idx, na.rm = TRUE)) sum(power[off_idx], na.rm = TRUE) else NA_real_
        
        CALL_SNR <- if (is.finite(CALL_power) && is.finite(OFF_power) && OFF_power > 0) {
          CALL_power / OFF_power
        } else NA_real_
      }
    }
    
    tibble(
      LF_power = LF_power,
      HF_power = HF_power,
      entropy  = entropy,
      CALL_power = CALL_power,
      CALL_SNR   = CALL_SNR,
      ok = TRUE,
      err = NA_character_
    )
  }, error = function(e) {
    tibble(
      LF_power = NA_real_, HF_power = NA_real_, entropy = NA_real_,
      CALL_power = NA_real_, CALL_SNR = NA_real_,
      ok = FALSE, err = as.character(e$message)
    )
  })
  
  out
}

# -------------------------
# 4) Compute metrics for each validated row
# -------------------------
val_metrics <- val2 %>%
  mutate(metrics = pmap(
    list(wav_path, start_time, end_time, call_low, call_high),
    ~ calc_metrics_3sec(..1, ..2, ..3,
                        LF_band = LF_band, HF_band = HF_band, entropy_band = entropy_band,
                        call_low = ..4, call_high = ..5)
  )) %>%
  unnest(metrics)

# Inspect failures
val_metrics %>%
  filter(!ok) %>%
  select(filename, start_time, end_time, err) %>%
  print(n = 50)

# Keep successful rows
dat <- val_metrics %>%
  filter(ok) %>%
  drop_na(LF_power, HF_power, entropy)

# -------------------------
# 5) Fit calibration GLM
#  main effects + one interaction you care about.
#  note to megan --> just doing one interaction for now because we probably need more validations
# -------------------------
# Manual z-scoring (so predict() behaves)
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

dat <- dat %>%
  mutate(
    LF_z = z(LF_power),
    HF_z = z(HF_power),
    ent_z = z(entropy),
    CALL_z = if (all(is.na(CALL_power))) NA_real_ else z(CALL_power),
    SNR_z  = if (all(is.na(CALL_SNR))) NA_real_ else z(CALL_SNR)
  )

m <- glm(tp ~ score + LF_z + HF_z + ent_z + score:LF_z,
         data = dat, family = binomial())

summary(m)


##  plots 


library(data.table)
library(tidyverse)
library(ggplot2)

# Tseng binning settings
min_conf <- 0.10
max_conf <- 1.00
step <- 0.05
threshold_grid <- seq(min_conf, 0.95, by = step)
target_precision <- 0.90

# Build N_i (how many detections occur in each 0.05 bin) from the FULL labels file
labels_dt <- as.data.table(labels)
labels_dt[, score := as.numeric(score)]
labels_dt <- labels_dt[!is.na(score)]
labels_dt <- labels_dt[score >= min_conf & score <= max_conf]

# Only species validated
validated_species <- unique(dat$species) %||% unique(dat$class_code)


labels_dt <- labels_dt[class_code %in% validated_species]

labels_dt[, bin_lower := min_conf + floor((score - min_conf)/step)*step]
labels_dt[bin_lower < min_conf, bin_lower := min_conf]
labels_dt[bin_lower > 0.95, bin_lower := 0.95]
labels_dt[, bin_mid := bin_lower + step/2]


Ni <- labels_dt[, .(N_i = .N), by = .(species = class_code, bin_lower, bin_mid)] %>%
  as_tibble()

Ni %>% arrange(species, bin_lower) %>% print(n = 50)

### HELLO - STOP RIGHT THERE! 
### THIS IS IMPORTANT FOR PLOTTING AND PICKING THRESHOLDS!!!!! 

# maybe some of these scenarios below would be based on different weaather? 

# Define LF conditions using validated-data, hold HF/entropy at medians aka constant
LF_lo <- quantile(dat$LF_z, 0.25, na.rm = TRUE)
LF_hi <- quantile(dat$LF_z, 0.75, na.rm = TRUE)
HF_med  <- median(dat$HF_z,  na.rm = TRUE)
ent_med <- median(dat$ent_z, na.rm = TRUE)

cond_lf <- tibble(
  LF_level = c("Low LF (25%)", "High LF (75%)"),
  LF_z  = c(LF_lo, LF_hi),
  HF_z  = c(HF_med, HF_med),
  ent_z = c(ent_med, ent_med)
)


#YOU COULD SWAP THE ABOVE AND PLAY WITH SETTINGS 
ent_lo <- quantile(dat$ent_z, 0.25, na.rm = TRUE)
ent_hi <- quantile(dat$ent_z, 0.75, na.rm = TRUE)
LF_med <- median(dat$LF_z, na.rm = TRUE)
HF_med <- median(dat$HF_z, na.rm = TRUE)

cond_ent <- tibble(
  level = c("Low entropy (25%)", "High entropy (75%)"),
  LF_z  = c(LF_med, LF_med),
  HF_z  = c(HF_med, HF_med),
  ent_z = c(ent_lo, ent_hi)
)

tseng_curve_one_condition <- function(Ni_one_species, cond_row) {
  
  # Explicit newdata with the exact variable names in the model
  newdata <- tibble(
    score = Ni_one_species$bin_mid,
    LF_z  = cond_row$LF_z,
    HF_z  = cond_row$HF_z,
    ent_z = cond_row$ent_z
  )
  
  TPR_i <- as.numeric(predict(m, newdata = newdata, type = "response"))
  
  pred <- Ni_one_species %>%
    mutate(TPR_i = TPR_i)
  
  purrr::map_dfr(threshold_grid, function(T) {
    kept <- pred %>% filter(bin_lower >= T)
    denom <- sum(kept$N_i)
    
    precision <- if (denom == 0) NA_real_ else sum(kept$TPR_i * kept$N_i) / denom
    prop_retained <- if (sum(pred$N_i) == 0) NA_real_ else sum(kept$N_i) / sum(pred$N_i)
    
    tibble(threshold = T, precision = precision, prop_retained = prop_retained)
  })
}
# Curves for each species and each   LF_level
curves_weighted_lf <- Ni %>%
  group_by(species) %>%
  group_modify(~{
    bind_rows(lapply(1:nrow(cond_lf), function(j) {
      cdf <- tseng_curve_one_condition(.x %>% select(bin_lower, bin_mid, N_i), cond_lf[j, ])
      cdf$LF_level <- cond_lf$LF_level[j]
      cdf
    }))
  }) %>%
  ungroup()

# Pick minimum threshold reaching 0.90 precision 
thresholds_weighted_lf <- curves_weighted_lf %>%
  group_by(species, LF_level) %>%
  summarise(
    thr = {tmp <- threshold[!is.na(precision) & precision >= target_precision];
    if (length(tmp)==0) NA_real_ else min(tmp)},
    precision_at_thr = ifelse(is.na(thr), NA_real_, precision[match(thr, threshold)]),
    prop_retained_at_thr = ifelse(is.na(thr), NA_real_, prop_retained[match(thr, threshold)]),
    .groups = "drop"
  )

thresholds_weighted_lf

ggplot(curves_weighted_lf, aes(threshold, precision, linetype = LF_level)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = target_precision, linetype = "dashed") +
  facet_wrap(~ species) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(
    x = "Confidence threshold (T)",
    y = "Precision of retained detections",
    linetype = "Low-frequency noise level",
    title = "Tseng-weighted precision vs threshold, conditional on LF noise",
    subtitle = "Weighted by full dataset score distribution"
  ) +
  theme_minimal(base_size = 13)+ scale_y_continuous(breaks = seq(0.5, 1, by = 0.1))

