# 🦅 HawkEars Validator

An R Shiny app for rapid visual and auditory validation of bird detections from [HawkEars](https://github.com/jhuus/HawkEars) (or similar CNN/embedding-based recognizers).

HawkEars Validator takes the model's output CSV — with filenames, start/end times, species, and scores — cuts 3-second audio snippets from long ARU recordings, and serves them up with synchronized spectrograms so you can efficiently **confirm**, **reject**, or **flag** predictions as unsure.

<img width="1076" height="869" alt="Screenshot of HawkEars Validator showing spectrogram with playhead bar and validation buttons" src="https://github.com/user-attachments/assets/abac8c97-5af9-4ffa-aa38-b7bb3e8bc841" />

---

## What It Does

### Reads a HawkEars-style CSV

Expects a CSV (e.g. `HawkEars_labels.csv`) with the following columns:

| Column | Description |
|---|---|
| `filename` | Full WAV file name |
| `start_time`, `end_time` | Clip bounds in seconds |
| `class_name`, `class_code` | Species names and codes |
| `score` | Classifier confidence |

### Extracts and displays audio clips

Locates matching WAVs on disk, then:

- Extracts **3-second snippets** on the fly using `tuneR`
- Renders a **high-contrast spectrogram** via `av::read_audio_fft`
- Displays a horizontal **playhead bar** that moves along the bottom of the spectrogram as audio plays

### Supports four validation labels

| Label | Meaning |
|---|---|
| ✅ **Yes** | Prediction correct |
| ❌ **No** | Prediction wrong |
| ❓ **Unsure** | Unclear — needs a second opinion |
| ⏭️ **Skip** | Skip for now |

### Writes results to CSV

Decisions are saved to `HawkEars_validation_results.csv` with: `clip_id`, file name, time bounds, species, score, `validator_id` (your initials), `label`, and `validated_at` timestamp.

---

## Queue Management

On startup (or when you click **Start / Refresh queue**), the app:

1. Reads the existing `HawkEars_validation_results.csv`
2. Removes any clips already marked **Yes** or **No**
3. Keeps **Unsure** and **Skip** clips in the pool for re-validation

**Optional cleanup:** toggle deletion of temporary 3-second snippets after each decision.

A running counter tracks clips labeled this session, and every 25 clips a modal pops up with a random boreal/grassland bird-themed compliment — because morale is a data quality issue. 🐦

---

## UI Overview

### Sidebar

- Enter **validator ID** (initials)
- Filter by **target species**
- Toggle **snippet deletion**
- View clip **metadata** and **queue progress**

### Main Panel

- **Predicted species + score**
- **Spectrogram** with moving playhead bar
- **Audio player** (play/pause, scrub, volume)
- Large, stacked buttons for **Yes / No / Unsure**, with **Skip** separate

---

## Requirements

**R packages:**

```r
install.packages(c("shiny", "dplyr", "readr", "tuneR", "av"))
```

**Data:**

- A HawkEars-style CSV and matching WAV files on disk

**Configuration** — set the following paths at the top of `app.R`:

```r
labels_file  <- "path/to/HawkEars_labels.csv"
audio_dir    <- "path/to/wav/files/"
results_file <- "path/to/HawkEars_validation_results.csv"
```

---

## Typical Workflow

1. Run `shiny::runApp("path/to/app")`
2. Enter your initials and (optionally) choose a species filter
3. Click **Start / Refresh queue** to load all clips not yet labeled Yes or No
4. For each clip: watch the spectrogram, play the audio, and choose **Yes**, **No**, **Unsure**, or **Skip**
5. Close the app whenever you're done — next time you start it, the app resumes from where you left off, serving only unvalidated clips
