
**
HawkEars Validator

HawkEars Validator is an R Shiny app for rapid visual and auditory validation of bird detections from HawkEars (or similar CNN/embedding-based recognizers). It takes the model’s output CSV (with filenames, start/end times, species, and scores), cuts 3-second audio snippets from long ARU recordings, and serves them up with synchronized spectrograms so you can efficiently confirm, reject, or flag predictions as unsure.**
<img width="1076" height="869" alt="image" src="https://github.com/user-attachments/assets/abac8c97-5af9-4ffa-aa38-b7bb3e8bc841" />

**What it does**

**Reads a HawkEars-style CSV (e.g. HawkEars_labels.csv) with:**

filename – full WAV file name

start_time, end_time – clip bounds in seconds

class_name, class_code – species names/codes

score – classifier confidence

Locates the corresponding WAVs on disk and:

Extracts 3-second snippets on the fly using tuneR

Renders a high-contrast spectrogram using av::read_audio_fft

Displays a horizontal “playhead” bar that moves along the bottom of the spectrogram as audio plays

**Lets the user label each clip as:**

Yes – prediction correct

No – prediction wrong

Unsure – unclear / needs a second opinion

Skip – skip for now

**Writes validation decisions to HawkEars_validation_results.csv with:**

clip_id, file name, time bounds, species, score

validator_id (your initials)

label (yes, no, unsure, skip)

validated_at timestamp

Smart behaviour

**On startup or when you click Start / Refresh queue, the app:**

Reads existing HawkEars_validation_results.csv

Removes any clips that have already been marked yes or no

Keeps unsure and skip clips in the pool for future re-validation

Optional cleanup: you can tell the app to delete the temporary 3-s snippets after each decision.

A running counter tracks how many clips you’ve labeled this session, and every 25 clips a modal pops up with a random boreal/grassland bird-themed compliment (because morale is a data quality issue).

**UI overview
**
Sidebar

Enter validator ID (initials)

Filter by target species

Toggle snippet deletion

View clip metadata and queue progress

Main panel

Predicted species + score

Spectrogram with moving playhead bar

Audio player (play/pause, scrub, volume)

Large, stacked buttons for Yes / No / Unsure, with Skip separate

**Requirements**

**R with the following packages:**

shiny, dplyr, readr, tuneR, av

A HawkEars-style CSV and matching WAV files on disk

**Paths to:**

labels_file (CSV)

audio_dir (WAVs)

results_file (output CSV)
are configured at the top of app.R.

**Typical workflow**

**Run shiny::runApp("path/to/app").**

Enter your initials and, optionally, choose a species filter.

Click Start / Refresh queue to load all clips not yet labeled yes or no.

**For each clip:**

Watch the spectrogram and play the audio.

Choose Yes, No, Unsure, or Skip.

Close the app whenever you’re done. Next time you start it, it will resume from where you left off, only serving unvalidated clips.
