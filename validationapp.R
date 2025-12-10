# app.R  --- HawkEars Validator ---------------------------------------------
# written by M Edgar 12/10/2025

# You normally run installs once in the console, not in the app:
# install.packages(c("shiny", "dplyr", "readr", "tuneR", "av"),
#                  repos = "https://cloud.r-project.org/")

library(shiny)
library(dplyr)
library(readr)
library(tuneR)
library(av)      # for spectrograms (read_audio_fft)

## ==================== CONFIG: EDIT THESE PATHS ============================

# HawkEars output file
labels_file <- "C:/Users/EdgarM/Desktop/Localization/hawkears_tags/cmd_all/HawkEars_labels.csv"

# Directory containing the resampled WAV files referenced in `filename`
audio_dir   <- "E:/BAR-LT_LocalizationProject/localizationtrim_new"

# Where to store validation decisions
results_file <- "C:/Users/EdgarM/Desktop/Localization/hawkears_tags/cmd_all/HawkEars_validation_results.csv"

# Directory to hold 3-second snippets (per session, in tempdir)
snippet_dir <- file.path(tempdir(), "hawkears_snips")
if (!dir.exists(snippet_dir)) dir.create(snippet_dir, recursive = TRUE)

# Serve snippet_dir to the browser under /snips
addResourcePath("snips", snippet_dir)

## ==================== READ HAWKEARS OUTPUT ================================

hawk_raw <- read_csv(labels_file, show_col_types = FALSE) %>%
  mutate(
    clip_id    = row_number(),
    start_time = as.numeric(start_time),
    end_time   = as.numeric(end_time)
  )

## ==================== HELPER: CREATE SNIPPET ==============================

make_snippet <- function(clip_row) {
  wav_path <- file.path(audio_dir, clip_row$filename)

  if (!file.exists(wav_path)) {
    warning("File not found: ", wav_path)
    return(NA_character_)
  }

  out_name <- paste0("clip_", clip_row$clip_id, ".wav")
  out_path <- file.path(snippet_dir, out_name)

  # Read just the needed section (seconds)
  snd <- readWave(
    filename = wav_path,
    from     = clip_row$start_time,
    to       = clip_row$end_time,
    units    = "seconds"
  )

  writeWave(snd, out_path)
  out_name
}

## ==================== UI ==================================================

ui <- fluidPage(
  # Comic Sans everywhere
  tags$head(
    tags$style(HTML("
      body, label, input, button, .btn, .navbar, .form-control,
      .shiny-text-output, .shiny-plot-output, .panel {
        font-family: 'Comic Sans MS', 'Comic Sans', cursive;
      }
      h1, h2, h3, h4, h5, h6 {
        font-family: 'Comic Sans MS', 'Comic Sans', cursive;
      }
    "))
  ),

  titlePanel("HawkEars Validator 🦤 "),

  sidebarLayout(
    sidebarPanel(
      textInput("validator_id", "Your initials / ID:", value = ""),
      selectInput(
        "species_filter", "Target species:",
        choices  = c("All species", sort(unique(hawk_raw$class_name))),
        selected = "All species"
      ),
      checkboxInput(
        "delete_after",
        "Delete snippet file after validation (cleans temp dir only)",
        value = FALSE
      ),
      actionButton("start_btn", "Start / Refresh queue"),
      hr(),
      h4("Current clip info"),
      verbatimTextOutput("clip_info"),
      hr(),
      textOutput("progress_text")
    ),

    mainPanel(
      h3(textOutput("current_species")),
      br(),
      # Spectrogram with overlaid progress bar (ABOVE audio)
      tags$div(
        id    = "spec_container",
        style = paste(
          "position: relative; width: 100%; max-width: 700px;",
          "background-color:black; padding:5px;"
        ),
        plotOutput("spec", height = "280px"),
        # Track at bottom of spectrogram
        tags$div(
          id    = "spec_progress_bg",
          style = paste(
            "position:absolute; left:8%; right:8%; bottom:10px;",
            "height:6px; background-color:rgba(255,255,255,0.2);",
            "border-radius:3px; border:1px solid #666;"
          )
        ),
        # Filled part that we update from JS
        tags$div(
          id    = "spec_progress_fg",
          style = paste(
            "position:absolute; left:8%; bottom:10px; height:6px;",
            "width:0%; background-color:#00ff80;",
            "border-radius:3px;"
          )
        )
      ),
      br(),
      # Audio player UNDER the spectrogram
      uiOutput("audio_player"),
      br(),
      fluidRow(
        column(
          width = 4,
          # Stack Yes / No / Unsure vertically
          div(
            actionButton(
              "yes_btn", "Yes – prediction correct",
              class = "btn btn-success btn-lg btn-block"
            ),
            br(), br(),
            actionButton(
              "no_btn", "No – prediction wrong",
              class = "btn btn-danger btn-lg btn-block"
            ),
            br(), br(),
            actionButton(
              "unsure_btn", "Unsure",
              class = "btn btn-warning btn-lg btn-block"
            )
          )
        ),
        column(
          width = 3,
          actionButton(
            "skip_btn", "Skip",
            class = "btn btn-secondary btn-lg btn-block"
          )
        )
      )
    )
  )
)

## ==================== SERVER ==============================================

server <- function(input, output, session) {

  # Bird compliments to choose from
   compliments <- c(
    # Grassland-themed
    "You’re nailing these IDs like a Sprague’s Pipit nailing a display flight.",
    "Your validation skills are sharper than a Baird’s Sparrow call on a calm morning.",
    "You’re more reliable than a Western Meadowlark on fencepost duty.",
    "You’re sorting these clips faster than a flock of Longspurs flushing from the prairie.",
    "Your spectrogram reading is as clean as fresh snow on a hayfield.",
    "Grassland sparrows wish they were as consistent as your labels.",
    "You’re as dialed-in as a Bobolink over a lush alfalfa field.",
    "Your focus is stronger than a Savannah Sparrow’s loyalty to that one fence line.",
    "You’re more thorough than a Le Conte’s Sparrow survey at dawn.",
    "These decisions are crisper than a Chestnut-collared Longspur’s flight song.",
    "You’re picking out calls like a Marsh Wren picking perfect cattails.",
    "You’re as steady as an Upland Sandpiper on a gatepost in a windstorm.",
    "Your accuracy could rival a sharp-eared meadowlark on territory.",
    "You’re turning noisy recordings into clean data like a magician of the mixed-grass prairie.",
    "Your ears are doing more work than a whole crew of point counters in June.",
    "You’re catching faint calls like a Nelson’s Sparrow in good headphones.",
    "You’re as dependable as Clay-colored Sparrows in a shelterbelt.",
    "You’re transforming chaos into science like the hero of a Breeding Bird Survey route.",
    
    # Boreal-themed
    "You’re more dialed-in than a Boreal Chickadee on a cold January morning.",
    "Your pattern recognition is smoother than a Winter Wren song in a mossy ravine.",
    "You’re picking Veeries out of the noise like a seasoned boreal bander at dawn.",
    "You’re as reliable as a Swainson’s Thrush on a good spruce ridge.",
    "Your focus is deeper than a black spruce bog in June.",
    "You’re catching faint calls like a distant Olive-sided Flycatcher on a snag.",
    "You’re sorting these clips cleaner than a Yellow-rumped Warbler sorting spruce budworms.",
    "Your consistency rivals a Hermit Thrush singing from the same branch every evening.",
    "You’re more committed than a Palm Warbler tail-bobbing through a peatland.",
    "You’re reading spectrograms like a Canada Jay reads campsite opportunities.",
    "You’re detecting boreal birds like a supercharged ARU grid in peak season.",
    "Your decisions are more solid than permafrost in a cold spring.",
    "You’re picking out warbler chips like a Blackpoll Warbler picking midges from the canopy.",
    "You’re as sharp as a Northern Goshawk on a forest edge.",
    "You’re handling these detections like a crossbill handles a cone crop.",
    "You’re more dependable than loons on a quiet shield lake at sunrise.",
    "You’re teasing signal from noise like a Rustic Blackbird researcher’s dream dataset.",
    "Your validation stamina is stronger than a Gray-cheeked Thrush migration.",
    "You’re as precise as a Bay-breasted Warbler on a spruce budworm outbreak.",
    "You’re turning messy boreal soundscapes into clean, usable data like an absolute pro."
  )
  rv <- reactiveValues(
    queue            = NULL,   # data.frame of clips to validate
    idx              = 1,      # current index
    validated_count  = 0       # how many clips you've labeled this session
  )

  # Read existing validation results if present
  read_results <- function() {
    if (file.exists(results_file)) {
      read_csv(results_file, show_col_types = FALSE)
    } else {
      tibble()
    }
  }

  # Build queue of clips that still need definitive validation
  observeEvent(input$start_btn, {
    res <- read_results()

    if (nrow(res) > 0) {
      # Only yes/no are treated as "validated" and removed from queue
      validated_ids <- res %>%
        filter(label %in% c("yes", "no")) %>%
        pull(clip_id)
    } else {
      validated_ids <- integer(0)
    }

    dat <- hawk_raw %>%
      filter(!clip_id %in% validated_ids)

    if (!is.null(input$species_filter) &&
        input$species_filter != "All species") {
      dat <- dat %>% filter(class_name == input$species_filter)
    }

    if (nrow(dat) == 0) {
      showNotification("Wow, you're done already? No clips left to validate under current filters.",
                       type = "message")
    }

    rv$queue <- dat
    rv$idx   <- 1
  }, ignoreNULL = FALSE)  # run once at app launch

  # Current clip
  current_clip <- reactive({
    req(rv$queue)
    if (nrow(rv$queue) == 0) return(NULL)
    if (rv$idx > nrow(rv$queue)) return(NULL)
    rv$queue[rv$idx, ]
  })

  # Audio player for current snippet + JS to drive on-spec progress bar
  output$audio_player <- renderUI({
    clip <- current_clip()
    req(clip)

    snip_name <- make_snippet(clip)
    req(!is.na(snip_name))

    rel_path <- file.path("snips", snip_name)

    tagList(
      tags$audio(
        id       = "audio_clip",
        controls = NA,
        src      = rel_path,
        type     = "audio/wav",
        style    = "width:100%;"
      ),
      # JS: update spec_progress_fg width as audio plays
      tags$script(HTML("
        (function() {
          var audio = document.getElementById('audio_clip');
          var bar   = document.getElementById('spec_progress_fg');
          if (!audio || !bar) return;

          bar.style.width = '0%';

          audio.addEventListener('timeupdate', function() {
            if (audio.duration) {
              var frac = audio.currentTime / audio.duration;
              bar.style.width = (frac * 100) + '%';
            }
          });

          audio.addEventListener('ended', function() {
            bar.style.width = '0%';
          });
        })();
      "))
    )
  })

  # Spectrogram for current snippet using av::read_audio_fft
  output$spec <- renderPlot({
    clip <- current_clip()
    req(clip)

    snip_name <- paste0("clip_", clip$clip_id, ".wav")
    snip_path <- file.path(snippet_dir, snip_name)

    if (!file.exists(snip_path)) {
      created <- make_snippet(clip)
      if (is.na(created)) return(NULL)
      snip_path <- file.path(snippet_dir, created)
    }

    # Compute FFT-based spectrogram (mono, via ffmpeg)
    fft_data <- av::read_audio_fft(snip_path)

    # High-contrast spec
    op <- par(bg = "black", mar = c(4, 4, 1, 1))
    on.exit(par(op), add = TRUE)

    plot(fft_data, main = "")
  })

  # Metadata printout
  output$clip_info <- renderPrint({
    clip <- current_clip()
    if (is.null(clip)) {
      cat("No clip loaded. Click 'Start / Refresh queue'.")
      return()
    }

    list(
      clip_id    = clip$clip_id,
      filename   = clip$filename,
      start_time = clip$start_time,
      end_time   = clip$end_time,
      class_name = clip$class_name,
      class_code = clip$class_code,
      score      = clip$score
    )
  })

  # Species text
  output$current_species <- renderText({
    clip <- current_clip()
    if (is.null(clip)) return("No clip")
    paste0(
      "Predicted species: ", clip$class_name,
      " (", clip$class_code, "), score = ", clip$score
    )
  })

  # Progress summary text
  output$progress_text <- renderText({
    q <- rv$queue
    if (is.null(q) || nrow(q) == 0) return("0 clips in queue.")
    paste0(
      "Clip ", rv$idx, " of ", nrow(q),
      " (", round(100 * rv$idx / nrow(q), 1), "%)"
    )
  })

  # Save a decision & go to next clip
  save_and_next <- function(label) {
    clip <- current_clip()
    if (is.null(clip)) return(NULL)

    if (is.null(input$validator_id) || input$validator_id == "") {
      showNotification("Please enter your initials / ID first.",
                       type = "error")
      return(NULL)
    }

    new_row <- tibble(
      clip_id      = clip$clip_id,
      filename     = clip$filename,
      start_time   = clip$start_time,
      end_time     = clip$end_time,
      class_name   = clip$class_name,
      class_code   = clip$class_code,
      score        = clip$score,
      validator_id = input$validator_id,
      label        = label,          # yes / no / unsure / skip
      validated_at = Sys.time()
    )

    if (!file.exists(results_file)) {
      write_csv(new_row, results_file)
    } else {
      write_csv(new_row, results_file, append = TRUE)
    }

    # Increment in-session validation counter (any label counts)
    rv$validated_count <- rv$validated_count + 1

    # Every 25 clips: bird compliment pop-up
    if (rv$validated_count %% 25 == 0) {
      showModal(
        modalDialog(
          title = "Nice work! 🥳 ",
          sample(compliments, 1),
          easyClose = TRUE,
          footer = modalButton("Back to the birds")
        )
      )
    }

    # Optionally delete snippet file from temp dir
    if (isTRUE(input$delete_after)) {
      snip_name <- paste0("clip_", clip$clip_id, ".wav")
      snip_path <- file.path(snippet_dir, snip_name)
      if (file.exists(snip_path)) unlink(snip_path)
    }

    # Advance index
    if (!is.null(rv$queue) && rv$idx < nrow(rv$queue)) {
      rv$idx <- rv$idx + 1
    } else {
      showNotification("Reached end of queue.", type = "message")
    }
  }

  observeEvent(input$yes_btn,    save_and_next("yes"))
  observeEvent(input$no_btn,     save_and_next("no"))
  observeEvent(input$unsure_btn, save_and_next("unsure"))
  observeEvent(input$skip_btn,   save_and_next("skip"))
}

shinyApp(ui, server)
