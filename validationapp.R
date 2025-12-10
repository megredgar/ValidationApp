# app.R ---------------------------------------------------------
# check where your libraries are
# written by M Edgar 12/10/2025 

install.packages(
  "promises",
  type = "source",
  repos = "https://cran.rstudio.com"
)

packageVersion("promises")
install.packages("shiny")
install.packages("seewave")
library(shiny)
library(dplyr)
library(readr)
library(tuneR)
library(seewave)

## ==== CONFIG: EDIT THESE IF PATHS CHANGE ======================

# HawkEars label file
labels_file <- "C:/Users/EdgarM/Desktop/Localization/hawkears_tags/cmd_all/HawkEars_labels.csv"

# Directory holding the resampled wavs referenced in `filename`
audio_dir   <- "E:/BAR-LT_LocalizationProject/localizationtrim_new"

# Where to store validation decisions
results_file <- "C:/Users/EdgarM/Desktop/Localization/hawkears_tags/cmd_all/HawkEars_validation_results.csv"

# Directory where we’ll write 3-second snippets for playback
snippet_dir <- file.path(tempdir(), "hawkears_snips")
dir.create(snippet_dir, showWarnings = FALSE)

# Serve snippet_dir under URL path /snips
addResourcePath("snips", snippet_dir)

## ==== READ HAWKEARS OUTPUT ===================================

hawk_raw <- read_csv(labels_file, show_col_types = FALSE) %>%
  mutate(
    clip_id    = row_number(),
    start_time = as.numeric(start_time),
    end_time   = as.numeric(end_time)
  )

## ==== HELPER: MAKE SNIPPET FOR CURRENT CLIP ===================

make_snippet <- function(clip_row) {
  wav_path <- file.path(audio_dir, clip_row$filename)

  if (!file.exists(wav_path)) {
    warning("File not found: ", wav_path)
    return(NA_character_)
  }

  # Output file name for this clip
  out_name <- paste0("clip_", clip_row$clip_id, ".wav")
  out_path <- file.path(snippet_dir, out_name)

  # Read only the needed section (in seconds)
  # If your files are stereo / 48kHz etc., tuneR handles it.
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
  titlePanel("HawkEars Validator"),

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
      fluidRow(
        column(
          width = 6,
          uiOutput("audio_player")
        ),
        column(
          width = 6,
          plotOutput("spec", height = "250px")
        )
      ),
      br(),
      fluidRow(
        column(
          width = 3,
          actionButton("yes_btn", "Yes – prediction correct",
                       class = "btn btn-success btn-lg")
        ),
        column(
          width = 3,
          actionButton("no_btn", "No – prediction wrong",
                       class = "btn btn-danger btn-lg")
        ),
        column(
          width = 3,
          actionButton("unsure_btn", "Unsure",
                       class = "btn btn-warning btn-lg")
        ),
        column(
          width = 3,
          actionButton("skip_btn", "Skip",
                       class = "btn btn-secondary btn-lg")
        )
      )
    )
  )
)

## ==================== SERVER ==============================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    queue = NULL,   # data.frame of clips to validate
    idx   = 1       # current index
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
      # Only yes/no are treated as "validated"
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
      showNotification("No clips left to validate under current filters.",
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

  # Audio player for current snippet
  output$audio_player <- renderUI({
    clip <- current_clip()
    req(clip)

    snip_name <- make_snippet(clip)
    req(!is.na(snip_name))

    rel_path <- file.path("snips", snip_name)

    tags$audio(
      controls = NA,
      src      = rel_path,
      type     = "audio/wav"
    )
  })

  # Spectrogram for current snippet
  output$spec <- renderPlot({
    clip <- current_clip()
    req(clip)

    snip_name <- paste0("clip_", clip$clip_id, ".wav")
    snip_path <- file.path(snippet_dir, snip_name)

    # Make snippet if missing (first time)
    if (!file.exists(snip_path)) {
      created <- make_snippet(clip)
      if (is.na(created)) return(NULL)
      snip_path <- file.path(snippet_dir, created)
    }

    w <- readWave(snip_path)

    seewave::spectro(
      w,
      flim  = c(0, 10),  # kHz; tweak if needed
      wl    = 512,
      ovlp  = 75,
      scale = FALSE,
      main  = ""
    )
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

  # Progress
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