# app.R  --- HawkEars Validator (Tseng-style bins) --------------------------
# written by M Edgar, updated with bin workflow + bird compliments

library(shiny)
library(dplyr)
library(readr)
library(tuneR)
library(av)

## ==================== CONFIG: EDIT THESE PATHS ============================

labels_file  <- "E:/BAR-LT_LocalizationProject/localization_05312025/hawkears_lowthresh/HawkEars_labels.csv"
audio_dir    <- "E:/BAR-LT_LocalizationProject/localization_05312025/localizationtrim_new"
results_file <- "E:/BAR-LT_LocalizationProject/localization_05312025/hawkears_lowthresh/HawkEars_validation_results.csv"

snippet_dir <- file.path(tempdir(), "hawkears_snips")
if (!dir.exists(snippet_dir)) dir.create(snippet_dir, recursive = TRUE)
addResourcePath("snips", snippet_dir)  # exposes snippet_dir at /snips

## ==================== TSENG-LOCKED SETTINGS ===============================

TSENG_MIN  <- 0.10
TSENG_MAX  <- 1.00
TSENG_STEP <- 0.05
TSENG_N    <- 10

# Fixed bin labels: 0.10–0.15, ..., 0.95–1.00 (18 bins)
tseng_lowers <- seq(TSENG_MIN, TSENG_MAX - TSENG_STEP, by = TSENG_STEP)
tseng_uppers <- tseng_lowers + TSENG_STEP
tseng_bins   <- paste0(sprintf("%.2f", tseng_lowers), "–", sprintf("%.2f", tseng_uppers))

# Breaks for cut(): left-closed [a,b); add tiny epsilon so 1.00 lands in last bin
tseng_breaks <- c(tseng_lowers, TSENG_MAX + 1e-9)

## ==================== READ HAWKEARS OUTPUT ================================

hawk_raw <- read_csv(labels_file, show_col_types = FALSE) %>%
  mutate(
    clip_id    = row_number(),
    start_time = as.numeric(start_time),
    end_time   = as.numeric(end_time),
    score      = as.numeric(score)
  )

## ==================== HELPERS =============================================

read_results <- function() {
  if (file.exists(results_file)) {
    read_csv(results_file, show_col_types = FALSE)
  } else {
    tibble()
  }
}

add_tseng_bins <- function(dat) {
  dat %>%
    mutate(score = as.numeric(score)) %>%
    filter(!is.na(score)) %>%
    filter(score >= TSENG_MIN, score <= TSENG_MAX) %>%
    mutate(
      score_bin = cut(
        score,
        breaks = tseng_breaks,
        right = FALSE,
        include.lowest = TRUE,
        labels = tseng_bins
      ),
      score_bin = as.character(score_bin),
      bin_lower = as.numeric(sub("^([0-9\\.]+)–.*$", "\\1", score_bin)),
      bin_upper = as.numeric(sub("^.*–([0-9\\.]+)$", "\\1", score_bin))
    )
}

make_snippet <- function(clip_row) {
  wav_path <- normalizePath(
    file.path(audio_dir, clip_row$filename),
    winslash = "/",
    mustWork = FALSE
  )
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
  
  titlePanel("HawkEars Validator 🦤"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("validator_id", "Your initials / ID:", value = ""),
      
      selectInput(
        "species_filter", "Target species:",
        choices  = sort(unique(hawk_raw$class_name))
      ),
      hr(),
      
      h4("Validate per confidence bin"),
      selectInput("bin_filter", "Bin:", choices = tseng_bins, selected = tseng_bins[1]),
      actionButton("next_bin_btn", "Next available bin ▶"),
      br(), br(),
      actionButton("start_btn", "Load this bin"),
      hr(),
      
      tags$div(
        tags$b("Locked settings: "),
        paste0("min=", TSENG_MIN, ", max=", TSENG_MAX,
               ", step=", TSENG_STEP, ", n/bin=", TSENG_N)
      ),
      numericInput("rng_seed", "Random seed (reproducible)", value = 1, min = 1, step = 1),
      
      hr(),
      h4("Bin availability"),
      verbatimTextOutput("bin_summary"),
      
      hr(),
      checkboxInput("delete_after", "Delete snippet file after validation (temp only)", value = FALSE),
      
      hr(),
      h4("Current clip info"),
      verbatimTextOutput("clip_info"),
      hr(),
      textOutput("progress_text")
    ),
    
    mainPanel(
      h3(textOutput("current_species")),
      br(),
      
      plotOutput("spec", height = "280px"),
      br(),
      
      uiOutput("audio_player"),
      br(),
      
      fluidRow(
        column(
          width = 4,
          actionButton("yes_btn", "Yes – prediction correct",
                       class = "btn btn-success btn-lg btn-block"),
          br(), br(),
          actionButton("no_btn", "No – prediction wrong",
                       class = "btn btn-danger btn-lg btn-block"),
          br(), br(),
          actionButton("unsure_btn", "Unsure",
                       class = "btn btn-warning btn-lg btn-block")
        ),
        column(
          width = 3,
          actionButton("skip_btn", "Skip",
                       class = "btn btn-secondary btn-lg btn-block")
        )
      )
    )
  )
)

## ==================== SERVER ==============================================

server <- function(input, output, session) {
  
  # Bird compliments to choose from (grassland + boreal)
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
    queue           = tibble(),
    idx             = 1,
    validated_count = 0,
    compliment_step = 0   # cycles through compliments
  )
  
  # Pool of available clips for current species, excluding already yes/no
  available_pool <- reactive({
    req(input$species_filter)
    
    res <- read_results()
    done_ids <- if (nrow(res) > 0) {
      res %>% filter(label %in% c("yes", "no")) %>% pull(clip_id)
    } else {
      integer(0)
    }
    
    hawk_raw %>%
      filter(class_name == input$species_filter) %>%
      filter(!clip_id %in% done_ids) %>%
      add_tseng_bins()
  })
  
  output$bin_summary <- renderPrint({
    dat <- available_pool()
    if (nrow(dat) == 0) {
      cat("No clips remaining for this species within Tseng bins.\n")
      return()
    }
    
    summ <- dat %>%
      count(score_bin, name = "available") %>%
      right_join(tibble(score_bin = tseng_bins), by = "score_bin") %>%
      mutate(available = ifelse(is.na(available), 0L, available)) %>%
      arrange(factor(score_bin, levels = tseng_bins))
    
    sel <- input$bin_filter
    cat("Selected bin:", sel, "\n")
    print(summ %>% filter(score_bin == sel), n = 50)
    
    cat("\nAll bins (available counts):\n")
    print(summ, n = 50)
  })
  
  observeEvent(input$next_bin_btn, {
    dat <- available_pool()
    if (nrow(dat) == 0) return()
    
    counts <- dat %>%
      count(score_bin, name = "available") %>%
      right_join(tibble(score_bin = tseng_bins), by = "score_bin") %>%
      mutate(available = ifelse(is.na(available), 0L, available)) %>%
      arrange(factor(score_bin, levels = tseng_bins))
    
    nonempty <- counts %>% filter(available > 0) %>% pull(score_bin)
    if (length(nonempty) == 0) {
      showNotification("No bins with clips remaining for this species.", type = "message")
      return()
    }
    
    cur <- input$bin_filter
    pos <- match(cur, nonempty)
    next_bin <- if (is.na(pos) || pos == length(nonempty)) nonempty[1] else nonempty[pos + 1]
    updateSelectInput(session, "bin_filter", selected = next_bin)
  })
  
  observeEvent(input$start_btn, {
    set.seed(input$rng_seed)
    
    dat <- available_pool() %>% filter(score_bin == input$bin_filter)
    
    if (nrow(dat) == 0) {
      rv$queue <- tibble()
      rv$idx <- 1
      showNotification("No clips remaining in that bin for this species.", type = "message")
      return()
    }
    
    dat <- dat %>%
      slice_sample(n = min(TSENG_N, nrow(dat))) %>%
      slice_sample(prop = 1)
    
    rv$queue <- dat
    rv$idx   <- 1
  }, ignoreNULL = FALSE)
  
  current_clip <- reactive({
    if (nrow(rv$queue) == 0) return(NULL)
    if (rv$idx > nrow(rv$queue)) return(NULL)
    rv$queue[rv$idx, ]
  })
  
  output$audio_player <- renderUI({
    clip <- current_clip()
    req(clip)
    
    snip_name <- make_snippet(clip)
    if (is.na(snip_name)) {
      return(tags$div("Audio snippet could not be created (missing file/path)."))
    }
    
    # Leading "/" + cache-buster so browser reloads per clip
    src <- paste0("/snips/", snip_name, "?v=", clip$clip_id)
    
    tags$audio(
      id       = "audio_clip",
      controls = NA,
      src      = src,
      type     = "audio/wav",
      style    = "width:100%;"
    )
  })
  
  output$spec <- renderPlot({
    clip <- current_clip()
    req(clip)
    
    snip_path <- file.path(snippet_dir, paste0("clip_", clip$clip_id, ".wav"))
    if (!file.exists(snip_path)) {
      created <- make_snippet(clip)
      if (is.na(created)) return(NULL)
      snip_path <- file.path(snippet_dir, created)
    }
    
    fft_data <- av::read_audio_fft(snip_path)
    op <- par(bg = "black", mar = c(4, 4, 1, 1))
    on.exit(par(op), add = TRUE)
    plot(fft_data, main = "")
  })
  
  output$clip_info <- renderPrint({
    clip <- current_clip()
    if (is.null(clip)) {
      cat("No clip loaded. Pick a bin and click 'Load this bin'.")
      return()
    }
    
    list(
      clip_id    = clip$clip_id,
      filename   = clip$filename,
      start_time = clip$start_time,
      end_time   = clip$end_time,
      class_name = clip$class_name,
      score      = clip$score,
      score_bin  = clip$score_bin,
      bin_lower  = clip$bin_lower,
      bin_upper  = clip$bin_upper
    )
  })
  
  output$current_species <- renderText({
    clip <- current_clip()
    if (is.null(clip)) {
      return(paste0("No clip | Species: ", input$species_filter,
                    " | Bin: ", input$bin_filter))
    }
    paste0(
      "Predicted species: ", clip$class_name,
      ", score = ", clip$score,
      " | bin = ", clip$score_bin
    )
  })
  
  output$progress_text <- renderText({
    if (nrow(rv$queue) == 0) return("0 clips loaded. Click 'Load this bin'.")
    paste0("Clip ", rv$idx, " of ", nrow(rv$queue))
  })
  
  save_and_next <- function(label) {
    clip <- current_clip()
    if (is.null(clip)) return()
    
    if (is.null(input$validator_id) || input$validator_id == "") {
      showNotification("Please enter your initials / ID first.", type = "error")
      return()
    }
    
    new_row <- tibble(
      clip_id      = clip$clip_id,
      filename     = clip$filename,
      start_time   = clip$start_time,
      end_time     = clip$end_time,
      class_name   = clip$class_name,
      class_code   = clip$class_code,
      score        = clip$score,
      score_bin    = clip$score_bin,
      bin_lower    = clip$bin_lower,
      bin_upper    = clip$bin_upper,
      tseng_min    = TSENG_MIN,
      tseng_max    = TSENG_MAX,
      tseng_step   = TSENG_STEP,
      tseng_n      = TSENG_N,
      rng_seed     = as.integer(input$rng_seed),
      validator_id = input$validator_id,
      label        = label,
      validated_at = Sys.time()
    )
    
    if (!file.exists(results_file)) {
      write_csv(new_row, results_file)
    } else {
      write_csv(new_row, results_file, append = TRUE)
    }
    
    # Increment in-session validation counter
    rv$validated_count <- rv$validated_count + 1
    
    # Every 25 clips: bird compliment pop-up, cycling through compliments
    if (rv$validated_count %% 25 == 0) {
      rv$compliment_step <- rv$compliment_step + 1
      idx <- ((rv$compliment_step - 1) %% length(compliments)) + 1
      
      showModal(
        modalDialog(
          title = "Nice work! 🥳",
          compliments[idx],
          easyClose = TRUE,
          footer = modalButton("Back to the birds")
        )
      )
    }
    
    # Optionally delete snippet file from temp dir
    if (isTRUE(input$delete_after)) {
      snip_path <- file.path(snippet_dir, paste0("clip_", clip$clip_id, ".wav"))
      if (file.exists(snip_path)) unlink(snip_path)
    }
    
    # Advance index within current bin queue
    if (rv$idx < nrow(rv$queue)) {
      rv$idx <- rv$idx + 1
    } else {
      showNotification(
        "End of queue for this bin. Click 'Load this bin' again or go to Next available bin.",
        type = "message"
      )
    }
  }
  
  observeEvent(input$yes_btn,    save_and_next("yes"))
  observeEvent(input$no_btn,     save_and_next("no"))
  observeEvent(input$unsure_btn, save_and_next("unsure"))
  observeEvent(input$skip_btn,   save_and_next("skip"))
}

shinyApp(ui, server)
