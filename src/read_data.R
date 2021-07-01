library(rjson)
require(tcltk)
library(tidyverse)
library(patchwork)
library(knitr)

print_line <- function(data_line, line_index) {
  if (nchar(data_line) == 0) {
    message(paste("Line", line_index, "empty"))
  } else if (nchar(data_line) < 60) {
    message(paste("Line", line_index, "contents:", data_line))
  } else {
    message(paste(
      "Line",
      line_index,
      "abbrev:",
      substring(data_line, 1, 25),
      "....",
      substring(data_line, nchar(data_line) - 25 + 1)
    ))
  }
}

parse_data_line <- function(data_line, line_index) {
  out <- tryCatch(
    {
      fromJSON(data_line)
    },
    error = function(cond) {
      print_line(data_line, line_index)
    },
    warning = function(cond) {
      print_line(data_line, line_index)
    }
  )
  return(out)
}


#' # Get data from a raw JSON file
#'
#' Prompts the user to select a file using a tk dialog, which should be a file
#' with each line representing a JSON object from the PrefSFT2020 code
#'
#' **returns** A list containing each participants data converted from JSON
get_data <- function(datafile) {
  connection <- file(datafile, open = "r")
  data_lines <- readLines(connection)
  raw_data <- c()
  for (line_index in seq_along(data_lines)) {
    line_data <- parse_data_line(data_lines[line_index], line_index)
    if (!is.null(line_data)) {
      raw_data <- append(raw_data, list(line_data))
    }
  }
  close(connection)
  raw_data
}


#' # Extract sampled values from raw data
#'
#' For a particular one line raw data list, extract the sampled price/rating
#' values and place into an appropriately formatted dataframe
#'
#' **param** raw_data: The data object returned from a call to getdata
#'
#' **returns** A data.frame with seeds, values, and price/rating attributes
getdf <- function(raw_data) {
  messages <- lapply(raw_data, "[[", "messages")
  jsonstr <- messages[[1]][which.max(unlist(lapply(messages[[1]], nchar)))]
  nested_lists <- fromJSON(jsonstr)
  seed_names <- paste0("seed_", seq(0, 1, .1))
  arr_psft <- array(NA_real_, dim = c(11, 900, 2))
  for (seed_idx in 1:11) {
    arr_psft[seed_idx, , ] <- t(array(unlist(nested_lists$psft[[seed_idx]]), dim = c(2, 900)))
  }
  dimnames(arr_psft) <- list("seeds" = seed_names, "values" = paste0("val_", 1:900), "attributes" = c("price", "rating"))
  psft_3d_cube <- as.tbl_cube(arr_psft)
  df_psft <- as_tibble(psft_3d_cube)
  df_psft$attributes <- factor(df_psft$attributes)
  df_psft$seeds <- factor(df_psft$seeds)
  df_psft
}


#' # Run the whole sampled values plotting checks
#'
#' Run the getdata/getdf functions to grab two datasets - and plot histograms
#' of price and rating in order to show the sampling process for each salience
#' level working
check_js_samples <- function() {
  seedvals <- seq(0, 1, 0.2)
  pthresh <- 100 + seedvals * 350
  rthresh <- 9.73 + -319.06 * (1 / pthresh)
  seeds <- paste0("seed_", seedvals)
  dfmu <- data.frame(seedvals, pthresh, rthresh, seeds)

  dfold <- getdf(getdata())
  dfold <- dfold[dfold$seeds %in% dfmu$seeds, ]
  dfold$seeds <- factor(dfold$seeds)
  dfnew <- getdf(getdata())
  dfnew <- dfnew[dfnew$seeds %in% dfmu$seeds, ]
  dfnew$seeds <- factor(dfnew$seeds)

  pdf("Price_Rating_Histograms.pdf")

  price_psft_old <- ggplot(data = dfold[dfold$attributes == "price", ], mapping = aes(x = arr_psft)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(. ~ seeds, scales = "free") +
    ggtitle("Histogram of Prices for each seed value - OldMethod") +
    xlab("Price") +
    geom_vline(data = dfmu, mapping = aes(xintercept = pthresh), color = "blue", linetype = "dashed")
  print(price_psft_old)

  price_psft_new <- ggplot(data = dfnew[dfnew$attributes == "price", ], mapping = aes(x = arr_psft)) +
    geom_histogram(binwidth = 1) +
    facet_wrap(. ~ seeds, scales = "free") +
    ggtitle("Histogram of Prices for each seed value - NewMethod") +
    xlab("Price") +
    geom_vline(data = dfmu, mapping = aes(xintercept = pthresh), color = "blue", linetype = "dashed")
  print(price_psft_new)

  rating_psft_old <- ggplot(data = dfold[dfold$attributes == "rating", ], mapping = aes(x = arr_psft)) +
    geom_histogram(binwidth = .1) +
    facet_wrap(. ~ seeds, scales = "free") +
    ggtitle("Histogram of Ratings for each seed value - Standard") +
    xlab("Rating") +
    geom_vline(data = dfmu, mapping = aes(xintercept = rthresh), color = "blue", linetype = "dashed")
  print(rating_psft_old)

  rating_psft_new <- ggplot(data = dfnew[dfnew$attributes == "rating", ], mapping = aes(x = arr_psft)) +
    geom_histogram(binwidth = .1) +
    facet_wrap(. ~ seeds, scales = "free") +
    ggtitle("Histogram of Ratings for each seed value - Standard") +
    xlab("Rating") +
    geom_vline(data = dfmu, mapping = aes(xintercept = rthresh), color = "blue", linetype = "dashed")
  print(rating_psft_new)

  dev.off()
}


extract_data <- function(raw_data) {
  part_data_list <- lapply(raw_data, FUN = function(data_line) {
    df <- read_csv(data_line$jsPsychData, col_types = cols(.default = "c"))
  })
  part_data_list <- bind_rows(lapply(
    part_data_list,
    function(dtt) {
      mutate_all(dtt, as.character)
    }
  ))
  bind_rows(part_data_list)
}

augment_data <- function(events) {
  # numeric columns
  num_cols <- c("trial_index", "time_elapsed", "price", "rating", "rt",
                "price_threshold", "rating_threshold", "key_press",
                "staircase_point")
  fac_cols <- c("type", "cell_name", "sonaID", "group", "hand", "attr_order",
                "staircase_result")
  log_cols <- c("accept")
  # Add computed columns
  min_price <- 100
  price_rng <- 450 - min_price
  min_rating <- 6.5
  rating_rng <- 9 - min_rating
  events %>%
    mutate_at(num_cols, as.numeric) %>%
    mutate_at(fac_cols, as.factor) %>%
    mutate_at(log_cols, as.logical) %>%
    rename(price_threshold_rescaled = price_threshold) %>%
    rename(rating_threshold_rescaled = rating_threshold) %>%
    mutate(price_thresh = min_price + price_threshold_rescaled * price_rng) %>%
    mutate(rating_thresh = min_rating + rating_threshold_rescaled * rating_rng)
}

filter_data <- function(events) {
  # Remove NA values in price threshold columns for slider responses
  complete <- events %>%
    filter(
      case_when(
        type == "slider" ~ !is.na(price_thresh),
        TRUE ~ TRUE
        )
      )

  common_columns <- c(
    "rt", "trial_index", "time_elapsed", "attr_order",
    "group", "hand", "sonaID", "price_threshold_rescaled",
    "rating_threshold_rescaled", "price_thresh", "rating_thresh"
  )
  pref_columns <- c("key_press", "accept", "price", "rating")
  sft_columns <- c("cell_name", pref_columns, common_columns)

  implicit_columns <- c("staircase_point", "staircase_result", pref_columns,
                         common_columns)
  explicit_columns <- c("response", common_columns)

  implicit <- complete %>%
    filter(type == "staircase") %>%
    select(all_of(implicit_columns)) %>%
    mutate(rt = as.numeric(rt) / 1000) %>%
    mutate(sonaID = as.factor(sonaID))

  explicit <- complete %>%
    filter(type == "slider") %>%
    select(all_of(explicit_columns)) %>%
    mutate(rt = as.numeric(rt) / 1000) %>%
    mutate(sonaID = as.factor(sonaID))

  sft <- complete %>%
    filter(type == "sft") %>%
    select(all_of(sft_columns)) %>%
    mutate(rt = as.numeric(rt) / 1000) %>%
    mutate(sonaID = as.factor(sonaID))

  list(
    complete = complete,
    explicit = explicit,
    implicit = implicit,
    sft = sft
  )
}

get_survey <- function(all_events) {
  # Messed up survey - no sonaID stored - so grab from previous trial
  survey_rows <- which(!is.na(all_events$responses))
  sona_ids <- all_events[survey_rows - 1, "sonaID"]
  survey_times <- all_events[
    survey_rows,
    c("rt", "time_elapsed")
  ]
  survey_responses <- pull(all_events[survey_rows, "responses"], responses)
  survey <- bind_rows(lapply(survey_responses, FUN = function(response) {
    as_tibble(fromJSON(response))
  }))
  bind_cols(sona_ids, survey, survey_times)
}



sft_rt_plot <- function(sft) {
  # Add mean lines
  sft_mu <- sft %>%
    group_by(sonaID) %>%
    summarize(grp.mean = mean(rt))

  sft_hist <- ggplot(sft, aes(x = rt, color = sonaID)) +
    geom_histogram(alpha = 0.2, fill = "white", position = "dodge", binwidth = 1) +
    geom_vline(
      data = sft_mu, aes(xintercept = grp.mean, color = sonaID),
      linetype = "dashed"
    ) +
    ggtitle("SFT trials") +
    theme(legend.position = "top")

  sft_elap <- ggplot(sft, aes(x = trial_index, y = time_elapsed / 60000, color = sonaID)) +
    geom_line() +
    ggtitle("Elapsed time") +
    ylab("Time in minutes") +
    theme(legend.position = "top")

  png("SFTPlots.png", width = 800)
  sft_hist / sft_elap
  print(sft_mu)
}

