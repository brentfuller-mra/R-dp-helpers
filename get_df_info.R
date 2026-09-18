get_df_info <- function(
    df,
    copy_to_clipboard = TRUE,
    max_freq = 30,
    max_unique = 50,
    output = c("none", "clipboard", "csv", "html"),
    output_path = NULL
) {
  output <- match.arg(output)
  
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or tibble", call. = FALSE)
  }
  
  # helper functions --------------------------------------------------------
  
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  get_label <- function(x) {
    attr(x, "label", exact = TRUE) %||% NA_character_
  }
  
  get_punches <- function(x, max_unique = 50) {
    labs <- attr(x, "labels", exact = TRUE)
    
    if (!is.null(labs) && length(labs) > 0) {
      # Has value labels (typical for SPSS)
      paste(paste(names(labs), labs, sep = " = "), collapse = "; ")
    } else {
      # No value labels → show unique values
      uniques <- unique(x[!is.na(x)])
      if (length(uniques) == 0) return(NA_character_)
      if (length(uniques) > max_unique) {
        return(paste0("Too many unique values (", length(uniques), ")"))
      }
      uniques <- sort(uniques)
      paste(uniques, collapse = "; ")
    }
  }
  
  get_value_label <- function(x, value) {
    labs <- attr(x, "labels", exact = TRUE)
    if (is.null(labs)) return(as.character(value))
    
    # value may be character; labs are usually named numeric
    # try matching on the underlying values
    value_num <- suppressWarnings(as.numeric(value))
    if (!is.na(value_num)) {
      idx <- match(value_num, labs)
    } else {
      # fallback: match on names (for character-coded labels)
      idx <- match(value, names(labs))
    }
    
    if (is.na(idx)) return(as.character(value))
    names(labs)[idx]
  }
  
  get_freq_table <- function(x, max_freq = 30) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) == 0) return("")
    if (length(unique(x_clean)) > max_freq) return("(too many unique values)")
    
    tbl <- table(x_clean, useNA = "no")
    pct <- round(prop.table(tbl) * 100, 1)
    
    # robust mapping: match by position, not direct [[val]] indexing
    out <- purrr::map2_chr(
      names(tbl),
      as.vector(tbl),
      function(val, count) {
        val_chr <- as.character(val)
        label_text <- get_value_label(x, val_chr)
        
        # match val to names(tbl) and use that index for pct
        idx <- match(as.character(val), as.character(names(tbl)))
        if (is.na(idx)) {
          percent <- NA_real_
        } else {
          percent <- pct[idx]
        }
        
        glue::glue("{label_text} ({val_chr}) : {count} ({percent}%)")
      }
    )
    
    paste(out, collapse = "\n")
  }
  
  # main table stats --------------------------------------------------------
  
  var_info <- tibble::tibble(
    index        = 1:ncol(df),
    variable     = names(df),
    label        = purrr::map_chr(df, get_label),
    punches      = purrr::map_chr(df, get_punches, max_unique = max_unique),
    format       = purrr::map_chr(df, ~ attr(.x, "format.spss", exact = TRUE) %||% NA_character_),
    type         = purrr::map_chr(df, ~ paste(class(.x), collapse = "/")),
    n_valid      = purrr::map_int(df, ~ sum(!is.na(.x))),
    n_missing    = purrr::map_int(df, ~ sum(is.na(.x))),
    pct_missing  = round(purrr::map_dbl(df, ~ mean(is.na(.x)) * 100), 1),
    
    min          = purrr::map_dbl(df, ~ { if (is.numeric(.x)) min(.x, na.rm = TRUE) else NA_real_ }),
    max          = purrr::map_dbl(df, ~ { if (is.numeric(.x)) max(.x, na.rm = TRUE) else NA_real_ }),
    mean         = purrr::map_dbl(df, ~ { if (is.numeric(.x)) mean(.x, na.rm = TRUE) else NA_real_ }),
    sd           = purrr::map_dbl(df, ~ { if (is.numeric(.x)) stats::sd(.x, na.rm = TRUE) else NA_real_ }),
    
    unique_count = purrr::map_int(df, ~ length(unique(.x[!is.na(.x)])))
  )
  
  var_info <- var_info %>%
    dplyr::mutate(
      freq_table = purrr::map_chr(df, ~ get_freq_table(.x, max_freq))
    )
  
  # output handling ---------------------------------------------------------
  
  message("Data frame has ", nrow(df), " rows and ", ncol(df), " columns\n")
  
  # clipboard (legacy flag + new output arg)
  if ((copy_to_clipboard || output == "clipboard") &&
      requireNamespace("clipr", quietly = TRUE)) {
    tryCatch({
      clipr::write_clip(var_info, row.names = FALSE)
      cat("✓ Summary copied to clipboard.\n")
    }, error = function(e) {})
  }
  
  if (!is.null(output_path)) {
    if (output == "csv") {
      readr::write_csv(var_info, output_path)
      cat("✓ Summary written to CSV: ", output_path, "\n")
    } else if (output == "html") {
      if (!requireNamespace("htmltools", quietly = TRUE) ||
          !requireNamespace("knitr", quietly = TRUE)) {
        warning("htmltools and knitr are required for HTML output.")
      } else {
        tbl_html <- knitr::kable(var_info, format = "html") %>%
          htmltools::HTML()
        
        page <- htmltools::tagList(
          htmltools::tags$html(
            htmltools::tags$head(
              htmltools::tags$title("Data frame summary")
            ),
            htmltools::tags$body(
              htmltools::tags$h1("Data frame summary"),
              tbl_html
            )
          )
        )
        
        htmltools::save_html(page, file = output_path)
        cat("✓ HTML report written to: ", output_path, "\n")
      }
    }
  }
  
  return(var_info)
}
