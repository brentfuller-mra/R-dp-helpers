library(tibble)
library(dplyr)
library(purrr)
library(rlang)
library(glue)

# summarize all variables in a df
# particularly those imported with haven::read_sav

get_df_info <- function(df, copy_to_clipboard = TRUE, max_freq = 30, max_unique = 50) {
  
  if (!is.data.frame(df)) {
    stop("Input must be a data frame or tibble", call. = FALSE)
  }
  
  # helper functions
  get_label <- function(x) {
    attr(x, "label", exact = TRUE) %||% NA_character_
  }
  
  # get_punches <- function(x) {
  #   labs <- attr(x, "labels", exact = TRUE)
  #   if (is.null(labs) || length(labs) == 0) return(NA_character_)
  #   paste(paste(names(labs), labs, sep = " = "), collapse = "; ")
  # }
  
  # get punch labels as well
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
  
  # help with value label for haven labelled vectors
  get_value_label <- function(x, value) {
    labs <- attr(x, "labels", exact = TRUE)
    if (is.null(labs)) return(as.character(value))
    
    # Correct indexing for haven::labelled
    match_idx <- match(value, labs)
    if (is.na(match_idx)) return(as.character(value))
    names(labs)[match_idx]
  }
  
  # freq table
  get_freq_table <- function(x, max_freq = 30) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) == 0) return("")
    if (length(unique(x_clean)) > max_freq) return("(too many unique values)")
    
    tbl <- table(x_clean, useNA = "no")
    pct <- round(prop.table(tbl) * 100, 1)
    
    map2_chr(names(tbl), as.vector(tbl), function(val, count) {
      val <- as.character(val)                    # important
      label_text <- get_value_label(x, val)
      
      percent <- pct[[val]]
      glue("{label_text} ({val}) : {count} ({percent}%)")
    }) %>% 
      paste(collapse = "\n")
  }
  
  # main table stats
  var_info <- tibble(
    index        = 1:ncol(df),
    variable     = names(df),
    label        = map_chr(df, get_label),
    punches      = map_chr(df, get_punches),
    format       = map_chr(df, ~ attr(.x, "format.spss", exact = TRUE) %||% NA_character_),
    type         = map_chr(df, ~ paste(class(.x), collapse = "/")),
    n_valid      = map_int(df, ~ sum(!is.na(.x))),
    n_missing    = map_int(df, ~ sum(is.na(.x))),
    pct_missing  = round(map_dbl(df, ~ mean(is.na(.x)) * 100), 1),
    
    min          = map_dbl(df, ~ { if (is.numeric(.x)) min(.x, na.rm = TRUE) else NA_real_ }),
    max          = map_dbl(df, ~ { if (is.numeric(.x)) max(.x, na.rm = TRUE) else NA_real_ }),
    mean         = map_dbl(df, ~ { if (is.numeric(.x)) mean(.x, na.rm = TRUE) else NA_real_ }),
    sd           = map_dbl(df, ~ { if (is.numeric(.x)) sd(.x, na.rm = TRUE) else NA_real_ }),
    
    unique_count = map_int(df, ~ length(unique(.x[!is.na(.x)])))
  )
  
  # add frequency columns
  var_info <- var_info %>%
    mutate(
      freq_table = map_chr(df, ~ get_freq_table(.x, max_freq)),
      # pct_table  = map_chr(df, ~ {
      #   x_clean <- .x[!is.na(.x)]
      #   if (length(x_clean) == 0 || length(unique(x_clean)) > max_freq) return("")
      #   tbl <- table(x_clean)
      #   pct <- round(prop.table(tbl) * 100, 1)
      #   paste0(names(tbl), ": ", pct, "%", collapse = ", ")
      # })
    )
  
  # clipboard option
  if (copy_to_clipboard && requireNamespace("clipr", quietly = TRUE)) {
    tryCatch({
      clipr::write_clip(var_info, row.names = FALSE)
      cat("✓ Summary copied to clipboard.\n")
    }, error = function(e) {})
  }
  
  cat("Data frame has", nrow(df), "rows and", ncol(df), "columns\n\n")
  #print(var_info, n = Inf, width = Inf)
  
  return(var_info)
}