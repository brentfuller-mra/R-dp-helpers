# =============================================================================
# R Data Processing & Stats Helpers for Survey / Market Research
# For documentation purposes only 
# Be careful - these functions are OLD and can be fragile
# There are more recent and robust solutions/packages to each of these
# for example: tidyverse, srvyr, gtsummary, gt, flextable, pollster,
# haven, labelled, rwa, relaimpo, emmeans, marginaleffects, collapse, 
# data.table, quarto, gt
# =============================================================================

library(tibble)
library(dplyr)
library(purrr)
library(rlang)
library(glue)
library(ggplot2)

#' Statistical mode(s)
#' @export
get_mode <- function(v, all = FALSE) {
  if (all(is.na(v))) return(NA)
  uniqv <- unique(v[!is.na(v)])
  tabv  <- tabulate(match(v, uniqv))
  modes <- uniqv[tabv == max(tabv)]
  if (all) sort(modes, decreasing = TRUE) else modes[1]
}

#' Get value labels from haven labelled variable
#' @export
spss_value_labels <- function(x) {
  labs <- attr(x, "labels", exact = TRUE)
  if (is.null(labs)) return(NULL)
  data.frame(value = unname(labs), label = names(labs), stringsAsFactors = FALSE)
}

#' Group descriptives with weighting option
#'
#' @param data Data frame containing variables to summarize
#' @param group Grouping variable (vector or column name)
#' @param weight Weighting variable (default = unweighted)
#' @param stat Statistic to compute: "mean", "sum", "median", "var", "min", "max"
#' @param total Include total row?
#' @param var_counts Include non-missing counts row?
#' @param f_stats Include ANOVA F-stat and p-value?
#' @export
group_means <- function(data,
                        group = rep(1, nrow(as.data.frame(data))),
                        weight = rep(1, nrow(as.data.frame(data))),
                        stat = "mean",
                        total = FALSE,
                        var_counts = FALSE,
                        f_stats = FALSE) {
  
  group <- as.matrix(group)
  the_vars <- names(as.data.frame(data))
  data <- data.frame(group = group, weight = weight, data, stringsAsFactors = FALSE)
  
  ns <- data %>% 
    dplyr::group_by(group) %>% 
    dplyr::count(wt = weight)
  
  # Statistic function dispatcher
  the_function <- switch(stat,
                         "mean"   = function(x, w) stats::weighted.mean(x, w, na.rm = TRUE),
                         "sum"    = function(x, w) sum(x * w, na.rm = TRUE),
                         "var"    = function(x, w) Hmisc::wtd.var(x, w, na.rm = TRUE),
                         "min"    = function(x, w) min(x, na.rm = TRUE),
                         "max"    = function(x, w) max(x, na.rm = TRUE),
                         "median" = function(x, w) matrixStats::weightedMedian(x, w, na.rm = TRUE),
                         stop("Unknown statistic. Use: mean, sum, median, var, min, max")
  )
  
  # Main calculations
  means <- data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_at(the_vars, ~ the_function(., weight)) %>%
    left_join(ns, by = "group") %>%
    mutate(group = as.character(group)) %>%
    select(group, n, everything())
  
  # Optional total row
  if (total) {
    total_row <- data %>%
      dplyr::group_by(group = "Total") %>%
      dplyr::summarise_at(the_vars, ~ the_function(., weight))
  }
  
  # Optional non-missing counts
  if (var_counts) {
    v_counts <- purrr::map_df(data[, -c(1:2)], 
                              ~ sum(as.numeric(!is.na(.)) * weight, na.rm = TRUE))
    v_counts$group <- "Counts (non-missing)"
  }
  
  # Optional F-stats
  if (f_stats) {
    Fstats <- purrr::map_dfr(data[, -c(1:2)], function(x) {
      mod <- stats::aov(x ~ as.factor(group), weight = weight)
      smod <- stats::anova(mod)
      fstat <- smod$`F value`[1]
      pval  <- stats::pf(fstat, smod$Df[1], smod$Df[2], lower.tail = FALSE)
      data.frame(fstat = fstat, pval = pval)
    }) %>%
      mutate(group = c("F Stat", "p value"))
  }
  
  output <- bind_rows(get0("v_counts"), get0("total_row"), get0("Fstats"), means)
  output %>% dplyr::select(group, one_of("n"), everything())
}

#' Sorted group means by standardized difference from total
#'
#' @param data Data frame of variables to analyze
#' @param group Grouping variable
#' @param weight Weighting variable
#' @export
group_means_sorted <- function(data,
                               group = NULL,
                               weight = rep(1, nrow(as.data.frame(data)))) {
  
  seg_name <- paste0(colnames(as.matrix(group)), "_")
  group <- as.character(as.matrix(group))
  the_vars <- names(data)
  data <- data.frame(group = group, weight = weight, data, stringsAsFactors = FALSE)
  
  # Total
  total_row <- data %>%
    group_by(group = "Total") %>%
    summarise_at(the_vars, ~ weighted.mean(., weight, na.rm = TRUE))
  
  # Group means
  means <- data %>%
    group_by(group) %>%
    summarise_at(the_vars, ~ weighted.mean(., weight, na.rm = TRUE))
  
  # Z-scores for sorting
  z_means <- data %>%
    mutate_at(vars(-c(1:2)), scale) %>%
    group_by(group) %>%
    summarise_at(the_vars, ~ weighted.mean(., weight, na.rm = TRUE))
  
  # Build sorted output
  sorted_list <- list()
  for (i in seq_along(unique(group))) {
    entry <- data.frame(
      Var   = names(means)[-1],
      Mean  = as.numeric(means[i, -1]),
      Diff  = as.numeric(means[i, -1] - total_row[-1]),
      Zmean = as.numeric(z_means[i, -1])
    ) %>%
      arrange(desc(Zmean))
    
    names(entry) <- c(paste0(seg_name, means$group[i]),
                      paste0("Mean_", means$group[i]),
                      paste0("Diff_", means$group[i]),
                      paste0("Zmean_", means$group[i]))
    
    sorted_list[[i]] <- entry
  }
  
  output <- bind_cols(sorted_list)
  clip_df(output, col.names = TRUE)
  output
}

#' Multiple comparisons (Tukey HSD) after ANOVA
#'
#' @param data Data frame
#' @param group Grouping variable
#' @param weight Weighting variable (note: TukeyHSD does not use weights directly)
#' @param conf_level Confidence level
#' @export
group_means_mult_comp <- function(data,
                                  group = rep(1, nrow(as.data.frame(data))),
                                  weight = rep(1, nrow(as.data.frame(data))),
                                  conf_level = 0.95) {
  
  group <- as.matrix(group)
  data <- data.frame(group = group, weight = weight, data, stringsAsFactors = FALSE)
  
  aov_fits <- purrr::map(data[, -c(1:2)], function(x) {
    aov(x ~ as.factor(group))
  })
  
  tukey_results <- purrr::map_dfr(aov_fits, function(model) {
    as.data.frame(TukeyHSD(model, conf.level = conf_level)[[1]]) %>%
      tibble::rownames_to_column("comparison")
  }, .id = "variable")
  
  tukey_results
}

#' Row center and/or standardize
#' @export
row_standardize <- function(data, center = TRUE, scale = TRUE) {
  data <- as.matrix(data)
  if (center) data <- data - rowMeans(data, na.rm = TRUE)
  if (scale)  data <- data / apply(data, 1, sd, na.rm = TRUE)
  as_tibble(data)
}

#' Top n Box
#' @export
make_tnb <- function(vec, scale = 1:7, n = 2) {
  threshold <- max(scale) - (n - 1)
  out <- tibble::as_tibble(ifelse(vec >= threshold, 1, 0))
  names(out) <- paste0(names(vec) %||% "var", "_T", n, "B")
  out
}

#' Bottom n Box
#' @export
make_bnb <- function(vec, scale = 1:7, n = 2) {
  threshold <- min(scale) + (n - 1)
  out <- tibble::as_tibble(ifelse(vec <= threshold, 1, 0))
  names(out) <- paste0(names(vec) %||% "var", "_B", n, "B")
  out
}

#' Top + Bottom n Box
#' @export
make_tbnb <- function(vec, scale = 1:7, top_n = 2, bottom_n = 2) {
  bind_cols(make_tnb(vec, scale, top_n), make_bnb(vec, scale, bottom_n))
}

#' Net Promoter Score
#' @export
nps <- function(x) mean(x >= 9, na.rm = TRUE)*100 - mean(x <= 6, na.rm = TRUE)*100

#' Create crosstab with counts and percentages
#'
#' @param data A data frame
#' @param row_var Row variable (unquoted)
#' @param col_var Column variable (unquoted)
#' @param weight Optional weighting variable
#' @param percent Type of percentage: "row" (default), "col", "total", or "none"
#' @param include_n Include raw counts column? Default TRUE
#'
#' @export
crosstab <- function(data, 
                     row_var, 
                     col_var, 
                     weight = NULL,
                     percent = c("row", "col", "total", "none"),
                     include_n = TRUE) {
  
  percent <- match.arg(percent)
  
  # Base count (weighted or unweighted)
  if (is.null(weight)) {
    tab <- data %>%
      count({{row_var}}, {{col_var}}, name = "n")
  } else {
    tab <- data %>%
      count({{row_var}}, {{col_var}}, wt = {{weight}}, name = "n")
  }
  
  # Add percentages
  if (percent == "row") {
    tab <- tab %>%
      group_by({{row_var}}) %>%
      mutate(pct = n / sum(n) * 100,
             pct = round(pct, 1))
  } else if (percent == "col") {
    tab <- tab %>%
      group_by({{col_var}}) %>%
      mutate(pct = n / sum(n) * 100,
             pct = round(pct, 1))
  } else if (percent == "total") {
    total_n <- sum(tab$n)
    tab <- tab %>%
      mutate(pct = n / total_n * 100,
             pct = round(pct, 1))
  } else {  # "none"
    tab <- tab %>% mutate(pct = NA_real_)
  }
  
  # Final arrangement
  tab <- tab %>%
    ungroup() %>%
    arrange({{row_var}}, {{col_var}})
  
  if (!include_n) {
    tab <- tab %>% select(-n)
  }
  
  tab
}

#' Quick bar plot (ggplot2)
#' @export
quick_bar <- function(data, x, y = NULL, title = NULL, fill = NULL) {
  p <- ggplot(data, aes(x = {{x}}, y = {{y}} %||% stat(count), fill = {{fill}})) +
    geom_bar(stat = if (is.null(y)) "count" else "identity") +
    theme_minimal() +
    labs(title = title)
  print(p)
  invisible(p)
}

#' Quick correlation heatmap
#' @export
quick_corrplot <- function(data, method = "pearson", title = "Correlation Matrix") {
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("Install corrplot package to use this function")
  }
  num_data <- data %>% select(where(is.numeric))
  corr_mat <- cor(num_data, use = "pairwise.complete.obs", method = method)
  corrplot::corrplot(corr_mat, method = "color", type = "upper", 
                     tl.cex = 0.8, title = title, mar = c(0,0,2,0))
}

#' Quick NPS-style stacked bar
#' @export
quick_nps_bar <- function(x, title = "NPS Breakdown") {
  df <- data.frame(score = factor(x, levels = 0:10)) %>%
    count(score) %>%
    mutate(pct = n / sum(n) * 100)
  
  ggplot(df, aes(x = score, y = pct, fill = score >= 9)) +
    geom_col() +
    scale_fill_manual(values = c("grey70", "#2ecc71")) +
    theme_minimal() +
    labs(title = title, y = "Percentage", x = "Likelihood to Recommend (0-10)") +
    theme(legend.position = "none")
}

#' nth largest / smallest
#' @export
nth_max <- function(x, N){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x). Setting N=length(x)')
    N <- length(x)
  }
  sort(x,partial=len-N+1)[len-N+1]
}

nth_min <- function(x, N){
  len <- length(x)
  if(N>len){
    warning('N greater than length(x). Setting N=length(x)')
    N <- length(x)
  }
  -sort(-x,partial=len-N+1)[len-N+1]
}

#' Fill NAs
#'
#' Fill NA's in a vector or factor with last non-NA value.
#' If firstBack is TRUE, it will fill in leading NA's with the first
#' non-NA value. If FALSE, it will not change leading NA's.
#' Taken from
#' http://www.cookbook-r.com/Manipulating_data/Filling_in_NAs_with_last_non-NA_value/
#' @export
fillNAgaps <- function(x, firstBack=FALSE) {
  
  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }
  
  goodIdx <- !is.na(x)
  
  # These are the non-NA values from x only
  # Add a leading NA or take the first good value, depending on firstBack
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])
  
  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1
  
  x <- goodVals[fillIdx]
  
  # If it was originally a factor, convert it back
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }
  
  x
}

#' NPS transformation on likelihood to purchase data
#'
#' Transform a vector containing 0-10 likelihood to recommend data into -100, 0,
#'  100 where 0 to 6 are assigned -100, 7 & 8 are assigned 0, and 9 & 10 are
#'  assigned 100.
#'
#' @param v A vector of likelihood to recommend data (0-10).
#' @param d The numeric threshold for detractors, default is 6 (inclusive).
#' @param p The numeric threshold for promotors, defualt is 9 (inclusive).
#'
#' @export
NPS_transformation <- function(v, d = 6, p = 9) {
  
  ifelse(v <= d, -100, ifelse(v < p, 0, 100))
}

#' NPS score (aggregate) from likelihood to purchase data
#'
#' Caclulate the NPS score given a vector of likelihood to recommend data
#' (0-10).
#'
#' @param v A vector of likelihood to recommend data (0-10).
#' @param weight A vector of weights that lines up with v respondent by
#' respondent
#' @param d The numeric threshold for detractors, default is 6 (inclusive).
#' @param p The numeric threshold for promotors, defualt is 9 (inclusive).
#'
#' @importFrom stats weighted.mean
#'
#' @export
NPS_calc <- function(v, weight = rep(1, length(v)), d = 6, p = 9) {
  # Compute the NPS statistic of a set of scores v
  stats::weighted.mean(NPS_transformation(v, d, p), weight, na.rm = T)
}

#' NPS variance
#'
#' Calculate the NPS variance given a vector of likelihood to recommend (0-10).
#' This replicates the variance calculations used in the old NPS projects.
#'
#' @param v A vector of likelihood to recommend data (0-10).
#' @param weight A vector of weights that lines up with v respondent by
#' respondent
#' @param d The numeric threshold for detractors, default is 6 (inclusive).
#' @param p The numeric threshold for promotors, defualt is 9 (inclusive).
#'
#' @export
NPS_var <- function(v, weight = rep(1, length(v)), d = 6, p = 9){
  
  W <- weight
  P <- (ifelse(v >= p, W, 0))
  D <- (ifelse(v <= d, W, 0))
  Pass <- (ifelse(v > d & v < p, W, 0))
  
  var <- 100^2 * (sum(P, na.rm = T)/sum(W, na.rm = T) *
                    (1-(sum(P, na.rm = T) /sum(W, na.rm = T))) +
                    sum(D, na.rm = T)/sum(W, na.rm = T) *
                    (1-(sum(D, na.rm = T)/sum(W, na.rm = T))) +
                    2*sum(P, na.rm = T)/sum(W, na.rm = T) *
                    sum(D, na.rm = T)/sum(W, na.rm = T) )
  return(var)
}

#' NPS confidence interval
#'
#' Calculate NPS confidence interval given a vector of likelhood to recommend
#' (1-10). The old original process used an alpha of 0.2. Here the defualt
#' is 0.05.
#'
#' @param v A vector of likelihood to recommend data (0-10).
#' @param weight A vector of weights that lines up with v respondent by
#' respondent
#' @param d The numeric threshold for detractors, default is 6 (inclusive).
#' @param p The numeric threshold for promotors, defualt is 9 (inclusive).
#' @param alpha The alpha level to be used in the calculation.
#'
#' @examples
#' # create 1000 fake likelihood to recommend answers with NPS = 10
#' set.seed(1234)
#' likelihood_to_recommend <- sample(c(sample(1:6, 30, TRUE),
#' sample(7:8, 30, TRUE), sample(9:10, 40, TRUE), TRUE), 1000, TRUE)
#' NPS_calc(likelihood_to_recommend)
#' NPS_ci(likelihood_to_recommend)
#'
#' @importFrom stats qnorm
#'
#' @export
NPS_ci <- function(v, weight = rep(1, length(v)), d = 6, p = 9, alpha = 0.05){
  
  W <- weight
  z <- -stats::qnorm(alpha / 2)
  
  VAR <- NPS_var(v, W, d, p)
  DE <- design_effect(W)
  CI_Range <- sqrt(VAR/sum(W, na.rm = T) * DE) * z
  
  NPS <- NPS_calc(v, weight, d, p)
  
  fin_out <- data.frame(NPS - CI_Range, NPS, NPS + CI_Range)
  names(fin_out) <- c("Lower", "Mean", "Upper")
  
  return(fin_out)
}

#' Implicit Association Tests (IAT)
#' Consider carefully before you actually do this, improved approaches are likely more desirable.
#' rescale on a -1 to 1 scale where -1 indicates fast and negative response
#' and 1 indicates fast and positive response, while 0s indicate slow responses.
#' shave off 5% on the tails.
#' Usually have a pre and post test, scale the pre and post testes separately.
#' A good sense check is that post tests are almost always quicker than pre tests.
#' In general, IAT is used when you are trying to get at
#' Type 1 (gut initial reaction) rather than Type 2 (logical reasoning) thinking.
#' seen it used effectively in either comparing different brands
#' (what is each brand associated with), or in comparing effects of ads (post-pre).
#' you want to know if your ad is conveying the gut reaction of excitement and adventure to a brand?
#' this is a way to measure its effectiveness.
#' you want to know your competitors gut reaction compares to yours in terms of honesty and trustworthiness?
#' this is a way to measure the difference.
#' the score is negative if they respond that it is not associated with the
#' attribute and positive if it is associated with the attribute.
#' if it takes people longer to associate or not associate the brand with the
#' attribute the score is weaker and closer to 0.
#' if they are able to make the determination quickly the score is closer
#' to 1 (associated) or -1 (not associated).
#' put in some quality control metrics to remove answers that were either
#' too fast of too slow before rescaling the scores between 1 and -1.
#'
#' A watered down version of Implicit Association Tests IAT.
#' For a more rigorous example of this methodology
#' see: https://implicit.harvard.edu/implicit/selectatest.html.
#' @param t A matrix of times (typically in milliseconds).
#' Each row is a respondent and each column is an item. 
#' The values indicate how long a respondent took to respond to that item.
#' @param pick A matrix of answers (1=yes, 2=no.)
#' @param negative Default is 2. This is typical to have punch 2 be no/disagree.
#' @param positive Default is 1. This is typical to have punch 1 be yes/agree.
#' @param trim Defualt is 0.05
#' @return Returns scores for each respondent.
#' @importFrom stats quantile
#' @importFrom tibble as_tibble
#' @export
iat <- function(t, pick, negative = 2, positive = 1, trim = 0.05) {
  q1 = stats::quantile(as.matrix(t), trim)
  q2 = stats::quantile(as.matrix(t), 1 - trim)
  rangeSize = q2 - q1
  scores <- apply(t, c(1, 2), function(x)
    ifelse(x >= q1 & x <= q2, (1 - ((x-q1)/rangeSize)), NA))
  
  pick[pick == negative] <- -1
  pick[pick == positive] <- 1
  
  scores <- tibble::as_tibble(scores * pick)
  names(scores) <- paste0(names(pick), "_score")
  cat("Score Means:", "\n")
  print(colMeans(scores, na.rm = T))
  return(scores)
}

#' Design effect
#' Calculate the design effect for a given vector of weights. Note that design
#' efficiency can then be manually calculated as 1 / design effect.
#' @param weight A vector of weights.
#' @seealso diagis::weight_plot
#' @export
design_effect <- function(weight){
  W <- weight
  DE <- length(W) * sum(W^2, na.rm = T)/(sum(W, na.rm = T)^2)
  return(DE)
}

#' Design efficiency
#' Calculate the design efficiency. This is 1 / design effect.
#' @param weight A vector of weights.
#' @seealso diagis::weight_plot
#' @export
design_efficiency <- function(weight){
  W <- weight
  DEF <- 1 / (length(W) * sum(W^2, na.rm = T)/(sum(W, na.rm = T)^2))
  return(DEF)
}

#' Effective sample size
#' @param weight A vector of weights.
#' @export
effective_sample_size <- function(weight){
  W <- weight
  DE <- (length(W) * sum(W^2, na.rm = T)/(sum(W, na.rm = T)^2))
  ESS <- length(weight) / DE
  return(ESS)
}

#' Johnson Relative Weights
#' Credit to Jacob Nelson
#' Takes the raw data (rows are respondents, columns are variables), with the
#'  dependent variable as the first variable in the data frame, and produces
#' @param data A data frame of variables with the dependent variable as the first.
#' @export
regression_jrw <- function(data) {
  #require(tibble)
  Labels <- names(data)[2:length(data)]
  numVar <- NCOL(data)
  Variables <-  names(data)[2:numVar]
  
  data <- cor(data, use = "pairwise.complete.obs")
  RXX <- data[2:numVar, 2:numVar]
  RXY <- data[2:numVar, 1]
  
  RXX.eigen <- eigen(RXX)
  D <- diag(RXX.eigen$val)
  delta <- sqrt(D)
  
  lambda <- RXX.eigen$vec %*% delta %*% t(RXX.eigen$vec)
  lambdasq <- lambda ^ 2
  beta <- solve(lambda) %*% RXY
  rsquare <- sum(beta ^ 2)
  
  RawWgt <- lambdasq %*% beta ^ 2
  import <- (RawWgt / rsquare)
  
  tibble(Variables, RelWeight = as.numeric(RawWgt), Importance = as.numeric(import))
}

#' Johnson Relative Weights - Logistic Models
#' @param data A data frame of variables with the dependent variable as the first.
#' @export
regression_log_jrw <- function(data){
  #require(tibble)
  numVar <- NCOL(data)
  Variables <- names(data)[2:numVar]
  
  Y <- data[,1]
  data.preds <- data[,2:length(data[1,])]
  X <- scale(data.preds)
  
  X.svd <- svd(X)
  Q <- X.svd$v
  P <- X.svd$u
  Z <- P %*% t(Q)
  
  Z.stand <- scale(Z)
  
  Lambda <- solve(t(Z.stand) %*% Z.stand) %*% t(Z.stand) %*% X
  
  logrfit <- glm(Y ~ Z.stand, family = binomial)
  unstCoefs <- coef(logrfit)
  b <- unstCoefs[2:length(unstCoefs)]
  LpredY <- predict(logrfit, newdata = data, type = "response")
  lYhat <- log(LpredY / (1 - LpredY))
  stdlYhat <- sd(lYhat)
  getting.Rsq <- lm(LpredY ~ Y)
  Rsq <- summary(getting.Rsq)$r.squared
  beta <- b * ((sqrt(Rsq)) / stdlYhat)
  
  epsilon <- Lambda ^ 2 %*% beta ^ 2
  R.sq <- sum(epsilon)
  PropWeights <- (epsilon / R.sq)
  
  tibble(Variables, RelWeight = as.numeric(epsilon), Importance = as.numeric(PropWeights))
}
