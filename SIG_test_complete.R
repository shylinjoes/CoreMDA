compare_variance_f_test <- function(core_ids, non_core_ids) {
  core_data <- as.vector(data.matrix(core_ids[, -1]))
  non_core_data <- as.vector(data.matrix(non_core_ids[, -1]))
  core_data <- core_data[!is.na(core_data)]
  non_core_data <- non_core_data[!is.na(non_core_data)]
  if (length(core_data) < 2 || length(non_core_data) < 2) {
    stop("Each group must have at least two non-NA observations for the F-test.")
  }
  ftest <- var.test(core_data, non_core_data, alternative = "two.sided")
  return(ftest)
}

compare_kruskal_wallis <- function(core_ids, non_core_ids) {
  core_data <- as.vector(data.matrix(core_ids[, -1]))
  non_core_data <- as.vector(data.matrix(non_core_ids[, -1]))
  core_data <- core_data[!is.na(core_data)]
  non_core_data <- non_core_data[!is.na(non_core_data)]
  if (length(core_data) < 1 || length(non_core_data) < 1) {
    stop("Each group must have at least one non-NA observation for the Kruskal-Wallis test.")
  }
  all_data <- c(core_data, non_core_data)
  groups <- factor(c(rep("core", length(core_data)), rep("non_core", length(non_core_data))))
  ktest <- kruskal.test(all_data ~ groups)
  return(ktest)
}

compare_variance_levene <- function(core_ids, non_core_ids) {
  core_data <- as.vector(data.matrix(core_ids[, -1]))
  non_core_data <- as.vector(data.matrix(non_core_ids[, -1]))
  core_data <- core_data[!is.na(core_data)]
  non_core_data <- non_core_data[!is.na(non_core_data)]
  if (length(core_data) < 2 || length(non_core_data) < 2) {
    stop("Each group must have at least two non-NA observations.")
  }
  all_data <- c(core_data, non_core_data)
  groups <- factor(c(rep("core", length(core_data)), rep("non_core", length(non_core_data))))
  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required. Install with: install.packages('car')")
  }
  levene_test <- car::leveneTest(all_data, groups, center = median)  # Non-parametric (median-based)
  return(levene_test)
}

compare_variance_fligner <- function(core_ids, non_core_ids) {
  core_data <- as.vector(data.matrix(core_ids[, -1]))
  non_core_data <- as.vector(data.matrix(non_core_ids[, -1]))
  core_data <- core_data[!is.na(core_data)]
  non_core_data <- non_core_data[!is.na(non_core_data)]
  if (length(core_data) < 2 || length(non_core_data) < 2) {
    stop("Each group must have at least two non-NA observations.")
  }
  all_data <- c(core_data, non_core_data)
  groups <- factor(c(rep("core", length(core_data)), rep("non_core", length(non_core_data))))
  fligner_test <- fligner.test(all_data ~ groups)
  
  return(fligner_test)
}

compute_otu_all_tests <- function(core_ids, non_core_ids) {
  results <- list()
  results$f_test <- tryCatch({
    compare_variance_f_test(core_ids, non_core_ids)
  }, error = function(e) list(error = e$message))
  results$levene_test <- tryCatch({
    compare_variance_levene(core_ids, non_core_ids)
  }, error = function(e) list(error = e$message))
  results$fligner_test <- tryCatch({
    compare_variance_fligner(core_ids, non_core_ids)
  }, error = function(e) list(error = e$message))
  results$kruskal_wallis_test <- tryCatch({
    compare_kruskal_wallis(core_ids, non_core_ids)
  }, error = function(e) list(error = e$message))
  results$summary_table <- data.frame(
    Test = c("F-test", "Levene", "Fligner-Killeen", "Kruskal-Wallis"),
    Statistic = c(
      ifelse(is.null(results$f_test$statistic), NA, results$f_test$statistic),
      ifelse(is.null(results$levene_test$`F value`), NA, results$levene_test$`F value`[1]),
      ifelse(is.null(results$fligner_test$statistic), NA, results$fligner_test$statistic),
      ifelse(is.null(results$kruskal_wallis_test$statistic), NA, results$kruskal_wallis_test$statistic)
    ),
    P_value = c(
      ifelse(is.null(results$f_test$p.value), NA, results$f_test$p.value),
      ifelse(is.null(results$levene_test$`Pr(>F)`), NA, results$levene_test$`Pr(>F)`[1]),
      ifelse(is.null(results$fligner_test$p.value), NA, results$fligner_test$p.value),
      ifelse(is.null(results$kruskal_wallis_test$p.value), NA, results$kruskal_wallis_test$p.value)
    ),
    CI_Lower = c(
      ifelse(is.list(results$f_test) && !is.na(results$f_test$conf.int[1]), results$f_test$conf.int[1], NA),
      NA, NA, NA
    ),
    CI_Upper = c(
      ifelse(is.list(results$f_test) && !is.na(results$f_test$conf.int[2]), results$f_test$conf.int[2], NA),
      NA, NA, NA
    ),
    stringsAsFactors = FALSE
  )
  results$summary_table$P_value <- sapply(results$summary_table$P_value, function(x) {
    if (is.na(x)) return(NA)
    if (is.character(x)) return(x)
    if (x < 0.0001) {
      return("<0.0001")
    }
    return(format(round(x, 4), nsmall = 4))
  })
  return(results)
}

