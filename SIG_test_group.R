compare_split_variance_f_test <- function(split_result_list) {
  group_names <- names(split_result_list)
  group_results <- list()
  all_core <- c()
  all_non_core <- c()
  for (group in group_names) {
    res <- split_result_list[[group]]
    core_vals <- as.numeric(unlist(res$core_otus))
    non_core_vals <- as.numeric(unlist(res$non_core_otus))
    core_vals <- core_vals[!is.na(core_vals)]
    non_core_vals <- non_core_vals[!is.na(non_core_vals)]
    all_core <- c(all_core, core_vals)
    all_non_core <- c(all_non_core, non_core_vals)
    if (length(core_vals) > 1 && length(non_core_vals) > 1) {
      ftest <- var.test(core_vals, non_core_vals, alternative = "two.sided")
      group_results[[group]] <- ftest
    } else {
      group_results[[group]] <- "Not enough observations for F-test"
    }
  }
  combined_result <- if (length(all_core) > 1 && length(all_non_core) > 1) {
    var.test(all_core, all_non_core, alternative = "two.sided")
  } else {
    "Not enough observations for combined F-test"
  }
  return(list(
    group_wise = group_results,
    combined = combined_result
  ))
}

compare_split_kruskal_wallis <- function(split_result_list) {
  group_names <- names(split_result_list)
  results <- list()
  all_values <- c()
  all_groups <- c()
  for (group in group_names) {
    res <- split_result_list[[group]]
    core_vals <- as.numeric(unlist(res$core_otus))
    non_core_vals <- as.numeric(unlist(res$non_core_otus))
    core_vals <- core_vals[!is.na(core_vals)]
    non_core_vals <- non_core_vals[!is.na(non_core_vals)]
    all_values <- c(all_values, core_vals, non_core_vals)
    all_groups <- c(all_groups, 
                    rep(paste0(group, "_core"), length(core_vals)),
                    rep(paste0(group, "_non_core"), length(non_core_vals)))
    if (length(core_vals) > 0 && length(non_core_vals) > 0) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      ktest <- kruskal.test(combined ~ groups)
      results[[group]] <- ktest
    } else {
      results[[group]] <- "Not enough observations for Kruskal-Wallis test"
    }
  }
  combined_result <- if (length(unique(all_groups)) > 1 && length(all_values) > 1) {
    kruskal.test(all_values ~ factor(all_groups))
  } else {
    "Not enough observations for combined Kruskal-Wallis test"
  }
  return(list(
    group_wise = results,
    combined = combined_result
  ))
}

compare_split_levene_test <- function(split_result_list) {
  group_names <- names(split_result_list)
  results <- list()
  all_values <- c()
  all_groups <- c()
  for (group in group_names) {
    res <- split_result_list[[group]]
    core_vals <- as.numeric(unlist(res$core_otus))
    non_core_vals <- as.numeric(unlist(res$non_core_otus))
    core_vals <- core_vals[!is.na(core_vals)]
    non_core_vals <- non_core_vals[!is.na(non_core_vals)]
    all_values <- c(all_values, core_vals, non_core_vals)
    all_groups <- c(all_groups, 
                    rep(paste0(group, "_core"), length(core_vals)),
                    rep(paste0(group, "_non_core"), length(non_core_vals)))
    if (length(core_vals) > 1 && length(non_core_vals) > 1) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      ltest <- car::leveneTest(combined ~ groups, center = median)
      results[[group]] <- ltest
    } else {
      results[[group]] <- "Not enough observations for Levene's test"
    }
  }
  combined_result <- if (length(unique(all_groups)) > 1 && length(all_values) > 1) {
    car::leveneTest(all_values ~ factor(all_groups), center = median)
  } else {
    "Not enough observations for combined Levene's test"
  }
  return(list(
    group_wise = results,
    combined = combined_result
  ))
}

compare_split_fligner_test <- function(split_result_list) {
  group_names <- names(split_result_list)
  results <- list()
  all_values <- c()
  all_groups <- c()
  for (group in group_names) {
    res <- split_result_list[[group]]
    core_vals <- as.numeric(unlist(res$core_otus))
    non_core_vals <- as.numeric(unlist(res$non_core_otus))
    core_vals <- core_vals[!is.na(core_vals)]
    non_core_vals <- non_core_vals[!is.na(non_core_vals)]
    all_values <- c(all_values, core_vals, non_core_vals)
    all_groups <- c(all_groups, 
                    rep(paste0(group, "_core"), length(core_vals)),
                    rep(paste0(group, "_non_core"), length(non_core_vals)))
    if (length(core_vals) > 1 && length(non_core_vals) > 1) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      ftest <- fligner.test(combined ~ groups)
      results[[group]] <- ftest
    } else {
      results[[group]] <- "Not enough observations for Fligner-Killeen test"
    }
  }
  combined_result <- if (length(unique(all_groups)) > 1 && length(all_values) > 1) {
    fligner.test(all_values ~ factor(all_groups))
  } else {
    "Not enough observations for combined Fligner-Killeen test"
  }
  return(list(
    group_wise = results,
    combined = combined_result
  ))
}

compute_split_all_results <- function(split_result_list) {
  car_available <- requireNamespace("car", quietly = TRUE)
  if (!car_available) {
    warning("Package 'car' not available. Levene's test will be skipped.")
  }
  group_names <- names(split_result_list)
  results <- list(
    f_test = list(group_wise = list(), combined = NULL),
    kruskal_wallis = list(group_wise = list(), combined = NULL),
    levene_test = list(group_wise = list(), combined = NULL),
    fligner_test = list(group_wise = list(), combined = NULL)
  )
  all_core <- numeric(0)
  all_non_core <- numeric(0)
  all_values <- numeric(0)
  all_groups <- character(0)
  for (group in group_names) {
    res <- split_result_list[[group]]
    core_vals <- as.numeric(unlist(res$core_otus))
    non_core_vals <- as.numeric(unlist(res$non_core_otus))
    core_vals <- core_vals[!is.na(core_vals)]
    non_core_vals <- non_core_vals[!is.na(non_core_vals)]
    all_core <- c(all_core, core_vals)
    all_non_core <- c(all_non_core, non_core_vals)
    all_values <- c(all_values, core_vals, non_core_vals)
    all_groups <- c(all_groups, 
                    rep(paste0(group, "_core"), length(core_vals)),
                    rep(paste0(group, "_non_core"), length(non_core_vals)))
    if (length(core_vals) > 1 && length(non_core_vals) > 1) {
      results$f_test$group_wise[[group]] <- var.test(core_vals, non_core_vals, alternative = "two.sided")
    } else {
      results$f_test$group_wise[[group]] <- "Not enough observations for F-test"
    }
    if (length(core_vals) > 0 && length(non_core_vals) > 0) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      results$kruskal_wallis$group_wise[[group]] <- kruskal.test(combined ~ groups)
    } else {
      results$kruskal_wallis$group_wise[[group]] <- "Not enough observations for Kruskal-Wallis test"
    }
    if (car_available && length(core_vals) > 1 && length(non_core_vals) > 1) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      results$levene_test$group_wise[[group]] <- car::leveneTest(combined ~ groups, center = median)
    } else if (!car_available) {
      results$levene_test$group_wise[[group]] <- "car package not available"
    } else {
      results$levene_test$group_wise[[group]] <- "Not enough observations for Levene's test"
    }
    if (length(core_vals) > 1 && length(non_core_vals) > 1) {
      combined <- c(core_vals, non_core_vals)
      groups <- factor(c(rep("core", length(core_vals)), rep("non_core", length(non_core_vals))))
      results$fligner_test$group_wise[[group]] <- fligner.test(combined ~ groups)
    } else {
      results$fligner_test$group_wise[[group]] <- "Not enough observations for Fligner-Killeen test"
    }
  }
  if (length(all_core) > 1 && length(all_non_core) > 1) {
    results$f_test$combined <- var.test(all_core, all_non_core, alternative = "two.sided")
  } else {
    results$f_test$combined <- "Not enough observations for combined F-test"
  }
  if (length(unique(all_groups)) > 1 && length(all_values) > 1) {
    results$kruskal_wallis$combined <- kruskal.test(all_values ~ factor(all_groups))
  } else {
    results$kruskal_wallis$combined <- "Not enough observations for combined Kruskal-Wallis test"
  }
  if (car_available && length(unique(all_groups)) > 1 && length(all_values) > 1) {
    results$levene_test$combined <- car::leveneTest(all_values ~ factor(all_groups), center = median)
  } else if (!car_available) {
    results$levene_test$combined <- "car package not available"
  } else {
    results$levene_test$combined <- "Not enough observations for combined Levene's test"
  }
  if (length(unique(all_groups)) > 1 && length(all_values) > 1) {
    results$fligner_test$combined <- fligner.test(all_values ~ factor(all_groups))
  } else {
    results$fligner_test$combined <- "Not enough observations for combined Fligner-Killeen test"
  }
  return(results)
}

create_split_summary_table <- function(all_results) {
  safe_extract <- function(test_result, stat_type = "p_value") {
    if (is.character(test_result)) {
      return(NA)
    }
    if (stat_type == "p_value") {
      if (is.list(test_result) && "p.value" %in% names(test_result)) {
        return(test_result$p.value)
      } else if (is.list(test_result) && "Pr(>F)" %in% names(test_result)) {
        return(test_result$`Pr(>F)`[1])
      } else if (is.list(test_result) && "Pr(>Chi)" %in% names(test_result)) {
        return(test_result$`Pr(>Chi)`)
      }
    } else if (stat_type == "statistic") {
      if (is.list(test_result) && "statistic" %in% names(test_result)) {
        return(as.numeric(test_result$statistic))
      } else if (is.list(test_result) && "F value" %in% names(test_result)) {
        return(test_result$`F value`[1])
      }
    }
    return(NA)
  }
  group_names <- names(all_results$f_test$group_wise)
  summary_tables <- list()
  for (group in group_names) {
    group_summary <- data.frame(
      Test = character(0),
      Statistic = numeric(0),
      P_value = numeric(0),
      CI_Lower = numeric(0),
      CI_Upper = numeric(0),
      stringsAsFactors = FALSE
    )
    f_result <- all_results$f_test$group_wise[[group]]
    f_stat <- safe_extract(f_result, "statistic")
    f_p <- safe_extract(f_result, "p_value")
    f_ci_lower <- if (is.list(f_result) && !is.na(f_result$conf.int[1])) f_result$conf.int[1] else NA
    f_ci_upper <- if (is.list(f_result) && !is.na(f_result$conf.int[2])) f_result$conf.int[2] else NA
    group_summary <- rbind(group_summary, data.frame(
      Test = "F-test",
      Statistic = f_stat,
      P_value = f_p,
      CI_Lower = f_ci_lower,
      CI_Upper = f_ci_upper
    ))
    kw_result <- all_results$kruskal_wallis$group_wise[[group]]
    kw_stat <- safe_extract(kw_result, "statistic")
    kw_p <- safe_extract(kw_result, "p_value")
    group_summary <- rbind(group_summary, data.frame(
      Test = "Kruskal-Wallis",
      Statistic = kw_stat,
      P_value = kw_p,
      CI_Lower = NA,
      CI_Upper = NA
    ))
    levene_result <- all_results$levene_test$group_wise[[group]]
    levene_stat <- safe_extract(levene_result, "statistic")
    levene_p <- safe_extract(levene_result, "p_value")
    group_summary <- rbind(group_summary, data.frame(
      Test = "Levene",
      Statistic = levene_stat,
      P_value = levene_p,
      CI_Lower = NA,
      CI_Upper = NA
    ))
    fligner_result <- all_results$fligner_test$group_wise[[group]]
    fligner_stat <- safe_extract(fligner_result, "statistic")
    fligner_p <- safe_extract(fligner_result, "p_value")
    group_summary <- rbind(group_summary, data.frame(
      Test = "Fligner-Killeen",
      Statistic = fligner_stat,
      P_value = fligner_p,
      CI_Lower = NA,
      CI_Upper = NA
    ))
    group_summary$Group <- group
    summary_tables[[group]] <- group_summary
  }
  final_summary <- do.call(rbind, summary_tables)
  rownames(final_summary) <- NULL
  final_summary <- final_summary[, c("Group", "Test", "Statistic", "P_value", "CI_Lower", "CI_Upper")]
  final_summary$P_value <- sapply(final_summary$P_value, function(p_val) {
    if (is.na(p_val) || is.character(p_val)) {
      return(p_val)
    }
    p_val <- as.numeric(p_val)
    if (p_val < 0.0001) {
      return("<0.0001")
    }
    return(format(round(p_val, 4), nsmall = 4))
  })
  
  return(final_summary)
}