filter_non_zero_otu_func <- function(otu_table, tax_table, md_table) {
  colnames(otu_table)[1] <- "OTU ID"
  colnames(tax_table)[1] <- "OTU ID"
  OTU_id.common <- fastmatch::fmatch(tax_table[[1]], otu_table[[1]])
  common_otu_ids_table1 <- otu_table[OTU_id.common, ]
  common_otu_ids_table <- as.data.frame(common_otu_ids_table1[, -1])
  arab_md <- as.data.frame(md_table[, 1])
  colnames(arab_md) <- "Sample"
  id.com.md <- fastmatch::fmatch(arab_md$Sample, colnames(common_otu_ids_table))
  id.col <- stats::na.omit(id.com.md)
  final_otu_table <- common_otu_ids_table[, id.col, drop = FALSE]
  final_otu_table_with_row <- cbind(common_otu_ids_table1[, 1, drop = FALSE], final_otu_table)
  colnames(final_otu_table_with_row)[1] <- "OTU ID"
  row_sums <- rowSums(final_otu_table_with_row[, -1])
  non_zero <- row_sums != 0
  final_sorted_otu_tab <- final_otu_table_with_row[non_zero, , drop = FALSE]
  return(final_sorted_otu_tab)
}

prepare_and_split_otu_data <- function(otu_table, tax_table, md_table, split_by_col = "Treatment") {
  colnames(otu_table)[1] <- "OTU_ID"
  colnames(tax_table)[1] <- "OTU_ID"
  colnames(md_table)[1] <- "Sample_ID"
  common_id <- fastmatch::fmatch(tax_table$OTU_ID, otu_table$OTU_ID)
  common_otu_table <- otu_table[common_id, ]
  common_otu_t <- as.data.frame(t(common_otu_table))
  common_otu_t <- cbind(RowNames = rownames(common_otu_t), common_otu_t)
  rownames(common_otu_t) <- NULL
  new_column_names <- as.character(common_otu_t[1, ])
  common_otu_t <- common_otu_t[-1, ]
  colnames(common_otu_t) <- new_column_names
  colnames(common_otu_t)[1] <- "Sample_ID"
  final_otu_table <- merge(md_table, common_otu_t, by = "Sample_ID")
  if (!split_by_col %in% colnames(final_otu_table)) {
    stop(paste("Column", split_by_col, "not found in the final merged table."))
  }
  split_data <- split(final_otu_table, final_otu_table[[split_by_col]])
  split_data_processed <- lapply(split_data, function(df) {
    sample_ids <- df$Sample_ID
    otu_df <- df[, !(colnames(df) %in% colnames(md_table)), drop = FALSE]
    rownames(otu_df) <- sample_ids
    otu_df <- otu_df[, colnames(otu_df) != "Sample_ID", drop = FALSE]
    transposed_df <- as.data.frame(t(otu_df))
    transposed_df <- cbind(OTU_ID = rownames(transposed_df), transposed_df)
    rownames(transposed_df) <- NULL
    return(transposed_df)
  })
  return(split_data_processed)
}

filter_otu_table <- function(filtered_otu_tab, filter_type, ...) {
  args <- list(...)
  selected_otus <- NULL
  
  abundance <- function(x, min_count, prop, min_total_count) {
    x <- as.numeric(x)  
    id1 <- which(x >= min_count)
    pr1 <- length(id1) / length(x)
    sel1 <- ifelse(pr1 >= prop, 1, 0)
    sel2 <- ifelse(sum(x) >= min_total_count, 1, 0)
    sel <- sel1 & sel2
    return(sel)
  }
  
  occupancy <- function(x, percent) {
    x <- as.numeric(x)  # Ensure numeric
    x[is.na(x)] <- 0    # Handle NA values safely
    nonzero_count <- sum(x != 0)
    total_samples <- length(x)
    prop_present <- nonzero_count / total_samples
    return(prop_present >= percent)
  }
  hard_cutoff <- function(x, cutoff = 25, sites = 5, taxa_as_rows = TRUE) {
    if (!taxa_as_rows) x <- t(x)
    x <- as.data.frame(x)
    names(x)[1] <- "X"
    otu_table_filtered <- x[rowSums(x[, -1] > cutoff) >= sites, ]
    return(otu_table_filtered$X)
  }
  membership_filter <- function(final_otu_table, threshold = 0.5, min_similarities = 10) {
    sorenson_index <- function(vec1, vec2) {
      a <- sum(vec1 == 1 & vec2 == 1)
      b <- sum(vec1 == 1 & vec2 == 0)
      c <- sum(vec1 == 0 & vec2 == 1)
      if ((2 * a + b + c) == 0) return(0)
      return((2 * a) / (2 * a + b + c))
    }
    presence_absence <- final_otu_table[, -1] > 0
    rownames(presence_absence) <- final_otu_table[[1]]
    n <- nrow(presence_absence)
    result <- matrix(0, n, n)
    rownames(result) <- rownames(presence_absence)
    colnames(result) <- rownames(presence_absence)
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        score <- sorenson_index(presence_absence[i, ], presence_absence[j, ])
        result[i, j] <- score
        result[j, i] <- score
      }
    }
    diag(result) <- 1
    similar_count <- apply(result, 1, function(x) sum(x >= threshold))
    rows_to_keep <- similar_count > min_similarities
    original_otu_matrix <- final_otu_table
    rownames(original_otu_matrix) <- final_otu_table[[1]]
    filtered_otu_result <- original_otu_matrix[rows_to_keep, ]
    return(filtered_otu_result)
  }
  selected_otus <- switch(
    filter_type,
    
    "abundance_filter" = {
      if (!all(c("min_count", "prop", "min_total_count") %in% names(args))) {
        stop("For 'abundance_filter', provide 'min_count', 'prop', and 'min_total_count'.")
      }
      sel.ind <- apply(filtered_otu_tab[, -1], 1, abundance,
                       min_count = args$min_count,
                       prop = args$prop,
                       min_total_count = args$min_total_count)
      filtered_otu_tab[sel.ind, ]
    },
    
    "occupancy_filter" = {
      if (!("percent" %in% names(args))) {
        stop("For 'occupancy_filter', provide 'percent'.")
      }
      
      sel.ind <- apply(filtered_otu_tab[, -1], 1, occupancy,
                       percent = args$percent)
      
      filtered_otu_tab[sel.ind, ]
    },
    
    "abu_occ_filter" = {
      if (!all(c("percent", "min_count", "prop", "min_total_count") %in% names(args))) {
        stop("Provide 'percent', 'min_count', 'prop', and 'min_total_count'.")
      }
      occ.ind <- apply(filtered_otu_tab[, -1], 1, occupancy,
                       percent = args$percent)
      selected_otus_occupancy <- filtered_otu_tab[occ.ind, ]
      sel.ind <- apply(selected_otus_occupancy[, -1], 1, abundance,
                       min_count = args$min_count,
                       prop = args$prop,
                       min_total_count = args$min_total_count)
      
      selected_otus_occupancy[sel.ind, ]
    },
    "hardcutoff_filter" = {
      if (!all(c("cutoff", "sites", "taxa_as_rows") %in% names(args))) {
        stop("For 'hardcutoff_filter', provide 'cutoff', 'sites', and 'taxa_as_rows'.")
      }
      otu_ids <- hard_cutoff(filtered_otu_tab,
                             cutoff = args$cutoff,
                             sites = args$sites,
                             taxa_as_rows = args$taxa_as_rows)
      filtered_otu_tab[filtered_otu_tab[, 1] %in% otu_ids, ]
    },
    
    "membership_filter" = {
      if (!all(c("threshold", "min_similarities") %in% names(args))) {
        stop("For 'membership_filter', provide 'threshold' and 'min_similarities'.")
      }
      membership_filter(filtered_otu_tab,
                        threshold = args$threshold,
                        min_similarities = args$min_similarities)
    },
    
    "none" = {
      filtered_otu_tab
    },
    stop("Invalid 'filter_type'. Use 'abundance_filter', 'occupancy_filter', 'abu_occ_filter', 'hardcutoff_filter', 'membership_filter', or 'none'.")
  )
  if (nrow(selected_otus) == 0) {
    message("All OTUs were filtered out using '", filter_type, "'. Adjust parameters.")
    return(NULL)
  }
  return(selected_otus)
}

process_and_filter <- function(subset_df, filter_type, ...) {
  numeric_part <- subset_df[, -1]  
  numeric_part <- as.data.frame(lapply(numeric_part, as.numeric))  
  filtered_otu_tab <- cbind(subset_df[, 1, drop = FALSE], numeric_part)  
  result <- filter_otu_table(filtered_otu_tab, filter_type, ...)
  return(result)
}

normalize_otus <- function(selected_otus, method = "default") {
  selected_otus[, -1] <- lapply(selected_otus[, -1], function(x) as.numeric(as.character(x)))
  selected_otus[is.na(selected_otus)] <- 0
  col_sums <- colSums(selected_otus[, -1], na.rm = TRUE)
  non_zero <- col_sums != 0
  selected_otus <- selected_otus[, c(TRUE, non_zero), drop = FALSE]
  if (ncol(selected_otus) <= 1 || all(colSums(selected_otus[, -1]) == 0)) {
    stop("Error: No valid samples with counts left after filtering. Cannot normalize.")
  }
  selected_otus_1 <- NULL
  switch(method,
         "rrarefy" = {
           library(vegan)
           rarefaction_depth <- min(colSums(selected_otus[, -1]))
           rarefied_table <- vegan::rrarefy(selected_otus[, -1], sample = rarefaction_depth)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], rarefied_table)
         },
         "srs" = {
           library(SRS)
           Cmin <- min(colSums(selected_otus[, -1]))
           SRS_output <- SRS::SRS(selected_otus[, -1], Cmin)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], SRS_output)
         },
         "css" = {
           cumulative_sums <- colSums(selected_otus[, -1])
           scaling_factors <- cumulative_sums / sum(cumulative_sums)
           css_output <- selected_otus[, -1] * scaling_factors
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], css_output)
         },
         "tmm" = {
           library(edgeR)
           d <- DGEList(counts = as.matrix(selected_otus[, -1]))
           d <- calcNormFactors(d, method = "TMM")
           tmm_output <- cpm(d, normalized.lib.sizes = TRUE, log = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], tmm_output)
         },
         "tmmwsp" = {
           library(edgeR)
           d <- DGEList(counts = as.matrix(selected_otus[, -1]))
           d <- calcNormFactors(d, method = "TMMwsp")
           tmmwsp_output <- cpm(d, normalized.lib.sizes = TRUE, log = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], tmmwsp_output)
         },
         "rle" = {
           library(edgeR)
           d <- DGEList(counts = as.matrix(selected_otus[, -1]))
           d <- calcNormFactors(d, method = "RLE")
           rle_output <- cpm(d, normalized.lib.sizes = FALSE, log = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], rle_output)
         },
         "upperquartile" = {
           library(edgeR)
           d <- DGEList(counts = as.matrix(selected_otus[, -1]))
           d <- calcNormFactors(d, method = "upperquartile")
           upperquartile_output <- cpm(d, normalized.lib.sizes = FALSE, log = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], upperquartile_output)
         },
         "gmpr" = {
           require(matrixStats)
           GMPR <- function (comm, intersect.no = 10, ct.min = 1, trace = FALSE) {
             comm[comm < ct.min] <- 0
             if (is.null(colnames(comm))) {
               colnames(comm) <- paste0('S', 1:ncol(comm))
             }
             comm.no <- numeric(ncol(comm))
             gmpr <- sapply(1:ncol(comm),  function(i) {
               x <- comm[, i]
               pr <- x / comm
               pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
               incl.no <- colSums(!is.na(pr))
               pr.median <- colMedians(pr, na.rm=TRUE)
               comm.no[i] <<- sum(incl.no >= intersect.no)
               if (comm.no[i] > 1) {
                 return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
               } else {
                 return(NA)
               }
             })
             names(gmpr) <- names(comm.no) <- colnames(comm)
             attr(gmpr, 'NSS') <- comm.no
             return(gmpr)
           }
           size_factors <- GMPR(t(selected_otus[, -1]))
           normalized_data <- t(t(selected_otus[, -1]) / size_factors)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], normalized_data)
         },
         "tss" = {
           otu_norm <- as.matrix(selected_otus[, -1])
           otu_TSS <- vegan::decostand(otu_norm, method = "total", MARGIN = 1, na.rm = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], as.data.frame(otu_TSS))
         },
         "clr" = {
           library(chemometrics)
           otu_norm <- as.matrix(selected_otus[, -1])
           otu_CLR <- chemometrics::clr(otu_norm + 1)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], as.data.frame(otu_CLR))
         },
         "log" = {
           otu_norm <- as.matrix(selected_otus[, -1])
           otu_norm[otu_norm == 0] <- min(otu_norm[otu_norm > 0]) / 2
           otu_LOG <- log2(otu_norm)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], as.data.frame(otu_LOG))
         },
         "invlogit" = {
           otu_norm <- as.matrix(selected_otus[, -1])
           otu_INVLOGIT <- stats::plogis(otu_norm)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], as.data.frame(otu_INVLOGIT))
         },
         "none" = {
           selected_otus_1 <- selected_otus  
         },
         "default" = {
           library(edgeR)
           d <- DGEList(counts = as.matrix(selected_otus[, -1]))
           d <- calcNormFactors(d, method = "TMM")
           tmm_output <- cpm(d, normalized.lib.sizes = TRUE)
           selected_otus_1 <- cbind(selected_otus[, 1, drop = FALSE], tmm_output)
         }
  )
  colnames(selected_otus_1) <- colnames(selected_otus)
  return(selected_otus_1)
}

process_and_normalize <- function(data_list, filter_type = "tmm", ...) {
  processed_data <- lapply(data_list, function(df) {
    if (ncol(df) <= 2) {
      warning("Skipping normalization: less than 2 samples in group.")
      return(df)
    }
    otu_ids <- df[[1]]
    abundance_data <- df
    abundance_data <- as.data.frame(lapply(abundance_data, function(x) as.numeric(as.character(x))))
    abundance_data[is.na(abundance_data)] <- 0
    col_sums <- colSums(abundance_data, na.rm = TRUE)
    abundance_data <- abundance_data[, col_sums != 0, drop = FALSE]
    df_clean <- cbind(OTU_ID = otu_ids, abundance_data)
    normalized_df <- tryCatch({
      normalize_otus(df_clean, method = filter_type, ...)
    }, error = function(e) {
      warning(paste("Normalization failed:", e$message))
      return(df_clean)
    })
    return(normalized_df)
  })
  return(processed_data)
}

calculate_alpha_diversity <- function(otu_table) {
  tryCatch({
    transposed_otu <- t(otu_table[, -1])
    richness <- vegan::specnumber(transposed_otu)
    shannon <- vegan::diversity(transposed_otu, index = "shannon")
    simpson <- vegan::diversity(transposed_otu, index = "simpson")
    evenness <- shannon / log(richness)
    
    format_to_3dp <- function(x) {
      x_formatted <- sprintf("%.3f", x)
      return(as.numeric(x_formatted))
    }
    shannon_formatted <- sapply(shannon, format_to_3dp)
    simpson_formatted <- sapply(simpson, format_to_3dp)
    evenness_formatted <- sapply(evenness, format_to_3dp)
    alpha_diversity <- rbind(richness, evenness_formatted, shannon_formatted, simpson_formatted)
    rownames(alpha_diversity) <- c("Richness", "Evenness", "Shannon", "Simpson")
    
    return(alpha_diversity)
  }, error = function(e) {
    message("Error in alpha diversity calculation: ", e$message)
    return(NULL)
  }, warning = function(w) {
    message("Warning in alpha diversity calculation: ", w$message)
    return(NULL)
  })
}

prepare_otu_table <- function(df) {
  df <- as.data.frame(df)
  if ("OTU_ID" %in% colnames(df)) {
    return(df %>% 
             select(OTU_ID, where(is.numeric)) %>%
             mutate(OTU_ID = as.character(OTU_ID)))
  }
  if (!is.null(rownames(df))) {
    return(df %>% 
             tibble::rownames_to_column("OTU_ID") %>% 
             select(OTU_ID, where(is.numeric)) %>%
             mutate(OTU_ID = as.character(OTU_ID)))
  }
  if (all(sapply(df, is.numeric))) {
    return(df %>% 
             t() %>% 
             as.data.frame() %>% 
             tibble::rownames_to_column("OTU_ID") %>% 
             select(OTU_ID, where(is.numeric)) %>%
             mutate(OTU_ID = as.character(OTU_ID)))
  }
  stop("Cannot determine OTU IDs. Data should have either:
       1. An 'OTU_ID' column, or
       2. OTUs as row names, or
       3. All numeric columns with samples as rows")
}

calculate_beta_distance <- function(group1, group2, method = "bray") {
  all_otus <- unique(c(group1$OTU_ID, group2$OTU_ID))
  otu_ref <- data.frame(OTU_ID = all_otus)
  join_with_ref <- function(df) {
    otu_ref %>% 
      left_join(df, by = "OTU_ID") %>% 
      mutate(across(-OTU_ID, ~replace_na(., 0))) %>% 
      column_to_rownames("OTU_ID") %>% 
      t()
  }
  mat1 <- join_with_ref(group1)
  mat2 <- join_with_ref(group2)
  combined <- rbind(mat1, mat2)
  dist_matrix <- as.matrix(vegdist(combined, method = method))
  dist_matrix[rownames(mat1), rownames(mat2)]
}

compute_all_beta_comparisons <- function(data_list, method = "bray") {
  prepared_data <- lapply(data_list, prepare_otu_table)
  group_names <- names(prepared_data)
  if (is.null(group_names)) {
    group_names <- paste0("Group", seq_along(prepared_data))
  }
  comparisons <- combn(length(prepared_data), 2, simplify = FALSE)
  results <- list()
  for (pair in comparisons) {
    i <- pair[1]
    j <- pair[2]
    comp_name <- paste(group_names[i], group_names[j], sep = "_vs_")
    results[[comp_name]] <- calculate_beta_distance(
      prepared_data[[i]],
      prepared_data[[j]],
      method = method
    )
  }
  return(results)
}

compute_alpha_split <- function(data_list) {
  alpha_results <- lapply(data_list, function(df) {
    df[, -1] <- lapply(df[, -1], function(x) as.numeric(as.character(x)))
    df[is.na(df)] <- 0
    calculate_alpha_diversity(df)
  })
  return(alpha_results)
}

create_groupwise_boxplots <- function(group_data, 
                                      log_transform = FALSE, 
                                      show_points = TRUE, 
                                      num_taxa = 10, 
                                      show_all_taxa = FALSE,
                                      color_scheme = "Default") {
  if(is.null(group_data) || ncol(group_data) <= 1) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = "No data available for this group") + 
             theme_void())
  }
  sample_ids <- group_data[, 1]
  taxa_names <- colnames(group_data)[-1]
  if(!show_all_taxa && length(taxa_names) > num_taxa) {
    taxa_means <- colMeans(group_data[, -1, drop = FALSE], na.rm = TRUE)
    top_taxa <- names(sort(taxa_means, decreasing = TRUE))[1:num_taxa]
    selected_cols <- c(1, which(colnames(group_data) %in% top_taxa))
    plot_data <- group_data[, selected_cols]
    taxa_names <- colnames(plot_data)[-1]
  } else {
    plot_data <- group_data
  }
  plot_data_long <- reshape2::melt(plot_data, id.vars = colnames(plot_data)[1],
                                   variable.name = "Sample", value.name = "Abundance")
  plot_data_long$Abundance[is.na(plot_data_long$Abundance)] <- 0
  if(log_transform) {
    # Add small value to avoid log(0)
    plot_data_long$Abundance <- log10(plot_data_long$Abundance + 1e-5)
    abundance_label <- "Log10 Abundance"
  } else {
    abundance_label <- "Abundance"
  }
  n_taxa <- length(unique(plot_data_long$Sample))
  color_values <- switch(color_scheme,
                         "Viridis" = viridis::viridis(n_taxa),
                         "Rainbow" = rainbow(n_taxa),
                         "Blues" = colorRampPalette(c("lightblue", "darkblue"))(n_taxa),
                         "Reds" = colorRampPalette(c("lightsalmon", "darkred"))(n_taxa),
                         # Default
                         colorRampPalette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", 
                                            "#FF7F00", "#FFFF33", "#A65628"))(n_taxa))

  p <- ggplot(plot_data_long, aes(x = Sample, y = Abundance, fill = Sample)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = color_values) +
    labs(
      x = "Sample",
      y = abundance_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none"
    )
  if(show_points) {
    p <- p + geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "black")
  }
  
  return(p)
}

plot_normalized_density <- function(normalized_data, title = "Density Plot of Normalized OTU Data") {
  if (!require(ggplot2, quietly = TRUE)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  values <- as.numeric(as.matrix(normalized_data))
  values <- values[is.finite(values)]
  if (length(values) == 0) stop("No numeric data to plot.")
  ggplot(data.frame(value = values), aes(x = value)) +
    geom_density(fill = "blue", alpha = 0.5, color = "black") +
    labs(title = title, x = "Normalized Abundance", y = "Density") +
    theme_minimal() +
    scale_x_continuous(trans = "log10",  
                       labels = scales::scientific)  
}

plot_normalized_density_split <- function(normalized_split_data, 
                                          title = "Density Plot of Normalized OTU Data",
                                          group_names = NULL) {
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(reshape2)) install.packages("reshape2") 
  library(ggplot2)
  library(reshape2)
  if (!is.list(normalized_split_data)) {
    stop("Input must be a list of normalized data frames/matrices.")
  }
  plot_data <- data.frame()
  for (i in seq_along(normalized_split_data)) {
    group_data <- normalized_split_data[[i]]
    values <- as.numeric(as.matrix(group_data))
    values <- values[is.finite(values)] 
    if (length(values) > 0) {
      group_name <- if (!is.null(group_names)) group_names[i] else paste("Group", i)
      plot_data <- rbind(plot_data, 
                         data.frame(value = values, 
                                    group = factor(group_name)))
    }
  }
  
  if (nrow(plot_data) == 0) stop("No valid numeric data to plot.")
  ggplot(plot_data, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5, color = NA) + 
    labs(title = title, 
         x = "Normalized Abundance (log10)", 
         y = "Density",
         fill = "Group") +
    theme_minimal() +
    scale_x_continuous(trans = "log10",
                       labels = scales::scientific) +
    scale_fill_viridis_d()
}

filter_top_otus_split <- function(split_selected_otus, top_percentage = 10) {
  lapply(split_selected_otus, function(selected_otus) {
    row_sums <- rowSums(selected_otus[, -1])
    otu_summary <- data.frame(OTU_ID = rownames(selected_otus), RowSum = row_sums)
    sorted_otu_summary <- otu_summary[order(otu_summary$RowSum, decreasing = TRUE), ]
    num_rows <- ceiling(nrow(sorted_otu_summary) * (top_percentage / 100))
    top_core_rows <- selected_otus[sorted_otu_summary$OTU_ID[1:num_rows], , drop = FALSE]
    non_core_rows <- selected_otus[sorted_otu_summary$OTU_ID[(num_rows + 1):nrow(sorted_otu_summary)], , drop = FALSE]
    return(list(core_otus = top_core_rows, non_core_otus = non_core_rows))
  })
}

filter_top_otus <- function(selected_otus, top_percentage = 10) {
  row_sums <- rowSums(selected_otus[, -1])
  otu_summary <- data.frame(OTU_ID = rownames(selected_otus), RowSum = row_sums)
  sorted_otu_summary <- otu_summary[order(otu_summary$RowSum, decreasing = TRUE), ]
  num_rows <- ceiling(nrow(sorted_otu_summary) * (top_percentage / 100))
  top_core_rows <- selected_otus[sorted_otu_summary$OTU_ID[1:num_rows], , drop = FALSE]
  non_core_rows <- selected_otus[sorted_otu_summary$OTU_ID[(num_rows + 1):nrow(sorted_otu_summary)], , drop = FALSE]
  return(list(core_otus = top_core_rows, non_core_otus = non_core_rows))
}

multidimensionality_plot <- function(normalized_data) {
  otu_matrix <- t(normalized_data[, -1]) 
  colnames(otu_matrix) <- normalized_data[, 1] 
  dist_matrix <- vegdist(otu_matrix, method = "bray")
  mds_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  mds_points <- as.data.frame(mds_result$points)
  colnames(mds_points) <- c("MDS1", "MDS2")
  mds_points$Sample <- rownames(mds_points)
  p <- ggplot(mds_points, aes(x = MDS1, y = MDS2, text = Sample)) +
    geom_point(size = 3, alpha = 0.7, color = "blue") +
    theme_minimal() +
    ggtitle("Interactive MDS Plot (Click on points to see sample names)") +
    xlab("MDS Axis 1") +
    ylab("MDS Axis 2")
  ggplotly(p, tooltip = "text")
}

multidimensionality_split_plot <- function(normalized_split_data) {
  normalized_split_data <- lapply(normalized_split_data, function(df) {
    if (!"OTU_ID" %in% colnames(df)) {
      colnames(df)[1] <- "OTU_ID"
    }
    df
  })

  merged_data <- Reduce(function(x, y) merge(x, y, by = "OTU_ID", all = TRUE), normalized_split_data)
  merged_data[is.na(merged_data)] <- 0  
  otu_matrix <- t(merged_data[, -1])  
  colnames(otu_matrix) <- merged_data$OTU_ID  
  sample_names <- rownames(otu_matrix)
  group_labels <- rep(names(normalized_split_data), sapply(normalized_split_data, function(x) ncol(x) - 1))
  dist_matrix <- vegdist(otu_matrix, method = "bray")
  mds_result <- cmdscale(dist_matrix, eig = TRUE, k = 2)
  mds_points <- as.data.frame(mds_result$points)
  colnames(mds_points) <- c("MDS1", "MDS2")
  mds_points$Sample <- sample_names
  mds_points$Group <- group_labels
  mds_points$HoverText <- paste("Sample:", mds_points$Sample, "<br>Group:", mds_points$Group)
  p <- ggplot(mds_points, aes(x = MDS1, y = MDS2, color = Group, text = HoverText)) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal() +
    ggtitle("MDS Plot of Samples by Group") +
    xlab("MDS Axis 1") +
    ylab("MDS Axis 2") +
    theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
  
  ggplotly(p, tooltip = "text")
}

multidimensionality_plot2 <- function(normalized_data, metadata) {
  if(ncol(normalized_data) < 3) {
    stop("Insufficient data: Need at least 2 samples (columns) for MDS analysis")
  }
  t_normalized <- as.data.frame(t(normalized_data[,-1])) 
  colnames(t_normalized) <- normalized_data[,1]  
  otu_matrix <- as.matrix(t_normalized)
  mode(otu_matrix) <- "numeric"
  if(any(is.na(otu_matrix))) {
    warning("NA values detected in OTU matrix - consider imputation or removal")
    otu_matrix[is.na(otu_matrix)] <- 0  
  }
  dist_matrix <- vegdist(otu_matrix, method = "bray")
  n_samples <- nrow(otu_matrix)
  k <- min(2, n_samples - 1)
  if(k < 1) {
    stop("Cannot perform MDS: Need at least 2 samples for dimensionality reduction")
  }
  mds_result <- cmdscale(dist_matrix, eig = TRUE, k = k)
  mds_points <- as.data.frame(mds_result$points)
  colnames(mds_points) <- paste0("MDS", 1:k)
  mds_points$Sample_ID <- rownames(mds_points)
  if(!"Sample_ID" %in% colnames(metadata)) {
    possible_id_cols <- colnames(metadata)[sapply(metadata, function(x) any(x %in% mds_points$Sample_ID))]
    if(length(possible_id_cols) > 0) {
      metadata <- metadata %>% rename(Sample_ID = !!possible_id_cols[1])
    } else {
      stop("Could not find matching sample IDs between data and metadata")
    }
  }
  if(!"Treatment" %in% colnames(metadata)) {
    if(length(colnames(metadata)) > 1) {
      treatment_col <- colnames(metadata)[!colnames(metadata) %in% "Sample_ID"][1]
      metadata <- metadata %>% rename(Treatment = !!treatment_col)
    } else {
      metadata$Treatment <- "Group 1" 
    }
  }
  if(!"Soil.Type" %in% colnames(metadata)) {
    metadata$Soil.Type <- metadata$Treatment  
  }
  plot_data <- merge(mds_points, metadata, by = "Sample_ID")
  p <- ggplot(plot_data, aes(x = MDS1, 
                             y = if(k > 1) MDS2 else 0,
                             color = Treatment,
                             text = paste0("Sample: ", Sample_ID,
                                           "\nGroup: ", Treatment))) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(title = "MDS Plot with Grouping",
         x = "MDS Axis 1",
         y = if(k > 1) "MDS Axis 2" else "") +
    scale_color_viridis_d()
  ggplotly(p, tooltip = "text") %>%
    layout(hoverlabel = list(bgcolor = "white"),
           legend = list(orientation = "h", x = 0.5, y = -0.1))
}

get_split_core_results <- function(split_result, abundance_split_data, tax_table) {
  split_core_otus <- lapply(split_result, function(res) res$core_otus)
  colnames(tax_table)[1] <- "OTU_ID"
  split_core_otu_counts <- mapply(function(core_otus, selected_otus) {
    colnames(core_otus)[1] <- "OTU_ID"
    match_core <- fastmatch::fmatch(core_otus$OTU_ID, selected_otus[, 1])
    selected_otus[match_core, , drop = FALSE]
  }, split_core_otus, abundance_split_data, SIMPLIFY = FALSE)
  split_common_otu <- lapply(split_core_otu_counts, function(core_otus) {
    colnames(core_otus)[2] <- "OTU_ID"
    tax_table[tax_table$OTU_ID %in% core_otus$OTU_ID, ]
  })
  return(list(
    Core_taxonomy = split_common_otu,
    Split_Core_OTU_counts = split_core_otu_counts
  ))
}

get_core_otu_data <- function(result, abundance_result, tax_table) {
  core_otus <- result$core_otus
  non_core_otus <- result$non_core_otus
  colnames(core_otus)[1] <- "OTU_ID"
  colnames(abundance_result)[1] <- "OTU_ID"
  colnames(tax_table)[1] <- "OTU_ID"
  match_core <- fastmatch::fmatch(core_otus$OTU_ID, abundance_result$OTU_ID)
  core_otu_counts <- abundance_result[match_core, , drop = FALSE]
  core_taxonomy <- tax_table[tax_table$OTU_ID %in% core_otus$OTU_ID, , drop = FALSE]
  return(list(
    Core_OTU_Counts = core_otu_counts,
    Core_Taxonomy = core_taxonomy,
    Core_OTUs = core_otus,
    Non_Core_OTUs = non_core_otus
  ))
}

compute_combined_beta_diversity <- function(data_list, method = "bray") {
  if (!is.list(data_list)) {
    stop("Input must be a list of data frames")
  }
  if (length(data_list) == 0) {
    stop("Input list is empty")
  }
  if (is.null(names(data_list))) {
    names(data_list) <- paste0("Group_", seq_along(data_list))
  }
  identify_id_column <- function(df) {
    if ("OTU_ID" %in% colnames(df)) return("OTU_ID")
    if ("Feature_ID" %in% colnames(df)) return("Feature_ID")
    if ("ASV_ID" %in% colnames(df)) return("ASV_ID")
    if (any(grepl("ID$", colnames(df)))) {
      return(colnames(df)[grepl("ID$", colnames(df))][1])
    }
    stop("No suitable ID column found. Expected column named OTU_ID, Feature_ID, or similar")
  }
  tryCatch({
    combined_data <- bind_rows(
      lapply(names(data_list), function(group_name) {
        df <- data_list[[group_name]]
        id_col <- identify_id_column(df)
        if (id_col != "OTU_ID") {
          message(paste("Using column", id_col, "as ID for group", group_name))
          df <- df %>% rename(OTU_ID = !!sym(id_col))
        }
        df[, -which(colnames(df) == "OTU_ID")] <- lapply(
          df[, -which(colnames(df) == "OTU_ID"), drop = FALSE], 
          function(x) as.numeric(as.character(x)))
        df[is.na(df)] <- 0
        df %>%
          pivot_longer(
            cols = -OTU_ID,
            names_to = "Sample",
            values_to = "Count"
          ) %>%
          mutate(Group = group_name)
      })
    )
    if (length(unique(combined_data$Sample)) < 2) {
      stop(paste("Insufficient samples (need at least 2). Found:", 
                 length(unique(combined_data$Sample))))
    }
    wide_data <- combined_data %>%
      select(-Group) %>%  
      pivot_wider(
        names_from = OTU_ID,
        values_from = Count,
        values_fill = 0
      ) %>%
      column_to_rownames("Sample") %>%
      as.matrix()
    if (!method %in% c("bray", "jaccard", "euclidean")) {
      method <- "bray"  
    }
    
    dist_matrix <- as.matrix(vegdist(wide_data, method = method))
    return(list(
      distance_matrix = dist_matrix,
      method = method,
      sample_groups = combined_data %>%
        select(Sample, Group) %>%
        distinct(),
      combined_data = combined_data
    ))
  }, error = function(e) {
    message("Detailed error in compute_combined_beta_diversity:")
    message(e$message)
    stop("Failed to calculate combined beta diversity. ", 
         "Please check: 1) Your data contains an ID column (OTU_ID, Feature_ID, etc.), ",
         "2) All count columns are numeric, ",
         "3) You have at least 2 samples across all groups")
  })
}