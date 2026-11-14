DA_fun <- function(filtered_data1, 
                   md_table, 
                   OTU_catVar, 
                   Variable, 
                   contrast_OTU,
                   OTU_norm.type, 
                   test, 
                   alpha, 
                   lfc, 
                   p.adj.method,
                   other_variable = NULL,
                   other_variable2 = NULL,
                   other_variable3 = NULL) {
  OTU_ID <- as.character(filtered_data1[[1]]) 
  filtered_data <- filtered_data1[, -1]
  rownames(filtered_data) <- OTU_ID 
  nf1 <- nrow(filtered_data1)
  metadata <- md_table[match(colnames(filtered_data), md_table$Sample), ]
  meta.col.names <- colnames(metadata)
  OTU_id1 <- which(meta.col.names == Variable)
  OTU_grp1 <- as.character(metadata[[OTU_id1]])
  OTU_grp.lev <- unique(OTU_grp1)
  OTU_group <- factor(OTU_grp1, levels = OTU_grp.lev)
  z <- DGEList(counts = filtered_data, group = OTU_group, genes = OTU_ID)
  nOTU <- nrow(z$counts)
  z <- calcNormFactors(z, method = OTU_norm.type)
  sample_names <- as.character(metadata[[which(meta.col.names == "Sample")]])
  ls <- colSums(filtered_data) / 1e6
  df.ls <- data.frame(SampleID = sample_names, LibrarySize = ls)
  df.sub <- df.ls[1:min(25, nrow(df.ls)), ]
  bp_OTU_ls <- ggplot(df.sub, aes(x = SampleID, y = LibrarySize, fill = SampleID)) +
    geom_bar(stat = "identity") +
    ylab("Library Size (in millions) after filtering") +
    ggtitle(paste("Bar plot of library Size (", contrast_OTU, ")")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90))
  mds_plot1 <- plotMDS(z, top = nOTU, labels = colnames(z$counts))
  mdsR_plotly <- function(plot_MDS, variableInterest, Groups) {
    data.frame(
      X = plot_MDS$x, 
      Y = plot_MDS$y, 
      CategoryLevel = Groups,
      Sample = colnames(plot_MDS$distance.matrix.squared)
    ) %>%
      plot_ly(
        x = ~X, y = ~Y, color = ~CategoryLevel,
        text = ~paste("Sample:", Sample, "<br>Group:", CategoryLevel),
        hoverinfo = "text", type = "scatter", mode = "markers",
        marker = list(size = 10, opacity = 0.8)
      ) %>%
      layout(
        title = paste("MDS plot (", contrast_OTU, ")"),
        xaxis = list(title = "Leading logFC dim1"),
        yaxis = list(title = "Leading logFC dim2"),
        legend = list(title = list(text = variableInterest))
      )
  }
  mds_plot_interactive <- mdsR_plotly(mds_plot1, Variable, OTU_group)
  OTU_id2 <- match(OTU_catVar, meta.col.names)
  metadata[OTU_id2] <- lapply(metadata[OTU_id2], function(x) factor(x))
  single_level_vars <- sapply(metadata[, OTU_id2, drop = FALSE], function(x) length(unique(x)) < 2)
  if (any(single_level_vars)) {
    warning("Dropping single-level variables: ", names(single_level_vars)[single_level_vars])
    OTU_catVar <- OTU_catVar[!OTU_catVar %in% names(single_level_vars)[single_level_vars]]
  }
  covariates <- unique(c(Variable, OTU_catVar, other_variable, other_variable2, other_variable3)) 
  create_stable_design <- function(metadata, covariates) {
    md_clean <- metadata
    for(cov in covariates) {
      if(n_distinct(md_clean[[cov]]) < 2) {
        warning("Removing ", cov, " - only one level present")
        covariates <- setdiff(covariates, cov)
        next
      }
      md_clean[[cov]] <- fct_drop(factor(md_clean[[cov]]))
    }
    design <- matrix(1, nrow = nrow(md_clean))
    colnames(design) <- "(Intercept)"
    for(cov in covariates) {
      temp_design <- tryCatch(
        model.matrix(as.formula(paste("~", cov)), data = md_clean),
        error = function(e) NULL
      )
      if(!is.null(temp_design)) {
        combined <- cbind(design, temp_design[,-1, drop=FALSE])
        if(qr(combined)$rank == ncol(combined)) {
          design <- combined
        }
      }
    }
    return(design)
  }
  OTU_design <- create_stable_design(metadata, covariates)
  z <- estimateDisp(z, OTU_design)
  A <- aveLogCPM(z$counts, offset = getOffset(z))
  disp <- sqrt(getDispersion(z))
  df <- data.frame(A, disp)
  if (!is.null(z$tagwise.dispersion)) df$tagwise <- sqrt(z$tagwise.dispersion)
  if (!is.null(z$common.dispersion)) df$common <- sqrt(z$common.dispersion)
  if (!is.null(z$trended.dispersion)) df$trended <- sqrt(z$trended.dispersion)
  OTU_bcv_plot <- ggplot(df, aes(A, disp)) +
    geom_point(alpha = 0.1) +
    geom_point(data = df, aes(A, tagwise), size = 1, alpha = 0.5) +
    geom_line(data = df, aes(A, common), color = "red") +
    geom_line(data = df, aes(A, trended), color = "blue") +
    xlab("Average log CPM") + ylab("Biological coefficient of variation") +
    ggtitle(paste("BCV plot (", contrast_OTU, ")")) + theme_bw(base_size = 10)
  OTU_pair <- strsplit(contrast_OTU, " vs. ")[[1]]
  OTU_det.id <- match(OTU_pair, OTU_grp.lev)
  OTU_coln1 <- paste(Variable, OTU_pair, sep = "")
  OTU_idc1 <- which(colnames(OTU_design) == OTU_coln1[1])
  OTU_idc2 <- which(colnames(OTU_design) == OTU_coln1[2])
  lfm <- 1
  if (test == "Exact") {
    de.test <- exactTest(z, pair = OTU_det.id)
  } else {
    contrast1 <- rep(0, ncol(OTU_design))
    if (length(OTU_idc1) > 0 && length(OTU_idc2) > 0) {
      contrast1[OTU_idc1] <- 1
      contrast1[OTU_idc2] <- -1
      cf1 <- NULL
    } else {
      cf1 <- ifelse(length(OTU_idc1) > 0, OTU_idc1, OTU_idc2)
      if (length(OTU_idc2) > 0) lfm <- -1
      contrast1 <- NULL
    }
    if (test == "glmQLF") {
      fit <- glmQLFit(z, design = OTU_design)
      de.test <- if (!is.null(cf1)) glmQLFTest(fit, coef = cf1) else glmQLFTest(fit, contrast = contrast1)
    } else if (test == "glmLRT") {
      fit <- glmFit(z, design = OTU_design)
      de.test <- if (!is.null(cf1)) glmLRT(fit, coef = cf1) else glmLRT(fit, contrast = contrast1)
    } else {
      stop("Invalid test method")
    }
  }
  res1 <- topTags(de.test, n = nOTU, adjust.method = p.adj.method)
  result1 <- res1$table
  colnames(result1)[1] <- "OTU ID"
  result1$logFC <- lfm * result1$logFC
  i.deg1 <- which(result1$PValue < alpha)
  i.deg2 <- which(result1$FDR < alpha)
  id.up1 <- which(result1$PValue < alpha & result1$logFC >= lfc)
  id.down1 <- which(result1$PValue < alpha & result1$logFC <= -lfc)
  id.up2 <- which(result1$FDR < alpha & result1$logFC >= lfc)
  id.down2 <- which(result1$FDR < alpha & result1$logFC <= -lfc)
  result.nm1 <- c("No. of OTU in the data", 
                  "No. of DA OTU", paste("Up-regulated OTU (logFC ≥", lfc, ")"),
                  paste("Down-regulated OTU (logFC ≤", -lfc, ")"),
                  paste("DA OTU between ±", lfc))
  
  result.summary <- cbind(c(nf1, length(i.deg1), length(id.up1), length(id.down1), sum(result1$PValue < alpha & abs(result1$logFC) < lfc)),
                          c(nf1, length(i.deg2), length(id.up2), length(id.down2), sum(result1$FDR < alpha & abs(result1$logFC) < lfc)))
  colnames(result.summary) <- c("Not adjusted", "Adjusted")
  rownames(result.summary) <- result.nm1
  
  inputs.summary <- data.frame(
    `Input parameters` = c("Categorical fixed effect(s)", 
                           "Categorical fixed effect variable", 
                           "Pairwise comparison", 
                           "Normalization type", 
                           "Method of testing", 
                           "Significance threshold", 
                           "Log fold change cutoff", 
                           "P-value adjustment"),
    `Inputs selected` = c(
      paste(OTU_catVar, collapse = ", "),
      Variable, contrast_OTU, OTU_norm.type, test,
      alpha, paste("±", lfc), p.adj.method
    )
  )
  rownames(inputs.summary) <- inputs.summary[, 1]
  
  list(inputs_selected = inputs.summary,
       barplot = bp_OTU_ls,
       mds_plot = mds_plot_interactive,
       bcv_plot = OTU_bcv_plot,
       summary_result = result.summary,
       table_result = result1)
}

md_plot_shiny <- function(id.nondeg, id.deg, id.up, id.down, c.lfc, lfc, av.lcpm, alpha, geneID, title) {
  fill.fr <- c("red", "blue", "black", "green")
  labels.fr = c("Up-regulated", "Down-regulated", "Non-significant", "Significant between logFC cut-offs")
  
  df.de <- cbind.data.frame(geneID, c.lfc, av.lcpm, labels.fr[4], fill.fr[4])
  colnames(df.de) <- c("Gene ID", "log fold change", "Average log CPM", "Status", "Fill")
  df.de$Status <- labels.fr[4]
  df.de$Fill <- fill.fr[4]
  if(length(id.up) > 0) {
    df.de$Status[id.up] <- labels.fr[1]
    df.de$Fill[id.up] <- fill.fr[1]
  }
  if(length(id.down) > 0) {
    df.de$Status[id.down] <- labels.fr[2]
    df.de$Fill[id.down] <- fill.fr[2]
  }
  if(length(id.nondeg) > 0) {
    df.de$Status[id.nondeg] <- labels.fr[3]
    df.de$Fill[id.nondeg] <- fill.fr[3]
  }
  Status.label <- unique(df.de$Status)
  Fill.label <- unique(df.de$Fill)
  df.de$Status <- factor(df.de$Status, levels = Status.label)
  df.de$Fill <- factor(df.de$Fill, levels = Fill.label)
  yl <- "log-fold-change"
  xl <- "Average log CPM"
  ggplot(df.de, aes(x = `Average log CPM`, y = `log fold change`, color = Status)) + 
    geom_point(size = 2, alpha = 0.8) + 
    scale_color_manual(values = Fill.label) + 
    geom_hline(yintercept = c(-lfc, lfc), colour = "orange", linetype = 2, lwd = 1) +
    ggtitle(paste("MD Plot (", title, ")")) +
    xlab(xl) + ylab(yl) +
    theme(legend.position = "bottom") +
    theme_bw(base_size = 10)
}
