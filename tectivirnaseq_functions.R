# Hierarchical clustering

num_to_letters <- function(n) {
  stopifnot(all(n >= 1))
  sapply(n, function(x) {
    res <- character(0)
    while (x > 0) {
      x <- x - 1
      res <- c(letters[(x %% 26) + 1], res)
      x <- x %/% 26
    }
    paste0(res, collapse = "")
  })
}

compute_stability <- function(clusters, strain, frac = 0.80, B = 100, seed = 0) {
  
  jaccard_mat <- matrix(0, nrow = length(clusters), ncol = B)

  set.seed(seed)
  
  for (b in 1:B) {
    
    keep <- sample(rownames(counts),
                   size = floor(frac * nrow(counts)),
                   replace = FALSE)
    
    counts_sub <- counts[keep, ]
    
    rhos_sub <- suppressMessages(propr(t(counts_sub), metric = "rho"))
    d_sub <- as.dist(1 - rhos_sub@matrix)
    hc_sub <- hclust(d_sub, method = "average")
    clusters_sub <- cutree(hc_sub, h = 0.25)
    
    for (i in seq_along(clusters)) {
      
      original_genes <- clusters_list[[i]]
      original_genes <- intersect(original_genes, keep)
      
      if (length(original_genes) < 5) next
      
      tab <- table(clusters_sub[original_genes])
      best_cluster <- names(which.max(tab))
      
      recovered_genes <- names(clusters_sub)[clusters_sub == best_cluster]
      
      inter <- length(intersect(original_genes, recovered_genes))
      union <- length(union(original_genes, recovered_genes))
      
      jaccard_mat[i, b] <- inter / union
    }
  }
  
  stability <- rowMeans(jaccard_mat, na.rm = TRUE)
  
  cluster_names <- sapply(clusters, function(id) {
    if (!is.na(cluster_to_name[as.character(id)])) {
      cluster_to_name[as.character(id)]
    } else {
      paste0(strain, "_SC", sprintf("%04d", id))
    }
  })
  
  stab_big <- data.frame(
    Cluster = cluster_names,
    Size = as.numeric(cluster_sizes[as.character(clusters)]),
    Stability = stability
  )
  
  stab_big <- stab_big[order(-stab_big$Stability), ]
  
  rownames(stab_big) <- NULL
  
  return(stab_big)
}

clusters_expression <- function(expression_matrix, metadata, time_points,
                                clusters = NULL, order = NULL,
                                de_genes = NULL, which_contrast = NULL) {

  logratio_long <- expression_matrix %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(cols = -Sample, names_to = "Gene", values_to = "CLR")

  meta_df <- metadata %>%
    dplyr::select(Sample, Time_point) %>%
    filter(Sample %in% rownames(expression_matrix)) %>%
    mutate(Time_point = as.numeric(as.character(Time_point))) %>%
    filter(Time_point %in% time_points)

  clu_map <- cluster_gene_table %>%
    dplyr::select(Gene, Cluster) %>%
    distinct()

  dat <- logratio_long %>%
    left_join(meta_df, by = "Sample") %>%
    left_join(clu_map, by = "Gene") %>%
    filter(!is.na(Cluster), !is.na(Time_point))

  if (!is.null(clusters)) {
    dat <- dat %>% filter(Cluster %in% clusters)
  }

  if (is.null(order)) order <- order_x

  dat <- dat %>% mutate(Cluster = factor(Cluster, levels = order))

  dat_gene <- dat %>%
    group_by(Gene, Cluster, Time_point) %>%
    summarise(gene_mean_CLR = mean(CLR), .groups = "drop")

  summary_dat <- dat_gene %>%
    group_by(Cluster, Time_point) %>%
    summarise(
      mean_CLR = mean(gene_mean_CLR),
      sd_CLR   = sd(gene_mean_CLR),
      n        = dplyr::n(),
      se       = sd_CLR / sqrt(n),
      t_crit   = qt(0.975, df = n - 1),
      ci95     = t_crit * se,
      .groups  = "drop"
    )

  cluster_colors <- NULL

  if (!is.null(de_genes) && !is.null(which_contrast)) {
    de_object <- readRDS(paste0(de_genes[which_contrast], "_de_genes.rds"))

    up_genes <- de_object$up_genes$Gene_ID
    down_genes <- de_object$down_genes$Gene_ID

    cluster_colors <- dat %>%
      distinct(Gene, Cluster) %>%
      group_by(Cluster) %>%
      summarise(
        n_total = dplyr::n(),
        n_up = sum(Gene %in% up_genes),
        n_down = sum(Gene %in% down_genes),
        de_score = (n_up - n_down) / n_total,
        .groups = "drop"
      )

    summary_dat <- summary_dat %>%
      left_join(cluster_colors, by = "Cluster")
  }

  label_dat <- summary_dat %>%
    group_by(Cluster) %>%
    summarise(
      n_up = unique(n_up),
      n_down = unique(n_down),
      de_score = unique(de_score),
      y_max = max(mean_CLR + ci95),
      y_min = min(mean_CLR - ci95),
      y_at_x = mean_CLR[Time_point == ((min(time_points) + max(time_points)) / 2)][1],
      .groups = "drop"
    ) %>%
    mutate(
      x_pos = (min(time_points) + max(time_points)) / 2,
      y_range = y_max - y_min,
      y_mid = (y_max + y_min) / 2,
      y_pos = ifelse(
        y_at_x > y_mid,
        y_min + 0.1 * y_range,
        y_max - 0.1 * y_range
      ),
      up_label = paste0("UP: ", n_up),
      down_label = paste0("DOWN: ", n_down)
    )

  label_clusters_all <- summary_dat %>%
    filter(abs(de_score) > 0.5) %>%
    group_by(Cluster) %>%
    filter(Time_point == max(Time_point)) %>%
    ungroup()
  
  all_clusters <- ggplot(summary_dat,
                         aes(Time_point, mean_CLR, group = Cluster, colour = de_score, fill = de_score)) +
    geom_ribbon(aes(ymin = mean_CLR - ci95, ymax = mean_CLR + ci95),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_text(
      data = label_clusters_all,
      aes(Time_point, mean_CLR, label = Cluster, colour = de_score),
      hjust = -0.1,
      size = 3,
      show.legend = FALSE
    ) +
    scale_colour_gradient2(
      name = "DE score",
      low = "blue", mid = "grey80", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    scale_fill_gradient2(
      name = "DE score",
      low = "blue", mid = "grey80", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    scale_x_continuous(
      breaks = time_points,
      limits = c(min(time_points) - 5, max(time_points) + 5),
      expand = c(0, 0)
    ) +
    labs(x = "Time points (minutes)", y = "CLR counts") +
    theme_bw(base_size = 11) +
    theme(legend.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"))

  if (!is.null(cluster_colors)) {
    facet_clusters <- ggplot(summary_dat,
                             aes(Time_point, mean_CLR, group = Cluster,
                                 colour = de_score, fill = de_score)) +
      geom_ribbon(aes(ymin = mean_CLR - ci95, ymax = mean_CLR + ci95),
                  alpha = 0.5, colour = NA) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.8) +
      geom_label(
        data = label_dat,
        aes(x = x_pos, y = y_pos, label = up_label),
        inherit.aes = FALSE,
        hjust = 0.5,
        vjust = 0,
        colour = "red",
        fill = "white",
        alpha = 0.5,
        size = 3,
        label.size = NA
      ) +
      geom_label(
        data = label_dat,
        aes(x = x_pos, y = y_pos, label = down_label),
        inherit.aes = FALSE,
        hjust = 0.5,
        vjust = 1,
        colour = "blue",
        fill = "white",
        alpha = 0.5,
        size = 3,
        label.size = NA
      ) +
      facet_wrap(~ Cluster) +
      scale_colour_gradient2(name = "DE score", low = "blue", mid = "grey80", high = "red",
                             midpoint = 0, limits = c(-1, 1)) +
      scale_fill_gradient2(name = "DE score", low = "blue", mid = "grey80", high = "red",
                           midpoint = 0, limits = c(-1, 1)) +
      scale_x_continuous(
        breaks = time_points,
        limits = c(min(time_points) - 5, max(time_points) + 5),
        expand = c(0, 0)
      ) +
      labs(x = "Time points (minutes)", y = "CLR counts") +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey92", colour = "grey70"),
            legend.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            axis.title = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"))
  } else {
    facet_clusters <- ggplot(summary_dat, aes(Time_point, mean_CLR, group = Cluster)) +
      geom_ribbon(aes(ymin = mean_CLR - ci95, ymax = mean_CLR + ci95),
                  fill = "steelblue", alpha = 0.18, colour = NA) +
      geom_line(colour = "steelblue4", linewidth = 0.8) +
      geom_point(colour = "steelblue4", size = 1.8) +
      facet_wrap(~ Cluster) +
      scale_x_continuous(
        breaks = time_points,
        limits = c(min(time_points) - 5, max(time_points) + 5),
        expand = c(0, 0)
      ) +
      labs(x = "Time points (minutes)", y = "CLR counts") +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "grey92", colour = "grey70"),
            legend.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold"))
  }

  list(
    all_clusters = all_clusters,
    facet_clusters = facet_clusters
  )
}

make_cluster_map <- function(lst, colname) {
  if (is.null(lst) || length(lst) == 0) {
    out <- tibble(Gene = character(), tmp = character())
  } else {
    d <- utils::stack(lst)
    out <- tibble(Gene = as.character(d$values), tmp = as.character(d$ind)) %>%
      distinct(Gene, .keep_all = TRUE)
  }
  names(out)[2] <- colname
  out
}

cluster_compare <- function(clusters1, clusters2) {
  mat <- matrix(0, nrow = length(clusters2) + 1, ncol = length(clusters1) + 1)
  rownames(mat) <- c("Total", names(clusters2))
  colnames(mat) <- c("Total", names(clusters1))
  mat[1, -1] <- sapply(clusters1, length)
  mat[-1, 1] <- sapply(clusters2, length)
  for (i in seq_along(clusters2)) {
    for (j in seq_along(clusters1)) {
      mat[i + 1, j + 1] <- length(intersect(clusters2[[i]], clusters1[[j]]))
    }
  }
  mat
}

heatmap_comparison <- function(comparison_matrix, rows, columns) {
  tr <- if (length(which(rownames(comparison_matrix) == "Total")) == 1) which(rownames(comparison_matrix) == "Total") else 1
  tc <- if (length(which(colnames(comparison_matrix) == "Total")) == 1) which(colnames(comparison_matrix) == "Total") else 1

  data <- comparison_matrix[-tr, -tc, drop = FALSE]

  row_tot <- comparison_matrix[-tr, tc, drop = TRUE]
  col_tot <- comparison_matrix[tr, -tc, drop = TRUE]

  norm_row <- sweep(data, 1, row_tot, "/")
  norm_col <- sweep(data, 2, col_tot, "/")
  
  p1 <- heatmaply(
    round(norm_row, 3),
    Rowv = FALSE,
    dendrogram = "none",
    xlab = paste(columns, "clusters"),
    ylab = paste(rows, "clusters"),
    label_names = c(paste0(rows, " cluster"), paste0(columns, " cluster"), "Gene ratio"),
    limits = c(0, 1),
    width = 1000,
    height = 800
  )

  p2 <- heatmaply(
    round(norm_col, 3),
    Rowv = FALSE,
    dendrogram = "none",
    xlab = paste(columns, "clusters"),
    ylab = paste(rows, "clusters"),
    label_names = c(paste0(rows, " cluster"), paste0(columns, " cluster"), "Gene ratio"),
    limits = c(0, 1),
    width = 1000,
    height = 800
  )

  invisible(list(per_column = norm_col, per_row = norm_row, heatmap_rows = p1, heatmap_columns = p2))

  return(list(heatmap_rows = p1, heatmap_columns = p2))
}

# Enrichment analysis

## COG

cog_enrichment <- function(genes_set, all_genes) {

  cog_annotation <- read.csv("cog_annotation.tsv", sep = "\t", header = TRUE)
  cog_names <- read.csv("cog_name.tsv", sep = "\t", header = TRUE)
  
  cog_annotation_sep <- cog_annotation %>%
    filter(!is.na(COG) & COG != "") %>%
    mutate(COG = gsub("[^A-Z]", "", COG)) %>%
    mutate(COG = strsplit(COG, "")) %>%
    unnest(COG)
  
  df_cog <- data.frame(ID = all_genes) %>%
    mutate(DE = ID %in% genes_set) %>%
    left_join(cog_annotation_sep, by = "ID") %>%
    mutate(COG = ifelse(is.na(COG), "Unassigned", COG))

  genes_by_cog <- df_cog %>%
    filter(COG != "Unassigned", DE) %>%
    group_by(COG) %>%
    summarise(Genes = paste(ID, collapse = "/"), .groups = "drop")

  tab_df <- df_cog %>%
    filter(COG != "Unassigned") %>%
    group_by(COG, DE) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = DE, values_from = n, values_fill = list(n = 0))    
    
  if (!"TRUE" %in% colnames(tab_df)) tab_df$`TRUE` <- 0
  if (!"FALSE" %in% colnames(tab_df)) tab_df$`FALSE` <- 0
    
  tab_df <- tab_df %>%
    select(COG, `TRUE`, `FALSE`) %>%
    rename(DE = `TRUE`, noDE = `FALSE`)

  tab_df <- tab_df %>%
      left_join(genes_by_cog, by = "COG") %>%
      relocate(Genes, .after = noDE)

  total_DE <- sum(tab_df$DE)
  total_noDE <- sum(tab_df$noDE)

  cog_enrichment_df <- tab_df %>%
    rowwise() %>%
    mutate(
      ft = list(
        fisher.test(
          matrix(c(DE, noDE, total_DE - DE, total_noDE - noDE), nrow = 2),
          alternative = "greater"
        )
      ),
      odds.ratio = unname(ft$estimate),
      p.value = ft$p.value
    ) %>%
    ungroup() %>%
    mutate(padj = p.adjust(p.value, method = "BH")) %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    dplyr::select(-ft) %>%
    dplyr::left_join(
      cog_names %>% dplyr::select(COG, Name),
      by = "COG"
    ) %>%
    dplyr::mutate(COG = paste0(Name, " (", COG, ")")) %>%
    dplyr::select(-Name)

  return(cog_enrichment_df)
}

cog_enrichment_by_cluster <- function(cluster_gene_table, all_genes) {

  cog_annotation <- read.csv("cog_annotation.tsv", sep = "\t", header = TRUE)
  cog_names <- read.csv("cog_name.tsv", sep = "\t", header = TRUE)
  
  cog_annotation_sep <- cog_annotation %>%
    filter(!is.na(COG) & COG != "") %>%
    mutate(COG = gsub("[^A-Z]", "", COG)) %>%
    mutate(COG = strsplit(COG, "")) %>%
    unnest(COG)
  
  results_list <- list()
  
  for (cl in levels(cluster_gene_table$Cluster)) {
    
    hc_genes <- cluster_gene_table %>%
      filter(Cluster == cl) %>%
      pull(Gene)
    
    df_cog <- data.frame(ID = all_genes) %>%
      mutate(IN = ID %in% hc_genes) %>%
      left_join(cog_annotation_sep, by = "ID") %>%
      mutate(COG = ifelse(is.na(COG), "Unassigned", COG))

    genes_by_cog <- df_cog %>%
      filter(COG != "Unassigned", IN) %>%
      group_by(COG) %>%
      summarise(Genes = paste(ID, collapse = "/"), .groups = "drop")
    
    tab_df <- df_cog %>%
      filter(COG != "Unassigned") %>%
      group_by(COG, IN) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = IN, values_from = n, values_fill = list(n = 0))      
    
    if (!"TRUE" %in% colnames(tab_df)) tab_df$`TRUE` <- 0
    if (!"FALSE" %in% colnames(tab_df)) tab_df$`FALSE` <- 0
    
    tab_df <- tab_df %>%
      dplyr::select(COG, `TRUE`, `FALSE`) %>%
      dplyr::rename(IN = `TRUE`, OUT = `FALSE`)

    tab_df <- tab_df %>%
      left_join(genes_by_cog, by = "COG") %>%
      relocate(Genes, .after = OUT)
    
    total_IN <- sum(tab_df$IN)
    total_OUT <- sum(tab_df$OUT)
    
    cog_enrichment_df <- tab_df %>%
      rowwise() %>%
      mutate(
        ft = list(
          fisher.test(
            matrix(c(IN, OUT, total_IN - IN, total_OUT - OUT), nrow = 2),
            alternative = "greater"
          )
        ),
        odds.ratio = unname(ft$estimate),
        p.value    = ft$p.value
      ) %>%
      ungroup() %>%
      mutate(padj = p.adjust(p.value, method = "fdr")) %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      dplyr::select(-ft) %>%
      mutate(Cluster = cl) %>%
      dplyr::left_join(
        cog_names %>% dplyr::select(COG, Name),
        by = "COG"
      ) %>%
      dplyr::mutate(COG = paste0(Name, " (", COG, ")")) %>%
      dplyr::select(-Name)
    
    results_list[[cl]] <- cog_enrichment_df
  }
  
  final_results <- bind_rows(results_list) %>% filter(padj < 0.05)
  
  return(final_results)
}

## KEGG

kegg_enrichment <- function(genes_set, all_genes) {

  kobas_map <- read.csv("kegg_annotation.tsv", header = TRUE, sep = "\t", col.names = c("original_id", "kegg_id"))

  mapped_genes <- dplyr::inner_join(
      data.frame(original_id = genes_set),
      kobas_map,
      by = "original_id"
    ) %>%
      dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
      dplyr::distinct(kegg_id) %>%
      dplyr::pull(kegg_id)

  universe_genes <- dplyr::inner_join(
      data.frame(original_id = all_genes),
      kobas_map,
      by = "original_id"
    ) %>%
      dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
      dplyr::distinct(kegg_id) %>%
      dplyr::pull(kegg_id)

  kk <- enrichKEGG(gene = mapped_genes,
                   organism = 'btn',
                   keyType = "kegg",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   universe = universe_genes)

  kk_df <- as.data.frame(kk)

  if (nrow(kk_df) == 0) {

      cols <- c(
        "category", "subcategory", "ID", "Description", "GeneRatio", "BgRatio",
        "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust",
        "qvalue", "geneID", "Count"
      )

      kk_df <- setNames(
        data.frame(matrix(ncol = length(cols), nrow = 0)),
        cols
      )

      return(kk_df)
  }

  if ("geneID" %in% colnames(kk_df)) {
      map_df <- kobas_map %>%
        dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
        dplyr::distinct(kegg_id, original_id) %>%
        dplyr::group_by(kegg_id) %>%
        dplyr::summarise(original_ids = paste(unique(original_id), collapse = "/"), .groups = "drop")

      key2orig <- stats::setNames(map_df$original_ids, map_df$kegg_id)

      kk_df$geneID <- vapply(strsplit(kk_df$geneID, "/"), function(x) {
        x <- trimws(x)
        m <- ifelse(x %in% names(key2orig), key2orig[x], x)
        paste(unique(unname(m)), collapse = "/")
      }, FUN.VALUE = character(1))
    }

  return(kk_df)
}

kegg_enrichment_by_cluster <- function(cluster_gene_table, all_genes) {

  kobas_map <- read.csv("kegg_annotation.tsv", header = TRUE, sep = "\t", col.names = c("original_id", "kegg_id"))

  mapped_genes <- dplyr::inner_join(
      data.frame(original_id = cluster_gene_table$Gene),
      kobas_map,
      by = "original_id"
    ) %>%
      dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
      dplyr::distinct(kegg_id) %>%
      dplyr::pull(kegg_id)

  universe_genes <- dplyr::inner_join(
      data.frame(original_id = all_genes),
      kobas_map,
      by = "original_id"
    ) %>%
      dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
      dplyr::distinct(kegg_id) %>%
      dplyr::pull(kegg_id)

  results_list <- list()

  for (cl in levels(cluster_gene_table$Cluster)) {

    hc_genes <- cluster_gene_table %>%
      dplyr::filter(Cluster == cl) %>%
      dplyr::pull(Gene)

    gene_kegg <- dplyr::inner_join(
        data.frame(original_id = hc_genes),
        kobas_map,
        by = "original_id"
      ) %>%
        dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
        dplyr::distinct(kegg_id) %>%
        dplyr::pull(kegg_id)

    kk <- enrichKEGG(
      gene         = gene_kegg,
      organism     = 'btn',
      keyType      = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      universe     = universe_genes
    )

    kk_df <- as.data.frame(kk) %>%
      dplyr::mutate(Cluster = cl)

    results_list[[as.character(cl)]] <- kk_df
  }

  final_results <- dplyr::bind_rows(results_list)

  map_df <- kobas_map %>%
    dplyr::filter(!is.na(kegg_id), kegg_id != "None") %>%
    dplyr::distinct(kegg_id, original_id) %>%
    dplyr::group_by(kegg_id) %>%
    dplyr::summarise(original_ids = paste(unique(original_id), collapse = "/"), .groups = "drop")

  key2orig <- stats::setNames(map_df$original_ids, map_df$kegg_id)

  final_results$geneID <- vapply(strsplit(final_results$geneID, "/"), function(x) {
    x <- trimws(x)
    m <- ifelse(x %in% names(key2orig), key2orig[x], x)
    paste(unique(unname(m)), collapse = "/")
  }, FUN.VALUE = character(1))

  return(final_results)
}

## GO

go_enrichment <- function(genes_set, all_genes, ont, plots = FALSE, showCategory = 10) {

  ego <- enrichGO(gene = genes_set,
                  OrgDb = org.Bthuringiensis.eg.db,
                  keyType = "GID",
                  ont = ont,
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  universe = all_genes)

  ego_df <- as.data.frame(ego)

  cols <- c(
    "ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
    "FoldEnrichment", "zScore", "pvalue", "p.adjust", "qvalue",
    "geneID", "Count"
  )

  if (nrow(ego_df) == 0) {
    ego_df <- setNames(data.frame(matrix(ncol = length(cols), nrow = 0)), cols)
  } else if (plots) {
    print(dotplot(ego, showCategory = showCategory))
    print(barplot(ego, showCategory = showCategory))
    print(cnetplot(ego, showCategory = showCategory))
  }

  return(ego_df)
}

go_enrichment_by_cluster <- function(cluster_gene_table, all_genes, ont) {

  results_list <- list()

  for (cl in levels(cluster_gene_table$Cluster)) {

    hc_genes <- cluster_gene_table %>%
      dplyr::filter(Cluster == cl) %>%
      dplyr::pull(Gene)

    ego <- enrichGO(
      gene          = hc_genes,
      OrgDb         = org.Bthuringiensis.eg.db,
      keyType       = "GID",
      ont           = ont,
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",      
      universe      = all_genes
    )

    ego_df <- as.data.frame(ego) %>%
      dplyr::mutate(Cluster = cl)

    results_list[[as.character(cl)]] <- ego_df
  }

  final_results <- dplyr::bind_rows(results_list) %>%
    dplyr::relocate(Cluster, .after = dplyr::last_col())

  return(final_results)
}

## EC

ec_enrichment <- function(genes_set, all_genes, number) {

  ec_annotation <- read.csv("ec_annotation.tsv", sep = "\t", header = TRUE)
  ec_names <- read.csv("ec_name.tsv", sep = "\t", header = TRUE)

  ec_annotation_sep <- ec_annotation %>%
    {
      if (number == 1) {
        mutate(., EC = sub("^([0-9]+)\\..*", "\\1", EC))
      } else if (number == 2) {        
        filter(., grepl("^[0-9]+\\.[0-9]+\\.", EC)) %>%
        mutate(EC = sub("^([0-9]+\\.[0-9]+)\\..*", "\\1", EC))
      } else if (number == 3) {        
        filter(., grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.", EC)) %>%
        mutate(EC = sub("^([0-9]+\\.[0-9]+\\.[0-9]+)\\..*", "\\1", EC))
      } else {
        stop("Number must be 1, 2 or 3.")
      }
    } %>%
    distinct()

  if (number == 1) {
    ec_annotation_sep <- ec_annotation_sep %>%
      left_join(ec_names %>% rename(Class = Name), by = "EC") %>%
      distinct(EC, Class)
  }

  if (number == 2) {
    ec_annotation_sep <- ec_annotation_sep %>%
      mutate(Class_EC = sub("\\..*", "", EC)) %>%
      left_join(ec_names %>% rename(Class = Name), by = c("Class_EC" = "EC")) %>%
      left_join(ec_names %>% rename(Subclass = Name), by = "EC") %>%
      select(-Class_EC) %>%
      distinct(EC, Class, Subclass)
  }

  if (number == 3) {
    ec_annotation_sep <- ec_annotation_sep %>%
      mutate(Class_EC = sub("\\..*", "", EC),
             Subclass_EC = sub("\\.[0-9]+$", "", EC)) %>%
      left_join(ec_names %>% rename(Class = Name), by = c("Class_EC" = "EC")) %>%
      left_join(ec_names %>% rename(Subclass = Name), by = c("Subclass_EC" = "EC")) %>%
      left_join(ec_names %>% rename(`Sub-subclass` = Name), by = "EC") %>%
      select(-Class_EC, -Subclass_EC) %>%
      distinct(EC, Class, Subclass, `Sub-subclass`)
  }

  df_ec <- data.frame(ID = all_genes) %>%
    mutate(DE = ID %in% genes_set) %>%
    left_join(ec_annotation, by = "ID") %>%
    mutate(EC =
      if (number == 1) sub("^([0-9]+)\\..*", "\\1", EC) else
      if (number == 2) sub("^([0-9]+\\.[0-9]+)\\..*", "\\1", EC) else
      if (number == 3) sub("^([0-9]+\\.[0-9]+\\.[0-9]+)\\..*", "\\1", EC) else EC) %>%
    mutate(EC = ifelse(is.na(EC), "Unassigned", EC))

  if (number == 1) {
    df_ec <- df_ec %>% filter(grepl("^[0-9]+$", EC) | EC == "Unassigned")
  }

  if (number == 2) {
    df_ec <- df_ec %>% filter(grepl("^[0-9]+\\.[0-9]+$", EC) | EC == "Unassigned")
  }

  if (number == 3) {
    df_ec <- df_ec %>% filter(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", EC) | EC == "Unassigned")
  }

  df_ec <- df_ec %>% distinct(ID, EC, .keep_all = TRUE)

  df_ec <- df_ec %>% left_join(ec_annotation_sep, by = "EC")

  genes_by_ec <- df_ec %>%
    filter(EC != "Unassigned", DE) %>%
    group_by(EC) %>%
    summarise(Genes = paste(ID, collapse = "/"), .groups = "drop")

  tab_df <- df_ec %>%
    filter(EC != "Unassigned") %>%
    group_by(EC, DE) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = DE, values_from = n, values_fill = 0)

  if (!"TRUE" %in% colnames(tab_df)) tab_df$`TRUE` <- 0
  if (!"FALSE" %in% colnames(tab_df)) tab_df$`FALSE` <- 0

  tab_df <- tab_df %>%
    select(EC, `TRUE`, `FALSE`) %>%
    rename(DE = `TRUE`, noDE = `FALSE`)

  if (number == 1) {
    tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
      relocate(Class, .after = EC)
  }

  if (number == 2) {
    tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
      relocate(Class, Subclass, .after = EC)
  }

  if (number == 3) {
    tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
      relocate(Class, Subclass, `Sub-subclass`, .after = EC)
  }

  tab_df <- tab_df %>%
    left_join(genes_by_ec, by = "EC") %>%
    relocate(Genes, .after = noDE)

  total_DE <- sum(tab_df$DE)
  total_noDE <- sum(tab_df$noDE)

  ec_enrichment_df <- tab_df %>%
    rowwise() %>%
    mutate(
      ft = list(
        fisher.test(
          matrix(c(DE, noDE, total_DE - DE, total_noDE - noDE), nrow = 2),
          alternative = "greater"
        )
      ),
      odds.ratio = unname(ft$estimate),
      p.value    = ft$p.value
    ) %>%
    ungroup() %>%
    mutate(padj = p.adjust(p.value, method = "fdr")) %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    select(-ft)

  return(ec_enrichment_df)
}

ec_enrichment_by_cluster <- function(cluster_gene_table, all_genes, number) {

  ec_annotation <- read.csv("ec_annotation.tsv", sep = "\t", header = TRUE)
  ec_names <- read.csv("ec_name.tsv", sep = "\t", header = TRUE)

  ec_annotation_sep <- ec_annotation %>%
    {
      if (number == 1) {
        mutate(., EC = sub("^([0-9]+)\\..*", "\\1", EC))
      } else if (number == 2) {        
        filter(., grepl("^[0-9]+\\.[0-9]+\\.", EC)) %>%
        mutate(EC = sub("^([0-9]+\\.[0-9]+)\\..*", "\\1", EC))
      } else if (number == 3) {        
        filter(., grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.", EC)) %>%
        mutate(EC = sub("^([0-9]+\\.[0-9]+\\.[0-9]+)\\..*", "\\1", EC))
      } else {
        stop("Number must be 1, 2 or 3.")
      }
    } %>%
    distinct()

  if (number == 1) {
    ec_annotation_sep <- ec_annotation_sep %>%
      left_join(ec_names %>% dplyr::rename(Class = Name), by = "EC") %>%
      distinct(EC, Class)
  }

  if (number == 2) {
    ec_annotation_sep <- ec_annotation_sep %>%
      mutate(Class_EC = sub("\\..*", "", EC)) %>%
      left_join(ec_names %>% dplyr::rename(Class = Name), by = c("Class_EC" = "EC")) %>%
      left_join(ec_names %>% dplyr::rename(Subclass = Name), by = "EC") %>%
      dplyr::select(-Class_EC) %>%
      distinct(EC, Class, Subclass)
  }

  if (number == 3) {
    ec_annotation_sep <- ec_annotation_sep %>%
      mutate(Class_EC = sub("\\..*", "", EC),
             Subclass_EC = sub("\\.[0-9]+$", "", EC)) %>%
      left_join(ec_names %>% dplyr::rename(Class = Name), by = c("Class_EC" = "EC")) %>%
      left_join(ec_names %>% dplyr::rename(Subclass = Name), by = c("Subclass_EC" = "EC")) %>%
      left_join(ec_names %>% dplyr::rename(`Sub-subclass` = Name), by = "EC") %>%
      dplyr::select(-Class_EC, -Subclass_EC) %>%
      distinct(EC, Class, Subclass, `Sub-subclass`)
  }

  results_list <- list()

  for (cl in levels(cluster_gene_table$Cluster)) {

    hc_genes <- cluster_gene_table %>%
      filter(Cluster == cl) %>%
      pull(Gene)

    df_ec <- data.frame(ID = all_genes) %>%
      mutate(IN = ID %in% hc_genes) %>%
      left_join(ec_annotation, by = "ID") %>%
      mutate(EC =
        if (number == 1) sub("^([0-9]+)\\..*", "\\1", EC) else
        if (number == 2) sub("^([0-9]+\\.[0-9]+)\\..*", "\\1", EC) else
        if (number == 3) sub("^([0-9]+\\.[0-9]+\\.[0-9]+)\\..*", "\\1", EC) else EC) %>%
      mutate(EC = ifelse(is.na(EC), "Unassigned", EC)) %>%
      distinct(ID, EC, .keep_all = TRUE)

    if (number == 1) {
      df_ec <- df_ec %>% filter(grepl("^[0-9]+$", EC) | EC == "Unassigned")
    }

    if (number == 2) {
      df_ec <- df_ec %>% filter(grepl("^[0-9]+\\.[0-9]+$", EC) | EC == "Unassigned")
    }

    if (number == 3) {
      df_ec <- df_ec %>% filter(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", EC) | EC == "Unassigned")
    }

    df_ec <- df_ec %>% left_join(ec_annotation_sep, by = "EC")

    genes_by_ec <- df_ec %>%
      filter(EC != "Unassigned", IN) %>%
      group_by(EC) %>%
      summarise(Genes = paste(ID, collapse = "/"), .groups = "drop")

    tab_df <- df_ec %>%
      filter(EC != "Unassigned") %>%
      group_by(EC, IN) %>%
      summarise(n = n(), .groups = "drop") %>%
      pivot_wider(names_from = IN, values_from = n, values_fill = 0)

    if (!"TRUE" %in% colnames(tab_df)) tab_df$`TRUE` <- 0
    if (!"FALSE" %in% colnames(tab_df)) tab_df$`FALSE` <- 0

    tab_df <- tab_df %>%
      dplyr::select(EC, `TRUE`, `FALSE`) %>%
      dplyr::rename(IN = `TRUE`, OUT = `FALSE`)

    if (number == 1) {
      tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
        relocate(Class, .after = EC)
    }

    if (number == 2) {
      tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
        relocate(Class, Subclass, .after = EC)
    }

    if (number == 3) {
      tab_df <- left_join(tab_df, ec_annotation_sep, by = "EC") %>%
        relocate(Class, Subclass, `Sub-subclass`, .after = EC)
    }

    tab_df <- tab_df %>%
      left_join(genes_by_ec, by = "EC") %>%
      relocate(Genes, .after = OUT)

    total_IN  <- sum(tab_df$IN)
    total_OUT <- sum(tab_df$OUT)

    ec_enrichment_df <- tab_df %>%
      rowwise() %>%
      mutate(
        ft = list(fisher.test(matrix(c(IN, OUT, total_IN - IN, total_OUT - OUT), nrow = 2), alternative = "greater")),
        odds.ratio = unname(ft$estimate),
        p.value    = ft$p.value
      ) %>%
      ungroup() %>%
      mutate(padj = p.adjust(p.value, method = "fdr")) %>%
      filter(padj < 0.05) %>%
      arrange(padj) %>%
      dplyr::select(-ft) %>%
      mutate(Cluster = cl)

    results_list[[as.character(cl)]] <- ec_enrichment_df
  }

  final_results <- bind_rows(results_list) %>%
    filter(padj < 0.05)

  return(final_results)
}

# Enrichment summary

enrichment_summary <- function(genes_set, categories) {
  collapse_vals <- function(x) {
    vals <- unique(stats::na.omit(x))
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = " / ")
  }

  pick_col <- function(df) {
    cands <- c("COG", "subcategory", "Description", "EC")
    hit <- cands[cands %in% names(df)]
    if (length(hit)) hit[[1]] else names(df)[1]
  }

  vals <- lapply(categories, function(df) {
    if (!is.data.frame(df) || nrow(df) == 0) return(NA_character_)
    col <- pick_col(df)
    collapse_vals(df[[col]])
  })

  res <- as.data.frame(as.list(vals), stringsAsFactors = FALSE)
  target_names <- c("COG", "KEGG", "GO-BP", "GO-MF", "GO-CC", "EC-1", "EC-2", "EC-3")
  names(res) <- target_names[seq_len(ncol(res))]
  res <- data.frame(Size = length(genes_set), res, check.names = FALSE, stringsAsFactors = FALSE)
  res
}

collapse_vals <- function(x) {
  vals <- unique(na.omit(x))
  if (length(vals) == 0) NA_character_ else paste(vals, collapse = " / ")
}

# Save results

create_summary <- function(contrast_name = "contrast", rds_name = "result"){
  if(!grepl("\\.rds$", rds_name, ignore.case = TRUE)){
    rds_name <- paste0(rds_name, ".rds")
  }

  prefix <- sub("\\.rds$", "", rds_name, ignore.case = TRUE)

  parts <- strsplit(prefix, "_")[[1]]

  is_time <- grepl("Time points", contrast_name, ignore.case = TRUE)

  if(is_time){
    group <- toupper(parts[1])
  } else {
    grp <- parts[3]
    if(grepl("^lysogeny$", grp, ignore.case = TRUE)){
      group <- "Lysogeny"
    } else {
      group <- grp
    }
  }

  get_n <- function(name, type){
    if(exists(name, inherits = TRUE)){
      obj <- get(name, inherits = TRUE)
      if(type == "genes") return(length(obj))
      return(nrow(obj))
    } else {
      return(0)
    }
  }

  contrast_summary <- data.frame(
    Size = c(get_n("up_genes", "genes"), get_n("down_genes", "genes")),
    COG = c(get_n("cog_up", "table"), get_n("cog_down", "table")),
    KEGG = c(get_n("kegg_up", "table"), get_n("kegg_down", "table")),
    `GO-BP` = c(get_n("go_bp_up", "table"), get_n("go_bp_down", "table")),
    `GO-MF` = c(get_n("go_mf_up", "table"), get_n("go_mf_down", "table")),
    `GO-CC` = c(get_n("go_cc_up", "table"), get_n("go_cc_down", "table")),
    `EC-1` = c(get_n("ec_1_up", "table"), get_n("ec_1_down", "table")),
    `EC-2` = c(get_n("ec_2_up", "table"), get_n("ec_2_down", "table")),
    `EC-3` = c(get_n("ec_3_up", "table"), get_n("ec_3_down", "table")),
    row.names = c("UP", "DOWN"),
    check.names = FALSE
  )

  attr(contrast_summary, "contrast") <- contrast_name
  attr(contrast_summary, "group") <- group
  saveRDS(contrast_summary, file = paste0(prefix, "_contrast_summary.rds"))

  if(exists("contrast_tectivirus", inherits = TRUE)){

    contrast_tectivirus_summary <- data.frame(
      Size = c(
        sum(contrast_tectivirus$Adjusted_p_value < 0.05 &
              contrast_tectivirus$LFC > log2(1.5), na.rm = TRUE),
        sum(contrast_tectivirus$Adjusted_p_value < 0.05 &
              contrast_tectivirus$LFC < -log2(1.5), na.rm = TRUE)
      ),
      row.names = c("UP", "DOWN"),
      check.names = FALSE
    )

    attr(contrast_tectivirus_summary, "contrast") <- contrast_name
    attr(contrast_tectivirus_summary, "group") <- group

    saveRDS(
      contrast_tectivirus_summary,
      file = paste0(prefix, "_contrast_tectivirus_summary.rds")
    )
  }

  if(exists("es_de", inherits = TRUE)){
    enrichment_summary <- get("es_de", inherits = TRUE)
  } else {
    enrichment_summary <- NULL
  }

  attr(enrichment_summary, "contrast") <- contrast_name
  attr(enrichment_summary, "group") <- group
  saveRDS(enrichment_summary, file = paste0(prefix,"_enrichment_summary.rds"))

  if(exists("up_genes", inherits = TRUE)){
    up_genes <- get("up_genes", inherits = TRUE)
  } else {
    up_genes <- NULL
  }

  if(exists("down_genes", inherits = TRUE)){
    down_genes <- get("down_genes", inherits = TRUE)
  } else {
    down_genes <- NULL
  }

  if(exists("contrast_df", inherits = TRUE)){
    contrast_df_obj <- get("contrast_df", inherits = TRUE)

    up_df <- contrast_df_obj[contrast_df_obj$Gene_ID %in% up_genes, c("Gene_ID", "LFC")]
    down_df <- contrast_df_obj[contrast_df_obj$Gene_ID %in% down_genes, c("Gene_ID", "LFC")]
  } else {
    up_df <- NULL
    down_df <- NULL
  }

  de_genes <- list(
    up_genes = up_df,
    down_genes = down_df
  )

  attr(de_genes,"contrast") <- contrast_name
  attr(de_genes,"group") <- group
  saveRDS(de_genes, file = paste0(prefix, "_de_genes.rds"))

  if(exists("contrast_tectivirus_df", inherits = TRUE) &&
    exists("contrast_tectivirus", inherits = TRUE)){

    contrast_tectivirus_df_obj <- get("contrast_tectivirus_df", inherits = TRUE)
    contrast_tectivirus_obj <- get("contrast_tectivirus", inherits = TRUE)

    up_tectivirus <- rownames(
      contrast_tectivirus_obj[
        contrast_tectivirus_obj$padj < 0.05 &
          contrast_tectivirus_obj$log2FoldChange > log2(1.5),
      ]
    )

    down_tectivirus <- rownames(
      contrast_tectivirus_obj[
        contrast_tectivirus_obj$padj < 0.05 &
          contrast_tectivirus_obj$log2FoldChange < -log2(1.5),
      ]
    )

    up_tectivirus_df <- contrast_tectivirus_df_obj[
      contrast_tectivirus_df_obj$Gene_symbol %in% up_tectivirus,
      c("Gene_symbol", "LFC")
    ]

    down_tectivirus_df <- contrast_tectivirus_df_obj[
      contrast_tectivirus_df_obj$Gene_symbol %in% down_tectivirus,
      c("Gene_symbol", "LFC")
    ]

    colnames(up_tectivirus_df)[1] <- "Gene_ID"
    colnames(down_tectivirus_df)[1] <- "Gene_ID"

    de_tectivirus_genes <- list(
      up_genes = up_tectivirus_df,
      down_genes = down_tectivirus_df
    )

    attr(de_tectivirus_genes, "contrast") <- contrast_name
    attr(de_tectivirus_genes, "group") <- group

    saveRDS(
      de_tectivirus_genes,
      file = paste0(prefix, "_de_tectivirus_genes.rds")
    )
  }

  make_cog_df <- function(x, direction){

    if(is.null(x) || nrow(x) == 0){
      return(NULL)
    }

    out <- data.frame(
      Direction = direction,
      COG = x$COG,
      DE = x$DE,
      No_DE = x$noDE,
      Total = x$DE + x$noDE,
      Ratio = x$DE / (x$DE + x$noDE),
      p = x$padj,
      check.names = FALSE
    )

    if(is_time){

      out <- cbind(
        Strain = group,
        Contrast = contrast_name,
        out
      )

    } else if(identical(group, "Lysogeny")){

      out <- cbind(
        Contrast = contrast_name,
        out
      )

    } else {

      out <- cbind(
        Contrast = contrast_name,
        Time_point = sub("^t", "", group),
        out
      )

    }

    out
  }

  cog_up_df <- if(exists("cog_up", inherits = TRUE)) {
    make_cog_df(get("cog_up", inherits = TRUE), "UP")
  } else {
    NULL
  }

  cog_down_df <- if(exists("cog_down", inherits = TRUE)) {
    make_cog_df(get("cog_down", inherits = TRUE), "DOWN")
  } else {
    NULL
  }

  cog_df <- do.call(
    rbind,
    Filter(Negate(is.null), list(cog_up_df, cog_down_df))
  )

  if(is.null(cog_df)){
    cog_df <- data.frame()
  }

  attr(cog_df, "contrast") <- contrast_name
  attr(cog_df, "group") <- group

  saveRDS(
    cog_df,
    file = paste0(prefix, "_cog.rds")
  )

  make_kegg_df <- function(x, direction){

    if(is.null(x) || nrow(x) == 0){
      return(NULL)
    }

    de <- as.numeric(sub("/.*", "", x$GeneRatio))
    total <- as.numeric(sub("/.*", "", x$BgRatio))

    out <- data.frame(
      Direction = direction,
      KEGG = x$subcategory,
      DE = de,
      No_DE = total - de,
      Total = total,
      Ratio = de / total,
      p = x$p.adjust,
      check.names = FALSE
    )

    if(is_time){

      out <- cbind(
        Strain = group,
        Contrast = contrast_name,
        out
      )

    } else if(identical(group, "Lysogeny")){

      out <- cbind(
        Contrast = contrast_name,
        out
      )

    } else {

      out <- cbind(
        Contrast = contrast_name,
        Time_point = sub("^t", "", group),
        out
      )

    }

    out
  }

  kegg_up_df <- if(exists("kegg_up", inherits = TRUE)) {
    make_kegg_df(get("kegg_up", inherits = TRUE), "UP")
  } else {
    NULL
  }

  kegg_down_df <- if(exists("kegg_down", inherits = TRUE)) {
    make_kegg_df(get("kegg_down", inherits = TRUE), "DOWN")
  } else {
    NULL
  }

  kegg_df <- do.call(
    rbind,
    Filter(Negate(is.null), list(kegg_up_df, kegg_down_df))
  )

  if(is.null(kegg_df)){
    kegg_df <- data.frame()
  }

  attr(kegg_df, "contrast") <- contrast_name
  attr(kegg_df, "group") <- group

  saveRDS(
    kegg_df,
    file = paste0(prefix, "_kegg.rds")
  )

  invisible(NULL)
}

# Networks analysis

network_hist <- function(x,
                         xlab = deparse(substitute(x)),
                         title = NULL,
                         n_intervals = 10,
                         y_pad = 0.05,
                         show_counts = TRUE,
                         fill = "lightgray",
                         color = "black",
                         rug = TRUE) {

  x <- x[!is.na(x)]
  
  brks <- pretty(range(x), n = n_intervals)
  if (length(brks) < 2) brks <- seq(min(x), max(x), length.out = n_intervals + 1)

  h <- hist(x, breaks = brks, plot = FALSE)
  max_count <- if (length(h$counts)) max(h$counts) else 0
  y_top <- if (max_count > 0) max_count * (1 + y_pad) else 1

  p <- ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(breaks = brks, fill = fill, color = color) +
    coord_cartesian(ylim = c(0, y_top)) +
    xlab(xlab) +
    ylab("Frequency") +
    ggtitle(if (is.null(title)) NULL else title) +
    theme_minimal(base_size = 14) +    
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"))

  if (show_counts && length(h$mids) > 0) {
    p <- p + 
      geom_text(
        data = data.frame(xmid = h$mids, count = h$counts),
        aes(x = xmid, y = count, label = count),
        vjust = -0.3, size = 3, inherit.aes = FALSE
      )
  }

  if (rug) p <- p + geom_rug(sides = "b", alpha = 0.25)

  p
}

apl <- function(g) {
  n <- vcount(g)

  apl_values <- sapply(seq_len(n), function(i) {
    q <- as.numeric(distances(g, v = V(g)[i]))
    finite <- is.finite(q)
    finite[i] <- FALSE
    num <- sum(q[finite], na.rm = TRUE)
    den <- sum(finite)
    if (den > 0) num / den else NA_real_
  })

  names(apl_values) <- if (!is.null(V(g)$name)) V(g)$name else as.character(seq_len(n))

  apl_values
}

build_nodes_data <- function(data) {
  
  result <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)

  for (col_name in names(data)) {
    vec <- data[[col_name]]
    col <- rep(NA_real_, length(all_genes))
    idx <- match(names(vec), all_genes)
    valid <- !is.na(idx)
    col[idx[valid]] <- unname(vec)[valid]
    result[[col_name]] <- col
  }

  rownames(result) <- NULL
  result
}

get_keystone_genes <- function(data, attributes, top) {
  k <- ceiling(top * nrow(data))

  ranks <- lapply(attributes, function(attr) rank(-data[[attr]], ties.method = "average"))
  ranks <- as.data.frame(ranks)
  names(ranks) <- attributes

  keep <- apply(ranks <= k, 1, all)

  out <- cbind(Gene = data$Gene[keep], ranks[keep, , drop = FALSE])
  out <- out[order(rowSums(out[attributes]), out$Gene), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Networks comparisons

load_networks_objects_rds <- function(
  dir = ".",
  pattern = "_network\\.rds$",
  envir = .GlobalEnv,
  overwrite = FALSE
) {
  files <- list.files(path = dir, pattern = pattern, full.names = TRUE)
  for (f in files) {
    suffix <- sub("_network\\.rds$", "", basename(f))
    x <- readRDS(f)
    for (nm in names(x)) {
      new_name <- paste0(nm, "_", suffix)
      if (!overwrite && exists(new_name, envir = envir, inherits = FALSE)) next
      assign(new_name, x[[nm]], envir = envir)
    }
  }
  invisible(NULL)
}

network_stats_tests <- function(data, names) {
  stopifnot(is.list(data), length(data) >= 2)
  stopifnot(length(names) == length(data))

  names(data) <- names

  common_genes <- Reduce(intersect, lapply(data, function(x) x$Gene))
  data <- lapply(data, function(x) x[x$Gene %in% common_genes, ])

  metrics <- setdiff(colnames(data[[1]]), "Gene")

  results <- lapply(metrics, function(m) {
    df <- do.call(
      cbind,
      lapply(data, function(x) x[[m]])
    )
    colnames(df) <- names(data)
    df <- df[complete.cases(df), , drop = FALSE]

    Friedman_p_value <- if (nrow(df) > 1) friedman.test(df)$p.value else NA_real_
    
    Kendalls_W <- if (nrow(df) > 1) {
      ft <- friedman.test(df)
      k <- ncol(df)
      n <- nrow(df)
      W <- ft$statistic / (n * (k - 1))
      as.numeric(W)
    } else {
      NA_real_
    }

    Pairwise_tests <- NULL
    if (!is.na(Friedman_p_value) && Friedman_p_value < 0.05) {
      combs <- combn(colnames(df), 2, simplify = FALSE)
      pw <- lapply(combs, function(cc) {
        p <- wilcox.test(df[, cc[1]], df[, cc[2]], paired = TRUE)$p.value
        data.frame(net1 = cc[1], net2 = cc[2], p_value = p)
      })
      Pairwise_tests <- do.call(rbind, pw)
      Pairwise_tests$p_adj <- p.adjust(Pairwise_tests$p_value, method = "holm")
    }

    list(
      Friedman_p_value = Friedman_p_value,
      Kendalls_W = Kendalls_W,
      Pairwise_tests = Pairwise_tests
    )
  })

  names(results) <- metrics
  results
}

networks_results_to_table <- function(results, alpha = 0.05) {
  metrics <- names(results)

  all_pairs <- unique(do.call(
    c,
    lapply(results, function(x) {
      if (is.null(x$Pairwise_tests)) return(NULL)
      paste(x$Pairwise_tests$net1, x$Pairwise_tests$net2, sep = "_vs_")
    })
  ))

  out <- lapply(metrics, function(m) {
    res_m <- results[[m]]
    row <- setNames(as.list(rep(NA_real_, length(all_pairs))), all_pairs)

    if (!is.null(res_m$Pairwise_tests) && res_m$Friedman_p_value < alpha) {
      pw <- res_m$Pairwise_tests
      pair_names <- paste(pw$net1, pw$net2, sep = "_vs_")
      row[pair_names] <- pw$p_adj
    }

    c(list(
      Attribute = m,
      Friedman_p_value = res_m$Friedman_p_value,
      Kendalls_W = res_m$Kendalls_W
    ), row)
  })

  df <- do.call(rbind, lapply(out, function(x) as.data.frame(x, optional = TRUE)))
  rownames(df) <- as.character(df$Attribute)
  df$Attribute <- NULL
  num_cols <- setdiff(colnames(df), "Attribute")
  df[num_cols] <- lapply(df[num_cols], function(x) suppressWarnings(as.numeric(x)))
  df
}

networks_attributes_differences <- function(data, names, gene_order, ranking_attributes = FALSE) {

  stopifnot(is.list(data), length(data) >= 2)
  stopifnot(length(names) == length(data))
  stopifnot(!missing(gene_order))
  names(data) <- names
  stopifnot(all(sapply(data, function(x) "Gene" %in% colnames(x))))

  nets <- names
  nets_lower <- tolower(nets)

  attributes_sets <- lapply(data, function(x) setdiff(colnames(x), "Gene"))
  all_attributes <- Reduce(intersect, attributes_sets)
  attr_labels <- gsub("_", " ", all_attributes)
  attributes <- all_attributes

  if (isFALSE(ranking_attributes)) {
    ranking_attributes <- all_attributes
    ranking_attr_labels <- attr_labels
  } else {
    ranking_attr_labels <- gsub("_", " ", ranking_attributes)
  }

  z_scores <- lapply(nets, function(net) {
    df <- data[[net]]
    out <- df
    for (m in attributes) {
      v <- df[[m]]
      out[[m]] <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    }
    out
  })
  names(z_scores) <- nets

  combs <- combn(nets, 2, simplify = FALSE)
  pair_levels <- vapply(combs, function(cc) paste(cc[2], cc[1], sep = " - "), character(1))

  diffs_long_list <- list()
  diffs_wide <- list()
  summary_list <- list()

  for (cc in combs) {
    net1 <- cc[1]; net2 <- cc[2]
    d1 <- data[[net1]]; d2 <- data[[net2]]

    gcommon <- intersect(d1$Gene, d2$Gene)
    if (length(gcommon) == 0) next

    m1 <- d1[match(gcommon, d1$Gene), c("Gene", attributes), drop = FALSE]
    m2 <- d2[match(gcommon, d2$Gene), c("Gene", attributes), drop = FALSE]

    z1 <- z_scores[[net1]][match(gcommon, z_scores[[net1]]$Gene), ]
    z2 <- z_scores[[net2]][match(gcommon, z_scores[[net2]]$Gene), ]

    df_wide <- data.frame(Gene_ID = gcommon)
    pair_key <- paste(net2, net1, sep = " - ")

    for (i in seq_along(attributes)) {
      m <- attributes[i]
      lbl <- attr_labels[i]

      val1 <- m1[[m]]
      val2 <- m2[[m]]
      diffv <- val2 - val1

      z1v <- z1[[m]]
      z2v <- z2[[m]]

      max_abs1 <- max(abs(z1v), na.rm = TRUE)
      max_abs2 <- max(abs(z2v), na.rm = TRUE)

      z1_norm <- if (is.finite(max_abs1) && max_abs1 > 0) z1v / max_abs1 else rep(0, length(z1v))
      z2_norm <- if (is.finite(max_abs2) && max_abs2 > 0) z2v / max_abs2 else rep(0, length(z2v))

      z_norm_diffv <- z2_norm - z1_norm

      diffs_long_list[[length(diffs_long_list) + 1]] <- data.frame(
        Gene_ID = gcommon,
        pair = pair_key,
        attribute = lbl,
        net_2 = val2,
        net_1 = val1,
        diff = diffv,
        z_norm_diff = z_norm_diffv
      )

      df_wide[[paste0(m, "_diff")]] <- diffv
    }

    diffs_wide[[pair_key]] <- df_wide

    summary_list[[pair_key]] <- do.call(
      rbind,
      lapply(seq_along(attributes), function(i) {
        m <- attributes[i]
        lbl <- attr_labels[i]
        v <- df_wide[[paste0(m, "_diff")]]
        v <- v[is.finite(v)]
        data.frame(
          pair = pair_key,
          attribute = lbl,
          n_genes = length(v),
          mean_abs = mean(abs(v), na.rm = TRUE),
          median_abs = stats::median(abs(v), na.rm = TRUE),
          sd_abs = stats::sd(v, na.rm = TRUE),
          mad_abs = stats::mad(v, constant = 1.4826, na.rm = TRUE),
          p95_abs = as.numeric(stats::quantile(abs(v), 0.95, na.rm = TRUE)),
          max_abs = max(abs(v), na.rm = TRUE)
        )
      })
    )
  }

  diffs_long <- if (length(diffs_long_list)) do.call(rbind, diffs_long_list) else data.frame()
  summary_by_pair_attribute <- if (length(summary_list)) do.call(rbind, summary_list) else data.frame()

  diffs_long$Gene_ID <- factor(diffs_long$Gene_ID, levels = gene_order)
  diffs_long$attribute <- factor(diffs_long$attribute, levels = attr_labels)
  diffs_long$pair <- factor(diffs_long$pair, levels = pair_levels)
  diffs_long <- diffs_long[order(diffs_long$Gene_ID, diffs_long$attribute, diffs_long$pair), ]

  diffs_wide <- if (length(diffs_wide)) {
    df <- do.call(rbind, lapply(names(diffs_wide), function(nm) {
      cbind(pair = nm, diffs_wide[[nm]])
    }))
    df <- df[, c("Gene_ID", "pair", setdiff(colnames(df), c("Gene_ID", "pair")))]
    df$Gene_ID <- factor(df$Gene_ID, levels = gene_order)
    df$pair <- factor(df$pair, levels = pair_levels)
    df[order(df$Gene_ID, df$pair), ]
  } else data.frame()

  summary_by_pair_attribute$attribute <- factor(summary_by_pair_attribute$attribute, levels = attr_labels)
  summary_by_pair_attribute$pair <- factor(summary_by_pair_attribute$pair, levels = pair_levels)
  summary_by_pair_attribute <- summary_by_pair_attribute[order(summary_by_pair_attribute$attribute, summary_by_pair_attribute$pair), ]

  sum_tmp <- diffs_long[diffs_long$attribute %in% ranking_attr_labels, ]

  sum_pair <- aggregate(z_norm_diff ~ pair + Gene_ID, data = sum_tmp, FUN = sum)
  colnames(sum_pair)[3] <- "sum_pair"

  sum_global <- aggregate(z_norm_diff ~ Gene_ID, data = sum_tmp, FUN = sum)
  colnames(sum_global)[2] <- "sum_global"

  sum_summary <- merge(sum_pair, sum_global[, "Gene_ID", drop = FALSE], by = "Gene_ID")
  sum_summary$pair <- factor(sum_summary$pair, levels = pair_levels)

  sum_summary$pair_chr <- as.character(sum_summary$pair)

  for (i in seq_along(nets)) {
    net <- nets[i]
    net_low <- nets_lower[i]
    included_rows <- grepl(net, sum_summary$pair_chr)

    tmp <- aggregate(sum_pair ~ Gene_ID, data = subset(sum_summary, included_rows), FUN = sum)
    colnames(tmp)[2] <- paste0("sum_regarding_", net_low)

    sum_summary <- merge(sum_summary, tmp, by = "Gene_ID")
  }

  sum_summary$pair_chr <- NULL
  first_regarding <- grep("^sum_regarding_", names(sum_summary), value = TRUE)[1]
  sum_summary <- sum_summary[order(-abs(sum_summary[[first_regarding]]), sum_summary$Gene_ID, sum_summary$pair), ]

  list(
    diffs_long = diffs_long,
    diffs_wide = diffs_wide,
    summary_by_pair_attribute = summary_by_pair_attribute,
    sum_summary = sum_summary
  )
}

# Results summaries

load_summaries <- function(path = "."){
  files <- list.files(
    path,
    pattern = ".*_(contrast|enrichment)_summary\\.rds$|.*_de_genes\\.rds$|.*_de_tectivirus_genes\\.rds$",
    full.names = TRUE
  )
  out <- list()
  for(f in files){
    obj <- readRDS(f)
    name <- sub("\\.rds$", "", basename(f))
    assign(name, obj, envir = .GlobalEnv)
    out[[name]] <- obj
  }
  invisible(out)
}

build_functional_enrichments <- function(
  dir = ".",
  envir = .GlobalEnv,
  overwrite = TRUE
){

  cog_files <- list.files(
    path = dir,
    pattern = "_cog\\.rds$",
    full.names = TRUE
  )

  kegg_files <- list.files(
    path = dir,
    pattern = "_kegg\\.rds$",
    full.names = TRUE
  )

  append_obj <- function(lst, obj){

    if(is.null(obj)) return(lst)

    if(!is.data.frame(obj)) return(lst)

    if(nrow(obj) == 0) return(lst)

    lst[[length(lst) + 1]] <- obj
    lst
  }

  lysogeny_cog_list <- list()
  intra_cog_list <- list()
  inter_cog_list <- list()

  for(f in cog_files){

    obj <- readRDS(f)
    grp <- attr(obj, "group")

    if(identical(grp, "Lysogeny")){

      lysogeny_cog_list <- append_obj(
        lysogeny_cog_list,
        obj
      )

    } else if(grp %in% c(
      "GBJ002",
      "GIL01",
      "GIL16"
    )){

      intra_cog_list <- append_obj(
        intra_cog_list,
        obj
      )

    } else if(grp %in% c(
      "t0",
      "t10",
      "t30",
      "t60",
      "t60-mock"
    )){

      inter_cog_list <- append_obj(
        inter_cog_list,
        obj
      )
    }
  }

  lysogeny_kegg_list <- list()
  intra_kegg_list <- list()
  inter_kegg_list <- list()

  for(f in kegg_files){

    obj <- readRDS(f)
    grp <- attr(obj, "group")

    if(identical(grp, "Lysogeny")){

      lysogeny_kegg_list <- append_obj(
        lysogeny_kegg_list,
        obj
      )

    } else if(grp %in% c(
      "GBJ002",
      "GIL01",
      "GIL16"
    )){

      intra_kegg_list <- append_obj(
        intra_kegg_list,
        obj
      )

    } else if(grp %in% c(
      "t0",
      "t10",
      "t30",
      "t60",
      "t60-mock"
    )){

      inter_kegg_list <- append_obj(
        inter_kegg_list,
        obj
      )
    }
  }

  lysogeny_cog <- if(length(lysogeny_cog_list)){
    do.call(rbind, lysogeny_cog_list)
  } else {
    data.frame()
  }

  intra_cog <- if(length(intra_cog_list)){
    do.call(rbind, intra_cog_list)
  } else {
    data.frame()
  }

  inter_cog <- if(length(inter_cog_list)){
    do.call(rbind, inter_cog_list)
  } else {
    data.frame()
  }

  lysogeny_kegg <- if(length(lysogeny_kegg_list)){
    do.call(rbind, lysogeny_kegg_list)
  } else {
    data.frame()
  }

  intra_kegg <- if(length(intra_kegg_list)){
    do.call(rbind, intra_kegg_list)
  } else {
    data.frame()
  }

  inter_kegg <- if(length(inter_kegg_list)){
    do.call(rbind, inter_kegg_list)
  } else {
    data.frame()
  }

  if(overwrite || !exists("lysogeny_cog", envir = envir, inherits = FALSE)){
    assign("lysogeny_cog", lysogeny_cog, envir = envir)
  }

  if(overwrite || !exists("intra_cog", envir = envir, inherits = FALSE)){
    assign("intra_cog", intra_cog, envir = envir)
  }

  if(overwrite || !exists("inter_cog", envir = envir, inherits = FALSE)){
    assign("inter_cog", inter_cog, envir = envir)
  }

  if(overwrite || !exists("lysogeny_kegg", envir = envir, inherits = FALSE)){
    assign("lysogeny_kegg", lysogeny_kegg, envir = envir)
  }

  if(overwrite || !exists("intra_kegg", envir = envir, inherits = FALSE)){
    assign("intra_kegg", intra_kegg, envir = envir)
  }

  if(overwrite || !exists("inter_kegg", envir = envir, inherits = FALSE)){
    assign("inter_kegg", inter_kegg, envir = envir)
  }

  invisible(NULL)
}

combine_summaries <- function(type = c("contrast", "enrichment"), group = c("strains", "time_points"), filter = FALSE){
  type <- match.arg(type)
  group <- match.arg(group)
  
  objs <- ls(.GlobalEnv, pattern = paste0("_", type, "_summary$"))
  
  obj_groups <- sapply(objs, function(nm) attr(get(nm, envir=.GlobalEnv), "group"))
  
  if(group == "strains"){
    valid_vals <- c("Lysogeny", "t0", "t10", "t30", "t60", "t60-mock")       
  } else if(group == "time_points"){
    valid_vals <- c("GBJ002", "GIL01", "GIL16") 
  } else {
    stop("Group must be 'strains' or 'time_points'.")
  }
  
  selected <- objs[obj_groups %in% valid_vals]
  
  if(!identical(filter, FALSE)){
    selected <- selected[
      sapply(obj_groups[selected], function(grp)
        any(sapply(filter, function(f) grepl(f, grp, ignore.case = TRUE)))
      )
    ]
  }
  
  dfs <- list()
  
  for(nm in selected){
    obj <- get(nm, envir = .GlobalEnv)
    grp <- attr(obj, "group")
    contrast_val <- attr(obj, "contrast")
    de_vals <- sub(".*\\.", "", rownames(obj))
    
    obj$Group <- grp
    obj$Contrast <- contrast_val
    obj$DE <- de_vals
    obj <- obj[, c("Group", "Contrast", "DE", setdiff(colnames(obj), c("Group", "Contrast", "DE")))]
    rownames(obj) <- NULL
    obj$Group <- factor(obj$Group)
    obj$DE <- factor(obj$DE, levels = unique(obj$DE))
    
    dfs[[nm]] <- obj
  }
  
  if(length(dfs) == 0) return(NULL)
  
  out <- do.call(rbind, dfs)
  rownames(out) <- NULL

  if(!identical(filter, FALSE) && length(filter) > 1){
    out$Group <- factor(out$Group, levels = filter)
  } else {
    out$Group <- factor(out$Group)
  }
  
  if(group == "time_points"){
    desired_order <- c(
      "Time points 10 vs 0",
      "Time points 30 vs 0",
      "Time points 60 vs 0",
      "Time points 60 vs 60-mock",
      "Time points 60-mock vs 0"
    )
  } else {
    desired_order <- c(
      "Strains GIL01 vs GBJ002",
      "Strains GIL16 vs GBJ002",
      "Strains GIL16 vs GIL01"
    )
  }
  
  valid_levels <- intersect(desired_order, unique(as.character(out$Contrast)))
  out$Contrast <- factor(out$Contrast, levels = valid_levels)

  out <- out[order(out$Group, out$Contrast, out$DE), ]
  rownames(out) <- NULL
  
  out
}

venn_contrasts <- function(group, de_genes = c("UP", "DOWN"), label_alpha) {
  de_genes <- match.arg(de_genes)
  gene_slot <- if (de_genes == "UP") "up_genes" else "down_genes"
  high_color <- if (de_genes == "UP") "red" else "blue"
  
  obj_names <- ls(envir = .GlobalEnv, pattern = "_de_genes$")
  
  objects <- lapply(obj_names, function(x) get(x, envir = .GlobalEnv))
  names(objects) <- obj_names
  
  objects <- Filter(function(x) {
    !is.null(attr(x, "group")) && attr(x, "group") == group
  }, objects)
  
  if (length(objects) < 2) {
    stop("Not enough contrasts found for the specified group.")
  }
  
  gene_sets <- lapply(objects, function(x) x[[gene_slot]]$Gene_ID)
  
  set_names <- sapply(objects, function(x) {
    lbl <- attr(x, "contrast")
    if (is.null(lbl)) NA_character_ else lbl
  })
  
  names(gene_sets) <- set_names
  
  ggVennDiagram(gene_sets, label_alpha = label_alpha) +
    scale_fill_gradient(low = "grey90", high = high_color, name = ifelse(de_genes == "UP", "UP genes", "DOWN genes")) +
    coord_cartesian(clip = "off") +
    theme(legend.title = element_text(face = "bold"), plot.margin = margin(30, 90, 30, 90))
}

upset_contrasts <- function(group, de_genes = c("UP", "DOWN"), order = NULL) {
  de_genes <- match.arg(de_genes)
  gene_slot <- if (de_genes == "UP") "up_genes" else "down_genes"
  
  obj_names <- ls(envir = .GlobalEnv, pattern = "_de_genes$")
  objects <- lapply(obj_names, function(x) get(x, envir = .GlobalEnv))
  names(objects) <- obj_names
  
  objects <- Filter(function(x) {
    !is.null(attr(x, "group")) && attr(x, "group") %in% group
  }, objects)
  
  if (length(objects) < 2) {
    stop("Not enough contrasts found for the specified group(s).")
  }
  
  gene_sets <- lapply(objects, function(x) x[[gene_slot]]$Gene_ID)
  
  if (length(group) == 1) {
    set_names <- vapply(objects, function(x) {
      attr(x, "contrast")
    }, character(1))
  } else {
    set_names <- vapply(objects, function(x) {
      grp <- attr(x, "group")
      grp_label <- sub("^t", "Time point ", grp)
      contrast <- attr(x, "contrast")
      paste(grp_label, contrast, sep = " - ")
    }, character(1))
  }
  
  names(gene_sets) <- set_names
  
  if (!is.null(order)) {
    if (!all(order %in% names(gene_sets))) {
      stop("Some values in 'order' are not valid contrast names.")
    }
    gene_sets <- gene_sets[rev(order)]
  }
  
  upset(
    fromList(gene_sets),
    sets = names(gene_sets),
    keep.order = TRUE,
    order.by = "freq",
    mainbar.y.label = "Shared genes",
    sets.x.label = "DE genes per contrast",
    main.bar.color = ifelse(de_genes == "UP", "red", "blue"),
    sets.bar.color = ifelse(de_genes == "UP", "red", "blue"),
    matrix.color = ifelse(de_genes == "UP", "red", "blue")
  )
}

create_de_genes_table <- function(group = c("strains", "time_points"), filter = FALSE){

  group <- match.arg(group)
  objs <- ls(.GlobalEnv, pattern = "_de_genes$")
  
  obj_groups <- sapply(objs, function(nm) attr(get(nm, envir = .GlobalEnv), "group"))
  
  if(group == "strains"){
    valid_vals <- c("Lysogeny", "t0", "t10", "t30", "t60", "t60-mock")
  } else if(group == "time_points"){
    valid_vals <- c("GBJ002", "GIL01", "GIL16")
  } else {
    stop("Group must be 'strains' or 'time_points'.")
  }
  
  selected <- objs[obj_groups %in% valid_vals]
  
  if(!identical(filter, FALSE)){
    selected <- selected[
      sapply(obj_groups[selected], function(grp)
        any(sapply(filter, function(f) grepl(f, grp, ignore.case = TRUE)))
      )
    ]
  }
  
  if(length(selected) == 0) return(NULL)
  
  dfs <- list()
  
  for(nm in selected){
    obj <- get(nm, envir = .GlobalEnv)
    grp <- attr(obj, "group")
    contrast <- attr(obj, "contrast")
    
    df_up <- obj$up_genes
    df_down <- obj$down_genes
    
    df <- rbind(df_up, df_down)
    
    if(nrow(df) == 0) next
    
    colnames(df) <- c("Gene_ID", "LFC")
    df$Contrast <- paste(grp, contrast, sep = " - ")
    
    dfs[[nm]] <- df
  }
  
  long_df <- do.call(rbind, dfs)
  
  wide_df <- reshape(long_df, timevar = "Contrast", idvar = "Gene_ID", direction = "wide")
  
  colnames(wide_df) <- sub("^LFC\\.", "", colnames(wide_df))
  
  if(exists("annotation", inherits = TRUE)){
    annot <- get("annotation", inherits = TRUE)
    annot_sel <- annot[, c("Gene_ID", "Gene_symbol", "Description")]
    
    final_df <- merge(annot_sel, wide_df, by = "Gene_ID", all.y = TRUE)
    order_idx <- match(final_df$Gene_ID, annot_sel$Gene_ID)
    final_df <- final_df[order(order_idx), ]
  } else {
    final_df <- wide_df
  }
  
  lfc_cols <- setdiff(colnames(final_df), c("Gene_ID", "Gene_symbol", "Description"))
  keep <- rowSums(!is.na(final_df[, lfc_cols])) > 0
  
  final_df <- final_df[keep, ]
  
  if(group == "time_points"){
    
    contrast_order <- c("GBJ002", "GIL01", "GIL16")
    
    desired_order <- c(
      "Time points 10 vs 0",
      "Time points 30 vs 0",
      "Time points 60 vs 0",
      "Time points 60 vs 60-mock",
      "Time points 60-mock vs 0"
    )
    
  } else {
    
    contrast_order <- if(!identical(filter, FALSE)) filter else unique(sub(" - .*", "",lfc_cols))
    
    desired_order <- c(
      "Strains GIL01 vs GBJ002",
      "Strains GIL16 vs GBJ002",
      "Strains GIL16 vs GIL01"
    )
  }
  
  lfc_cols <- setdiff(
    colnames(final_df),
    c("Gene_ID", "Gene_symbol", "Description")
  )
  
  group_tag <- sub(" - .*", "", lfc_cols)
  base_names <- sub(".* - ", "", lfc_cols)
  
  ord <- order(
    match(group_tag, contrast_order),
    match(base_names, desired_order),
    lfc_cols
  )
  
  final_df <- final_df[, c(
    "Gene_ID", "Gene_symbol", "Description",
    lfc_cols[ord]
  )]
  
  rownames(final_df) <- NULL
  
  final_df
}

create_de_tectivirus_genes_table <- function(){

  counts <- readRDS(file = "lysogeny_counts.rds")
  
  valid_vals <- c(
    "Lysogeny",
    "t0",
    "t10",
    "t30",
    "t60",
    "t60-mock"
  )

  tectivirus_annotation <- combine_tectivirus_genes()$tectivirus_annotation

  objs <- ls(.GlobalEnv, pattern = "_de_tectivirus_genes$")

  obj_groups <- sapply(
    objs,
    function(nm) attr(get(nm, envir = .GlobalEnv), "group")
  )

  selected <- objs[obj_groups %in% valid_vals]

  if(length(selected) == 0) return(NULL)

  dfs <- list()

  for(nm in selected){

    obj <- get(nm, envir = .GlobalEnv)
    grp <- attr(obj, "group")

    df_up <- obj$up_genes
    df_down <- obj$down_genes

    df <- rbind(df_up, df_down)

    if(nrow(df) == 0) next

    colnames(df) <- c("Gene_symbol", "LFC")
    df$Group <- grp

    dfs[[nm]] <- df
  }

  if(length(dfs) == 0) return(NULL)

  long_df <- do.call(rbind, dfs)

  wide_df <- reshape(
    long_df,
    timevar = "Group",
    idvar = "Gene_symbol",
    direction = "wide"
  )

  colnames(wide_df) <- sub("^LFC\\.", "", colnames(wide_df))

  final_df <- merge(
    tectivirus_annotation[, c("Gene_symbol", "Description")],
    wide_df,
    by = "Gene_symbol",
    all.x = TRUE
  )

  order_idx <- match(
    final_df$Gene_symbol,
    tectivirus_annotation$Gene_symbol
  )

  final_df <- final_df[order(order_idx), ]

  lfc_cols <- intersect(valid_vals, colnames(final_df))

  keep <- rowSums(!is.na(final_df[, lfc_cols, drop = FALSE])) > 0

  final_df <- final_df[keep, ]

  missing_cols <- setdiff(valid_vals, colnames(final_df))

  for(col in missing_cols){
    final_df[[col]] <- NA_real_
  }

  final_df <- final_df[, c(
    "Gene_symbol",
    "Description",
    valid_vals
  )]

  rownames(final_df) <- NULL

  final_df
}

# Other functions

mark_pBtic235_genes <- function(genes, annotation) {
  sapply(genes, function(g) {
    chr <- annotation$Element[annotation$Gene_ID == g]
    if (length(chr) > 0 && chr == "pBtic235") {
      paste0(g, "*")
    } else {
      g
    }
  })
}

volcano_plot <- function(contrast_df, limits = NULL, tectivirus = FALSE) {

  if (is.null(limits)) {
    max_lfc <- ceiling(max(abs(contrast_df$LFC), na.rm = TRUE))
    limits <- c(-max_lfc, max_lfc)
  }

  if (tectivirus) {
    contrast_df$Gene_ID <- contrast_df$Gene_symbol
  }

  contrast_df$gene <- contrast_df$Gene_ID

  contrast_df$direction <- "N. S."
  contrast_df$direction[contrast_df$Adjusted_p_value < 0.05 & contrast_df$LFC > log2(1.5)] <- "UP"
  contrast_df$direction[contrast_df$Adjusted_p_value < 0.05 & contrast_df$LFC < -log2(1.5)] <- "DOWN"

  contrast_df$direction <- factor(contrast_df$direction, levels = c("UP", "DOWN", "N. S."))

  contrast_df <- contrast_df %>%
    dplyr::mutate(
      Gene_symbol = dplyr::na_if(Gene_symbol, ""),
      Name_display = dplyr::case_when(
        !is.na(Gene_symbol) & !is.na(Description) ~ paste0(Description, " (", Gene_symbol, ")"),
        !is.na(Gene_symbol) & is.na(Description) ~ Gene_symbol,
        TRUE ~ Description
      )
    )

  if (!tectivirus) {
    contrast_df <- contrast_df %>%
      dplyr::left_join(
        annotation %>% dplyr::select(Gene_ID, Element),
        by = c("gene" = "Gene_ID")
      ) %>%
      dplyr::mutate(
        shape_type = dplyr::case_when(
          Element == "chromosome" ~ "chromosome",
          Element == "pBtic235" ~ "pBtic235",
          Element == "GIL01" ~ "GIL01",
          Element == "GIL16" ~ "GIL16",
          TRUE ~ NA_character_
        ),
        shape_type = factor(
          shape_type,
          levels = c("chromosome", "pBtic235", "GIL01", "GIL16")
        )
      )
  }

  contrast_df <- contrast_df %>%
    dplyr::mutate(
      hover_text = if (tectivirus) {
        paste0(
          ifelse(!is.na(Name_display), Name_display, ""),
          "<br>LFC: ", round(LFC, 2),
          "<br>padj: ", signif(Adjusted_p_value, 3)
        )
      } else {
        paste0(
          gene,
          ifelse(!is.na(Name_display), paste0("<br>", Name_display), ""),
          "<br>LFC: ", round(LFC, 2),
          "<br>padj: ", signif(Adjusted_p_value, 3)
        )
      }
    )

  if (tectivirus) {

    p <- ggplot(
      data = contrast_df,
      mapping = aes(
        x = LFC,
        y = -log10(Adjusted_p_value),
        color = direction,
        text = hover_text
      )
    ) +
      geom_point(size = 2, alpha = 0.5) +
      scale_color_manual(
        values = c(
          "UP" = "#E41A1C",
          "DOWN" = "#377EB8",
          "N. S." = "lightgrey"
        )
      ) +
      labs(
        x = "Fold change (log₂)",
        y = "Adjusted p value (-log₁₀)",
        color = NULL
      ) +
      scale_x_continuous(limits = limits) +
      theme_linedraw(base_size = 14) +
      theme(panel.grid = element_blank())

  } else {

    p <- ggplot(
      data = contrast_df,
      mapping = aes(
        x = LFC,
        y = -log10(Adjusted_p_value),
        color = direction,
        shape = shape_type,
        text = hover_text
      )
    ) +
      geom_point(size = 2, alpha = 0.5) +
      scale_color_manual(
        values = c(
          "UP" = "#E41A1C",
          "DOWN" = "#377EB8",
          "N. S." = "lightgrey"
        )
      ) +
      scale_shape_manual(
        values = c(
          "chromosome" = 19,
          "pBtic235" = 15,
          "GIL01" = 17,
          "GIL16" = 17
        )
      ) +
      labs(
        x = "Fold change (log₂)",
        y = "Adjusted p value (-log₁₀)",
        color = NULL,
        shape = NULL
      ) +
      scale_x_continuous(limits = limits) +
      theme_linedraw(base_size = 14) +
      theme(panel.grid = element_blank())
  }

  ggplotly(p, tooltip = "text")
}

combine_tectivirus_genes <- function() {

  gil01_annotation <- annotation %>%
    dplyr::filter(
      grepl("^GIL01_", Gene_ID),
      !is.na(Classification),
      Classification != ""
    ) %>%
    dplyr::select(Gene_ID, Classification)

  gil16_annotation <- annotation %>%
    dplyr::filter(
      grepl("^GIL16_", Gene_ID),
      !is.na(Classification),
      Classification != ""
    ) %>%
    dplyr::select(Gene_ID, Classification)

  homology_map <- gil01_annotation %>%
    dplyr::left_join(
      gil16_annotation,
      by = "Classification",
      suffix = c("_GIL01", "_GIL16")
    ) %>%
    dplyr::filter(
      !is.na(Gene_ID_GIL01),
      !is.na(Gene_ID_GIL16),
      Gene_ID_GIL01 %in% rownames(counts),
      Gene_ID_GIL16 %in% rownames(counts)
    )

  dnt_row <- data.frame(
    Classification = "dnt",
    Gene_ID_GIL01 = annotation$Gene_ID[grepl("^GIL01-DNT_", annotation$Gene_ID)][1],
    Gene_ID_GIL16 = annotation$Gene_ID[grepl("^GIL16-DNT_", annotation$Gene_ID)][1]
  )

  dnt_pos <- which(grepl("^GIL01-DNT_", annotation$Gene_ID))[1]

  insert_after <- sum(
    match(homology_map$Gene_ID_GIL01, annotation$Gene_ID) < dnt_pos
  )

  homology_map <- dplyr::bind_rows(
    homology_map[seq_len(insert_after), ],
    dnt_row,
    homology_map[(insert_after + 1):nrow(homology_map), ]
  )

  homology_map$gene_name <- ifelse(
    homology_map$Classification == "dnt",
    "dnt",
    paste0(
      "g",
      sprintf("%02d", as.integer(homology_map$Classification))
    )
  )

  tectivirus_annotation <- homology_map %>%
    dplyr::filter(Classification != "dnt") %>%
    dplyr::left_join(
      annotation %>%
        dplyr::select(
          Gene_ID,
          Description,
          Type,
          Subtype,
          Source
        ),
      by = c("Gene_ID_GIL01" = "Gene_ID")
    ) %>%
    dplyr::transmute(
      Gene_symbol = gene_name,
      Description,
      Type,
      Subtype,
      Source
    )

  tectivirus_annotation <- homology_map %>%
    dplyr::left_join(
      annotation %>%
        dplyr::select(
          Gene_ID,
          Description,
          Type,
          Subtype,
          Source
        ),
      by = c("Gene_ID_GIL01" = "Gene_ID")
    ) %>%
    dplyr::transmute(
      Gene_symbol = gene_name,
      Description,
      Type,
      Subtype,
      Source
    )

  tectivirus_counts <- do.call(
    rbind,
    lapply(seq_len(nrow(homology_map)), function(i) {

      x <- rep(0, ncol(counts))
      names(x) <- colnames(counts)

      gene_gil01 <- homology_map$Gene_ID_GIL01[i]
      gene_gil16 <- homology_map$Gene_ID_GIL16[i]

      cols_gil01 <- grep("^GIL01_", colnames(counts))
      cols_gil16 <- grep("^GIL16_", colnames(counts))

      x[cols_gil01] <- as.numeric(counts[gene_gil01, cols_gil01])
      x[cols_gil16] <- as.numeric(counts[gene_gil16, cols_gil16])

      x
    })
  )

  rownames(tectivirus_counts) <- homology_map$gene_name
  tectivirus_counts <- as.data.frame(tectivirus_counts)

  genes_to_replace <- unique(c(
    homology_map$Gene_ID_GIL01,
    homology_map$Gene_ID_GIL16
  ))

  background_counts <- counts[
    !rownames(counts) %in% genes_to_replace,
    ,
    drop = FALSE
  ]

  tectivirus_counts <- rbind(
    background_counts,
    tectivirus_counts
  )

  list(
    tectivirus_counts = tectivirus_counts,
    tectivirus_annotation = tectivirus_annotation
  )
}

concat_or_replace <- function(old, new) {
  ifelse(is.na(old), new,
         ifelse(is.na(new), old,
                paste(old, new, sep = " / ")))
}

get_matrix <- function(matrix,
                       genes,
                       name,
                       neighbours = FALSE,
                       neighbours_edges = "all",
                       threshold = NULL) {

  if (!is.null(threshold)) {
    matrix[matrix > threshold] <- 0
  }

  genes <- intersect(genes, rownames(matrix))

  if (length(genes) == 0) {
    stop("None of the genes are present in the matrix.")
  }

  all_genes <- genes

  if (neighbours) {
    neighbour_genes <- colnames(matrix)[
      colSums(matrix[genes, , drop = FALSE] != 0) > 0
    ]
    all_genes <- unique(c(genes, neighbour_genes))
  }

  submat <- matrix[all_genes, all_genes, drop = FALSE]

  if (neighbours_edges != "all") {

    neighbour_set <- setdiff(all_genes, genes)

    if (neighbours_edges == "genes") {
      submat[neighbour_set, neighbour_set] <- 0
    }

    if (neighbours_edges == "neighbours") {

      keep <- matrix(FALSE, nrow(submat), ncol(submat),
                     dimnames = dimnames(submat))

      keep[genes, genes] <- TRUE
      keep[genes, neighbour_set] <- TRUE
      keep[neighbour_set, genes] <- TRUE

      for (g in genes) {
        neigh_g <- rownames(matrix)[matrix[g, ] != 0]
        neigh_g <- intersect(neigh_g, all_genes)
        if (length(neigh_g) > 1) {
          keep[neigh_g, neigh_g] <- TRUE
        }
      }

      submat[!keep] <- 0
    }
  }

  edges <- which(submat != 0, arr.ind = TRUE)
  edges <- edges[edges[,1] < edges[,2], ]

  edge_list <- data.frame(
    source = rownames(submat)[edges[,1]],
    target = colnames(submat)[edges[,2]],
    weight = submat[edges]
  )

  ordered_levels <- c(genes, setdiff(all_genes, genes))

  edge_list$source <- factor(edge_list$source, levels = ordered_levels)
  edge_list$target <- factor(edge_list$target, levels = ordered_levels)

  edge_list <- edge_list[order(edge_list$source, edge_list$target), ]

  output_file <- paste0(name, ".tsv")

  write.table(edge_list,
              file = output_file,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)

  message("File successfully saved as ", output_file)
}