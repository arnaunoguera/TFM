library(ComplexHeatmap)
library(phyloseq)
library(randomcoloR)
library(tidyverse)
library(reshape2)
library(ggtext)
library(Hmisc)
library(corrplot)
library(circlize)


plot_rpkm_heatmap <- function(phy, cluster_rows = F, colside_vars = NULL, colside_cols = NULL, rowside_vars = NULL, rowside_cols = NULL,
                              col_title = NA, row_title = NA, axis_name = "log RPKM\nabundance", rsize = 0, csize = 0,
                              cluster_columns = T, row_dend_width = 0.1, clustering_distance_columns = "manhattan",
                              clustering_method_columns = "complete", ...){
  # Function to plot a heatmap from a phyloseq object for plasmids
  plot_matrix <- as.matrix.data.frame(data.frame(otu_table(phy))) %>% t()
  # Prepare rowside variables to plot
  n_colors <- 0
  r_ha <- c()
  c_ha <- c()
  if (! is.null(rowside_vars)) {
    # Get how many rowside variables we have
    n_row_vars <- length(rowside_vars)
    for (n in 1:n_row_vars) {
      # For each, get how many colors there are
      values <- sample_data(phy) %>% data.frame() %>% pull(rowside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(rowside_vars)[n]
      if (is_null(var_name) || var_name == "") {
        var_name <- rowside_vars[n]
      }
      # Prepare palette
      if (is.null(rowside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- rowside_cols[[var_name]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      r_ha <- append(r_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        sample_data(phy)[, rowside_vars[n]] %>% data.frame() %>% rename(setNames(object = rowside_vars[n], nm = var_name))
          %>% replace_na(replace = setNames(object = list("NA"), nm = c(rowside_vars[n]))), 
        list(col = cols, which = "row")
      )))
    }
  }
  # Do the same for the colside variables 
  if (! is.null(colside_vars)) {
    # Get how many colside variables we have
    n_col_vars <- length(colside_vars)
    for (n in 1:n_col_vars) {
      # For each, get how many colors there are
      values <- tax_table(phy) %>% data.frame() %>% pull(colside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(colside_vars)[n]
      if (var_name == "" || is.null(var_name)) {
        var_name <- colside_vars[n]
      }
      # Prepare palette
      if (is.null(colside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- colside_cols[[var_name]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      c_ha <- append(c_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        tax_table(phy)[, colside_vars[n]] %>% as.data.frame() %>% rename(setNames(object = colside_vars[n], nm = var_name))
          %>% replace_na(replace = setNames(object = list("NA"), nm = c(colside_vars[n]))),
        list(col = cols, which = "col")
      )))
    }
  }
  # Plot the heatmap 
  Heatmap(plot_matrix,
          # cluster rows will either be F or the dendrogram
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          row_names_gp = gpar(fontsize = rsize),
          column_names_gp = gpar(fontsize = csize),
          #clustering_distance_rows = braycurtis,
          #clustering_method_rows = "average",
          clustering_distance_columns = clustering_distance_columns,
          clustering_method_columns = clustering_method_columns,
          right_annotation = r_ha, 
          top_annotation = c_ha,
          name = axis_name, column_title = col_title,
          row_title = row_title, 
          use_raster = F,
          row_dend_width = unit(row_dend_width, units = "npc"),
          col = colorRamp2(c(-10, -6, -3, 0, 2, 5, 7), c("black", viridis_pal(begin = 0, end = 1, direction = 1)(6))),
          #col = colorRamp2(c(0, 1000, 2000, 3000, 4000, 5000, 8000), c("black", "blue", "cyan", "green", "yellow","orange", "red")),
          ...
  )
}


plot_coverage_heatmap <- function(phy, cov_matrix, cluster_rows = F, colside_vars = NULL, colside_cols = NULL, rowside_vars = NULL, rowside_cols = NULL,
                              col_title = NA, row_title = NA, axis_name = "Detection", rsize = 0, csize = 0,
                              cluster_columns = T, row_dend_width = 0.1, clustering_distance_columns = "manhattan",
                              colorpalette = colorRamp2(c(0, 50, 100), c("white", "white", "green4")),
                              clustering_method_columns = "complete", ...){
  # Function to plot a heatmap from a coverage matrix for plasmids
  plot_matrix <- cov_matrix %>% t()
  # Prepare rowside variables to plot
  n_colors <- 0
  r_ha <- c()
  c_ha <- c()
  if (! is.null(rowside_vars)) {
    # Get how many rowside variables we have
    n_row_vars <- length(rowside_vars)
    for (n in 1:n_row_vars) {
      # For each, get how many colors there are
      values <- sample_data(phy) %>% data.frame() %>% pull(rowside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(rowside_vars)[n]
      if (var_name == "") {
        var_name <- rowside_vars[n]
      }
      # Prepare palette
      if (is.null(rowside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- rowside_cols[[n]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      r_ha <- append(r_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        sample_data(phy)[, rowside_vars[n]] %>% as.data.frame() %>% rename(setNames(object = rowside_vars[n], nm = var_name))
          %>% replace_na(replace = setNames(object = list("NA"), nm = c(rowside_vars[n]))),
        list(col = cols, which = "row")
      )))
    }
  }
  # Do the same for the colside variables 
  if (! is.null(colside_vars)) {
    # Get how many colside variables we have
    n_col_vars <- length(colside_vars)
    for (n in 1:n_col_vars) {
      # For each, get how many colors there are
      values <- tax_table(phy) %>% data.frame() %>% pull(colside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(colside_vars)[n]
      if (var_name == "" || is.null(var_name)) {
        var_name <- colside_vars[n]
      }
      # Prepare palette
      if (is.null(colside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- colside_cols[[n]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      c_ha <- append(c_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        tax_table(phy)[, colside_vars[n]] %>% as.data.frame() %>% rename(setNames(object = colside_vars[n], nm = var_name))
        %>% replace_na(replace = setNames(object = list("NA"), nm = c(colside_vars[n]))),
        list(col = cols, which = "col")
      )))
    }
  }
  # Plot the heatmap 
  Heatmap(plot_matrix,
          # cluster rows will either be F or the dendrogram
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          row_names_gp = gpar(fontsize = rsize),
          column_names_gp = gpar(fontsize = csize),
          #clustering_distance_rows = braycurtis,
          #clustering_method_rows = "average",
          clustering_distance_columns = clustering_distance_columns,
          clustering_method_columns = clustering_method_columns,
          right_annotation = r_ha, 
          top_annotation = c_ha,
          name = axis_name, column_title = col_title,
          row_title = row_title, 
          use_raster = F,
          #row_dend_width = unit(row_dend_width, units = "npc"),
          col = colorpalette, 
          ...
  )
}

plot_taxa_abundances_rel <- function(phy_obj, facet_condition = F, taxonomic_rank = "phylum", plotlegend = T, 
                                     n = NULL, tax_glom_done = F, legend_label = "Taxa", italics = F, palette = NULL, 
                                     facet_space = "fixed", top_func = sum, drop_none = F, non_unique = F, clean_names = T, 
                                     return_abundance = F){
  #extract data from physeq object
  if (tax_glom_done == T) {
    phy_subset <- phy_obj
    if (drop_none) {
      phy_subset <- subset_taxa(phy_subset, get(taxonomic_rank) != "none")
    }
  } else {
    phy_subset <- tax_glom(phy_obj, taxrank = taxonomic_rank)
  }
  current_abundance <- data.frame(otu_table(phy_subset))
  current_metadata <- data.frame(sample_data(phy_subset))
  
  if (non_unique){
    #rename taxa
    rownames(current_abundance) <- paste0(1:nrow(current_abundance), tax_table(phy_subset)[rownames(current_abundance), taxonomic_rank])
  } else {
    #rename taxa
    rownames(current_abundance) <- tax_table(phy_subset)[rownames(current_abundance), taxonomic_rank]
  }
  relative_abundance <- apply(current_abundance, 2, function(column) column / sum(column))
  
  if (!is_null(n)) {
    # if n is specified, keep only the top n taxa with the most global abundance and group the rest at "Others"
    if (n < nrow(relative_abundance)) {
      # Calculate plasmid systems with the most abundance and keep the top 20
      topn_abundance <- apply(relative_abundance, 1, function(row) top_func(row)) %>% sort(decreasing = T) %>% head(n = n) %>% names()
      
      # Get the abundance of the n selected systems
      selected_abundance <- relative_abundance[topn_abundance, ]
      # Get the others, to sum 
      other <- relative_abundance[! rownames(relative_abundance) %in% topn_abundance, ] %>% 
        apply(2, function(column) sum(column))
      current_abundance <- data.frame(rbind(selected_abundance, other))
    } else {
      topn_abundance <- apply(relative_abundance, 1, function(row) top_func(row)) %>% sort(decreasing = T) %>% names()
      current_abundance <- relative_abundance[topn_abundance, ]
    }
    
  } else {
    topn_abundance <- apply(relative_abundance, 1, function(row) top_func(row)) %>% sort(decreasing = T) %>% names()
    current_abundance <- relative_abundance[topn_abundance, ]
  }
  
  #convert to long df format for plotting purposes
  current_abundance <- data.frame(t(current_abundance))
  current_abundance$sample <- rownames(current_abundance)
  if (facet_condition != F) {
    current_abundance$condition <- current_metadata[rownames(current_abundance), facet_condition]
    current_abundance <- current_abundance %>% arrange(desc(get(make.names(topn_abundance[1])))) 
    current_abundance_long <- melt(current_abundance, id.vars = c("sample", "condition"))
    colnames(current_abundance_long) <- c("sample", "condition", "taxa", "abundance")
    current_abundance_long$sample <- factor(current_abundance_long$sample, levels = unique(current_abundance_long$sample))
  } else {
    current_abundance <- current_abundance %>% arrange(desc(get(make.names(topn_abundance[1])))) 
    current_abundance_long <- melt(current_abundance, id.vars = c("sample"))
    colnames(current_abundance_long) <- c("sample", "taxa", "abundance")
    current_abundance_long$sample <- factor(current_abundance_long$sample, levels = unique(current_abundance_long$sample))
  }
  
  if (clean_names){
    current_abundance_long$taxa <- gsub(pattern = "_", replacement = " ", x = current_abundance_long$taxa, fixed = T)
    current_abundance_long$taxa <- gsub(pattern = ".", replacement = "-", x = current_abundance_long$taxa, fixed = T)
  }
  
  if (italics){
    current_abundance_long <- current_abundance_long %>% mutate(taxa = case_when(taxa == "other" ~ "other",
                                                                                 TRUE ~ paste0("*", taxa, "*"))) 
  } 
  
  current_abundance_long$taxa <- factor(current_abundance_long$taxa, levels = unique(current_abundance_long$taxa))
  
  if (is_null(palette)) {
    #set the palette to use depending on number of taxa
    library(RColorBrewer)
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    final_palette <- getPalette(length(unique(current_abundance_long$taxa)))
  } else {
    #set the palette to use depending on number of taxa
    library(RColorBrewer)
    getPalette <- colorRampPalette(palette)
    final_palette <- getPalette(length(unique(current_abundance_long$taxa)))
  }
 
  
  #plot
  p <- ggplot(current_abundance_long, aes(x  = sample, fill = taxa, y = abundance)) +
    geom_bar(stat = "identity", width = 1) +
    theme_classic() +
    # theme(axis.text.x = element_text(angle = 90),
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.line.x = element_blank(),
          axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"),
          # legend.text = element_text(size = 12, face = "bold", colour = "black"),
          axis.text.y = element_text(colour = "black", size = 8, face = "bold"),
          strip.text.x = element_text(size = 16)) +
    # scale_y_continuous(expand = c(0,0)) +
    labs(x = "", y = "Relative Abundance", fill = legend_label) +
    scale_fill_manual(values = final_palette) +
    scale_y_continuous(expand = c(0,0))
  
  if (facet_condition != F) {
    p <- p + facet_grid(~condition, scales = "free", space = facet_space )
  }
  
  if (plotlegend == F) {
    p <- p + theme(legend.position = "none")
  }
  
  if (italics) {
    p <- p + theme(legend.text = element_markdown())
  }
  
  if (return_abundance){
    return(list(plot = p, table = current_abundance_long))
  } else {
    return(p)  
  }
  
}


scree_plot <- function(pcoa_results, num = 10, color = "orange"){
  results <- pcoa_results$values
  results$n <- as.numeric(rownames(results))
  p <- ggplot(results %>% head(n = num), mapping = aes(x = n, y = Relative_eig)) + geom_col(fill = color) +
    theme_minimal() + xlab("Dimensions") + ylab("Percentage of variance explained") + labs(title = "Scree plot") +
    scale_x_continuous(breaks = results$n) + # Ensure a tick for each x-axis value
    scale_y_continuous(labels = percent_format()) 
  return(p)
}


plot_pcoa <- function(pcoa_table, pcoa_results, colorvar, shapevar, palette, ellipse = F){
  # Plot from pcoa results
  p <- ggplot(pcoa_table, mapping = aes(x = Axis.1, y = Axis.2, color = get(colorvar), shape = get(shapevar))) +
    geom_point(alpha = 0.9) +
    xlab(paste0('Axis 1 (', round(pcoa_results$values$Relative_eig[1] * 100, 2), '%)')) +
    ylab(paste0('Axis 2 (', round(pcoa_results$values$Relative_eig[2] * 100, 2), '%)')) +
    theme_minimal() +
    scale_color_manual(values = palette) +
    guides(color = guide_legend(title = colorvar), shape = guide_legend(title = shapevar), )
  if (ellipse){
    p <- p + stat_ellipse(data=pcoa_table, mapping = aes(color=get(colorvar), fill=get(colorvar), x=Axis.1, y=Axis.2, shape = NULL),
                          type="norm", level = 0.9, geom = "polygon", alpha=0.04) +
      guides(fill = guide_legend(title = colorvar))
  }
  return(p)
}

plot_umap <- function(phy_rpkm, n_neighbors, colorvar, shapevar, palette, nn_method = "nndescent", metric = 'braycurtis', 
                      n_components = 2, min_dist = 0.5, ellipse = F){
  umap_matrix <- otu_table(phy_rpkm) %>% data.frame() %>% as.matrix.data.frame() %>% t()
  umap_result <- uwot::umap(umap_matrix, nn_method = nn_method, metric = metric, n_neighbors = n_neighbors,
                            n_components = n_components, min_dist = min_dist)
  umap_plot_table <- cbind(umap_result, sample_data(phy_rpkm)[rownames(umap_result)])
  p <- ggplot(umap_plot_table, mapping = aes(x = `1`, y = `2`, color = get(colorvar), shape = get(shapevar))) +
    geom_point(alpha = 0.9) + theme_minimal() +
    xlab('UMAP1') + ylab('UMAP2') + labs(title = paste0("UMAP (",nrow(tax_table(phy_rpkm)) ," plasmids)")) +
    guides(color = guide_legend(title = colorvar), shape = guide_legend(title = shapevar)) + 
    scale_color_manual(values = palette)
  if (ellipse){
    p <- p + stat_ellipse(data=umap_plot_table, mapping = aes(color=get(colorvar), fill=get(colorvar), x = `1`, y = `2`, shape = NULL),
                          type="norm", level = 0.9, geom = "polygon", alpha=0.04) +
      guides(fill = guide_legend(title = colorvar))
  }
  return(p)
}


corrplot_from_physeq_complexHeatmap_arnau <- function(phy_obj, interest_metadata, is_taxa = T, corrmethod = "spearman", qval = 0.05, 
                                                      adjust_method = "BH", legend = T, symbol_fontsize = 10, row_fontsize = 12, 
                                                      column_fontsize = 12, AncomBC_results = NULL, colnames_table = NULL,
                                                      col_fun = colorRamp2(c(-0.3, 0, 0.5), c("darkgoldenrod2", "white", "blue3")),
                                                      row_label = "Plasmid", cols_label = "Personal data", n_cols = F, n_rows = F
                                                      ) { #interest_metadata is a vector with the colnamaes of the metadata variables to correlate
  
  current_metadata <- data.frame(sample_data(phy_obj))
  current_metadata <- current_metadata[,interest_metadata]
  current_metadata <- as.matrix(current_metadata)
  current_species <- as.matrix(data.frame(otu_table(phy_obj)))
  current_taxtab <- tax_table(phy_obj) %>% data.frame
  
  a <- rcorr(current_metadata, t(current_species), type = corrmethod)
  a$r[lower.tri(a$r)] <- NA
  a$P[lower.tri(a$P)] <- NA
  a$n[lower.tri(a$n)] <- NA
  
  melted_pvals <- reshape2::melt(a$P)
  melted_pvals <- melted_pvals[! is.na(melted_pvals$value),]
  melted_pvals <- melted_pvals[melted_pvals$Var1 %in% interest_metadata,]
  melted_pvals <- melted_pvals[!(melted_pvals$Var2 %in% interest_metadata),]
  melted_pvals$FDR <- p.adjust(melted_pvals$value, method = adjust_method)
  
  melted_Rsquared <- reshape2::melt(a$r)
  melted_Rsquared <- melted_Rsquared[! is.na(melted_Rsquared$value),]
  melted_Rsquared <- melted_Rsquared[melted_Rsquared$Var1 %in% interest_metadata,]
  melted_Rsquared <- melted_Rsquared[!(melted_Rsquared$Var2 %in% interest_metadata),]
  melted_pvals$R <- melted_Rsquared$value
  #create a column with R coefficient only if correlation is significant
  melted_pvals$r_if_sig <- melted_pvals$R
  melted_pvals$r_if_sig[melted_pvals$FDR > qval] <- NA
  
  
  rm(melted_Rsquared)
  
  #melted_pvals$Var2  <- current_taxtab[melted_pvals$Var2 %>% as.character(), 7]
  # corr_food_species_table <- corr_food_species_table[abs(corr_food_species_table$r_if_sig) >= 0.4, ]
  
  
  all_corrs <- melted_pvals
  significant_species <-  all_corrs$Var2[all_corrs$FDR <= qval] %>%  unique
  significant_metadata <-  all_corrs$Var1[all_corrs$FDR <= qval] %>%  unique
  
  all_corrs_dcast <- all_corrs[all_corrs$Var2 %in% significant_species,]
  all_corrs_dcast <- all_corrs_dcast[all_corrs_dcast$Var1 %in% significant_metadata,]
  # all_corrs_dcast <- all_corrs[abs(all_corrs$R) >= 0.4, ]
  
  #create the values matrix
  all_corrs_dcast$value <- all_corrs_dcast$r_if_sig
  value_mat <- dcast(all_corrs_dcast, Var2 ~ Var1) %>% column_to_rownames("Var2")
  value_mat <- ifelse(value_mat < 0, "-", "+")
  value_mat[is.na(value_mat)] <- ""
  
  
  all_corrs_dcast$value <- all_corrs_dcast$R
  
  all_corrs_dcast <- dcast(all_corrs_dcast, Var2 ~ Var1) %>% column_to_rownames("Var2")
  
  # all_corrs_dcast[is.na(all_corrs_dcast)] <- 0
  
  
  
  
  if (!is.null(AncomBC_results)) {
    rowcolors_notation <- data.frame(Plasmid_ID = rownames(all_corrs_dcast))
    rowcolors_notation <- left_join(rowcolors_notation, AncomBC_results, by = "Plasmid_ID") %>% select(Plasmid_ID, label) %>% rename("Abundance" = label)
    rownames(rowcolors_notation) <- rowcolors_notation$Plasmid_ID
    rowcolors <- rowcolors_notation[match(rownames(all_corrs_dcast), rowcolors_notation$Plasmid_ID), ] %>% select(Abundance)
    row_anno <- HeatmapAnnotation(df = rowcolors, col = list(Abundance = c("CD+" = "#ED6925FF","CD-" = "#414487FF", "UC+" = "#CF4446FF", "UC-" = "#FDE725FF",
                                                                           "CD-,UC-" = "#4B0C6BFF", "CD+,UC+" = "#B4DE2CFF")), 
                                  na_col = "white", which = "row", name = "Abundance"
    )
  } else {
    row_anno <- c()
    rowcolors_notation <- NA
  }
  
  if (!is_null(colnames_table)){
    colnames(all_corrs_dcast) <- colnames_table[colnames(all_corrs_dcast), "names"]
  }
  
  if (n_cols){
    cols_label <- paste0(ncol(all_corrs_dcast), " ", cols_label)
  }
  if (n_rows){
    row_label <- paste0(nrow(all_corrs_dcast), " ", row_label)
  }
  
  plot_result <- Heatmap(all_corrs_dcast,
                         name = "Spearman's Rho", 
                         col = col_fun, na_col = "white", 
                         rect_gp = gpar(col = "grey95", lwd = 1), 
                         show_heatmap_legend = legend,
                         column_names_gp = grid::gpar(fontsize = column_fontsize),
                         row_names_gp = grid::gpar(fontsize = row_fontsize),
                         column_title = cols_label,  
                         row_title = row_label,
                         right_annotation = row_anno,
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           if(value_mat[i,j] == "+"){
                             grid.text(sprintf("%s", value_mat[i, j]), x, y, gp = gpar(fontsize = symbol_fontsize, col = "black", fontface = 2))
                           } else if(value_mat[i,j] == "-") {
                             grid.text(sprintf("%s", value_mat[i, j]), x, y, gp = gpar(fontsize = symbol_fontsize, col = "white", fontface = 2))
                           }
                         }
  )
  
  return(list(table = melted_pvals, plot = plot_result, ANCOMBC_match = rowcolors_notation))
}


ancombc_heatmap <- function(physeq, ancombc_results, disease, color_ramp = NULL, plot_title = NULL, control = "Healthy", fdr_threshold = 0.05, 
                            fold_change_threshold = log(2), tax_level = "Plasmid_ID", addtitle = "",
                            legendname = NA, col_fontsize = 7, row_fontisize = 0, struc_zeros = NULL, only_struc_zeroes = F, n_struc_zeroes = NULL, 
                            whatisdifferential = "Plasmids", rowside_vars = NULL, rowside_cols = NULL, numeric_rowsidevars = NULL, 
                            numeric_rowside_cols = NULL, 
                            ...){
  
  temp_taxa_table <- tax_table(physeq) %>% data.frame()
  temp_metadata <- sample_data(physeq) %>% data.frame()
  
  ancombc_results <- ancombc_results  %>% filter(Disease == disease)
  
  
  if (!only_struc_zeroes){
    #relative abundance
    df <- data.frame(otu_table(physeq)[ancombc_results$taxon[ancombc_results$`q` <= fdr_threshold & 
                                                               abs(ancombc_results$lfc) >= fold_change_threshold & 
                                                               ancombc_results$passed_ss == T],])
    if (!is.null(struc_zeros)){
      struc_zeros <- struc_zeros %>% rename(sz_disease = paste0("sz_", disease))
      rownames(struc_zeros) <- struc_zeros$taxon
      if (is.null(n_struc_zeroes)){
        # Get the otu information for the structural zero plasmids 
        zeros_df <- data.frame(otu_table(physeq)[struc_zeros %>% filter(sz_healthy != sz_disease) %>% pull(taxon), colnames(df)])
      } else {
        # get a table witht the median of each structural zero plasmid in healthy and disease samples
        pre_zeros_df <- struc_zeros %>% filter(sz_healthy != sz_disease)
        pre_zeros_df$median_healthy <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == control, "SubjectID"]] %>%
          apply(MARGIN = 1, FUN = median)
        pre_zeros_df$median_disease <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == disease, "SubjectID"]] %>%
          apply(MARGIN = 1, FUN = median)
        pre_zeros_df$mean_healthy <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == control, "SubjectID"]] %>%
          apply(MARGIN = 1, FUN = mean)
        pre_zeros_df$mean_disease <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == disease, "SubjectID"]] %>%
          apply(MARGIN = 1, FUN = mean)
        pre_zeros_df <- pre_zeros_df %>% mutate(mean_diff = abs(mean_healthy - mean_disease), 
                                                median_diff = abs(median_healthy - median_disease)) %>% 
          arrange(desc(median_diff), desc(mean_diff)) %>% 
          head(n = n_struc_zeroes)
        zeros_df <- data.frame(otu_table(physeq)[pre_zeros_df$taxon, colnames(df)])
      }
      df <- rbind(df, zeros_df)
    } else {
      zeros_df <- NULL
    }
  } else {
    struc_zeros <- struc_zeros %>% rename(sz_disease = paste0("sz_", disease))
    rownames(struc_zeros) <- struc_zeros$taxon
    if (is.null(n_struc_zeroes)){
      # Get the otu information for the structural zero plasmids 
      zeros_df <- data.frame(otu_table(physeq)[struc_zeros %>% filter(sz_healthy != sz_disease) %>% pull(taxon), ])
    } else {
      # get a table witht the median of each structural zero plasmid in healthy and disease samples
      pre_zeros_df <- struc_zeros %>% filter(sz_healthy != sz_disease)
      pre_zeros_df$median_healthy <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == control, "SubjectID"]] %>%
        apply(MARGIN = 1, FUN = median)
      pre_zeros_df$median_disease <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == disease, "SubjectID"]] %>%
        apply(MARGIN = 1, FUN = median)
      pre_zeros_df$mean_healthy <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == control, "SubjectID"]] %>%
        apply(MARGIN = 1, FUN = mean)
      pre_zeros_df$mean_disease <- otu_table(physeq)[pre_zeros_df$taxon, temp_metadata[temp_metadata$Disease == disease, "SubjectID"]] %>%
        apply(MARGIN = 1, FUN = mean)
      pre_zeros_df <- pre_zeros_df %>% mutate(mean_diff = abs(mean_healthy - mean_disease), 
                                              median_diff = abs(median_healthy - median_disease)) %>% 
        arrange(desc(median_diff), desc(mean_diff)) %>% 
        head(n = n_struc_zeroes)
      zeros_df <- data.frame(otu_table(physeq)[pre_zeros_df$taxon, ])
    }
    df <- zeros_df
  }
  
  rownames(df) <- temp_taxa_table[rownames(df), tax_level]  #we use the taxonomy_SGB because with tax_glom it is possible to take different taxonID for same species name
  
  # df <- filter_prevalence_abundance_dataframe_no_loop(df, abundance = 0.001, prevalence = 0.1)[1]  %>% as.data.frame()#filter low abundance, 0.1% of abundance in 10% of samples
  df <- data.frame(t(df))
  # df <- df[,colSums(df) >= 5] #rm low abundance samples
  
  #create the logFC dataframe
  lfolds <- data.frame(taxon = ancombc_results$taxon, logFC = ancombc_results$lfc)
  if (!is_null(zeros_df)){
    zeros_lfolds <- data.frame(taxon = colnames(df)[!colnames(df) %in% lfolds$taxon], 
                               sz_healthy = struc_zeros[colnames(df)[!colnames(df) %in% lfolds$taxon], "sz_healthy"], 
                               sz_disease = struc_zeros[colnames(df)[!colnames(df) %in% lfolds$taxon], "sz_disease"]) %>% 
      mutate(logFC = case_when(sz_disease == T ~ min(lfolds$logFC) - 0.5, sz_disease == F ~ max(lfolds$logFC) + 0.5)) %>%
      # assign max value + sth for the logFC value of structural zeroes
      select(taxon, logFC)
    lfolds <- rbind(lfolds, zeros_lfolds)
  }
  rownames(lfolds) <- lfolds$taxon
  lfolds <- lfolds[colnames(df),]
  lfolds <- lfolds[order(lfolds$logFC),]
  df <- df[,rownames(lfolds)]
  
  #df <- cbind(Enriched = temp_metadata[rownames(df),]$Enriched, group = temp_metadata[rownames(df),]$group, df)
  #df <- df[df$group != "CTRL",]
  df <- cbind(Disease = temp_metadata[rownames(df),]$Disease, df)
  
  
  df <- df[order(df$Disease),]
  # heatmap(df[,c(-1, -2)] %>% as.matrix(), RowSideColors = df$disease, Rowv = NA)
  
  # #col <- list(Enriched = c("Enriched" = "tomato2", "Nonenriched" = "dodgerblue2"), group = c("CTRL" = "blue", "INT" = "red"))
  # col <- list(Disease = colors[["disease"]])
  # 
  # 
  # ha <- HeatmapAnnotation(
  #   Disease = df$Disease,
  #   col = col,
  #   which = "row"
  # )
  # 
  colmat <- lfolds$logFC
  # abs_max <- quantile(abs(colmat - 0.5), 0.95, na.rm = TRUE)
  col_fun_logFC = colorRamp2(c(min(colmat), 0, abs(0.5 + max(colmat))), c("blue4", "white", "red"))
  
  logFC_annotation <- HeatmapAnnotation(
    logFC = lfolds$logFC,
    col = list( logFC = col_fun_logFC)
  )
  
  # Prepare rowside variables to plot
  n_colors <- 0
  r_ha <- c()
  if (! is.null(rowside_vars)) {
    # Get how many rowside variables we have
    n_row_vars <- length(rowside_vars)
    for (n in 1:n_row_vars) {
      # For each, get how many colors there are
      values <- sample_data(physeq) %>% data.frame() %>% pull(rowside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(rowside_vars)[n]
      if (is_null(var_name) || var_name == "") {
        var_name <- rowside_vars[n]
      }
      # Prepare palette
      if (is.null(rowside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- rowside_cols[[var_name]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      r_ha <- append(r_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        sample_data(physeq)[rownames(df), rowside_vars[n]] %>% data.frame() %>% rename(setNames(object = rowside_vars[n], nm = var_name))
        %>% replace_na(replace = setNames(object = list("NA"), nm = c(rowside_vars[n]))), 
        list(col = cols, which = "row")
      )))
    }
  }
  
  # FOR NUMERIC VARIABLES
  
  if (!is_null(numeric_rowsidevars)) {
    # Get how many rowside variables we have
    n_row_vars <- length(numeric_rowsidevars)
    for (n in 1:n_row_vars) {
      var_name <- names(numeric_rowsidevars)[n]
      if (is_null(var_name) || var_name == "") {
        var_name <- numeric_rowsidevars[n]
      }
      # Prepare palette
      if (is.null(numeric_rowside_cols)) {
        col_fun = colorRamp2(c(min(sample_data(physeq)[, numeric_rowsidevars[n]], na.rm = T), 
                               max(sample_data(physeq)[, numeric_rowsidevars[n]], na.rm = T)), 
                                randomColor(count = 2))
      } else {
        col_fun <- numeric_rowside_cols[[var_name]]
      }
      cols <- list()
      cols[[var_name]] <- col_fun
      r_ha <- append(r_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        sample_data(physeq)[rownames(df), numeric_rowsidevars[n]] %>% data.frame() %>% rename(setNames(object = numeric_rowsidevars[n], 
                                                                                                       nm = var_name)),
        #          %>% replace_na(replace = setNames(object = list(0), nm = c(numeric_rowsidevars[n]))), 
        list(col = cols, which = "row")
      )))
    }
  }
  
  if (is_null(plot_title)){
    plot_title <- paste0(addtitle, disease, " vs. ", control, " (FC >= ", exp(fold_change_threshold) , "): ", ncol(df) - 1, 
                         " Differentially Abundant ", whatisdifferential)
  }
  
  hm_plot <- Heatmap(as.matrix(df[,c(-1)]),
                     cluster_rows = T,
                     cluster_columns = F,
                     split = df$Disease,
                     row_names_gp = gpar(fontsize = row_fontisize),
                     column_names_gp = gpar(fontsize = col_fontsize),
                     clustering_distance_rows = "manhattan",
                     clustering_method_rows = "average",
                     right_annotation = r_ha,
                     top_annotation = logFC_annotation,
                     name = legendname, 
                     column_title = plot_title, 
                     column_title_side = "top",
                     column_title_gp = gpar(fontsize = 16),
                     use_raster = F, 
                     #heatmap_legend_param = list(title = "name"),
                     col = color_ramp, 
                     ...
  )
  
  return(list(plot = hm_plot, FC = lfolds, profiling = df))  
}


corrplot_from_physeq_metadata_complexHeatmap_arnau <- function(phy_obj, interest_metadata1, interest_metadata2, corrmethod = "spearman", 
                                                               qval = 0.05, adjust_method = "BH", legend = T, symbol_fontsize = 10, 
                                                               row_fontsize = 12, column_fontsize = 12, keep1 = NA,
                                                               colnames_table = NULL, 
                                                               col_fun = colorRamp2(c(-0.3, 0, 0.5), c("darkgoldenrod2", "white", "blue3"))
                                                               ) { #interest_metadata is a vector with the colnamaes of the metadata variables to correlate
  library(Hmisc)
  library(corrplot)
  current_metadata <- data.frame(sample_data(phy_obj))
  current_metadata <- current_metadata[,append(interest_metadata1, interest_metadata2)]
  current_metadata <- as.matrix(current_metadata)
  
  a <- rcorr(current_metadata, type = corrmethod)
  a$r[lower.tri(a$r)] <- NA
  a$P[lower.tri(a$P)] <- NA
  a$n[lower.tri(a$n)] <- NA
  
  melted_pvals <- reshape2::melt(a$P)
  melted_pvals <- melted_pvals[! is.na(melted_pvals$value),]
  melted_pvals <- melted_pvals[melted_pvals$Var1 %in% interest_metadata1,]
  melted_pvals <- melted_pvals[!(melted_pvals$Var2 %in% interest_metadata1),]
  melted_pvals$FDR <- p.adjust(melted_pvals$value, method = adjust_method)
  
  melted_Rsquared <- reshape2::melt(a$r)
  melted_Rsquared <- melted_Rsquared[! is.na(melted_Rsquared$value),]
  melted_Rsquared <- melted_Rsquared[melted_Rsquared$Var1 %in% interest_metadata1,]
  melted_Rsquared <- melted_Rsquared[!(melted_Rsquared$Var2 %in% interest_metadata1),]
  melted_pvals$R <- melted_Rsquared$value
  #create a column with R coefficient only if correlation is significant
  melted_pvals$r_if_sig <- melted_pvals$R
  melted_pvals$r_if_sig[melted_pvals$FDR > qval] <- NA
  
  rm(melted_Rsquared)
  
  to_plot <- melted_pvals
  
  # if (! plotnonsig) {
  #   #Subset only species that have at least one significant correlation in variable1
  #   significantvar1 <- unique(melted_pvals$Var1[melted_pvals$FDR <= qval]) %>% as.vector
  #   to_plot <- to_plot[to_plot$Var1 %in% significantvar1,]
  #   
  #   #Subset only species that have at least one significant correlation in variable2
  #   significantvar2 <- unique(melted_pvals$Var2[melted_pvals$FDR <= qval]) %>% as.vector
  #   to_plot <- to_plot[to_plot$Var2 %in% significantvar2,]
  # }
  
  all_corrs <- to_plot
  significant_species <-  all_corrs$Var2[all_corrs$FDR <= qval] %>%  unique
  significant_metadata <-  all_corrs$Var1[all_corrs$FDR <= qval] %>%  unique
  if (!any(is.na(keep1))){
    significant_metadata <- factor(unique(c(as.character(significant_metadata), keep1)))
  }
  
  all_corrs_dcast <- all_corrs[all_corrs$Var2 %in% significant_species,]
  all_corrs_dcast <- all_corrs_dcast[all_corrs_dcast$Var1 %in% significant_metadata,]
  
  #create the values matrix
  all_corrs_dcast$value <- all_corrs_dcast$r_if_sig
  value_mat <- dcast(all_corrs_dcast, Var2 ~ Var1) %>% column_to_rownames("Var2")
  value_mat <- ifelse(value_mat < 0, "-", "+")
  value_mat[is.na(value_mat)] <- ""
  
  
  all_corrs_dcast$value <- all_corrs_dcast$R
  
  all_corrs_dcast <- dcast(all_corrs_dcast, Var2 ~ Var1) %>% column_to_rownames("Var2")
  
  # all_corrs_dcast[is.na(all_corrs_dcast)] <- 0
  
  
  library(circlize)
  
  
  if (!is_null(colnames_table)){
    rownames(all_corrs_dcast) <- colnames_table[rownames(all_corrs_dcast), "names"]
  }
  
  plot_result <- Heatmap(all_corrs_dcast,
                         name = "Spearman's Rho", 
                         col = col_fun, na_col = "white", 
                         rect_gp = gpar(col = "grey95", lwd = 1), 
                         show_heatmap_legend = legend,
                         column_names_gp = grid::gpar(fontsize = column_fontsize),
                         row_names_gp = grid::gpar(fontsize = row_fontsize),
                         # column_title = "Species",  
                         # row_title = "Personal data",
                         cell_fun = function(j, i, x, y, width, height, fill) {
                           if(value_mat[i,j] == "+"){
                             grid.text(sprintf("%s", value_mat[i, j]), x, y, gp = gpar(fontsize = symbol_fontsize, col = "black", fontface = 2))
                           } else if(value_mat[i,j] == "-") {
                             grid.text(sprintf("%s", value_mat[i, j]), x, y, gp = gpar(fontsize = symbol_fontsize, col = "white", fontface = 2))
                           }
                         }
  )
  return(list(table = melted_pvals, plot = plot_result))
}




plot_rpkm_heatmap_longitudinal <- function(phy, cluster_rows = F, colside_vars = NULL, colside_cols = NULL, rowside_vars = NULL, rowside_cols = NULL,
                              col_title = NA, row_title = NA, axis_name = "log RPKM\nabundance", rsize = 0, csize = 0,
                              cluster_columns = T, row_dend_width = 0.1, clustering_distance_columns = "manhattan",
                              clustering_method_columns = "complete", ...){
  # Function to plot a heatmap from a phyloseq object for plasmids
  plot_matrix <- as.matrix.data.frame(data.frame(otu_table(phy))) %>% t()
  group_annotation <- phy@sam_data$Subject
  # Prepare rowside variables to plot
  n_colors <- 0
  r_ha <- c()
  c_ha <- c()
  if (! is.null(rowside_vars)) {
    # Get how many rowside variables we have
    n_row_vars <- length(rowside_vars)
    for (n in 1:n_row_vars) {
      # For each, get how many colors there are
      values <- sample_data(phy) %>% data.frame() %>% pull(rowside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(rowside_vars)[n]
      if (is_null(var_name) || var_name == "") {
        var_name <- rowside_vars[n]
      }
      # Prepare palette
      if (is.null(rowside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- rowside_cols[[var_name]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      r_ha <- append(r_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        sample_data(phy)[, rowside_vars[n]] %>% data.frame() %>% rename(setNames(object = rowside_vars[n], nm = var_name))
        %>% replace_na(replace = setNames(object = list("NA"), nm = c(rowside_vars[n]))), 
        list(col = cols, which = "row")
      )))
    }
  }
  # Do the same for the colside variables 
  if (! is.null(colside_vars)) {
    # Get how many colside variables we have
    n_col_vars <- length(colside_vars)
    for (n in 1:n_col_vars) {
      # For each, get how many colors there are
      values <- tax_table(phy) %>% data.frame() %>% pull(colside_vars[n]) %>% unique()
      # Raplace the NAs for "NA" if there are any 
      values <- replace(x = values, list = is.na(values), values = "NA")
      diff_values <- length(values)
      # Get the names from the vector
      var_name <- names(colside_vars)[n]
      if (var_name == "" || is.null(var_name)) {
        var_name <- colside_vars[n]
      }
      # Prepare palette
      if (is.null(colside_cols)) {
        new_n_colors <- diff_values + n_colors
        colors <- distinctColorPalette(new_n_colors)
        palette <- colors[(n_colors+1):new_n_colors]
        n_colors <- new_n_colors
      } else {
        palette <- colside_cols[[var_name]][1:diff_values]
      }
      if (is.null(names(palette))) {
        names(palette) <- values
      }
      cols <- list()
      cols[[var_name]] <- palette
      # prepare heatmap annotation
      c_ha <- append(c_ha, do.call(HeatmapAnnotation, c(
        # Also replace NA values by "NA"
        tax_table(phy)[, colside_vars[n]] %>% as.data.frame() %>% rename(setNames(object = colside_vars[n], nm = var_name))
        %>% replace_na(replace = setNames(object = list("NA"), nm = c(colside_vars[n]))),
        list(col = cols, which = "col")
      )))
    }
  }
  # Plot the heatmap 
  Heatmap(plot_matrix,
          # cluster rows will either be F or the dendrogram
          cluster_rows = cluster_rows,
          cluster_columns = cluster_columns,
          row_names_gp = gpar(fontsize = rsize),
          column_names_gp = gpar(fontsize = csize),
          #clustering_distance_rows = braycurtis,
          #clustering_method_rows = "average",
          clustering_distance_columns = clustering_distance_columns,
          clustering_method_columns = clustering_method_columns,
          right_annotation = r_ha, 
          top_annotation = c_ha,
          name = axis_name, column_title = col_title,
          row_title = row_title, 
          use_raster = F,
          split = group_annotation,
          row_dend_width = unit(row_dend_width, units = "npc"),
          col = colorRamp2(c(-10, -6, -3, 0, 2, 5, 7), c("black", viridis_pal(begin = 0, end = 1, direction = 1)(6))),
          #col = colorRamp2(c(0, 1000, 2000, 3000, 4000, 5000, 8000), c("black", "blue", "cyan", "green", "yellow","orange", "red")),
          ...
  )
}

plot_beta_div_plasmids <- function(physeq, dist, group, show_centroids = TRUE, ellipse = TRUE,
                          theme = c('pubr', 'prism', 'classic'), title = NULL, pairwise_adonis = FALSE, colors = NULL,
                          add_pval = FALSE){
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PLOTS BETA DIVERSITY METRICS
  #
  # Parameters
  # ----------
  #   physeq: phyloseq object
  #   dist: distance object (distances between samples)
  #   group: grouping variable
  #   show_centroids: if True, plot centroids
  #   ellipse: if True, show group ellipses
  #   theme: ggplot theme
  #   title: plot title
  #
  # Returns: ggplot plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  df_comp <- NULL
  
  
  metadata <- sample_data(physeq)
  theme = match.arg(theme)
  
  
  
  # filter distance and metadata to keep only matching sample ids
  
  ids <- intersect(rownames(metadata), rownames(as.matrix(dist)))
  
  metadata <- metadata[ids, ]
  dist <- dist(as.matrix(dist)[ids, ids])
  
  adonis_res <- adonis2(dist ~ metadata[[group]],
                        na.action = na.omit,
                        by = "terms",
                        permutations = 1999
  )
  
  
  # Compute PCOA
  
  mod <- betadisper(dist, metadata[[group]])
  d <- permutest(mod)
  
  
  centroids <- data.frame(group=rownames(mod$centroids), data.frame(mod$centroids))
  vectors <- data.frame(group = mod$group, data.frame(mod$vectors))
  
  # class(vectors$group) <- class(metadata[[group]])
  
  # Plot PCoA
  p <- ggplot() +
    geom_point(data = vectors, aes(color = group, x=PCoA1, y= PCoA2)) +
    guides(color = guide_legend(title = group))
  
  # theme
  p <- switch(
    theme,
    "prism"= p + theme_prism(),
    "pubr"= p + theme_pubr(),
    "classic"= p + theme_classic()
  )
  
  if (show_centroids){
    p <- p + geom_point(data = centroids, aes(color = group, x=PCoA1, y=PCoA2, size = 5)) +
      guides(size = FALSE)
  }
  
  if (ellipse){
    p <- p + stat_ellipse(data=vectors, mapping = aes(color=group, fill=group, x=PCoA1, y=PCoA2),
                          type="norm", level = 0.9, geom = "polygon", alpha=0.04) +
      guides(fill = guide_legend(title = group))
  }
  
  if (!is.null(title)){
    p <- p + ggtitle(title)
  }
  
  if (!is.null(colors)){
    p <- p + scale_color_manual(values = colors)
  }
  
  
  
  if (pairwise_adonis){
    # all comparisons
    comp <- combn(levels(as.factor(metadata[[group]])), 2)
    
    df_comp <- setNames(data.frame(matrix(ncol = 4, nrow=0)),
                        c('group1', 'group2', 'adonis_pval', 'permutest_pval'))
    
    for (i in 1:ncol(comp)){
      
      pair <- comp[, i]
      
      # physeq_filt <- subset_samples(physeq, get(group) %in% pair)
      
      # m <- sample_data(physeq_filt)
      
      m <- metadata[metadata[[group]] %in% pair, ]
      
      # filter distance and metadata to keep only matching sample ids
      
      ids <- intersect(rownames(m), rownames(as.matrix(dist)))
      
      m <- m[ids, ]
      dist_filt <- dist(as.matrix(dist)[ids, ids])
      
      adonis_pair_res <- adonis2(dist_filt ~ m[[group]],
                                 na.action = na.omit,
                                 by = "terms",
                                 permutations = 1999
      )
      pval <- adonis_pair_res$`Pr(>F)`[1]
      
      mod_pair <- betadisper(dist_filt, m[[group]])
      d_pair <- permutest(mod_pair)
      pempval <- d_pair$tab$`Pr(>F)`[1]
      
      df_comp[nrow(df_comp) + 1,] <- c(pair, pval, pempval)
      
    }
    
    df_comp$adonis_p.adj <- p.adjust(df_comp$adonis_pval, method = 'BH')
    
    df_comp <- df_comp %>%
      mutate(adonis_pval = as.numeric(adonis_pval)) %>%
      mutate(adonis_p.signif = case_when(
        adonis_pval <=0.001 ~ '***',
        adonis_pval >0.001 & adonis_pval <=0.01 ~ '**',
        adonis_pval > 0.01 & adonis_pval <=0.05 ~ '*',
        adonis_pval > 0.05 & adonis_pval <=0.1 ~ '.',
        adonis_pval > 0.1 ~ 'ns'
      )) %>%
      mutate(permutest_pval = as.numeric(permutest_pval)) %>%
      mutate(permutest_p.signif = case_when(
        permutest_pval <=0.001 ~ '***',
        permutest_pval >0.001 & permutest_pval <=0.01 ~ '**',
        permutest_pval > 0.01 & permutest_pval <=0.05 ~ '*',
        permutest_pval > 0.05 & permutest_pval <=0.1 ~ '.',
        permutest_pval > 0.1 ~ 'ns'
      ))  %>%
      mutate(adonis_p.adj = as.numeric(adonis_p.adj)) %>%
      mutate(adonis_p.adj.signif = case_when(
        adonis_p.adj <=0.001 ~ '***',
        adonis_p.adj >0.001 & adonis_p.adj <=0.01 ~ '**',
        adonis_p.adj > 0.01 & adonis_p.adj <=0.05 ~ '*',
        adonis_p.adj > 0.05 & adonis_p.adj <=0.1 ~ '.',
        adonis_p.adj > 0.1 ~ 'ns'
      ))
    
  }
  
  
  if (add_pval){
    
    pval_text <- paste0('PERMANOVA, R2=', round(adonis_res$R2[1], 2), ', p = ', round(adonis_res$`Pr(>F)`[1], 3))
    
    p <- p + annotate(geom = 'text',
                      label = pval_text,
                      x = Inf, y = -Inf, hjust = 1, vjust = -1)
    
    if (pairwise_adonis){
      p <- plot_grid(p,
                     plot_grid(
                       tableGrob(get_pairwise_matrix(df_comp, values_from = 'adonis_p.adj'), rows = NULL),
                       tableGrob(get_pairwise_matrix(df_comp, values_from = 'adonis_p.adj.signif'), rows = NULL),
                       # tableGrob(get_pairwise_matrix(df_comp, values_from = 'permutest_p.signif'), rows = NULL),
                       nrow=2,
                       labels = c('adonis - adjusted p-values', 'adonis - significance'), hjust = c(0,0), label_fontface = "plain"),
                     
                     ncol = 2, rel_widths = c(2, 1)
                     
      )
    }
    
  }
  
  
  return(list('plot'=p, 'adonis' = adonis_res, 'permutest' = d, 'pairwise' = df_comp))
  
  
  
}

