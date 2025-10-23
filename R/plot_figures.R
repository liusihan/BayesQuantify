#' Visualize the distribution of variants
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param classification_col The column name for variant classification (str)
#' @param consequence_col The column name for the annotation results of variant consequences(str)
#' @param gene_col The column name for the gene where the variant is located(str)
#'
#'
#' @import stringr
#' @import reshape2
#' @import ggpie
#' @importFrom gridExtra tableGrob
#' @import patchwork
#' @import ggplot2
#' @importFrom dplyr `%>%` mutate group_by summarize full_join n arrange desc
#' @importFrom utils globalVariables
#'
#'
#' @return Figures
#' @export
#'
#' @examples
#' data("ClinGen_dataset")
#' ClinGen_dataset <- add_info(ClinGen_dataset, "Assertion")
#' ClinGen_dataset <- VUS_classify(ClinGen_dataset, "Assertion", "Applied Evidence Codes (Met)")
#' multi_plot(ClinGen_dataset, "Assertion", "HGNC Gene Symbol")
#'
multi_plot <- function(data, classification_col, gene_col, consequence_col = NULL) {
  classification_num <- which(colnames(data) == classification_col)
  gene_num <- which(colnames(data) == gene_col)
  consequence_num <- which(colnames(data) == consequence_col)
  if ("Pathogenic" %in% names(table(data[, classification_num]))) {
    data[, classification_num] <- factor(data[, classification_num], levels = c("Pathogenic", "Likely Pathogenic", "Uncertain Significance", "Likely Benign", "Benign"))
  } else {
    data[, classification_num] <- factor(data[, classification_num], levels = c("P", "LP", "VUS", "LB", "B"))
  }
  f1 <- ggdonut(
    data = data, group_key = colnames(data)[classification_num], count_type = "full",
    label_info = "all", label_type = "horizon",
    label_size = 6, label_pos = "out"
  ) + scale_fill_manual(values = c("#780000", "#c1121f", "#bfdbf7", "#cad2c5", "#84a98c")) + labs(fill = "Variant classification")
  VUS <- data.frame(table(data[data$VUS_class != "", ]$VUS_class))
  colnames(VUS) <- c("VUS_temperature_scale", "Number of variants")
  VUS_data <- data.frame(matrix(, nrow = 6, ncol = 0))
  VUS_data$Posterior_probability <- c("[81.2%-90)", "[67.5%-81.2%)", "[50%-67.5%)", "[32.5%-50%)", "[18.8%-32.5%)", "[10%-18.8%)")
  VUS_data$VUS_temperature_scale <- c("Hot", "Warm", "Tepid", "Cool", "Cold", "IceCold")
  VUS_data <- full_join(VUS_data, VUS)
  rownames(VUS_data) <- VUS_data$VUS_temperature_scale
  f2 <- gridExtra::tableGrob(VUS_data[, c(1, 3)])
  count_data <- data.frame(data %>%
    group_by(get(colnames(data)[classification_num]), get(colnames(data)[gene_num])) %>%
    summarize(n = n()))
  colnames(count_data) <- c("Classification", "Gene", "n")
  count_data <- arrange(count_data, desc(n), Gene)
  count_data <- count_data[count_data$Classification == "P" | count_data$Classification == "LP" | count_data$Classification == "Pathogenic" | count_data$Classification == "Likely Pathogenic", ]
  gene <- data.frame(count_data %>%
    group_by(Gene) %>%
    summarize(sum(n)))
  gene <- arrange(gene, desc(sum.n.))
  count_data$Gene <- factor(count_data$Gene, levels = as.character(unique(gene$Gene)))
  if ("Pathogenic" %in% names(table(count_data$Classification))) {
    count_data$Classification <- factor(count_data$Classification, levels = c("Pathogenic", "Likely Pathogenic"))
  } else {
    count_data$Classification <- factor(count_data$Classification, levels = c("P", "LP"))
  }
  f3 <- ggplot(data = count_data[count_data$Gene %in% gene$Gene[1:20], ], mapping = aes(x = Gene, y = n, fill = Classification)) +
    geom_bar(stat = "identity", colour = "black", width = 1, lwd = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "italic")) +
    labs(fill = "", x = "", y = "Number of P/LP variants") +
    scale_fill_manual(values = c("#780000", "#c1121f")) +
    guides(fill = FALSE)
  if ((!(is.null(consequence_col)))) {
    cover <- data.frame(data %>%
      group_by(get(colnames(data)[classification_num]), Classification_3, get(colnames(data)[consequence_num])) %>%
      summarize(n = n()))
    colnames(cover)[1] <- "Classification"
    colnames(cover)[3] <- "Consequence"
    cover <- arrange(cover, desc(n), Consequence)
    count <- data.frame(data %>%
      group_by(Classification_3) %>%
      summarize(n = n()))
    cover <- full_join(cover, count, by = c("Classification_3"))
    cover <- cover %>% mutate(percent = n.x / n.y)
    cover2 <- data.frame(data %>%
      group_by(Classification_3, get(colnames(data)[consequence_num])) %>%
      summarize(n = n()))
    colnames(cover2)[2] <- "Consequence"
    cover2 <- arrange(cover2, desc(n), Consequence)
    count <- data.frame(data %>%
      group_by(Classification_3) %>%
      summarize(n = n()))
    cover2 <- full_join(cover2, count, by = c("Classification_3"))
    cover2 <- cover2 %>% mutate(percent = n.x / n.y)
    cover$Consequence <- factor(cover$Consequence, levels = as.character(unique(c(cover2[cover2$Classification_3 == "P/LP", ]$Consequence, cover2[cover2$Classification_3 == "VUS", ]$Consequence, cover2[cover2$Classification_3 == "B/LB", ]$Consequence))))
    cover$Classification_3 <- factor(cover$Classification_3, levels = c("P/LP", "VUS", "B/LB"))
    if ("Pathogenic" %in% names(table(cover$Classification))) {
      cover$Classification <- factor(cover$Classification, levels = c("Pathogenic", "Likely Pathogenic", "Uncertain Significance", "Likely Benign", "Benign"))
    } else {
      cover$Classification <- factor(cover$Classification, levels = c("P", "LP", "VUS", "LB", "B"))
    }
    cover3 <- data.frame(data %>%
      group_by(Classification_3, get(colnames(data)[consequence_num])) %>%
      summarize(n = n()))
    colnames(cover3)[2] <- "Consequence"
    cover3 <- arrange(cover3, desc(n), Consequence)
    cover3 <- full_join(cover3, count, by = c("Classification_3"))
    cover3 <- cover3 %>% mutate(percent = n.x / n.y)
    cover3$Classification_3 <- factor(cover3$Classification_3, levels = c("P/LP", "VUS", "B/LB"))
    cover3$Consequence <- factor(cover3$Consequence, levels = as.character(unique(c(cover3[cover3$Classification_3 == "P/LP", ]$Consequence, cover3[cover3$Classification_3 == "VUS", ]$Consequence, cover3[cover3$Classification_3 == "B/LB", ]$Consequence))))
    cover4 <- full_join(cover, cover3, by = c("Classification_3", "Consequence"))
    cover4$number <- ifelse(duplicated(paste(cover4$Classification_3, cover4$Consequence, sep = "_")), "", cover4$n.x.y)
    f4 <- ggplot(cover4, aes(x = Consequence, y = percent.x, fill = Classification)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      scale_fill_manual(values = c("#780000", "#c1121f", "#bfdbf7", "#cad2c5", "#84a98c")) +
      guides(fill = guide_legend(title = NULL)) +
      theme(axis.title = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
      facet_grid(Classification_3 ~ .) +
      ylim(0, 1) +
      guides(fill = FALSE) +
      geom_text(aes(x = Consequence, y = percent.y, label = number), size = 3, vjust = 0, nudge_y = 0.01)
    fall <- (f1 + f2) / (f3 + f4) + plot_annotation(tag_levels = "A")
  } else {
    fall <- (f1 + f2) / f3 + plot_annotation(tag_levels = "A")
  }
  return(fall)
}

#' Visualize the results of LR+ for each evaluated cutoff
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param direction The direction of evidence pathogenic (Pathogenic or Benign)
#' @param op_list A list of odds path corresponding to each level of evidence strength
#'
#' @return Figures
#' @import ComplexHeatmap
#' @import grid
#' @import scales
#' @import circlize
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' data("LR_result")
#' op_list <- c(2.08, 4.33, 18.70, 350)
#' heatmap_LR(LR_result, "Pathogenic", op_list)
#'
heatmap_LR <- function(data, direction, op_list) {
  if (direction != "Pathogenic" && direction != "Benign") {
    return(message("Error,the direction of evidence pathogenic must be Pathogenic or Benign"))
  }
  fea_plot <- data[, c(1:7, 13, 14, 12, 15, 16)]
  if (direction == "Benign") {
    fea_plot <- data[, c(1:7, 13, 14, 12, 18, 19)]
  }
  p1 <- Heatmap(fea_plot[, c(1:4)], cluster_rows = FALSE, name = "Number of variants", row_names_side = "left", cluster_columns = FALSE, column_names_rot = 45, rect_gp = gpar(col = "black", lwd = 1.5), col = circlize::colorRamp2(unname(quantile(as.vector(t(fea_plot[, c(1:4)])), c(0.05, 0.25, 0.5, 0.75, 0.95))), c("white", "#CBCDE0", "#A4ABD6", "#ACA0D2", "#8076A3"), transparency = 0.2), column_dend_height = unit(2, "cm"), row_names_gp = gpar(fontsize = 10, col = "black"), column_names_gp = gpar(fontsize = 10, col = c(rep("#8076A3", 4))), cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.0f", fea_plot[, c(1:4)][i, j]), x, y, gp = gpar(fontsize = 14))
  })
  p2 <- Heatmap(fea_plot[, c(5:10)], cluster_rows = FALSE, name = "Evaluation metrics", row_names_side = "left", cluster_columns = FALSE, column_names_rot = 45, rect_gp = gpar(col = "black", lwd = 1.5), col = circlize::colorRamp2(unname(quantile(as.vector(t(fea_plot[, c(5:10)])), c(0.05, 0.25, 0.5, 0.75, 0.95))), c("white", "#aed4e5", "#81b5d5", "#5795c7", "#3371b3"), transparency = 0.2), column_dend_height = unit(2, "cm"), row_names_gp = gpar(fontsize = 10, col = "black"), column_names_gp = gpar(fontsize = 10, col = c(rep("#3371b3", 6))), cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", fea_plot[, c(5:10)][i, j]), x, y, gp = gpar(fontsize = 14))
  })
  value <- rep(0, nrow(fea_plot) * 2)
  cols <- structure(rep("#fef9ef", nrow(fea_plot) * 2), names = c(fea_plot[, c(11)], fea_plot[, c(12)]))
  if (direction == "Pathogenic") {
    strength_cols <- structure(c("#e5989b", "#b5838d", "#bf4342", "#780000"), names = c(1, 2, 3, 4))
    anno_cols <- structure(c("#fef9ef", "#e5989b", "#b5838d", "#bf4342", "#780000"), names = c(0, op_list))
  }
  if (direction == "Benign") {
    strength_cols <- structure(c("#eff6e0", "#aec3b0", "#84a98c", "#52796f"), names = c(1, 2, 3, 4))
    anno_cols <- structure(c("#fef9ef", "#eff6e0", "#aec3b0", "#84a98c", "#52796f"), names = c(0, op_list))
  }
  if (direction == "Pathogenic") {
    for (i in 1:length(fea_plot[, 11])) {
      for (k in 1:length(op_list)) {
        if (fea_plot[, 11][i] >= op_list[k]) {
          value[i] <- op_list[k]
          cols[i] <- strength_cols[k]
        }
      }
    }
    for (i in 1:length(fea_plot[, 12])) {
      for (k in 1:length(op_list)) {
        if (fea_plot[, 12][i] >= op_list[k]) {
          value[i + nrow(fea_plot)] <- op_list[k]
          cols[i + nrow(fea_plot)] <- strength_cols[k]
        }
      }
    }
  }
  if (direction == "Benign") {
    for (i in 1:length(fea_plot[, 11])) {
      for (k in 1:length(op_list)) {
        if (fea_plot[, 11][i] <= op_list[k]) {
          value[i] <- op_list[k]
          cols[i] <- strength_cols[k]
        }
      }
    }
    for (i in 1:length(fea_plot[, 12])) {
      for (k in 1:length(op_list)) {
        if (fea_plot[, 12][i] <= op_list[k]) {
          value[i + nrow(fea_plot)] <- op_list[k]
          cols[i + nrow(fea_plot)] <- strength_cols[k]
        }
      }
    }
  }
  strength <- factor(value[(nrow(fea_plot) + 1):(nrow(fea_plot) * 2)])
  cols_list <- anno_cols[which(names(anno_cols) %in% levels(strength))]
  anno_lable <- c("NA", "Supporting", "Moderate", "Strong", "Very Strong")
  p3 <- Heatmap(fea_plot[, c(11:12)], cluster_rows = FALSE, name = "Strength of evidence", row_names_side = "left", column_names_rot = 45, cluster_columns = FALSE, rect_gp = gpar(col = "black", lwd = 1.5), col = cols, column_dend_height = unit(2, "cm"), row_names_gp = gpar(fontsize = 10, col = "#1f78b4"), column_names_gp = gpar(fontsize = 10, col = "#780000"), cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", fea_plot[, c(11:12)][i, j]), x, y, gp = gpar(fontsize = 14))
  }, show_heatmap_legend = FALSE, right_annotation = HeatmapAnnotation(which = "row", Strength = strength, col = list(Strength = cols_list), annotation_name_rot = 45, annotation_legend_param = list(Strength = list(at = levels(strength), labels = anno_lable[which(names(anno_cols) %in% levels(strength))]))))
  Figure_heatmap <- p1 + p2 + p3
  return(Figure_heatmap)
}

#' Generate plots depicting the results of lr+ for each tested cutoff
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param postp_list A list of posterior probability corresponding to each level of evidence strength
#' @param direction The direction of evidence pathogenic (Pathogenic or Benign)
#' @import ggplot2
#' @importFrom utils globalVariables
#' @return Figures
#' @export
#'
#' @examples
#' data("lr_CI_result")
#' # data <- add_info(ClinVar_2019_dataset, "clnsig")
#' # local_bootstrapped_lr(data, "PrimateAI_score", 0.0441, 10000, 100, 0.01, "test_dir")
#' postp_list <- c(0.100, 0.211, 0.608, 0.981)
#' # lr_CI_result <- lr_CI(30, "test_dir")
#' plot_lr(lr_CI_result, "Pathogenic", postp_list)
#'
plot_lr <- function(data, direction, postp_list) {
  if (direction != "Pathogenic" && direction != "Benign") {
    return(message("Error,the direction of evidence pathogenic must be Pathogenic or Benign"))
  }
  if (direction == "Pathogenic") {
    Figurelr <- ggplot() +
      geom_line(data = data, mapping = aes(x = test_cutoff, y = Posterior), colour = "black", size = 1) +
      geom_line(data = data, mapping = aes(x = test_cutoff, y = Posterior1), colour = "grey", size = 1) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      geom_hline(yintercept = postp_list[4], colour = "red", linetype = 1, size = 1) +
      geom_hline(yintercept = postp_list[3], colour = "red", linetype = 2, size = 1) +
      geom_hline(yintercept = postp_list[2], colour = "red", linetype = 3, size = 1) +
      geom_hline(yintercept = postp_list[1], colour = "red", linetype = 4, size = 1) +
      ylab("Posterior probability")
  }
  if (direction == "Benign") {
    Figurelr <- ggplot() +
      geom_line(data = data, mapping = aes(x = test_cutoff, y = Posterior), colour = "black", size = 1) +
      geom_line(data = data, mapping = aes(x = test_cutoff, y = Posterior1), colour = "grey", size = 1) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid")) +
      geom_hline(yintercept = postp_list[4], colour = "#0971BA", linetype = 1, size = 1) +
      geom_hline(yintercept = postp_list[3], colour = "#0971BA", linetype = 2, size = 1) +
      geom_hline(yintercept = postp_list[2], colour = "#0971BA", linetype = 3, size = 1) +
      geom_hline(yintercept = postp_list[1], colour = "#0971BA", linetype = 4, size = 1) +
      ylab("Posterior probability") +
      ylim(0.9, 1)
  }
  return(Figurelr)
}

#' Visualize the correlation of ACMG/AMP evidence
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param evidence_col The column name for the ACMG/AMP evidence (str)
#'
#'
#' @import stringr
#' @import reshape2
#' @importFrom dplyr mutate
#' @import ggplot2
#' @import corrplot
#' @importFrom Hmisc rcorr
#' @import RColorBrewer
#'
#' @return Figures
#' @export
#'
#' @examples
#' data("ClinGen_dataset")
#' evidence_corplot(ClinGen_dataset, "Applied Evidence Codes (Met)")
#'
evidence_corplot <- function(data, evidence_col) {
  evidence_num <- which(colnames(data) == evidence_col)
  P_evidence <- c("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM2", "PM3", "PM4", "PM5", "PM6", "PP1", "PP2", "PP3", "PP4", "PP5")
  B_evidence <- c("BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7")
  data_temp <- data
  for (i in 1:length(P_evidence))
  {
    data_temp <- mutate(data_temp, test = ifelse(str_detect(data_temp[, evidence_num], P_evidence[i]) & (!grepl(paste(P_evidence[i], "", sep = "_"), data_temp[, evidence_num])), ifelse(grepl("PVS", P_evidence[i]), 8, ifelse(grepl("PS", P_evidence[i]), 4, ifelse(grepl("PM", P_evidence[i]), 2, ifelse(grepl("PP", P_evidence[i]), 1, 0)))), ifelse(str_detect(data_temp[, evidence_num], paste(P_evidence[i], "_VeryStrong", sep = "")), 8, ifelse(str_detect(data_temp[, evidence_num], paste(P_evidence[i], "_Strong", sep = "")), 4, ifelse(str_detect(data_temp[, evidence_num], paste(P_evidence[i], "_Moderate", sep = "")), 2, ifelse(str_detect(data_temp[, evidence_num], paste(P_evidence[i], "_Supporting", sep = "")), 1, 0)))))) %>% plyr::rename(c("test" = P_evidence[i]))
  }
  for (i in 1:length(B_evidence))
  {
    data_temp <- mutate(data_temp, test = ifelse(str_detect(data_temp[, evidence_num], B_evidence[i]) & (!grepl(paste(B_evidence[i], "", sep = "_"), data_temp[, evidence_num])), ifelse(grepl("BS", P_evidence[i]), 4, ifelse(grepl("BP", P_evidence[i]), 1, 0)), ifelse(str_detect(data_temp[, evidence_num], paste(B_evidence[i], "_VeryStrong", sep = "")), 8, ifelse(str_detect(data_temp[, evidence_num], paste(B_evidence[i], "_Strong", sep = "")), 4, ifelse(str_detect(data_temp[, evidence_num], paste(B_evidence[i], "_Moderate", sep = "")), 2, ifelse(str_detect(data_temp[, evidence_num], paste(B_evidence[i], "_Supporting", sep = "")), 1, 0)))))) %>% plyr::rename(c("test" = B_evidence[i]))
  }
  Evidence <- data_temp[, (ncol(data_temp) - 26):ncol(data_temp)]
  Evidence <- Evidence[, which(colSums(Evidence) != 0)]
  cor_result <- rcorr(as.matrix(Evidence), type = "spearman")
  col <- colorRampPalette(c("#4B5AA1", "white", "#D2352C"))(5)
  cor_plot <- corrplot(
    cor_result$r,
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.pos = "td",
    tl.cex = .8,
    tl.srt = 45,
    pch.cex = 0.8,
    is.corr = F,
    sig.level = 0.05,
    insig = "label_sig",
    col = rev(COL2("PiYG", 10)),
    col.lim = c(-1, 1)
  )
  return(cor_plot)
}
