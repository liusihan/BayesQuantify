#' Classifying variants into five distinct categories according to the 2015 ACMG/AMP guidelines
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param evidence_col The column name for ACMG evidence(str)
#' @param public_p DataFrame contains reported P/LP variants
#'
#' @return A new DataFrame that incorporates the input data and the results of variant classification
#' @import stringr
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinGen_dataset")
#' ACMG_Classification(ClinGen_dataset, "Applied Evidence Codes (Met)")
#' }
#'
ACMG_Classification <- function(data, evidence_col, public_p = NULL) {
  evidence_num <- which(colnames(data) == evidence_col)
  if (!"BayesQuantify_ACMG" %in% colnames(data)) {
    data$BayesQuantify_ACMG <- ""
  }
  if ((!(is.null(public_p)))){
    required_public_cols <- c("chr", "pos", "ref", "alt", "classification")
    if (!all(required_public_cols %in% tolower(colnames(public_p)))) {
      stop("The following columns must present in the reported P/LP variants datasets: chr, pos, ref, alt, classification")
    }
    required_data_cols <- c("chr", "pos", "ref", "alt")
    if (!all(required_data_cols %in% tolower(colnames(data)))) {
      stop("The following columns must present in the input data framae: chr, pos, ref, alt")
    }
    public_keys <- paste(public_p$chr, public_p$pos, public_p$ref, public_p$alt, sep = "_")
    class_map <- setNames(public_p$classification, public_keys)
    data_keys <- paste(data$chr, data$pos, data$ref, data$alt, sep = "_")
    matched_indices <- which(data_keys %in% public_keys)
    data$BayesQuantify_ACMG[matched_indices] <- class_map[data_keys[matched_indices]]
  }
  for (k in 1:nrow(data)) {
    if (data$BayesQuantify_ACMG[k] != ""){
      next
    }
    if (is.na(data[k, evidence_num])) {
      data$BayesQuantify_ACMG[k] <- "VUS"
    }
    if (!is.na(data[k, evidence_num])) {
      # remove evidence for benign
      p_evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      # remove evidence for pathogenic
      b_evidence <- gsub(pattern = "[P][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      p_VeryStrong <- str_count(string = p_evidence, pattern = "_Very Strong") + str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      p_Strong <- str_count(string = p_evidence, pattern = "_Strong") + str_count(string = p_evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      p_Moderate <- str_count(string = p_evidence, pattern = "_Moderate") + str_count(string = p_evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      p_Supporting <- str_count(string = p_evidence, pattern = "_Supporting") + str_count(string = p_evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
      Stand_alone <- str_count(string = b_evidence, pattern = "BA1[^_]|BA1$") + str_count(string = b_evidence, pattern = "_Stand Alone")
      b_Strong <- str_count(string = b_evidence, pattern = "_Strong") + str_count(string = b_evidence, pattern = "(BS[0-9][^_])|(BS[0-9]$)")
      b_Moderate <- str_count(string = b_evidence, pattern = "_Moderate")
      b_Supporting <- str_count(string = b_evidence, pattern = "_Supporting") + str_count(string = b_evidence, pattern = "(BP[0-9][^_])|(BP[0-9]$)")
      PVS1 <- str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      PM2 <- str_count(string = p_evidence, pattern = "PM2_Supporting")
      p_class <- ""
      b_class <- ""
      if (p_VeryStrong >= 1 && p_Moderate == 1) {
        p_class <- "LP"
      }
      if (p_Strong == 1 && (p_Moderate == 1 || p_Moderate == 2)) {
        p_class <- "LP"
      }
      if (p_Strong == 1 && p_Supporting >= 2) {
        p_class <- "LP"
      }
      if (p_Moderate >= 3) {
        p_class <- "LP"
      }
      if (p_Moderate == 2 && p_Supporting >= 2) {
        p_class <- "LP"
      }
      if (p_Moderate == 1 && p_Supporting >= 4) {
        p_class <- "LP"
      }
      if (PVS1 == 1 && PM2 == 1) {
        p_class <- "LP"
      }
      if (p_VeryStrong >= 1) {
        if (p_Strong >= 1 || p_Moderate >= 2 || p_Supporting >= 2 || (p_Moderate == 1 && p_Supporting == 1)) {
          p_class <- "P"
        }
      }
      if (p_VeryStrong >= 2) {
        p_class <- "P"
      }
      if (p_Strong >= 2) {
        p_class <- "P"
      }
      if (p_Strong == 1) {
        if (p_Moderate >= 3 || (p_Moderate == 2 && p_Supporting >= 2) || (p_Moderate == 1 && p_Supporting >= 4)) {
          p_class <- "P"
        }
      }
      if ((b_Strong == 1 && b_Supporting == 1) || (b_Supporting >= 2) || (b_Strong == 1 && b_Moderate == 1) || ((b_Supporting == 1 && b_Moderate >= 1))) {
        b_class <- "LB"
      }
      if (Stand_alone >= 1 || b_Strong >= 2) {
        b_class <- "B"
      }
      if (b_class == "" && p_class == "") {
        data$BayesQuantify_ACMG[k] <- "VUS"
      } else if (b_class != "" && p_class != "") {
        data$BayesQuantify_ACMG[k] <- "VUS"
      } else if (b_class == "" && p_class != "") {
        data$BayesQuantify_ACMG[k] <- p_class
      } else {
        data$BayesQuantify_ACMG[k] <- b_class
      }
    }
  }
  return(data)
}

#' Classifying variants into five distinct categories according to the Bayesian classification framework
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param evidence_col The column name for ACMG evidence(str)
#' @param prior_p The prior probability of pathogenicity (proportion of P/LP variants in a set of variants)
#' @param op_vs Odds of pathogenicity (OP) of "Very String"
#' @param public_p DataFrame contains reported P/LP variants
#'
#' @import stringr
#' @return A new DataFrame that incorporates the input data and the results of variant classification
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinGen_dataset")
#' BCF(ClinGen_dataset, "Applied Evidence Codes (Met)", 0.1, 350)
#' }
#'
BCF <- function(data, evidence_col, prior_p, op_vs, public_p = NULL) {
  evidence_num <- which(colnames(data) == evidence_col)
  data$BCF_postp <- ""
  data$BCF_classification <- ""
  op_s <- op_vs^(2^-1)
  op_m <- op_vs^(2^-2)
  op_p <- op_vs^(2^-3)
  if ((!(is.null(public_p)))){
    required_public_cols <- c("chr", "pos", "ref", "alt", "classification")
    if (!all(required_public_cols %in% tolower(colnames(public_p)))) {
      stop("The following columns must present in the reported P/LP variants datasets: chr, pos, ref, alt, classification")
    }
    required_data_cols <- c("chr", "pos", "ref", "alt")
    if (!all(required_data_cols %in% tolower(colnames(data)))) {
      stop("The following columns must present in the input data framae: chr, pos, ref, alt")
    }
    public_keys <- paste(public_p$chr, public_p$pos, public_p$ref, public_p$alt, sep = "_")
    class_map <- setNames(public_p$classification, public_keys)
    data_keys <- paste(data$chr, data$pos, data$ref, data$alt, sep = "_")
    matched_indices <- which(data_keys %in% public_keys)
    data$BCF_classification[matched_indices] <- class_map[data_keys[matched_indices]]
  }
  for (k in 1:nrow(data)) {
    if (data$BCF_classification[k] != ""){
      # remove evidence for benign
      p_evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      # remove evidence for pathogenic
      b_evidence <- gsub(pattern = "[P][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      p_VeryStrong <- str_count(string = p_evidence, pattern = "_Very Strong") + str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      p_Strong <- str_count(string = p_evidence, pattern = "_Strong") + str_count(string = p_evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      p_Moderate <- str_count(string = p_evidence, pattern = "_Moderate") + str_count(string = p_evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      p_Supporting <- str_count(string = p_evidence, pattern = "_Supporting") + str_count(string = p_evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
      Stand_alone <- str_count(string = b_evidence, pattern = "BA1[^_]|BA1$")+ str_count(string = b_evidence, pattern = "_Stand Alone")
      b_VeryStrong <- str_count(string = b_evidence, pattern = "_Very Strong")
      b_Strong <- str_count(string = b_evidence, pattern = "_Strong") + str_count(string = b_evidence, pattern = "(BS[0-9][^_])|(BS[0-9]$)")
      b_Moderate <- str_count(string = b_evidence, pattern = "_Moderate")
      b_Supporting <- str_count(string = b_evidence, pattern = "_Supporting") + str_count(string = b_evidence, pattern = "(BP[0-9][^_])|(BP[0-9]$)")
      op <- op_vs^(p_Supporting / 8 + p_Moderate / 4 + p_Strong / 2 + p_VeryStrong / 1 - b_Supporting / 8 - b_Moderate / 4 - b_Strong / 2 - b_VeryStrong / 1)
      post_p <- round((op * prior_p) / ((op - 1) * prior_p + 1),3)
      data$BCF_postp[k] <- post_p
      next
    }
    if (is.na(data[k, evidence_num])) {
      data$BCF_postp[k] <- prior_p
      data$BCF_classification[k] <- "VUS"
    }
    if (!is.na(data[k, evidence_num])) {
      # remove evidence for benign
      p_evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      # remove evidence for pathogenic
      b_evidence <- gsub(pattern = "[P][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      p_VeryStrong <- str_count(string = p_evidence, pattern = "_Very Strong") + str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      p_Strong <- str_count(string = p_evidence, pattern = "_Strong") + str_count(string = p_evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      p_Moderate <- str_count(string = p_evidence, pattern = "_Moderate") + str_count(string = p_evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      p_Supporting <- str_count(string = p_evidence, pattern = "_Supporting") + str_count(string = p_evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
      Stand_alone <- str_count(string = b_evidence, pattern = "BA1[^_]|BA1$")+ str_count(string = b_evidence, pattern = "_Stand Alone")
      b_VeryStrong <- str_count(string = b_evidence, pattern = "_Very Strong")
      b_Strong <- str_count(string = b_evidence, pattern = "_Strong") + str_count(string = b_evidence, pattern = "(BS[0-9][^_])|(BS[0-9]$)")
      b_Moderate <- str_count(string = b_evidence, pattern = "_Moderate")
      b_Supporting <- str_count(string = b_evidence, pattern = "_Supporting") + str_count(string = b_evidence, pattern = "(BP[0-9][^_])|(BP[0-9]$)")
      op <- op_vs^(p_Supporting / 8 + p_Moderate / 4 + p_Strong / 2 + p_VeryStrong / 1 - b_Supporting / 8 - b_Moderate / 4 - b_Strong / 2 - b_VeryStrong / 1)
      post_p <- round((op * prior_p) / ((op - 1) * prior_p + 1),3)
      data$BCF_postp[k] <- post_p
      if (post_p < 0.001) {
        data$BCF_classification[k] <- "B"
      } else if (post_p < 0.1) {
        data$BCF_classification[k] <- "LB"
      } else if (post_p > 0.99) {
        data$BCF_classification[k] <- "P"
      } else if (post_p >= 0.9) {
        data$BCF_classification[k] <- "LP"
      } else {
        data$BCF_classification[k] <- "VUS"
      }
      if (Stand_alone == 1) {
        data$BCF_classification[k] <- "B"
        data$BCF_postp[k] <- ""
      }
    }
  }
  return(data)
}

#' Classifying variants into five distinct categories according to the scaled point system.
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param evidence_col The column name for ACMG evidence(str)
#' @param public_p DataFrame contains reported P/LP variants
#'
#' @import stringr
#' @return A new DataFrame that incorporates the input data and the results of variant classification
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinGen_dataset")
#' Point_Classification(ClinGen_dataset, "Applied Evidence Codes (Met)")
#' }
#'
Point_Classification <- function(data, evidence_col, public_p = NULL) {
  evidence_num <- which(colnames(data) == evidence_col)
  data$scaled_point <- ""
  data$point_classification <- ""
  if ((!(is.null(public_p)))){
    required_public_cols <- c("chr", "pos", "ref", "alt", "classification")
    if (!all(required_public_cols %in% tolower(colnames(public_p)))) {
      stop("The following columns must present in the reported P/LP variants datasets: chr, pos, ref, alt, classification")
    }
    required_data_cols <- c("chr", "pos", "ref", "alt")
    if (!all(required_data_cols %in% tolower(colnames(data)))) {
      stop("The following columns must present in the input data framae: chr, pos, ref, alt")
    }
    public_keys <- paste(public_p$chr, public_p$pos, public_p$ref, public_p$alt, sep = "_")
    class_map <- setNames(public_p$classification, public_keys)
    data_keys <- paste(data$chr, data$pos, data$ref, data$alt, sep = "_")
    matched_indices <- which(data_keys %in% public_keys)
    data$point_classification[matched_indices] <- class_map[data_keys[matched_indices]]
  }
  for (k in 1:nrow(data)) {
    if (data$point_classification[k] != ""){
      # remove evidence for benign
      p_evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      # remove evidence for pathogenic
      b_evidence <- gsub(pattern = "[P][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      p_VeryStrong <- str_count(string = p_evidence, pattern = "_Very Strong") + str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      p_Strong <- str_count(string = p_evidence, pattern = "_Strong") + str_count(string = p_evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      p_Moderate <- str_count(string = p_evidence, pattern = "_Moderate") + str_count(string = p_evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      p_Supporting <- str_count(string = p_evidence, pattern = "_Supporting") + str_count(string = p_evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
      Stand_alone <- str_count(string = b_evidence, pattern = "BA1[^_]|BA1$")+ str_count(string = b_evidence, pattern = "_Stand Alone")
      b_VeryStrong <- str_count(string = b_evidence, pattern = "_Very Strong")
      b_Strong <- str_count(string = b_evidence, pattern = "_Strong") + str_count(string = b_evidence, pattern = "(BS[0-9][^_])|(BS[0-9]$)")
      b_Moderate <- str_count(string = b_evidence, pattern = "_Moderate")
      b_Supporting <- str_count(string = b_evidence, pattern = "_Supporting") + str_count(string = b_evidence, pattern = "(BP[0-9][^_])|(BP[0-9]$)")
      point_score <- 1 * p_Supporting + 2 * p_Moderate + 4 * p_Strong + 8 * p_VeryStrong - 1 * b_Supporting - 2 * b_Moderate - 4 * b_Strong - 8 * b_VeryStrong
      data$scaled_point[k] <- point_score
      next
    }
    if (is.na(data[k, evidence_num])) {
      data$point_classification[k] <- "VUS"
      data$scaled_point[k] <- 0
    }
    if (!is.na(data[k, evidence_num])) {
      # remove evidence for benign
      p_evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      # remove evidence for pathogenic
      b_evidence <- gsub(pattern = "[P][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      p_VeryStrong <- str_count(string = p_evidence, pattern = "_Very Strong") + str_count(string = p_evidence, pattern = "PVS1[^_]|PVS1$")
      p_Strong <- str_count(string = p_evidence, pattern = "_Strong") + str_count(string = p_evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      p_Moderate <- str_count(string = p_evidence, pattern = "_Moderate") + str_count(string = p_evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      p_Supporting <- str_count(string = p_evidence, pattern = "_Supporting") + str_count(string = p_evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
      Stand_alone <- str_count(string = b_evidence, pattern = "BA1[^_]|BA1$")+ str_count(string = b_evidence, pattern = "_Stand Alone")
      b_VeryStrong <- str_count(string = b_evidence, pattern = "_Very Strong")
      b_Strong <- str_count(string = b_evidence, pattern = "_Strong") + str_count(string = b_evidence, pattern = "(BS[0-9][^_])|(BS[0-9]$)")
      b_Moderate <- str_count(string = b_evidence, pattern = "_Moderate")
      b_Supporting <- str_count(string = b_evidence, pattern = "_Supporting") + str_count(string = b_evidence, pattern = "(BP[0-9][^_])|(BP[0-9]$)")
      point_score <- 1 * p_Supporting + 2 * p_Moderate + 4 * p_Strong + 8 * p_VeryStrong - 1 * b_Supporting - 2 * b_Moderate - 4 * b_Strong - 8 * b_VeryStrong
      data$scaled_point[k] <- point_score
      if (point_score >= 10) {
        data$point_classification[k] <- "P"
      } else if (point_score >= 6 && point_score <= 9) {
        data$point_classification[k] <- "LP"
      } else if (point_score >= 0 && point_score <= 5) {
        data$point_classification[k] <- "VUS"
      } else if (point_score >= -6 && point_score <= -1) {
        data$point_classification[k] <- "LB"
      } else if (point_score <= -7) {
        data$point_classification[k] <- "B"
      } else {
        data$point_classification[k] <- "NA"
      }
      if (Stand_alone == 1) {
        data$point_classification[k] <- "B"
        data$scaled_point[k] <- ""
      }
    }
  }
  return(data)
}
