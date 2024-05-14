#' Count the number of "supporting", "moderate", "strong" and "very strong" strengths of evidence for pathogenicity
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param classification_col The column name for variant classification (str). Variants should be classified into five distinct categories: "P," "LP," "B," "LB," and "VUS."
#'
#'
#' @return A new DataFrame that includes the input data and four new columns, these four columns count the number of different pathogenic evidence strengths for each variant, which can be used for further categorization
#' @export
#'
#' @examples
#' data("ClinGen_dataset")
#' ClinGen_dataset <- add_info(ClinGen_dataset, "Assertion")
#'
add_info <- function(data, classification_col) {
  classification_num <- which(colnames(data) == classification_col)
  data$Classification_P <- ""
  data$Classification_3 <- ""
  for (i in 1:nrow(data)) {
    if (data[i, classification_num] == "P" || data[i, classification_num] == "LP" || data[i, classification_num] == "Pathogenic" || data[i, classification_num] == "Likely Pathogenic" || data[i, classification_num] == "Likely_pathogenic" || data[i, classification_num] == "Pathogenic/Likely_pathogenic") {
      data$Classification_P[i] <- "P"
    } else {
      data$Classification_P[i] <- "NonP"
    }
  }
  for (i in 1:nrow(data)) {
    if (data[i, classification_num] == "P" || data[i, classification_num] == "LP" || data[i, classification_num] == "Pathogenic" || data[i, classification_num] == "Likely Pathogenic" || data[i, classification_num] == "Likely_pathogenic" || data[i, classification_num] == "Pathogenic/Likely_pathogenic") {
      data$Classification_3[i] <- "P/LP"
    } else if (data[i, classification_num] == "VUS" || data[i, classification_num] == "Uncertain Significance") {
      data$Classification_3[i] <- "VUS"
    } else {
      data$Classification_3[i] <- "B/LB"
    }
  }
  return(data)
}

#' Variants of uncertain significance (VUS) were categorized into six levels (hot, warm, tepid, cool, cold, and ice cold), according to the Association for Clinical Genomic Science (ACGS) Best Practice Guidelines.
#' hot: 1 very strong or 1 strong + 1 supporting or 2 moderate + 1 supporting or 1 moderate + 3 supporting evidence;
#' warm: 1 strong or 2 moderate or 1 moderate + 2 supporting or 4 supporting evidence;
#' tepid: 1 moderate + 1 supporting or 3 supporting evidence;
#' cool: 1 moderate or 2 supporting evidence;
#' cold: 1 supporting evidence;
#' ice cold: no supporting evidence (https://www.acgs.uk.com/quality/best-practice-guidelines/#VariantGuidelines).
#' Variants classified as cool, cold, or ice cold were considered as benign-leaning VUS, unlikely to be disease-causing.
#' Variants classified as hot, warm, or tepid were considered to be pathogenic-leaning VUS.
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param classification_col The column name for variant classification (str). Variants should be classified into five distinct categories: "P," "LP," "B," "LB," and "VUS."
#' @param evidence_col The column name for ACMG evidence(str). The content of this column should be composed of evidence names and their strengths, connected by semicolons or comma, such as "PM2_Supporting;PM5;BP4"
#'
#' @import stringr
#' @return A new DataFrame that includes the input data and VUS classification
#' @export
#'
#' @examples
#' data("ClinGen_dataset")
#' ClinGen_dataset <- VUS_classify(ClinGen_dataset, "Assertion", "Applied Evidence Codes (Met)")
#'
VUS_classify <- function(data, classification_col, evidence_col) {
  classification_num <- which(colnames(data) == classification_col)
  evidence_num <- which(colnames(data) == evidence_col)
  data$VeryStrong <- 0
  data$Strong <- 0
  data$Moderate <- 0
  data$Supporting <- 0
  data$VUS_class <- ""
  for (k in 1:nrow(data)) {
    # remove evidence for benign
    if (!is.na(data[k, evidence_num])) {
      evidence <- gsub(pattern = "[B][A-Z]+\\d[:_:]?[A-Z]*", replacement = "", x = data[k, evidence_num])
      data$VeryStrong[k] <- str_count(string = evidence, pattern = "_Very Strong") + str_count(string = evidence, pattern = "PVS1[^_]|PVS1$")
      data$Strong[k] <- str_count(string = evidence, pattern = "_Strong") + str_count(string = evidence, pattern = "(PS[0-9][^_])|(PS[0-9]$)")
      data$Moderate[k] <- str_count(string = evidence, pattern = "_Moderate") + str_count(string = evidence, pattern = "(PM[0-9][^_])|(PM[0-9]$)")
      data$Supporting[k] <- str_count(string = evidence, pattern = "_Supporting") + str_count(string = evidence, pattern = "(PP[0-9][^_])|(PP[0-9]$)")
    }
  }
  for (i in 1:nrow(data)) {
    if (data[i, classification_num] == "VUS" || data[i, classification_num] == "Uncertain Significance") {
      if ((data$Strong[i] == 1 && data$Supporting[i] == 1) || (data$Moderate[i] == 2 && data$Supporting[i] == 1) || (data$Moderate[i] == 1 && data$Supporting[i] == 3)) {
        data$VUS_class[i] <- "Hot"
      } else if ((data$Strong[i] == 1) || (data$Moderate[i] == 2) || (data$Moderate[i] == 1 && data$Supporting[i] == 2) || (data$Supporting[i] == 4)) {
        data$VUS_class[i] <- "Warm"
      } else if ((data$Moderate[i] == 1 && data$Supporting[i] == 1) || (data$Supporting[i] == 3)) {
        data$VUS_class[i] <- "Tepid"
      } else if ((data$Moderate[i] == 1) || (data$Supporting[i] == 2)) {
        data$VUS_class[i] <- "Cool"
      } else if (data$Supporting[i] == 1) {
        data$VUS_class[i] <- "Cold"
      } else if ((data$Supporting[i] == 0) && (data$Moderate[i] == 0) && (data$Strong[i] == 0) && (data$VeryStrong[i] == 0)) {
        data$VUS_class[i] <- "IceCold"
      } else {
        data$VUS_class[i] <- "NA"
      }
    }
  }
  return(data)
}
