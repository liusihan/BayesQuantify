#' Introducing columns to assess if the observed value is above (1) or below (0) a tested cutoff. A value of 1 indicates being above the tested cutoff, while 0 indicates being below the tested cutoff
#'
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param feature The column name that requires testing for optimizing the thresholds
#' @param range Evaluated intervals
#' @param criteria ACMG/AMP guidelines criteria (str)
#'
#' @importFrom plyr rename
#' @importFrom dplyr `%>%` mutate
#'
#' @return A fresh DataFrame incorporating the input data with additional column
#' @export
#'
#' @examples
#' data("VCI_data")
#' discrete_cutoff(VCI_data, "Applied Evidence Codes (Met)", criteria = "PM2")
#'
discrete_cutoff <- function(data, feature, range = NULL, criteria = NULL) {
  feature_col <- which(colnames(data) == feature)
  if ((!(is.null(range))) && is.null(criteria)) {
    for (k in range) {
      data <- data %>%
        mutate(test = ifelse(eval(parse(text = feature)) > k, 1, 0)) %>%
        plyr::rename(c("test" = paste(feature, k, sep = "_")))
    }
    return(data)
  } else if (is.null(range) && (!(is.null(criteria)))) {
    data <- mutate(data, test = ifelse(str_detect(data[, feature_col], criteria), 1, 0)) %>% plyr::rename(c("test" = paste(criteria, "", sep = "")))
    return(data)
  } else {
    return(message("Error: Please use either the range or criteria parameters."))
  }
}

#' Calculating positive likelihood ratio (LR) for each tested cutoff (for discrete cutoffs)
#' For each cutoff, true positive (TP, the number of P/LP variants above a tested cutoff),
#' false positive (FP, the number of BL-VUS/B/LB variants above a tested cutoff),
#' true negative (TN, the number of BL-VUS/B/LB variants below a tested cutoff),
#' and false negative (FN, the number of P/LP variants below a tested cutoff) were estimated.
#' Subsequently, LR+, overall accuracy, true positive rate (sensitivity), true negative rate (specificity), positive predictive value (PPV), negative predictive value (NPV), and F1 score were calculated.
#' Estimates of 95% CI of LR+ were generated using bootstrapping in the R package, bootLR.
#'
#'
#' @param data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param start The beginning column index of the evaluated cutoffs
#' @param end The concluding column index of evaluated cutoffs
#'
#' @import bootLR
#' @importFrom dplyr filter
#' @return A DataFrame comprising the evaluation metrics for each assessed cutoff
#' @export
#'
#' @examples
#' \dontrun{
#' data("VCI_data")
#' data <- add_info(VCI_data, "Assertion")
#' data <- VUS_classify(data, "Assertion", "Applied Evidence Codes (Met)")
#' library(dplyr)
#' truth_set <- filter(data,VUS_class %in% c("IceCold","Cold","Cool",""))
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PM2")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PP3")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PP1")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PVS1")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PS1")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PS2")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PM3")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PM4")
#' truth_set <- discrete_cutoff(truth_set, "Applied Evidence Codes (Met)", criteria = "PM5")
#' LR(truth_set, 28, 36)
#' }
#'
LR <- function(data, start, end) {
  output <- as.data.frame(matrix(nrow = 0, ncol = 21))
  for (i in start:end) {
    if (length(table(data$Classification_P, data[, i])) != 2) {
      if (table(data$Classification_P, data[, i])[2] != 0 && table(data$Classification_P, data[, i])[3] != 0 && table(data$Classification_P, data[, i])[4] != 0 && table(data$Classification_P, data[, i])[1] != 0) {
        a <- BayesianLR.test(table(data$Classification_P, data[, i])[4], table(data$Classification_P, data[, i])[4] + table(data$Classification_P, data[, i])[2], table(data$Classification_P, data[, i])[1], table(data$Classification_P, data[, i])[1] + table(data$Classification_P, data[, i])[3], R = 10000)
        TP <- table(data$Classification_P, data[, i])[4]
        FN <- table(data$Classification_P, data[, i])[2]
        FP <- table(data$Classification_P, data[, i])[3]
        TN <- table(data$Classification_P, data[, i])[1]
        total <- TP + FN + TN + FP
        # confu_final<-cbind(a$truePos,a$totalDzPos,a$trueNeg,a$totalDzNeg,a[1],a[2],a[3],a[4])
        Prevalence <- (TP + FN) / total
        Accuracy <- (TP + TN) / total
        PPV <- TP / (TP + FP)
        NPV <- TN / (FN + TN)
        F1_score <- 2 * (a$statistics[1] * PPV / (a$statistics[1] + PPV))
        FNR <- FN / (FN + TP)
        FPR <- FP / (FP + TN)
        FOR <- FN / (FN + TN)
        FDR <- FP / (FP + TP)
        confu_final <- cbind(colnames(data)[i], TP, FN, FP, TN, Accuracy, PPV, NPV, FNR, FPR, FOR, FDR, F1_score, a$statistics[1], a$statistics[2], a$posLR, a$posLR.ci[1], a$posLR.ci[2], a$negLR, a$negLR.ci[1], a$negLR.ci[2])
        output <- rbind(output, confu_final)
      }
    }
  }
  colnames(output) <- c("name", "TP", "FN", "FP", "TN", "Accuracy", "PPV", "NPV", "FNR", "FPR", "FOR", "FDR", "F1", "Sensitivity", "Specificity", "posLR", "posLR_LB", "posLR_UB", "negLR", "negLR_LB", "negLR_UB")
  return(output)
}
