#' Calculating the local positive likelihood ratio (lr+) value, which is applicable to continuous evidence proposed by Pejaver et al.
#' First, all unique tested cutoff values were sorted, then each value was positioned at the center of a sliding window.
#' The posterior probability was calculated for each tested cutoff value within the interval, considering a minimum of selected pathogenic and non-pathogenic variants.
#'
#'
#' @param input_data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param feature The column name that requires testing for optimizing the thresholds
#' @param alpha Prior probability
#' @param minpoints The number of at least pathogenic and non-pathogenic variants
#' @param increment Sliding window
#'
#' @return The posterior probability value for each tested cutoff
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar2020_AJHG_Pejaver_data")
#' data <- add_info(ClinVar2020_AJHG_Pejaver_data, "clnsig")
#' local_lr(data, "PrimateAI_score", 0.1, 100, 0.1)
#' }
#'
local_lr <- function(input_data, feature, alpha, minpoints, increment) {
  input_data <- input_data[!is.na(input_data[[feature]]), ]
  w <- (1 - alpha) * (sum(input_data$Classification_P == "P")) / (alpha * (sum(input_data$Classification_P == "NonP")))
  thrs <- sort(unique(c(input_data[[feature]], floor(min(input_data[[feature]])), ceiling(max(input_data[[feature]])))))
  post <- rep(0, times = length(thrs))
  finalwindow <- rep(0, times = length(thrs))
  maxthrs <- ceiling(max(input_data[[feature]]))
  minthrs <- floor(min(input_data[[feature]]))
  output <- as.data.frame(matrix(nrow = 0, ncol = 3))
  for (i in c(1:length(thrs))) {
    halfwindow <- 0
    lo <- thrs[i] - halfwindow
    hi <- thrs[i] + halfwindow
    while (1) {
      pos <- length(input_data[input_data[[feature]] > lo & input_data[[feature]] < hi & input_data$Classification_P == "P", ][[feature]])
      neg <- length(input_data[input_data[[feature]] > lo & input_data[[feature]] < hi & input_data$Classification_P == "NonP", ][[feature]])
      if (hi > maxthrs) {
        c <- (maxthrs - lo) / (hi - lo)
      } else if (lo < minthrs) {
        c <- (hi - minthrs) / (hi - lo)
      } else {
        c <- 1
      }
      if (c <= 0) {
        break
      }
      if (pos + neg < c * minpoints) {
        halfwindow <- halfwindow + increment
        lo <- thrs[i] - halfwindow
        hi <- thrs[i] + halfwindow
      }
      if (pos + neg >= c * minpoints) {
        break
      }
    }
    post[i] <- pos / (pos + w * neg)
    finalwindow[i] <- halfwindow
    temp <- cbind(thrs[i], post[i], finalwindow[i])
    output <- rbind(output, temp)
  }
  colnames(output) <- c("thrs", "post_p", "finalwindow")
  return(output)
}

#'  The one-sided 95% confidence bound for each estimated lr+  was determined through bootstrapping iterations, enabling the assessment of evidence strength.
#'
#' @param input_data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param feature The column name that requires testing for optimizing the thresholds
#' @param alpha Prior probability
#' @param bootstrap The number of bootstrapping iterations
#' @param minpoints The number of at least pathogenic and non-pathogenic variants
#' @param increment Sliding window
#' @param output_dir Output directory
#'
#'
#' @importFrom utils write.table
#' @return The posterior probability values for each bootstrap iteration
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar2020_AJHG_Pejaver_data")
#' data <- add_info(ClinVar2020_AJHG_Pejaver_data, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score", 0.1, 10, 100, 0.1, "test_dir")
#' }
#'
local_bootstrapped_lr <- function(input_data, feature, alpha, bootstrap, minpoints, increment, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  previousWorkPath <- getwd()
  WorkPath <- setwd(output_dir)
  lr_result <- local_lr(input_data, feature, alpha, minpoints, increment)
  write.table(lr_result, file = "local_lr.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  for (i in c(1:bootstrap)) {
    idx <- sample(1:nrow(input_data), replace = T)
    sample_data <- input_data[idx, ]
    bootstrap_i <- local_lr(sample_data, feature, alpha, minpoints, increment)
    write.table(bootstrap_i, file = paste("bootstrap_", i, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
  }
  setwd(previousWorkPath)
}

#' Establish the thresholds for each level of evidence strength
#'
#' @param postp_list A list of posterior probability corresponding to each level of evidence strength
#' @param discountonesided The one-sided confidence intervals
#' @param bootstrap The number of bootstrapping iterations
#' @param dir The directory containing the results of bootstrapping
#'
#' @return A list of optimized thresholds
#' @importFrom utils read.table
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar2020_AJHG_Pejaver_data")
#' data <- add_info(ClinVar2020_AJHG_Pejaver_data, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score", 0.1, 10, 100, 0.1, "test_dir")
#' postp_list <- c(0.1778350, 0.3129676, 0.6689245, 0.9754584)
#' get_lr_threshold(postp_list, 0.05, 10, "test_dir"
#' }
#'
get_lr_threshold <- function(postp_list, discountonesided, bootstrap, dir) {
  thresh <- as.data.frame(matrix(nrow = bootstrap, ncol = 4))
  DiscountedThreshold <- c(0, 0, 0, 0)
  if (!dir.exists(dir)) {
    return(message("Error: The input directory does not exist!"))
  } else {
    previousWorkPath <- getwd()
    WorkPath <- setwd(dir)
    for (i in c(1:bootstrap)) {
      input <- read.table(paste("bootstrap_", i, ".txt", sep = ""), sep = "\t", header = T)
      for (j in c(1:length(postp_list))) {
        ind <- min(which(input$post_p >= postp_list[j]))
        if (ind == 1) {
          thresh[i, j] <- input$thrs[2]
        } else if (ind > 1) {
          thresh[i, j] <- input$thrs[ind]
        } else {
          (thresh[i, j] <- NA)
        }
      }
    }
    for (j in c(1:length(postp_list))) {
      invalids <- sum(is.na(thresh[, j]))
      if (invalids > (discountonesided * bootstrap)) {
        DiscountedThreshold[j] <- NA
      } else {
        t <- sort(thresh[which(!is.na(thresh[, j])), j], decreasing = TRUE)
        DiscountedThreshold[j] <- t[floor(discountonesided * (bootstrap - invalids)) + 1]
      }
    }
    setwd(previousWorkPath)
    return(DiscountedThreshold)
  }
}

#' Merging the results from bootstrap
#'
#' @param bootstrap The number of bootstrapping iterations
#' @param dir The directory containing the results of bootstrapping
#' @importFrom stats t.test
#' @importFrom utils read.table
#' @importFrom dplyr left_join
#'
#'
#' @return A DataFrame containing posterior probabilities and the 95% confidence interval lower bounds of posterior probabilities for each cutoff
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar2020_AJHG_Pejaver_data")
#' data <- add_info(ClinVar2020_AJHG_Pejaver_data, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score", 0.1, 30, 100, 0.1, "test_dir")
#' lr_CI_result <- lr_CI(30, "test_dir"
#' }
#'
lr_CI <- function(bootstrap, dir) {
  if (!dir.exists(dir)) {
    return(message("Error: The input directory does not exist!"))
  } else {
    previousWorkPath <- getwd()
    WorkPath <- setwd(dir)
    postp <- read.table("local_lr.txt", sep = "\t", header = T)
    postp_matrix <- postp[, c(1:2)]
    for (i in 1:bootstrap) {
      input <- read.table(paste("bootstrap_", i, ".txt", sep = ""), sep = "\t", header = T)
      postp_matrix <- left_join(postp_matrix, input[, c(1:2)], by = "thrs")
      colnames(postp_matrix)[i + 2] <- paste("postp_", i, sep = "")
    }
    output <- as.data.frame(matrix(nrow = nrow(postp_matrix), ncol = 3))
    colnames(output) <- c("test_cutoff", "Posterior", "Posterior1")
    output[, 1] <- postp_matrix[, 1]
    output[, 2] <- postp_matrix[, 2]
    for (j in 1:nrow(postp_matrix)) {
      a <- t.test(postp_matrix[j, c(3:ncol(postp_matrix))])
      output[j, 3] <- a$conf.int[1]
    }
    setwd(previousWorkPath)
    return(output)
  }
}
