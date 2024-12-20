#' Calculating the local positive likelihood ratio (lr+) value, which is applicable to continuous evidence proposed by Pejaver et al.
#' First, all unique tested cutoff values were sorted, then each value was positioned at the center of a sliding window.
#' The posterior probability was calculated for each tested cutoff value within the interval, considering a minimum of selected pathogenic and non-pathogenic variants.
#'
#'
#' @param input_data DataFrame comprising fundamental variant information, evidence labeling, and classification details
#' @param feature The column name that requires testing for optimizing the thresholds
#' @param direction The direction of evidence pathogenic (Pathogenic or Benign)
#' @param alpha Prior probability
#' @param minpoints The number of at least pathogenic and non-pathogenic variants
#' @param increment Sliding window
#'
#' @return The posterior probability value for each tested cutoff
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar_2019_dataset")
#' data <- add_info(ClinVar_2019_dataset, "clnsig")
#' local_lr(data, "PrimateAI_score", "Pathogenic",0.0441, 100, 0.01)
#' }
#'
local_lr <- function(input_data, feature, direction,alpha, minpoints, increment) {
  if(direction!="Pathogenic" && direction!="Benign"){
    return(message("Error,the direction of evidence pathogenic must be Pathogenic or Benign"))
  }
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
    if(direction=="Pathogenic"){
      post[i] <- pos / (pos + w * neg)
    }
    else if(direction=="Benign"){
      post[i] <- (w * neg) / (pos + w * neg)
    }
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
#' @param direction The direction of evidence pathogenic (Pathogenic or Benign)
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
#' data("ClinVar_2019_dataset")
#' data <- add_info(ClinVar_2019_dataset, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score","Pathogenic", 0.0441, 10000, 100, 0.01, "test_dir")
#' }
#'
local_bootstrapped_lr <- function(input_data, feature, direction,alpha, bootstrap, minpoints, increment, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  if(direction!="Pathogenic" && direction!="Benign"){
    return(message("Error,the direction of evidence pathogenic must be Pathogenic or Benign"))
  }
  previousWorkPath <- getwd()
  WorkPath <- setwd(output_dir)
  lr_result <- local_lr(input_data, feature, direction,alpha, minpoints, increment)
  write.table(lr_result, file = "local_lr.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  for (i in c(1:bootstrap)) {
    idx <- sample(1:nrow(input_data), replace = T)
    sample_data <- input_data[idx, ]
    bootstrap_i <- local_lr(sample_data, feature, direction,alpha, minpoints, increment)
    write.table(bootstrap_i, file = paste("bootstrap_", i, ".txt", sep = ""), sep = "\t", col.names = T, row.names = F, quote = F)
  }
  setwd(previousWorkPath)
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
#' data("ClinVar_2019_dataset")
#' data <- add_info(ClinVar_2019_dataset, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score", "Pathogenic",0.0441, 10000, 100, 0.01, "test_dir")
#' lr_CI_result <- lr_CI(10000, "test_dir")
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
      if(all(na.omit(postp_matrix[j, c(3:ncol(postp_matrix))])) == postp_matrix[j, 2]){output[j, 3] <- postp_matrix[j, 2]}
      else{
        a <- t.test(postp_matrix[j, c(3:ncol(postp_matrix))])
        output[j, 3] <- a$conf.int[1]
      }
    }
    setwd(previousWorkPath)
    return(output)
  }
}

#' Establish the thresholds for each level of evidence strength
#'
#' @param data The results of bootstrapping
#' @param postp_list A list of posterior probability corresponding to each level of evidence strength
#' @param direction The direction of evidence pathogenic (Pathogenic or Benign)
#'
#' @return A list of optimized thresholds
#' @importFrom utils read.table
#' @importFrom stats na.omit
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("ClinVar_2019_dataset")
#' data <- add_info(ClinVar_2019_dataset, "clnsig")
#' local_bootstrapped_lr(data, "PrimateAI_score", "Pathogenic",0.0441, 10000, 100, 0.01, "test_dir")
#' postp_list <- c(0.100, 0.211, 0.608, 0.981)
#' lr_CI_result <- lr_CI(10000, "test_dir")
#' get_lr_threshold(lr_CI_result, postp_list, "Pathogenic")
#' }
#'
get_lr_threshold <- function(data, postp_list, direction) {
  if(direction!="Pathogenic" && direction!="Benign"){
    return(message("Error,the direction of evidence pathogenic must be Pathogenic or Benign"))
  }
  thresh <- c(0, 0, 0, 0)
  if (!exists("data")) {
    return(message("Error: The input data does not exist!"))
  } else {
    data<-na.omit(data)
    for (j in c(1:length(postp_list))) {
      valid_thrs <- data$test_cutoff[data$Posterior1 > postp_list[j]]
      thresh[j] <- NA
      while (length(valid_thrs) > 0) {
        if(direction=="Pathogenic"){
          candidate_min_thrs <- min(valid_thrs)
          if (all(data$Posterior1[data$test_cutoff >= candidate_min_thrs] > postp_list[j])) {
            thresh[j] <- candidate_min_thrs
            break
          } else{
            valid_thrs <- valid_thrs[valid_thrs != candidate_min_thrs]
          }
        }
        else if(direction=="Benign"){
          candidate_max_thrs <- max(valid_thrs)
          if (all(data$Posterior1[data$test_cutoff <= candidate_max_thrs] > postp_list[j])) {
            thresh[j] <- candidate_max_thrs
            if(candidate_max_thrs==0){thresh[j] <- NA}
            break
          } else{
            valid_thrs <- valid_thrs[valid_thrs != candidate_max_thrs]
          }
        }
      }
    }
    return(thresh)
  }
}
