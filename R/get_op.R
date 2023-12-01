#' Calculate the corresponding combined odds_path and posterior probability of
#' 17 combination rules for a given prior_probability and odds_path of pathogenicity
#'
#' @param prior_probability The prior probability of pathogenicity (proportion of P/LP variants in a set of variants)
#' @param op_vs Odds path of "Very String"
#'
#' @return Combined odds_path and posterior probability of 17 combination rules outlined by avtigian et al.(2018)
#' @export
#'
#' @examples
#' op_postp(0.1, 350)
#'
op_postp <- function(prior_probability, op_vs) {
  op_s <- op_vs^(2^-1)
  op_m <- op_vs^(2^-2)
  op_p <- op_vs^(2^-3)
  postp_vs <- op_vs * prior_probability / ((op_vs - 1) * prior_probability + 1)
  postp_s <- op_s * prior_probability / ((op_s - 1) * prior_probability + 1)
  postp_m <- op_m * prior_probability / ((op_m - 1) * prior_probability + 1)
  postp_p <- op_p * prior_probability / ((op_p - 1) * prior_probability + 1)
  print(c(op_vs, op_s, op_m, op_p, postp_vs, postp_s, postp_m, postp_p))
  output <- matrix(nrow = 17, ncol = 2)
  colnames(output) <- c("Combined_Odds_Path", "Post_P")
  rownames(output) <- c("LP(i)", "LP(ii)", "LP(iii)", "LP(iv)", "LP(v)", "LP(vi)", "P(ia)", "P(ib)", "P(ic)", "P(id)", "P(ii)", "P(iiia)", "P(iiib)", "P(iiic)", "LB(i)", "LB(ii)", "B(i)")
  Combined_Odds_Path <- op_vs * op_m
  output[1, 1] <- Combined_Odds_Path
  output[1, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_m
  output[2, 1] <- Combined_Odds_Path
  output[2, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_p * op_p
  output[3, 1] <- Combined_Odds_Path
  output[3, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_m * op_m * op_m
  output[4, 1] <- Combined_Odds_Path
  output[4, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_m * op_m * op_p * op_p
  output[5, 1] <- Combined_Odds_Path
  output[5, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_m * op_p * op_p * op_p * op_p
  output[6, 1] <- Combined_Odds_Path
  output[6, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_vs * op_s
  output[7, 1] <- Combined_Odds_Path
  output[7, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_vs * op_m * op_m
  output[8, 1] <- Combined_Odds_Path
  output[8, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_vs * op_m * op_p
  output[9, 1] <- Combined_Odds_Path
  output[9, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_vs * op_p * op_p
  output[10, 1] <- Combined_Odds_Path
  output[10, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_s
  output[11, 1] <- Combined_Odds_Path
  output[11, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_m * op_m * op_m
  output[12, 1] <- Combined_Odds_Path
  output[12, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_m * op_m * op_p * op_p
  output[13, 1] <- Combined_Odds_Path
  output[13, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- op_s * op_m * op_p * op_p * op_p * op_p
  output[14, 1] <- Combined_Odds_Path
  output[14, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- (op_s^-1) * (op_p^-1)
  output[15, 1] <- Combined_Odds_Path
  output[15, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- (op_p^-2)
  output[16, 1] <- Combined_Odds_Path
  output[16, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  Combined_Odds_Path <- (op_s^-2)
  output[17, 1] <- Combined_Odds_Path
  output[17, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
  print(output)
}

#' Automatic definition of posterior probability and positive likelihood ratio values
#' for different strengths of evidence
#'
#' @param prior_probability The prior probability of pathogenicity (proportion of P/LP variants in a set of variants)
#'
#' @return Selected odds path for each evidence level, combined odds path and posterior probability of 17 combination rules
#' @export
#'
#' @examples
#' auto_select_postp(0.1)
#'
auto_select_postp <- function(prior_probability) {
  select_op <- 0
  select_matrix <- NULL
  for (op_vs in 1:10000) {
    op_s <- op_vs^(2^-1)
    op_m <- op_vs^(2^-2)
    op_p <- op_vs^(2^-3)
    postp_vs <- op_vs * prior_probability / ((op_vs - 1) * prior_probability + 1)
    postp_s <- op_s * prior_probability / ((op_s - 1) * prior_probability + 1)
    postp_m <- op_m * prior_probability / ((op_m - 1) * prior_probability + 1)
    postp_p <- op_p * prior_probability / ((op_p - 1) * prior_probability + 1)
    output <- matrix(nrow = 17, ncol = 2)
    colnames(output) <- c("Combined_Odds_Path", "Post_P")
    rownames(output) <- c("LP(i)", "LP(ii)", "LP(iii)", "LP(iv)", "LP(v)", "LP(vi)", "P(ia)", "P(ib)", "P(ic)", "P(id)", "P(ii)", "P(iiia)", "P(iiib)", "P(iiic)", "LB(i)", "LB(ii)", "B(i)")
    Combined_Odds_Path <- op_vs * op_m
    output[1, 1] <- Combined_Odds_Path
    output[1, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_m
    output[2, 1] <- Combined_Odds_Path
    output[2, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_p * op_p
    output[3, 1] <- Combined_Odds_Path
    output[3, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_m * op_m * op_m
    output[4, 1] <- Combined_Odds_Path
    output[4, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_m * op_m * op_p * op_p
    output[5, 1] <- Combined_Odds_Path
    output[5, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_m * op_p * op_p * op_p * op_p
    output[6, 1] <- Combined_Odds_Path
    output[6, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_vs * op_s
    output[7, 1] <- Combined_Odds_Path
    output[7, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_vs * op_m * op_m
    output[8, 1] <- Combined_Odds_Path
    output[8, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_vs * op_m * op_p
    output[9, 1] <- Combined_Odds_Path
    output[9, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_vs * op_p * op_p
    output[10, 1] <- Combined_Odds_Path
    output[10, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_s
    output[11, 1] <- Combined_Odds_Path
    output[11, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_m * op_m * op_m
    output[12, 1] <- Combined_Odds_Path
    output[12, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_m * op_m * op_p * op_p
    output[13, 1] <- Combined_Odds_Path
    output[13, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- op_s * op_m * op_p * op_p * op_p * op_p
    output[14, 1] <- Combined_Odds_Path
    output[14, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- (op_s^-1) * (op_p^-1)
    output[15, 1] <- Combined_Odds_Path
    output[15, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- (op_p^-2)
    output[16, 1] <- Combined_Odds_Path
    output[16, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    Combined_Odds_Path <- (op_s^-2)
    output[17, 1] <- Combined_Odds_Path
    output[17, 2] <- (Combined_Odds_Path * prior_probability) / ((Combined_Odds_Path - 1) * prior_probability + 1)
    if (output[6, 2] > 0.9 && output[14, 2] > 0.99 && output[16, 2] < 0.1 && output[17, 2] < 0.001) {
      print(output)
      print(c(op_vs, op_s, op_m, op_p, postp_vs, postp_s, postp_m, postp_p))
      select_op <- 1
      break
    } else {
      n <- 0
      for (i in 1:6) {
        if (output[i, 2] > 0.9 && output[i, 2] <= 0.99) {
          n <- n + 1
        }
      }
      for (i in 7:14) {
        if (output[i, 2] > 0.99) {
          n <- n + 1
        }
      }
      for (i in 15:16) {
        if (output[i, 2] >= 0.001 && output[i, 2] < 0.1) {
          n <- n + 1
        }
      }
      if (output[17, 2] < 0.001) {
        n <- n + 1
      }
      select_matrix[op_vs] <- n
      post_p_matrix <- paste("post_p", op_vs, sep = "_")
      assign(post_p_matrix, output)
      next
    }
  }
  if (select_op == 0) {
    select_op <- min(which(select_matrix >= max(select_matrix)))
    print(eval(parse(text = paste("post_p", select_op, sep = "_"))))
    op_vs <- select_op
    op_s <- op_vs^(2^-1)
    op_m <- op_vs^(2^-2)
    op_p <- op_vs^(2^-3)
    postp_vs <- op_vs * prior_probability / ((op_vs - 1) * prior_probability + 1)
    postp_s <- op_s * prior_probability / ((op_s - 1) * prior_probability + 1)
    postp_m <- op_m * prior_probability / ((op_m - 1) * prior_probability + 1)
    postp_p <- op_p * prior_probability / ((op_p - 1) * prior_probability + 1)
    print(c(op_vs, op_s, op_m, op_p, postp_vs, postp_s, postp_m, postp_p))
  }
}
