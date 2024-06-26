% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_op.R
\name{op_postp}
\alias{op_postp}
\title{Calculate the corresponding combined odds_path and posterior probability of
17 combination rules for a given prior_probability and odds_path of pathogenicity}
\usage{
op_postp(prior_probability, op_vs)
}
\arguments{
\item{prior_probability}{The prior probability of pathogenicity (proportion of P/LP variants in a set of variants)}

\item{op_vs}{Odds of pathogenicity (OP) of "Very String"}
}
\value{
Prior_probability, OP for each evidence level and Combined odds_path and posterior probability of 17 combination rules outlined by avtigian et al.(2018)
}
\description{
Calculate the corresponding combined odds_path and posterior probability of
17 combination rules for a given prior_probability and odds_path of pathogenicity
}
\examples{
op_postp(0.1, 350)

}
