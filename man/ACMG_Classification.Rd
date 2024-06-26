% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_classification.R
\name{ACMG_Classification}
\alias{ACMG_Classification}
\title{Classifying variants into five distinct categories according to the 2015 ACMG/AMP guidelines}
\usage{
ACMG_Classification(data, evidence_col)
}
\arguments{
\item{data}{DataFrame comprising fundamental variant information, evidence labeling, and classification details}

\item{evidence_col}{The column name for ACMG evidence(str)}
}
\value{
A new DataFrame that incorporates the input data and the results of variant classification
}
\description{
Classifying variants into five distinct categories according to the 2015 ACMG/AMP guidelines
}
\examples{
\dontrun{
data("ClinGen_dataset")
ACMG_Classification(ClinGen_dataset, "Applied Evidence Codes (Met)")
}

}
