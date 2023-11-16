#' Clinical variant classification in ClinGen Variant Curation Interface
#'
#' A dataset containing the curated 5724 variants by ClinGen.
#'
#' @format A data frame with 5724 rows and 20 variables:
#' \describe{
#'   \item{#Variation}{Variation, in HGVSc}
#'   \item{ClinVar Variation Id}{ClinVar Variation ID}
#'   \item{Allele Registry Id}{ClinGen Allele Registry ID}
#'   \item{HGVS Expressions}{HGVS Expressions in ClinVar}
#'   \item{HGNC Gene Symbol}{Gene Symbol}
#'   \item{Disease}{Variant related disease}
#'   \item{Mondo Id}{Mondo Disease Ontology ID}
#'   \item{Mode of Inheritance}{Genetic Inheritance pattern}
#'   \item{Assertion}{Variant Classification}
#'   \item{Applied Evidence Codes (Met)}{Criteria, represent following the SVI's recommendations}
#'   \item{Applied Evidence Codes (Not Met)}{Criteria not used}
#'   \item{Summary of interpretation}{Detailed information for each applied criteria}
#'   \item{PubMed Articles}{PubMed ID}
#'   \item{Expert Panel}{The name of variant curation expert panel}
#'   \item{Guideline}{Links of specific guidelines}
#'   \item{Approval Date}{Approval Date}
#'   \item{Published Date}{Published Date}
#'   \item{Retracted}{Retracted, in logical}
#'   \item{Evidence Repo Link}{Evidence Repo Link}
#'   \item{Uuid}{ID}
#'   ...
#' }
#' @source \url{https://erepo.clinicalgenome.org/evrepo/}
"VCI_data"


#' Dataset in the paper "Calibration of computational tools for missense variant pathogenicity classification and ClinGen recommendations for PP3/BP4 criteria".
#'
#' A dataset containing the ClinVar 2020 set to validate the calibration procedure proposed by Pejaver et al (2022).
#'
#' @format A data frame with 9114 rows and 29 variables:
#' \describe{
#'   \item{hg19_chr}{Chromosome}
#'   \item{hg19_pos.1.based.}{Position}
#'   \item{ref}{Reference allele}
#'   \item{alt}{Alternative allele}
#'   \item{rs_dbSNP151}{rsID}
#'   \item{genename}{Gene name}
#'   \item{Ensembl_geneid}{GeneID}
#'   \item{Ensembl_transcriptid}{TranscriptID}
#'   \item{Ensembl_proteinid}{ProteinID}
#'   \item{Uniprot_acc}{Uniprot Accession}
#'   \item{Uniprot_entry}{UniProt entry name}
#'   \item{aavar}{AA change}
#'   \item{clnsig}{ClinVar Significance}
#'   \item{MAF}{Minor allele frequency}
#'   \item{SIFT_score}{SIFT score}
#'   \item{FATHMM_score}{FATHMM score}
#'   \item{VEST4_score}{VEST4 score}
#'   \item{REVEL_score}{REVEL score}
#'   \item{GERP.._RS}{GERP++ score}
#'   \item{phyloP100way_vertebrate}{phyloP score}
#'   \item{EA_1.0}{EA score}
#'   \item{BayesDel_nsfp33a_noAF}{BayesDel score}
#'   \item{MutPred2.0_score}{MutPred score}
#'   \item{CADDv1.6_PHRED}{CADD score}
#'   \item{pph2_prob}{pph2 score}
#'   \item{MPC_score}{MPC score}
#'   \item{PrimateAI_score}{PrimateAI score}
#'   ...
#' }
#' @source \url{https://zenodo.org/records/8347415}
"ClinVar2020_AJHG_Pejaver_data"


#' LR results for VCI_data
#'
#'
#' @format A data frame with 8 rows and 20 variables:
#' \describe{
#'   \item{TP}{True positive}
#'   \item{FN}{False negative}
#'   \item{FP}{False positive}
#'   \item{TN}{True negative}
#'   \item{Accuracy}{(TP+TN)/Total}
#'   \item{PPV}{Positive predictive values}
#'   \item{NPV}{Negative predictive values}
#'   \item{FNR}{False negative rate}
#'   \item{FPR}{False positive rate}
#'   \item{FOR}{False omission rate}
#'   \item{FDR}{False discovery rate}
#'   \item{F1}{F1 score}
#'   \item{Sensitivity}{True positive rate}
#'   \item{Specificity}{True negative rate}
#'   \item{posLR}{Positive likelihood ratio}
#'   \item{posLR_LB}{The 95% CI lower boundry of posLR}
#'   \item{posLR_UB}{The 95% CI upper boundry of posLR}
#'   \item{negLR}{Negative likelihood ratio}
#'   \item{negLR_LB}{The 95% CI lower boundry of negLR}
#'   \item{negLR_UB}{The 95% CI upper boundry of negLR}
#'   ...
#' }
#' @source VCI_data
"LR_result"

#' locallr results for PrimateAI_score in ClinVar2020_AJHG_Pejaver_data
#'
#'
#' @format A data frame with 8586 rows and 3 variables:
#' \describe{
#'   \item{test_cutoff}{Each PrimateAI score}
#'   \item{Posterior}{Posterior probability}
#'   \item{Posterior1}{The 95% CI lower boundry of posterior probability}
#'   ...
#' }
#' @source ClinVar2020_AJHG_Pejaver_data
"lr_CI_result"
