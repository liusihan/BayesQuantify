#' Define the global Variables
#'
#' @param libname lib name
#' @param pkgname package name
#'
#' @return global variables
#' @export
#'
#' @examples
#' #null
#'
.onLoad <- function (libname, pkgname)
{
  # set global variables in order to avoid CHECK notes
  utils::globalVariables (c("Gene","sum.n.","Classification","Classification_3","Consequence","n.x","n.y","percent.x","percent.y","test_cutoff","Posterior","Posterior1"))
  invisible ()
}
