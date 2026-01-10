#' Example GWAS Data
#'
#' A simulated dataset containing GWAS results for demonstration purposes.
#'
#' @format A data frame with 99226 rows and 4 variables:
#' \describe{
#'   \item{SNP}{SNP identifier}
#'   \item{CHR}{Chromosome number (1-25)}
#'   \item{BP}{Base pair position}
#'   \item{P}{P-value from association test}
#' }
#'
#' @source Simulated data for package demonstration.
#'
#' @examples
#' data(gwas_example)
#' head(gwas_example)
#' ggManhattan(gwas_example)
"gwas_example"
