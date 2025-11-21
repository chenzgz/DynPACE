#' PBC Example Dataset
#'
#' A processed subset of the PBC (Primary Biliary Cholangitis) data.
#'
#' @format A data frame with 1945 rows and 8 variables:
#' \describe{
#'   \item{id}{patients identifier; in total there are 312 patients.}
#'   \item{time}{number of years between registration and the earlier of death or study analysis time.}
#'   \item{rtime}{number of years between enrollment and this visit date.}
#'   \item{status}{Event status: 0 = alive, 1 = dead.}
#'   \item{age}{at registration in years.}
#'   \item{serBilir}{serum bilirubin in mg/dl.}
#'   \item{prothrombin}{prothrombin time in seconds.}
#'   \item{albumin}{albumin in gm/dl.}
#' }
#' @source Modified from the pbc2 dataset in the JMbayes package
#' @examples
#' data(pbc2_example)
#' head(pbc2_example)
#' summary(pbc2_example)
"pbc2_example"
