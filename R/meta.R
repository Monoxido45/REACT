#' Meta analysis data
#'
#' Data collect from the meta analysis originally conducted by da Silva Teixeira et al (2022).
#' The study evaluates patients' adherence to tobacco cessation protocols comparing
#' traditional approaches combined to computer assisted health technologies to traditional
#' approaches themselves. The treatments are compared by two main outcomes: the adherence
#' to the follow-up period of treatments without any drug ("Follow-up") and the adherence
#' to the pharmacotherapy, on studies that used any drug besides nicotine reposition ("Pharmacothetapy").
#'
#' @format Table with 9 columns and 11 rows:
#' \describe{
#'   \item{studlab}{Study label}
#'   \item{study_class}{Study outcome}
#'   \item{ne}{Number of patients under new treatment}
#'   \item{evente}{Numbers of patients adherent to the new treatment}
#'   \item{nc}{Number of patients under the traditional treatment (control)}
#'   \item{eventc}{Number of patients adherent to the traditional treatment}
#'   \item{te}{Relative risk between new and traditional treatment}
#'   \item{lowerte}{Lower bound of each study 95% relative risk-based confidence interval}
#'   \item{upperte}{Upper bound of each study 95% relative risk-based confidence interval}
#' }
#' @source
#' da Silva Teixeira, R., I. F. Nazareth, L. C. de Paula, G. P. do Nascimento Duque, and F. A. B. Colugnati,
#' 2022: Adherence to computational technologies for the treatment of smoking cessation: Systematic review and
#' meta-analysis. International Journal of Mental Health and Addiction
"meta"