#' Meta analysis data
#'
#' Data from the meta analysis conducted by da Silva Teixeira et al (2022).
#' The study evaluates patients' adherence to tobacco cessation protocols
#' from two main outcomes: adherence to the follow-up period of treatments
#' without any drug (“Follow-up”) and adherence to the pharmacotherapy on
#' studies that used any drug besides nicotine replacement (“Pharmacothetapy”).
#'
#' @format Table with 9 columns and 11 rows:
#' \describe{
#'   \item{studlab}{Study label}
#'   \item{study_class}{Study outcome class}
#'   \item{ne}{Number of patients under new treatment}
#'   \item{evente}{Numbers of patients adherent to the new treatment}
#'   \item{nc}{Number of patients under the control group}
#'   \item{eventc}{Number of patients adherent in the control group}
#'   \item{te}{Relative risk between new and traditional treatment}
#'   \item{lowerte}{Lower bound of the 95% confidence interval for the study's relative risk}
#'   \item{upperte}{Upper bound of the 95% confidence interval for the study's relative risk}
#' }
#' @source
#' R. da Silva Teixeira, I.F. Nazareth, L.C. de Paula, G.P. do Nascimento Duque,
#' and F.A.B. Colugnati, 2022: Adherence to computational technologies for the
#' treatment of smoking cessation: Systematic review and meta-analysis.
#' International Journal of Mental Health and Addiction
"meta"
