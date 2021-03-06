#' Apgar score of 1755 babies
#'
#' A dataset containing the apgar score after 5 minutes, the condition of the babies, the birthweight of the babies, the gestational age and the smoking status of the mothers   
#'
#' @format  A data frame with 1755 rows and 7 variables:
#' \describe{
#'   \item{apgar5}{Apgar score after 5 minutes}
#'   \item{babycon}{Condition of the baby}
#'   \item{birthwt}{Birthweight of the babies}
#'   \item{gest}{gestational age}
#'   \item{smoke}{smoking status of the mother (0 non-smoker, 1 smoker)}
#'   \item{apgar_10}{10 minus apgar5, so the apgarscore can be seen as gamma distributed}
#'   \item{smoke2}{Allowing the smoking status to have three characteristics (0 non-smoker, 1 smoker, 2 no information)}
#'   \item{momid}{ID-number of the mother}
#' }
"bwsmokecompl"