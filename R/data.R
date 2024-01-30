
# DOCUMENTATION OF DATASETS


# Dataset bookshop_orders -------------------------------------------------------
#' @title Antiquarian Bookshop Orders
#'
#' @description
#' Individual orders of a Czech antiquarian bookshop from June 8, 2018 to December 20, 2018.
#' This dataset is analyzed in Tomanová and Holý (2021).
#'
#' @format
#' A data frame with columns:
#' \describe{
#'   \item{id}{ID of the order.}
#'   \item{datetime}{Date and time of the order.}
#'   \item{quantity}{Number of purchased books.}
#'   \item{duration}{Number of minutes since the last order.}
#'   \item{duration_adj}{Duration since the last order adjusted for diurnal pattern.}
#' }
#'
#' @references
#' Tomanová, P. and Holý, V. (2021). Clustering of Arrivals in Queueing Systems: Autoregressive Conditional Duration Approach. \emph{Central European Journal of Operations Research}, \strong{29}(3), 859–874. \doi{10.1007/s10100-021-00744-7}.
#'
"bookshop_orders"
# ------------------------------------------------------------------------------


# Dataset ice_hockey_championships ---------------------------------------------
#' @title Results of the Ice Hockey World Championships
#'
#' @description
#' The dataset contains the results of the annual men's Ice Hockey World Championships from 1998 to 2023.
#' In 1998, the International Ice Hockey Federation set the number of teams participating in the championships at 16.
#' Since 1998, a total of 24 teams have qualified for the championship division.
#' This dataset is analyzed in Holý and Zouhar (2022).
#'
#' @source
#' International Ice Hockey Federation (\href{https://www.iihf.com/}{www.iihf.com}).
#'
#' @format
#' A list with components:
#' \describe{
#'   \item{rankings}{A matrix of final rankings. Rows correspond to years, columns to teams. Value \code{Inf} means that the team did not advance to the championship. Value \code{NA} means that the championship did not take place.}
#'   \item{hosts}{A matrix of dummy variables indicating whether the team hosted the championship. Rows correspond to years, columns to teams. Multiple hosts of one championship is possible. Value \code{NA} means that the championship did not take place.}
#' }
#'
#' @references
#' Holý, V. and Zouhar, J. (2022). Modelling Time-Varying Rankings with Autoregressive and Score-Driven Dynamics. Journal of the Royal Statistical Society: Series C (Applied Statistics), \strong{71}(5). \doi{10.1111/rssc.12584}.
#'
"ice_hockey_championships"
# ------------------------------------------------------------------------------


# Dataset bookshop_orders -------------------------------------------------------
#' @title Daily Toilet Paper Sales
#'
#' @description
#' The daily number of toilet paper packs sold in a European store in 2001 and 2002.
#' The \code{promo} variable indicates whether the product was promoted in a campaign.
#' Missing values correspond to the days when the store was closed.
#'
#' @format
#' A data frame with columns:
#' \describe{
#'   \item{date}{Date.}
#'   \item{monday}{Dummy variable indicating whether it is Monday.}
#'   \item{tuesday}{Dummy variable indicating whether it is Tuesday.}
#'   \item{wednesday}{Dummy variable indicating whether it is Wednesday.}
#'   \item{thursday}{Dummy variable indicating whether it is Thursday.}
#'   \item{friday}{Dummy variable indicating whether it is Friday.}
#'   \item{saturday}{Dummy variable indicating whether it is Saturday.}
#'   \item{sunday}{Dummy variable indicating whether it is Sunday.}
#'   \item{promo}{Dummy variable indicating whether the product is promoted.}
#'   \item{quantity}{Number of packs sold.}
#' }
#'
"toilet_paper_sales"
# ------------------------------------------------------------------------------


