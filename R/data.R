
# DOCUMENTATION OF DATASETS


# Dataset bookshop_sales -------------------------------------------------------
#' @title Antiquarian Bookshop Sales
#'
#' @description
#' Individual orders of a Czech antiquarian bookshop from June 8, 2018 to December 20, 2018.
#' This dataset is analyzed in Tomanová and Holý (2021).
#' Details on the bookshop can be also found in Tomanová and Černý (2022).
#'
#' @source
#' Petra Tomanová (\email{petra.tomanova@@vse.cz}).
#'
#' @format
#' A data frame with columns:
#' \describe{
#'   \item{order}{ID of the order.}
#'   \item{time}{Time of the order.}
#'   \item{quantity}{Number of purchased books. Zero value means the order was canceled.}
#' }
#'
#' @references
#' Tomanová, P. and Černý, M. (2022). Efficiency of Antiquarian Bookshops in Informationally Complete Markets. \emph{Central European Journal of Operations Research}, \strong{30}(2), 573–593. \doi{10.1007/s10100-021-00780-3}.
#'
#' Tomanová, P. and Holý, V. (2021). Clustering of Arrivals in Queueing Systems: Autoregressive Conditional Duration Approach. \emph{Central European Journal of Operations Research}, \strong{29}(3), 859–874. \doi{10.1007/s10100-021-00744-7}.
#'
"bookshop_sales"
# ------------------------------------------------------------------------------


# Dataset german_car_market_cap ------------------------------------------------
#' @title Market Capitalization of German Car Manufacturers
#'
#' @description
#' Market capitalization of the "Germany's Big Three" automobile manufacturers – Volkswagen Group, Mercedes-Benz Group, and BMW.
#' Market capitalization is reported in billions of euros and covers the period 1994–2021.
#'
#' @source
#' Thomson Reuters (\href{https://www.thomsonreuters.com/en.html}{www.thomsonreuters.com}).
#'
#' @format
#' A data frame with columns:
#' \describe{
#'   \item{year}{Year.}
#'   \item{car_manufacturer}{Car manufacturer.}
#'   \item{market_cap}{Market capitalization in billions of euros.}
#' }
#'
"german_car_market_cap"
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


# Dataset sp500_daily ----------------------------------------------------------
#' @title Daily S&P 500 Prices
#'
#' @description
#' Daily opening, highest, lowest, and closing prices of the Standard and Poor's 500 stock market index (SPX) from 2013.
#'
#' @source
#' Nasdaq (\href{https://www.nasdaq.com/market-activity/index/spx}{www.nasdaq.com/market-activity/index/spx}).
#'
#' @format
#' A data frame with columns:
#' \describe{
#'   \item{date}{Trading day.}
#'   \item{open}{Opening price of the day.}
#'   \item{high}{Highest price of the day.}
#'   \item{low}{Lowest price of the day.}
#'   \item{close}{Closing price of the day.}
#' }
#'
"sp500_daily"
# ------------------------------------------------------------------------------


