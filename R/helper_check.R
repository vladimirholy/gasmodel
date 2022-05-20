
# VARIABLE CHECKS FUNCTIONS


# Check Method -----------------------------------------------------------------
check_my_method <- function(method, values) {
  if (is.vector(method) && !is.list(method) && length(method) == 1L && !is.na(method)) {
    method <- method
  } else {
    stop("Invalid type of argument method.")
  }
  if (!is.null(values) && !(method %in% values)) {
    stop("Invalid method.")
  }
  return(method)
}
# ------------------------------------------------------------------------------


# Check Time Series ------------------------------------------------------------
check_my_y <- function(y = NULL, t = NULL, n = NULL, dim = NULL, type = NULL) {
  if (!is.null(type) && type == "categorical" && is.null(n) && is.vector(y) && !is.list(y) && length(y) > 0L && all(as.numeric(y) == as.integer(y), na.rm = TRUE) && min(as.integer(y), na.rm = TRUE) >= 1L && max(as.integer(y), na.rm = TRUE) >= 2L) {
    y <- t(sapply(y, function(i) { yy <- rep(ifelse(is.na(i), NA, 0), max(as.integer(y), na.rm = TRUE)); yy[i] <- 1; yy }))
  } else if (is.vector(y) && !is.list(y) && length(y) > 0L && !is.null(dim) && dim == "uni") {
    y <- matrix(as.numeric(y), ncol = 1)
  } else if (is.vector(y) && !is.list(y) && length(y) > 0L && !is.null(dim) && dim == "multi") {
    y <- matrix(as.numeric(y), nrow = 1)
  } else if (is.matrix(y) && nrow(y) > 0L && ncol(y) > 0L) {
    y <- matrix(as.numeric(y), nrow = nrow(y), ncol = ncol(y))
  } else {
    stop("Invalid type of argument y.")
  }
  y[rowSums(is.na(y)) > 0L, ] <- NA_real_
  if (!is.null(type) && type == "binary" && all(y == 0L | y == 1L, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "categorical" && all(y == 0L | y == 1L, na.rm = TRUE) && all(rowSums(y) == 1L, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "ranking" && all(y > 0, na.rm = TRUE) && all(apply(y, 1, function(e) { all(is.na(e)) || all(1:sum(is.finite(e)) %in% e) }))) {
    y <- y
  } else if (!is.null(type) && type == "count" && all(!is.infinite(y), na.rm = TRUE) && all(y >= 0, na.rm = TRUE) && all(as.vector(y) == as.integer(y), na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "integer" && all(!is.infinite(y), na.rm = TRUE) && all(as.vector(y) == as.integer(y), na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "directional" && all(!is.infinite(y), na.rm = TRUE) && all(y >= 0 | y <= 2 * pi, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "interval" && all(!is.infinite(y), na.rm = TRUE) && all(y >= 0 | y <= 1, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "compositional" && all(!is.infinite(y), na.rm = TRUE) && all(y >= 0, na.rm = TRUE) && all(abs(rowSums(y) - 1) < 1e-6, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "duration" && all(!is.infinite(y), na.rm = TRUE) && all(y >= 0, na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type) && type == "real" && all(!is.infinite(y), na.rm = TRUE)) {
    y <- y
  } else if (!is.null(type)) {
    stop("Invalid values of the time series y.")
  }
  if (!is.null(t) && nrow(y) != t) {
    stop("Invalid length of the time series y.")
  }
  if (!is.null(n) && ncol(y) != n) {
    stop("Invalid dimension of the time series y.")
  }
  if (!is.null(dim) && ((dim == "uni" && ncol(y) != 1L) || (dim == "multi" && ncol(y) <= 1L))) {
    stop("Invalid dimension of the time series y.")
  }
  if (all(rowSums(is.na(y)) > 0L)) {
    stop("All values of the time series y are missing.")
  }
  return(y)
}
# ------------------------------------------------------------------------------


# Check Exogeneous Variables ---------------------------------------------------
check_my_x <- function(x = NULL, t = NULL, m = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(x) && !is.null(t) && !is.null(par_num)) {
    x <- rep(list(matrix(NA_real_, nrow = t, ncol = 0L)), times = par_num)
  } else if (is.vector(x) && !is.list(x) && length(x) == 0L && !is.null(t) && !is.null(par_num)) {
    x <- rep(list(matrix(NA_real_, nrow = t, ncol = 0L)), times = par_num)
  } else if (is.vector(x) && !is.list(x) && length(x) > 0L && !is.null(par_num)) {
    x <- rep(list(matrix(x, ncol = 1L)), times = par_num)
  } else if (is.matrix(x) && !is.null(par_num)) {
    x <- rep(list(x), times = par_num)
  } else if (is.list(x) && !is.null(group_num) && !is.null(par_in_group_num) && length(x) == group_num && all(sapply(x, FUN = is.vector) | sapply(x, FUN = is.matrix)) && all(sapply(x, FUN = length) > 0L)) {
    x <- sapply(1:group_num, FUN = function(i) { rep(list(as.matrix(x[[i]])), times = par_in_group_num[i]) })
  } else if (is.list(x) && !is.null(group_num) && !is.null(par_in_group_num) && length(x) == group_num && all(sapply(x, FUN = is.vector) | sapply(x, FUN = is.matrix)) && !is.null(t)) {
    x <- sapply(1:group_num, FUN = function(i) { if (length(x[[i]]) == 0L) { rep(list(matrix(NA_real_, nrow = t, ncol = 0L)), times = par_in_group_num[i]) } else { rep(list(as.matrix(x[[i]])), times = par_in_group_num[i]) } })
  } else if (is.list(x) && length(x) > 0L && all(sapply(x, FUN = is.vector) | sapply(x, FUN = is.matrix)) && all(sapply(x, FUN = length) > 0L)) {
    x <- lapply(x, FUN = function(e) { as.matrix(e) })
  } else if (is.list(x) && length(x) > 0L && all(sapply(x, FUN = is.vector) | sapply(x, FUN = is.matrix)) && !is.null(t)) {
    x <- lapply(x, FUN = function(e) { if (length(e) == 0L) { matrix(NA_real_, nrow = t, ncol = 0L) } else { as.matrix(e) } })
  } else {
    stop("Invalid type of argument x.")
  }
  if (length(unique(sapply(x, FUN = nrow))) > 1L) {
    stop("Different lengths of exogeneous variables x.")
  }
  if (!is.null(t) && any(sapply(x, FUN = nrow) != t)) {
    stop("Invalid length of exogeneous variables x.")
  }
  if (!is.null(par_num) && length(x) != par_num) {
    stop("Ivalid number of exogeneous variables x.")
  }
  if (!is.null(m) && length(x) != length(m)) {
    stop("Ivalid number of exogeneous variables x.")
  }
  if (!is.null(m) && any(sapply(x, FUN = ncol) != m)) {
    stop("Ivalid number of exogeneous variables x.")
  }
  if (any(!sapply(x, is.numeric))) {
    stop("Exogeneous variables x must be numeric.")
  }
  if (nrow(x[[1]]) > 0L && all(rowSums(as.matrix(sapply(x, function(e) { rowSums(is.na(e)) }))) > 0L)) {
    stop("All periods have at least one exogeneous variable with missing value.")
  }
  return(x)
}
# ------------------------------------------------------------------------------


# Check Parameters -------------------------------------------------------------
check_my_f <- function(f = NULL, y = NULL, t = NULL, par_num = NULL) {
  if (is.vector(f) && !is.list(f) && length(f) > 0L && !is.null(t)) {
    f <- matrix(f, nrow = t, ncol = length(f), byrow = TRUE)
  } else if (is.vector(f) && !is.list(f) && length(f) > 0L && !is.null(y)) {
    f <- matrix(f, nrow = nrow(y), ncol = length(f), byrow = TRUE)
  } else if (is.vector(f) && !is.list(f) && length(f) > 0L) {
    f <- matrix(f, nrow = 1L)
  } else if (is.matrix(f) && nrow(f) == 1L && ncol(f) > 0L && !is.null(t)) {
    f <- matrix(as.vector(f), nrow = t, ncol = ncol(f), byrow = TRUE)
  } else if (is.matrix(f) && nrow(f) == 1L && ncol(f) > 0L && !is.null(y)) {
    f <- matrix(as.vector(f), nrow = nrow(y), ncol = ncol(f), byrow = TRUE)
  } else if (is.matrix(f) && nrow(f) > 0L && ncol(f) > 0L) {
    f <- f
  } else {
    stop("Invalid type of argument f.")
  }
  if (!is.null(y) && nrow(f) != nrow(y)) {
    stop("Argument f must have the same number of rows as argument y.")
  }
  if (!is.null(t) && nrow(f) != t) {
    stop("Argument f must have the number of rows equal to argument t.")
  }
  if (!is.null(par_num) && ncol(f) != par_num) {
    stop("Invalid number of columns of argument f.")
  }
  if (!is.numeric(f)) {
    stop("Argument f must be numeric.")
  }
  return(f)
}
# ------------------------------------------------------------------------------


# Check Distribution -----------------------------------------------------------
check_my_distr <- function(distr = NULL) {
  if (is.vector(distr) && !is.list(distr) && length(distr) == 1L && nrow(distr(filter_distr = distr)) > 0L) {
    distr <- distr
  } else {
    stop("Unknown distribution given by argument distr.")
  }
  return(distr)
}
# ------------------------------------------------------------------------------


# Check Parametrization --------------------------------------------------------
check_my_param <- function(param = NULL, distr = NULL) {
  if (!is.null(distr) && is.null(param)) {
    param <- distr(filter_distr = distr, filter_param = param, filter_default = TRUE)$param[1]
  } else if (!is.null(distr) && is.vector(param) && !is.list(param) && length(param) == 1L && nrow(distr(filter_distr = distr, filter_param = param)) > 0L) {
    param <- param
  } else {
    stop("Unknown parametrization given by argument param.")
  }
  return(param)
}
# ------------------------------------------------------------------------------


# Check Scaling ----------------------------------------------------------------
check_my_scaling <- function(scaling = NULL) {
  if (is.null(scaling)) {
    scaling <- "unit"
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "unit") {
    scaling <- "unit"
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && (scaling == "fisher_inv")) {
    scaling <- "fisher_inv"
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && (scaling == "fisher_inv_sqrt")) {
    scaling <- "fisher_inv_sqrt"
  } else {
    stop("Unknown scaling given by argument scaling.")
  }
  return(scaling)
}
# ------------------------------------------------------------------------------


# Check Equation ---------------------------------------------------------------
check_my_spec <- function(spec = NULL) {
  if (is.null(spec)) {
    scaling <- "joint"
  } else if (is.vector(spec) && !is.list(spec) && length(spec) == 1L && spec == "joint") {
    spec <- "joint"
  } else if (is.vector(spec) && !is.list(spec) && length(spec) == 1L && (spec == "reg_err")) {
    scaling <- "reg_err"
  } else {
    stop("Unknown specification of the dynamic equation given by argument spec.")
  }
  return(scaling)
}
# ------------------------------------------------------------------------------


# Check Length of Time Series --------------------------------------------------
check_my_t <- function(t = NULL, y = NULL, x = NULL, f = NULL, positive = TRUE) {
  if (is.null(t) && !is.null(y)) {
    t <- nrow(y)
  } else if (is.null(t) && !is.null(x)) {
    t <- nrow(x[[1]])
  } else if (is.null(t) && !is.null(f)) {
    t <- nrow(f)
  } else if (is.null(t) && positive) {
    t <- 1L
  } else if (is.null(t) && !positive) {
    t <- 0L
  } else if (is.vector(t) && !is.list(t) && length(t) == 1L && is.numeric(t) && !is.na(t) && is.finite(t) && t >= as.integer(positive)) {
    t <- as.integer(t)
  } else {
    stop("Invalid type of argument t.")
  }
  if (!is.null(y) && t != nrow(y)) {
    stop("Length t must be equal to the number of rows of argument y.")
  }
  if (!is.null(x) && t != nrow(x[[1]])) {
    stop("Length t must be equal to the number of rows of elements of argument x.")
  }
  return(t)
}
# ------------------------------------------------------------------------------


# Check Dimension Number -------------------------------------------------------
check_my_n <- function(n = NULL, y = NULL, f = NULL, distr = NULL, param = NULL, dim = NULL) {
  if (is.null(n) && !is.null(dim) && dim == "uni") {
    n <- 1L
  } else if (is.null(n) && !is.null(f) && !is.null(distr) && !is.null(param) && dim == "multi") {
    n <- NA
    for (i in 2:max(2, ncol(f))) {
      if (length(do.call(paste("distr", distr, param, "parameters", sep = "_"), args = list(n = i))$par_names) == ncol(f)) {
        n <- i
        break()
      }
    }
  } else if (is.null(n) && !is.null(y)) {
    n <- ncol(y)
  } else if (is.vector(n) && !is.list(n) && length(n) == 1L && is.numeric(n) && !is.na(n) && is.finite(n) && n >= 1L) {
    n <- as.integer(n)
  } else {
    stop("Invalid type of argument n.")
  }
  if (!is.null(y) && ncol(y) != n) {
    stop("Dimension n must be equal to the number of columns of argument y.")
  }
  if (!is.null(dim) && dim == "uni" && n != 1L) {
    stop("Dimension n must be equal to 1 for univariate time series.")
  }
  if (!is.null(dim) && dim == "multi" && n < 2L) {
    stop("Dimension n must be greater than or equal to 2 for multivariate time series.")
  }
  if (!is.numeric(n)) {
    stop("Invalid dimension n.")
  }
  if (is.na(n)) {
    stop("Invalid dimension n.")
  }
  return(n)
}
# ------------------------------------------------------------------------------


# Check Number of Exogeneous Variables -----------------------------------------
check_my_m <- function(m = NULL, x = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(m) && !is.null(par_num)) {
    m <- rep(0L, times = par_num)
  } else if (is.null(m) && !is.null(x)) {
    m <- sapply(x, FUN = ncol)
  } else if (is.vector(m) && !is.list(m) && !is.null(par_num) && length(m) == 1L && is.numeric(m) && !is.na(m) && is.finite(m) && m >= 0L) {
    m <- rep(as.integer(m), times = par_num)
  } else if(is.vector(m) && !is.list(m) && !is.null(group_num) && !is.null(par_in_group_num) && length(m) == group_num && is.numeric(m) && all(!is.na(m)) && all(is.finite(m)) && all(m >= 0L)) {
    m <- unlist(lapply(1:group_num, function(i) { rep(as.integer(m)[i], times = par_in_group_num[i])}))
  } else if (is.vector(m) && !is.list(m) && length(m) > 0L && is.numeric(m) && all(!is.na(m)) && all(is.finite(m)) && all(m >= 0L)) {
    m <- as.integer(m)
  } else {
    stop("Invalid type of argument m.")
  }
  if (!is.null(par_num) && length(m) != par_num) {
    stop("Ivalid number of elements in argument m.")
  }
  if (!is.null(x) && length(m) != length(x)) {
    stop("Ivalid number of elements in argument m.")
  }
  if (!is.null(x) && any(m != sapply(x, ncol))) {
    stop("Ivalid number of elements in argument m.")
  }
  return(m)
}
# ------------------------------------------------------------------------------


# Check Score Order ------------------------------------------------------------
check_my_p <- function(p = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(p) && !is.null(par_num)) {
    p <- rep(0L, times = par_num)
  } else if (is.vector(p) && !is.list(p) && !is.null(par_num) && length(p) == 1L && is.numeric(p) && !is.na(p) && is.finite(p) && p >= 0L) {
    p <- rep(as.integer(p), times = par_num)
  } else if(is.vector(p) && !is.list(p) && !is.null(group_num) && !is.null(par_in_group_num) && length(p) == group_num && is.numeric(p) && all(!is.na(p)) && all(is.finite(p)) && all(p >= 0L)) {
    p <- unlist(lapply(1:group_num, function(i) { rep(as.integer(p)[i], times = par_in_group_num[i])}))
  } else if (is.vector(p) && !is.list(p) && length(p) > 0L && is.numeric(p) && all(!is.na(p)) && all(is.finite(p)) && all(p >= 0L)) {
    p <- as.integer(p)
  } else {
    stop("Invalid type of argument p.")
  }
  if (!is.null(par_num) && length(p) != par_num) {
    stop("Ivalid number of elements in argument p.")
  }
  return(p)
}
# ------------------------------------------------------------------------------


# Check Autoregressive Order ---------------------------------------------------
check_my_q <- function(q = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(q) && !is.null(par_num)) {
    q <- rep(0L, times = par_num)
  } else if (is.vector(q) && !is.list(q) && !is.null(par_num) && length(q) == 1L && is.numeric(q) && !is.na(q) && is.finite(q) && q >= 0L) {
    q <- rep(as.integer(q), times = par_num)
  } else if(is.vector(q) && !is.list(q) && !is.null(group_num) && !is.null(par_in_group_num) && length(q) == group_num && is.numeric(q) && all(!is.na(q)) && all(is.finite(q)) && all(q >= 0L)) {
    q <- unlist(lapply(1:group_num, function(i) { rep(as.integer(q)[i], times = par_in_group_num[i])}))
  } else if (is.vector(q) && !is.list(q) && length(q) > 0L && is.numeric(q) && all(!is.na(q)) && all(is.finite(q)) && all(q >= 0L)) {
    q <- as.integer(q)
  } else {
    stop("Invalid type of argument q.")
  }
  if (!is.null(par_num) && length(q) != par_num) {
    stop("Ivalid number of elements in argument q.")
  }
  return(q)
}
# ------------------------------------------------------------------------------


# Check Constant Parameters ----------------------------------------------------
check_my_par_static <- function(par_static = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(par_static) && !is.null(par_num)) {
    par_static <- rep(FALSE, par_num)
  } else if (is.vector(par_static) && !is.list(par_static) && !is.null(par_num) && length(par_static) == 1L && is.logical(par_static) && !is.na(par_static)) {
    par_static <- rep(par_static, times = par_num)
  } else if (is.vector(par_static) && !is.list(par_static) && !is.null(group_num) && !is.null(par_in_group_num) && length(par_static) == group_num && is.logical(par_static) && all(!is.na(par_static))) {
    par_static <- unlist(lapply(1:group_num, function(i) { rep(par_static[i], times = par_in_group_num[i])}))
  } else if (is.vector(par_static) && !is.list(par_static) && length(par_static) > 0L && is.logical(par_static) && all(!is.na(par_static))) {
    par_static <- par_static
  } else {
    stop("Invalid type of argument par_static.")
  }
  if (!is.null(par_num) && length(par_static) != par_num) {
    stop("Ivalid number of elements in argument par_static.")
  }
  return(par_static)
}
# ------------------------------------------------------------------------------


# Check Links of Parameters ----------------------------------------------------
check_my_par_link <- function(par_link = NULL, par_static = NULL, par_num = NULL, group_num = NULL, par_in_group_num = NULL) {
  if (is.null(par_link) && !is.null(par_static)) {
    par_link <- !par_static
  } else if (is.null(par_link) && is.null(par_static) && !is.null(par_num)) {
    par_link <- rep(FALSE, times = par_num)
  } else if (is.vector(par_link) && !is.list(par_link) && !is.null(par_num) && length(par_link) == 1L && is.logical(par_link) && !is.na(par_link)) {
    par_link <- rep(par_link, times = par_num)
  } else if (is.vector(par_link) && !is.list(par_link) && !is.null(group_num) && !is.null(par_in_group_num) && length(par_link) == group_num && is.logical(par_link) && all(!is.na(par_link))) {
    par_link <- unlist(lapply(1:group_num, function(i) { rep(par_link[i], times = par_in_group_num[i])}))
  } else if (is.vector(par_link) && !is.list(par_link) && length(par_link) > 0L && is.logical(par_link) && all(!is.na(par_link))) {
    par_link <- par_link
  } else {
    stop("Invalid type of argument par_link.")
  }
  if (!is.null(par_num) && length(par_link) != par_num) {
    stop("Ivalid number of elements in argument par_link")
  }
  return(par_link)
}
# ------------------------------------------------------------------------------


# Check Initial Values of Parameters -------------------------------------------
check_my_par_init <- function(par_init = NULL, par_num = NULL) {
  if (is.null(par_init) && !is.null(par_num)) {
    par_init <- rep(NA_real_, times = par_num)
  } else if (is.vector(par_init) && !is.list(par_init) && length(par_init) > 0L && all(is.na(par_init))) {
    par_init <- rep(NA_real_, times = length(par_init))
  } else if (is.vector(par_init) && !is.list(par_init) && length(par_init) > 0L) {
    par_init <- par_init
  } else {
    stop("Invalid type of argument par_init.")
  }
  if (!is.null(par_num) && length(par_init) != par_num) {
    stop("Invalid length of argument par_init.")
  }
  if (!is.numeric(par_init)) {
    stop("Argument par_init must be numeric.")
  }
  return(par_init)
}
# ------------------------------------------------------------------------------


# Check the Number of Skipped Observations -------------------------------------
check_my_lik_skip <- function(lik_skip = NULL, t = NULL, p = NULL, q = NULL) {
  if (is.null(lik_skip) && !is.null(p) && !is.null(q)) {
    lik_skip <- max(p, q)
  } else if (is.null(lik_skip)) {
    lik_skip <- 0L
  } else if (is.vector(lik_skip) && !is.list(lik_skip) && length(lik_skip) == 1L && is.numeric(lik_skip) && !is.na(lik_skip) && is.finite(lik_skip) && lik_skip >= 0L) {
    lik_skip <- as.integer(lik_skip)
  } else {
    stop("Invalid type of argument lik_skip.")
  }
  if (!is.null(t) && lik_skip >= t) {
    stop("Argument lik_skip must be lower than the length of the time series t.")
  }
  return(lik_skip)
}
# ------------------------------------------------------------------------------


# Check the Number of Burned Observations in Simulations -----------------------
check_my_init_skip <- function(init_skip = NULL, p = NULL, q = NULL) {
  if (is.null(init_skip) && !is.null(p) && !is.null(q)) {
    init_skip <- max(p, q)
  } else if (is.null(init_skip)) {
    init_skip <- 0L
  } else if (is.vector(init_skip) && !is.list(init_skip) && length(init_skip) == 1L && is.numeric(init_skip) && !is.na(init_skip) && is.finite(init_skip) && init_skip >= 0L) {
    init_skip <- as.integer(init_skip)
  } else {
    stop("Invalid type of argument init_skip.")
  }
  return(init_skip)
}
# ------------------------------------------------------------------------------


# Check Coefficients -----------------------------------------------------------
check_my_coef_est <- function(coef_est = NULL, coef_num = NULL) {
  if (is.vector(coef_est) && !is.list(coef_est) && length(coef_est) > 0L && is.numeric(coef_est) && all(is.finite(coef_est))) {
    coef_est <- coef_est
  } else {
    stop("Invalid type of argument coef_est.")
  }
  if (!is.null(coef_num) && length(coef_est) != coef_num) {
    stop("Invalid length of argument coef_est.")
  }
  return(coef_est)
}
# ------------------------------------------------------------------------------


# Check Variance-Covariance Matrix ---------------------------------------------
check_my_coef_vcov <- function(coef_vcov = NULL, coef_num = NULL) {
  if (is.matrix(coef_vcov) && nrow(coef_vcov) > 0L && ncol(coef_vcov) > 0L && is.numeric(coef_vcov) && all(is.finite(coef_vcov))) {
    coef_vcov <- coef_vcov
  } else {
    stop("Invalid type of argument coef_vcov.")
  }
  if (!is.null(coef_num) && (nrow(coef_vcov) != coef_num || ncol(coef_vcov) != coef_num)) {
    stop("Invalid dimensions of argument coef_vcov.")
  }
  return(coef_vcov)
}
# ------------------------------------------------------------------------------


# Check Set of Coefficients ----------------------------------------------------
check_my_coef_set <- function(coef_set = NULL, coef_num = NULL) {
  if (is.vector(coef_set) && !is.list(coef_set) && length(coef_set) > 0L && is.numeric(coef_set) && all(is.finite(coef_set))) {
    coef_set <- matrix(coef_set, nrow = 1L)
  } else if (is.matrix(coef_set) && nrow(coef_set) > 0L && ncol(coef_set) > 0L && is.numeric(coef_set) && all(is.finite(coef_set))) {
    coef_set <- coef_set
  } else {
    stop("Invalid type of argument coef_set.")
  }
  if (!is.null(coef_num) && ncol(coef_set) != coef_num) {
    stop("Invalid length of argument coef_set.")
  }
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Check Fixed Values of Coefficients -------------------------------------------
check_my_coef_fix_value <- function(coef_fix_value = NULL, coef_num = NULL) {
  if (is.null(coef_fix_value) && !is.null(coef_num)) {
    coef_fix_value <- rep(NA_real_, coef_num)
  } else if (is.vector(coef_fix_value) && !is.list(coef_fix_value) && length(coef_fix_value) > 0L && is.numeric(coef_fix_value) && any(is.na(coef_fix_value))) {
    coef_fix_value <- as.numeric(coef_fix_value)
  } else {
    stop("Invalid type of argument coef_fix_value.")
  }
  if (!is.null(coef_num) && length(coef_fix_value) != coef_num) {
    stop("Invalid length of argument coef_fix_value.")
  }
  return(coef_fix_value)
}
# ------------------------------------------------------------------------------


# Check Parameters Fixed on Other Parameters -----------------------------------
check_my_coef_fix_other <- function(coef_fix_other = NULL, coef_fix_value = NULL, coef_num = NULL) {
  if (is.null(coef_fix_other) && !is.null(coef_num)) {
    coef_fix_other <- matrix(0, nrow = coef_num, ncol = coef_num)
  } else if (is.matrix(coef_fix_other) && nrow(coef_fix_other) == ncol(coef_fix_other) && nrow(coef_fix_other) > 0L && is.numeric(coef_fix_other)) {
    coef_fix_other <- coef_fix_other
  } else {
    stop("Invalid type of argument coef_fix_other.")
  }
  if (!is.null(coef_fix_value) && any(!is.finite(coef_fix_other[!is.na(coef_fix_value), is.na(coef_fix_value)]))) {
    stop("A subset of argument coef_fix_other given by coef_fix_value must have finite elements.")
  }
  return(coef_fix_other)
}
# ------------------------------------------------------------------------------


# Check Special Structure for Fixed Coefficients -------------------------------
check_my_coef_fix_special <- function(coef_fix_special = NULL) {
  if (is.null(coef_fix_special)) {
    coef_fix_special <- character()
  } else if (is.vector(coef_fix_special) && !is.list(coef_fix_special)) {
    coef_fix_special <- unname(unique(coef_fix_special))
  } else {
    stop("Invalid type of argument coef_fix_special.")
  }
  if (any(!(coef_fix_special %in% c("panel_structure", "zero_sum_intercept", "random_walk")))) {
    stop("Invalid value of argument coef_fix_special.")
  }
  return(coef_fix_special)
}
# ------------------------------------------------------------------------------


# Check Lower bound on Coefficients --------------------------------------------
check_my_coef_bound_lower <- function(coef_bound_lower = NULL, coef_bound_upper = NULL, par_static = NULL, par_support = NULL, par_num = NULL, coef_in_par_num = NULL, coef_num = NULL) {
  if (is.null(coef_bound_lower) && !is.null(par_static) && !is.null(par_support) && !is.null(par_num) && !is.null(coef_in_par_num)) {
    coef_bound_lower <- unlist(lapply(1:par_num, function(i) { if (!par_static[i]) { rep(-Inf, times = coef_in_par_num[i]) } else if (par_support[i] == "positive" || par_support[i] == "probability") { 1e-12 } else { -Inf } }))
  } else if (is.vector(coef_bound_lower) && !is.list(coef_bound_lower) && length(coef_bound_lower) > 0L && is.numeric(coef_bound_lower) && !is.null(par_static) && !is.null(par_support) && !is.null(par_num) && !is.null(coef_in_par_num)) {
    coef_bound_lower <- pmax(coef_bound_lower, unlist(lapply(1:par_num, function(i) { if (!par_static[i]) { rep(-Inf, times = coef_in_par_num[i]) } else if (par_support[i] == "positive" || par_support[i] == "probability") { 1e-12 } else { -Inf } })), na.rm = TRUE)
  } else if (is.vector(coef_bound_lower) && !is.list(coef_bound_lower) && length(coef_bound_lower) > 0L && is.numeric(coef_bound_lower) && all(!is.na(coef_bound_lower))) {
  } else {
    stop("Invalid argument coef_bound_lower.")
  }
  if (!is.null(coef_num) && length(coef_bound_lower) != coef_num) {
    stop("Invalid length of coef_bound_lower.")
  }
  if (!is.null(coef_bound_upper) && any(coef_bound_upper <= coef_bound_lower)) {
    stop("The lower bound must be lower than the upper bound.")
  }
  return(coef_bound_lower)
}
# ------------------------------------------------------------------------------


# Check Upper Bound on Coefficients --------------------------------------------
check_my_coef_bound_upper <- function(coef_bound_upper = NULL, coef_bound_lower = NULL, par_static = NULL, par_support = NULL, par_num = NULL, coef_in_par_num = NULL, coef_num = NULL) {
  if (is.null(coef_bound_upper) && !is.null(par_static) && !is.null(par_support) && !is.null(par_num) && !is.null(coef_in_par_num)) {
    coef_bound_upper <- unlist(lapply(1:par_num, function(i) { if (!par_static[i]) { rep(Inf, times = coef_in_par_num[i]) } else if (par_support[i] == "probability") { 1 - 1e-12 } else { Inf } }))
  } else if (is.vector(coef_bound_upper) && !is.list(coef_bound_upper) && length(coef_bound_upper) > 0L && is.numeric(coef_bound_upper) && !is.null(par_static) && !is.null(par_support) && !is.null(par_num) && !is.null(coef_in_par_num)) {
    coef_bound_upper <- pmin(coef_bound_upper, unlist(lapply(1:par_num, function(i) { if (!par_static[i]) { rep(Inf, times = coef_in_par_num[i]) } else if (par_support[i] == "probability") { 1 - 1e-12 } else { Inf } })), na.rm = TRUE)
  } else if (is.vector(coef_bound_upper) && !is.list(coef_bound_upper) && length(coef_bound_upper) > 0L && is.numeric(coef_bound_upper) && all(!is.na(coef_bound_upper))) {
  } else {
    stop("Invalid argument coef_bound_upper.")
  }
  if (!is.null(coef_num) && length(coef_bound_upper) != coef_num) {
    stop("Invalid length of coef_bound_upper.")
  }
  if (!is.null(coef_bound_lower) && any(coef_bound_upper <= coef_bound_lower)) {
    stop("The upper bound must be higher than the lower bound.")
  }
  return(coef_bound_upper)
}
# ------------------------------------------------------------------------------


# Check Starting Values of Coefficients  ---------------------------------------
check_my_coef_start <- function(coef_start = NULL, coef_bound_lower = NULL, coef_bound_upper = NULL, coef_num = NULL) {
  if (is.null(coef_start) && !is.null(coef_num)) {
    coef_start <- rep(NA_real_, times = coef_num)
  } else if (is.vector(coef_start) && !is.list(coef_start) && length(coef_start) > 0L && all(is.na(coef_start))) {
    coef_start <- rep(NA_real_, times = length(coef_start))
  } else if (is.vector(coef_start) && !is.list(coef_start) && length(coef_start) > 0L && is.numeric(coef_start)) {
    coef_start <- coef_start
  } else {
    stop("Invalid type of argument coef_start.")
  }
  if (!is.null(coef_num) && length(coef_start) != coef_num) {
    stop("Invalid length of coef_start.")
  }
  if (!is.null(coef_bound_lower) && any(coef_start <= coef_bound_lower, na.rm = TRUE)) {
    stop("Starting coefficients must be higher than the lower bound.")
  }
  if (!is.null(coef_bound_upper) && any(coef_start >= coef_bound_upper, na.rm = TRUE)) {
    stop("Starting coefficients must be lower than the upper bound.")
  }
  return(coef_start)
}
# ------------------------------------------------------------------------------


# Check for Generic Logical Scalar ---------------------------------------------
check_generic_logical_scalar <- function(arg, arg_name) {
  if (is.vector(arg) && !is.list(arg) && length(arg) == 1L && is.logical(arg) && !is.na(arg)) {
    arg <- arg
  } else {
    stop("Invalid type of argument ", arg_name)
  }
  return(arg)
}
# ------------------------------------------------------------------------------


# Check for Generic Positive Integer Scalar ------------------------------------
check_generic_positive_integer_scalar <- function(arg, arg_name) {
  if (is.vector(arg) && !is.list(arg) && length(arg) == 1L && is.numeric(arg) && !is.na(arg) && arg >= 1L) {
    arg <- as.integer(arg)
  } else {
    stop("Invalid type of argument ", arg_name)
  }
  return(arg)
}
# ------------------------------------------------------------------------------


# Check for Generic Positive Integer Scalar ------------------------------------
check_generic_probability_vector <- function(arg, arg_name) {
  if (is.vector(arg) && !is.list(arg) && length(arg) > 0L && is.numeric(arg) && all(!is.na(arg)) && all(arg >= 0) && all(arg <= 1)) {
    arg <- arg
  } else {
    stop("Invalid type of argument ", arg_name)
  }
  return(arg)
}
# ------------------------------------------------------------------------------


# Check for Generic List -------------------------------------------------------
check_generic_list <- function(arg, arg_name) {
  if (is.list(arg)) {
    arg <- arg
  } else {
    stop("Invalid type of argument ", arg_name)
  }
  return(arg)
}
# ------------------------------------------------------------------------------


# Check for Generic Function ---------------------------------------------------
check_generic_function <- function(arg, arg_name) {
  if (is.function(arg)) {
    arg <- arg
  } else {
    stop("Invalid type of argument ", arg_name)
  }
  return(arg)
}
# ------------------------------------------------------------------------------


