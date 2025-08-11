
# SILENT FUNCTIONS

# Message but Theta ------------------------------------------------------------
but_theta <- function(msg) {
  txt <- conditionMessage(msg)
  if (startsWith(txt, "Theta")) {
    message(txt)
  }
  invokeRestart("muffleMessage")
}
# ------------------------------------------------------------------------------


# No Warnings and Messages -----------------------------------------------------
be_silent <- function(expr) {
  suppressMessages(suppressWarnings(expr))
}
# ------------------------------------------------------------------------------


# No Warnings and Messages Except Theta ----------------------------------------
be_silent_but_theta <- function(expr) {
    withCallingHandlers(suppressWarnings(expr), message = but_theta)
}
# ------------------------------------------------------------------------------


# Try Without Warnings and Messages --------------------------------------------
try_and_be_silent <- function(expr) {
  suppressMessages(suppressWarnings(try(expr, silent = TRUE)))
}
# ------------------------------------------------------------------------------


# Try Without Warnings and Messages  Except Theta ------------------------------
try_and_be_silent_but_theta <- function(expr) {
  withCallingHandlers(suppressWarnings(try(expr, silent = TRUE)), message = but_theta)
}
# ------------------------------------------------------------------------------


