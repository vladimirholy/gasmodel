
# SILENT FUNCTIONS

# No Warnings and Messages -----------------------------------------------------
be_silent <- function(obj) {
  suppressMessages(suppressWarnings(obj))
}
# ------------------------------------------------------------------------------


# Try Without Warnings and Messages --------------------------------------------
try_and_be_silent <- function(obj) {
  suppressMessages(suppressWarnings(try(obj, silent = TRUE)))
}
# ------------------------------------------------------------------------------


