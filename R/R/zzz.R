#' @useDynLib phylotriplet, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Explicitly load the shared library
  # Use tryCatch to avoid errors during devtools::document()
  tryCatch({
    library.dynam("phylotriplet", pkgname, libname)
  }, error = function(e) {
    # During documentation generation, the package may not be installed yet
    # This is expected and can be safely ignored
    if (!interactive()) {
      # Silent fail during non-interactive use (like documentation)
      return(invisible(NULL))
    } else {
      # In interactive use, show the error
      warning("Could not load phylotriplet shared library: ", e$message)
    }
  })
}

.onUnload <- function(libpath) {
  # Unload the shared library
  tryCatch({
    library.dynam.unload("phylotriplet", libpath)
  }, error = function(e) {
    # Silently ignore errors during unload
    invisible(NULL)
  })
}
