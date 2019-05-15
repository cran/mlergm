#' Print summary of a \code{gof_mlergm} object. 
#'
#' Prints a formatted summary output for \code{gof_mlergm} object which was produced by \code{\link{gof.mlergm}}. 
#' 
#' @param x An object of class \code{gof_mlergm}, probably produced by \code{\link{gof.mlergm}}. 
#' @param \dots Additional arguments to be passed if necessary. 
#'
#' @seealso \code{\link{gof.mlergm}}
#' @export
print.gof_mlergm <- function(x, ...) { 
  
  if (!is.gof_mlergm(x)) { 
    stop("Argument must be a 'gof_mlergm' object. See 'help(mlergm)' and 'help(gof.mlergm)' for details.\n",
         call. = FALSE)
  }  

  cat("\n") 
  cat("gof_mlergm object.") 
  cat("\n\n    Contains the following GOF statistics:\n") 
  for (stat_name in x$stat_names) { 
    cat("      ")
    cat(paste0("\n", stat_name))
  }
  cat("\n\n")
}

