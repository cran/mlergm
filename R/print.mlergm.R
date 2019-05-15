#' @describeIn mlergm
#'
#' Print method for objects of class \code{\link{mlergm}}. Indicates whether the model was succesfully estimated, as well as the model formula provided. 
#'
#' @param x An object of class \code{mlergm}, probably produced by \code{\link{mlergm}}. 
#'
#' @export
print.mlergm <- function(x, ...)  { 
  
  if (!is.mlergm(x)) { 
    stop("Argument must be an 'mlegm' object. See 'help(mlergm)' for details.\n", call. = FALSE)
  }

  if (x$estimation_status == "success") { 
    cat("\n")
    cat("mlergm model.")
    cat("\n\nEstimated model:\n\n")
    cat("  ")
    form_to_cat <- format_form_for_cat(x$formula, 2)
    cat(form_to_cat) 
    cat("\n")
  } else { 
    cat("\n") 
    cat("mlergm model.")
    cat("\n\nEstimation failed for:\n\n")
    cat("  ")
    form_to_cat <- format_form_for_cat(x$formula, 2)  
    cat(form_to_cat)
    cat("\n")
  }
}

