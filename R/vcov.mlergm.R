#' @describeIn mlergm
#'
#' Extracts the estimated Fisher information matrix of the model from an estimated \code{mlergm} model. When the model is a canonical exponential family (i.e., no curved terms such as gwesp or gwdegree), this is equal to the variance-covariance matrix of the canonical statistics. In the case that specified model is a curved exponential family, the Fisher information matrix is a transformation of the variance-covariance matrix of the exponential family and is given by Equation (3.2) of Hunter & Handcock (Journal of Computational and Graphical Statistics, 2006, DOI: 10.1198/106186006X133069). 
#'
#' @param object An object of class \code{mlergm}, probably produced by \code{\link{mlergm}}.
#' @param \dots Additional arguments to be passed if necessary.  
#'
#' @references
#'
#' Hunter, D. R., and Handcock, M. S. (2006).
#' Inference in curved exponential family models for networks.
#' Journal of Computational and Graphical Statistics, 15(3), 565-583.
#'
#' @export 
vcov.mlergm <- function(object, ...) { 

  if (!is.mlergm(object)) { 
    stop("Argument must be an 'mlergm' object. See 'help(mlergm)' for details.\n", call. = FALSE)
  }

  if (object$estimation_status == "success") { 
    return(object$information_matrix)
  } else { 
    cat("\n")
    cat("Model estimation unsuccesful. Information matrix not estimated.\n")
    return(NULL)
  }
}


