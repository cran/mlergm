#' @describeIn mlergm
#'
#' Constructs a confidence interval at a desird significance level for an estimated parameter vector. 
#'
#' @param object An object of class \code{mlergm}, probably produced by \code{\link{mlergm}}.
#' @param alpha Desired significance level for the confidence interval. 
#'
#' @importFrom stats qnorm
#' @export 
compute_CI <- function(object, alpha) { 

  if (!is.mlergm(object)) { 
    stop("Argument must be an 'mlergm' object. See 'help(mlergm)' for details.\n", call. = FALSE)
  }

  if (missing(alpha)) {
    stop("Argument 'alpha' not provided. A significance level must be provided.\n", call. = FALSE)
  } else if ((alpha <= 0) | (alpha >= 1)) { 
    stop("Argument 'alpha' must be a number in the interval (0, 1).\n", call. = FALSE)
  }

  if (object$estimation_status == "success") { 
    param_vec <- object$theta 
    se_ <- object$se 
    CI_interval <- matrix(0, nrow = length(se_), ncol = 2)
    rownames(CI_interval) <- names(se_)
    colnames(CI_interval) <- c("Lower", "Upper")
    for (p in 1:length(se_)) { 
      CI_interval[p, ] <- c(param_vec[p] - qnorm(1 - alpha / 2) * se_[p], param_vec[p] + qnorm(1 - alpha / 2) * se_[p])
    }
    return(CI_interval)
  } else { 
    cat("\n")
    cat("Model estimation unsuccesful. Parameter vector not estimated and confidence interval not available.\n")
    return(NULL)
  }
}


