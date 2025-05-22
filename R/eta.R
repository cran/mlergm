#'
#' Extracts the natural parameters from an estimated \code{mlergm} model. When the model is a canonical exponential family, this is the the identity mapping from the parameter vector. However, when the model is a curved exponential family this returns the resulting canonical parameters of the curved exponential family. 
#'
#' @param object An object of class \code{mlergm}, probably produced by \code{\link{mlergm}}.
#' @param block Block identifier which should be an element of node_memb in the network. 
#'
#' @references
#'
#' Hunter, D. R., and Handcock, M. S. (2006).
#' Inference in curved exponential family models for networks.
#' Journal of Computational and Graphical Statistics, 15(3), 565-583.
#'
#' @export 
eta <- function(object, block = NULL) { 

  if (!is.mlergm(object)) { 
    stop("Argument must be an 'mlergm' object. See 'help(mlergm)' for details.\n", call. = FALSE)
  }

  if (object$estimation_status == "success") { 
    eta_vec <- ergm.eta(object$theta, object$etamap)
    names(eta_vec) <- colnames(object$mcmc_chain) 
    if (object$parameterization == "offset") { 
      if (is.null(block)) { 
        stop("When parameterization type is equal to 'offset', argument 'block' cannot be NULL.\n", call. = FALSE)
      }
      where_edge <- which("edges" == names(eta_vec))
      where_mutual <- which("mutual" == names(eta_vec))
      if (length(where_edge) > 0) { 
        eta_vec[where_edge] <- eta_vec[where_edge] - log(sum(object$node_memb == block))
      }
      if (length(where_mutual) > 0) { 
        eta_vec[where_mutual] <- eta_vec[where_mutual] + log(sum(object$node_memb == block))
      }
    }
    if (object$parameterization == "size") { 
      if (is.null(block)) { 
        stop("When parameterization type is equal to 'size', argumebt 'block' cannot be NULL.\n", call. = FALSE)
      }
      eta_vec <- eta_vec * log(sum(object$node_memb == block))
    }
    return(eta_vec)
  } else { 
    cat("\n")
    cat("Model estimation unsuccesful and estimated parameters are not available.\n")
    return(NULL)
  }
}


