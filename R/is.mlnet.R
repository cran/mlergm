#' Check if object is of class \code{mlnet} 
#'
#' Function checks if a provided object is of class \code{mlnet} (see \code{\link[mlergm]{mlnet}} for details). 
#' 
#' @param x An object to be checked. 
#'
#' @return \code{TRUE} if the provided object \code{x} is of class \code{mlnet}, \code{FALSE} otherwise. 
#' @export 
#' @seealso \code{\link[mlergm]{mlnet}}  
is.mlnet <- function(x) { 
  if (is(x)[1] == "mlnet") { 
    return(TRUE)
  } else { 
    return(FALSE)
  }
}
