#' Check if object is of class \code{mlnet} 
#'
#' Function checks if a provided object is of class \code{mlnet} (see \code{\link{mlnet}} for details). 
#' 
#' @param x An object to be checked. 
#'
#' @return \code{TRUE} if the provided object \code{x} is of class \code{mlnet}, \code{FALSE} otherwise. 
#' @export 
#' @seealso \code{\link{mlnet}}  
is.mlnet <- function(x) { 
  if (is(x)[1] == "mlnet") { 
    return(TRUE)
  } else { 
    return(FALSE)
  }
}
