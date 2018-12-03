#' Check if the object is of class \code{mlergm} 
#'
#' Function checks if a provided object is of class \code{mlergm} (see \code{\link{mlergm}} for details). 
#'
#' @param x An objected to be checked. 
#'
#' @return \code{TRUE} if the provided object \code{x} is of class \code{mlergm}, \code{FALSE} otherwise. 
#' @export  
#' @seealso \code{\link{mlergm}}
is.mlergm <- function(x) { 
  if (is(x) == "mlergm") { 
    return(TRUE) 
  } else { 
    return(FALSE)
  }
}
