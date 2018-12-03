#' Check if object is of class \code{gof_mlergm}
#'
#' Function checks if a provided object is of class \code{gof_mlergm} (see \code{\link{gof.mlergm}} for details). 
#'
#' @param x An object to be checked. 
#'
#' @return \code{TRUE} if the provided object \code{x} is of class \code{gof_mlergm}, \code{FALSE} otherwise. 
#' @export 
#' @seealso \code{\link{mlergm}}, \code{\link{gof.mlergm}}
is.gof_mlergm <- function(x) {
  if (is(x)[1] == "gof_mlergm") { 
    return(TRUE)
  } else { 
    return(FALSE)
  }
}
