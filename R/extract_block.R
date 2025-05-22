#' @describeIn mlnet
#'
#' Extracts a specified block subgraph from a network object of type \code{mlnet}. 
#'
#' @param net An object of class \code{mlnet}, possibly produced by \code{\link{mlnet}} or \code{\link{simulate_mlnet}}.
#' @param block Block identifier which should be an element of node_memb in the network.  
#'
#' @export 
extract_block <- function(net, block) { 

  if (!is.mlnet(net)) { 
    stop("Argument must be an 'mlnet' object. See 'help(mlnet)' for details.\n", call. = FALSE)
  }

  # Find out how many blocks there are in the network 
  sub_net <- get.inducedSubgraph(net, v = which(get.vertex.attribute(net, "node_memb") == block))

  return(sub_net)

}




