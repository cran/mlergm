#' Multilevel Network 
#'
#' Function creates a multilevel network object of class \code{mlnet}. The object inherits the \code{\link{network}} class, with additional information concerning the multilevel structure. 
#'
#' The \code{mlnet} function creates an object of class \code{mlnet} which is used to access methods designed specifically for multilevel networks, including visualization methods as well as direct interface with some of the main functions, such as \code{\link{mlergm}}. Presently, the \code{mlnet} function and object class cover multilevel structure where the set of nodes is nested within known block structure.  
#' 
#' @param network Either a \code{\link{network}} object, an adjacency matrix, or an edge list. 
#' @param node_memb Vector (length equal to the number of nodes in the network) indicating to which block or group the nodes belong.
#' @param directed (\code{TRUE} or \code{FALSE}) Indicates whether the supplied network is directed or undirected. Default is \code{FALSE}.
#' @return 
#' \code{mlnet} returns an object of class \code{mlnet} which inherits the \code{\link{network}} class, with the additional vector attribute \code{node_memb}, which encodes the block membership of the multilevel netwrok. 
#' @export 
#' @examples
#'  # Show how the sampson dataset can be turned into an mlnet object 
#'  data(sampson)
#'  net <- mlnet(samplike, get.vertex.attribute(samplike, "group"))
mlnet <- function(network, node_memb, directed = FALSE) {

  # Test if the network provided can be convereted to a 'network' object  
  can_make_network <- tryCatch(as.network(network, directed = directed),
                               error = function(err) { 
                                 return(err)
                               })
  if (any(grepl("Error", can_make_network))) {
    msg1 <- "The network provided cannot be converted to an object of class 'network'."
    msg2 <- "Compatibility can be checked using as.network() and network() functions." 
    cat("\n")
    stop(paste(msg1, msg2))
  } else {
    network <- can_make_network 
  } 

  # Assign node memberships
  set.vertex.attribute(network, "node_memb", node_memb)

  # Class object as 'mlnet' to enable multilevel plotting and method features
  class(network) <- c("mlnet", "network")

  return(network)
}





