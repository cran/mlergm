#' Simulate a multilevel network 
#'
#' Function simulates a multilevel network by specifying a network size, node block memberships, and within-block and between-block models. The function currently only suppports block-models where between-block edges are dyad-independent. 
#'
#' Simulation of multilevel block networks is done with a Monte-Carlo Markov chain (MCMC) and can be done in parallel where \code{\link{set_options}} can be used to adjust the simulation settings (such as \code{burnin}, \code{interval}, and \code{sample_size}). Each within-block subgraph is given its own Markov chain, and so these settings are the settings to be used for each within-block chain. 
#'
#' @param form A \code{\link{formula}} object of the form \code{network ~ model terms} which specifies how the within-block subgraphs are modeled.
#' @param node_memb Vector of node block memberships.
#' @param theta A vector of model parameters (coefficients) for the ERGM governing the within-subgraph edges. 
#' @param parameterization Parameterization options include 'standard', 'offset', or 'size'.
#' \itemize{
#' \item 'standard' : Does not adjust the individual block parameters for size.
#' \item 'offset' : The offset parameterization uses edge and mutual offsets along the lines of Krivitsky, Handcock, and Morris (2011) and Krivitsky and Kolaczyk (2015). The edge parameter is offset by \eqn{-log n(k)} and the mutual parameter is offset by \eqn{+log n(k)}, where \eqn{n(k)} is the size of the kth block.
#' \item 'size' : Multiplies the block parameters by \eqn{log n(k)}, where \eqn{n(k)} is the size of the kth block.
#' }
#' @param seed Seed to be provided for reproducibility.  
#' @param between_form A \code{\link{formula}} object of the form \code{~ model terms} which specifies how the within-block subgraphs are modeled. 
#' @param between_theta A vector of model parameters (coefficients) for the ERGM governing the between-subgraph edges.
#' @param between_prob A probability which specifies how edges between blocks are governerd. An ERGM (\code{between_form} and \code{between_theta}) cannot be specified together with \code{between_prob}.  
#' @param options Use \code{\link{set_options}} to change the simulation options. Note that some options are only valid for estimation using \code{\link{mlergm}}.
#' @return \code{simulate_mlnet} returns an objects of class \code{\link{mlnet}}.
#' @export
#' @importFrom stats rbinom 
#' @keywords simulation
#' @examples 
#' \donttest{
#' # Create a K = 2 block network with edge + gwesp term 
#' net <- simulate_mlnet(form = network.initialize(30, directed = FALSE) ~ edges + gwesp, 
#'                       node_memb = c(rep(1, 15), rep(2, 15)),
#'                       theta = c(-3, 0.5, 1.0), 
#'                       between_prob = 0.01,
#'                       options = set_options(number_cores = 2, burnin = 2000))
#'
#' # Simulate a K = 2 block directed network, specifying a formula for between edges
#' net <- simulate_mlnet(form = network.initialize(30, directed = TRUE) ~ edges + gwesp,
#'                       node_memb = c(rep(1, 15), rep(2, 15)),
#'                       theta = c(-3, 0.5, 1.0),
#'                       between_form = ~ edges + mutual, 
#'                       between_theta = c(-4, 2),
#'                       options = set_options(number_cores = 2, burnin = 2000))
#' }
simulate_mlnet <- function(form, node_memb, theta, parameterization = "standard",
                           seed = NULL, between_form = NULL, between_theta = NULL, between_prob = NULL,
                           options = set_options()) { 

  # Check required arguments are provided correctly 
  if (missing(form)) { 
    cat("\n")
    stop("\nArgument 'form' not provided. A formula object must be provided.\n", call. = FALSE)
  }
  if (missing(theta)) { 
    cat("\n")
    stop("\nArgument 'theta' not provided. The within-block model parameters must be provided.",
         call. = FALSE)
  }

  base_net <- get_network_from_formula(form)
  if (missing(node_memb)) {
    if (any(grepl("mlnet", is(base_net)))) {
      node_memb <- get.vertex.attribute(base_net, "node_memb")
    } else {  
      cat("\n")
      stop("\nArgument 'node_memb' not provided. The node memberships must be provided.\n",
           call. = FALSE)
    }
  }

  # Check entries 
  if (!(is.null(between_form) == is.null(between_theta))) { 
    cat("\n")
    stop("\nThe arguments 'between_form' and 'between_theta' must be both specified.\n",
         call. = FALSE)
  }
  if (!is.null(between_form) & !is.null(between_prob)) { 
    cat("\n")
    stop("\nArguments 'between_form' and 'between_prob' cannot both be specified.\n",
         call. = FALSE)
  }


  # If a seed is provided, set it
  if (!is.null(seed)) {
    check_integer(seed, "seed");
    set.seed(seed, "L'Ecuyer");
  }

  # Extract the network properties and create an mlnet 
  base_net <- get_network_from_formula(form)
  directed <- is.directed(base_net)
  net <- mlnet(base_net, node_memb, directed = directed)
  model_temp <- ergm_model(form, base_net)
  rm(base_net); clean_mem() 
 
  memb_list <- check_and_convert_memb(node_memb)
  memb_labels <- memb_list$memb_labels
  memb_internal <- memb_list$memb_internal
  rm(memb_list); clean_mem()
  net_list <- make_net_list(net, memb_internal)
 
  model_temp <- ergm_model(form, net);
  within_terms <- get_terms_from_formula(form, net);
  param_list <- check_parameterization_type(net_list, within_terms, parameterization, model_temp);
  edge_loc <- param_list$edge_loc
  mutual_loc <- param_list$mutual_loc
  rm(param_list); 
  rm(model_temp); clean_mem() 

  obj <- initialize_simulation_object(options = options, 
                                      net_list = net_list,
                                      edge_loc = edge_loc, 
                                      mutual_loc = mutual_loc) 

 
  # Get within model terms 
  within_terms <- get_terms_from_formula(form, obj$mlnet)
  
  # Simulate the within block
    if (.Platform$OS.type == "unix") { 
    within_nets <- mclapply(net_list, 
                            function(sub_net, parameterization, theta, within_terms, burnin, interval, edge_loc, mutual_loc) {
                            sub_form <- as.formula(paste0("sub_net ~ ", within_terms))
                              if (parameterization == "offset") { 
                                if (is.numeric(edge_loc)) { 
                                  theta[edge_loc] <- theta[edge_loc] - log_fun(network.size(sub_net))
                                }
                                if (is.numeric(mutual_loc)) {
                                  theta[mutual_loc] <- theta[mutual_loc] + log_fun(network.size(sub_net))
                                }
                              } else if (parameterization == "size") { 
                                model <- ergm_model(sub_form, sub_net)
                                which_canonical <- which(model$etamap$canonical != 0)
                                theta[which_canonical] <- theta[which_canonical] * log_fun(network.size(sub_net))
                                if (sum(model$etamap$canonical == 0) > 0) {
                                  which_ <- which(model$etamap$canonical == 0)
                                  if (length(which_) > 2) {
                                    for (ii in seq(1, length(which_), by = 2)) {
                                      theta[which_[ii]] <- theta[which_[ii]] * log_fun(network.size(sub_net))
                                    }
                                  } else { 
                                    theta[which_[1]] <- theta[which_[1]] * log_fun(network.size(sub_net))
                                  }
                                }
                              }
                              sub_net <- simulate(sub_form, 
                                                  nsim = 1,
                                                  coef = theta, 
                                                  control = control.simulate(MCMC.burnin = burnin, 
                                                                             MCMC.interval = interval))
                              return(sub_net)

                            }, 
                            parameterization = parameterization, 
                            theta = theta, 
                            within_terms = within_terms, 
                            burnin = obj$burnin, 
                            interval = obj$interval, 
                            edge_loc = obj$edge_loc, 
                            mutual_loc = obj$mutual_loc, 
                            mc.cores = obj$par_n_cores)
  } else { 
    cl <- makeCluster(obj$par_n_cores)
    clusterEvalQ(cl, library(mlergm))
    within_nets <- parLapply(cl,
                             net_list, 
                             function(sub_net, parameterization, theta, within_terms, burnin, interval, edge_loc, mutual_loc) {
                               sub_form <- as.formula(paste0("sub_net ~ ", within_terms))
                               if (parameterization == "offset") {
                                if (is.numeric(edge_loc)) {
                                  theta[edge_loc] <- theta[edge_loc] - log_fun(network.size(sub_net))
                                }
                                if (is.numeric(mutual_loc)) {
                                  theta[mutual_loc] <- theta[mutual_loc] + log_fun(network.size(sub_net))
                                }
                              } else if (parameterization == "size") { 
                                model <- ergm_model(sub_form, sub_net)
                                which_canonical <- which(model$etamap$canonical != 0)
                                theta[which_canonical] <- theta[which_canonical] * log_fun(network.size(sub_net))
                                if (sum(model$etamap$canonical == 0) > 0) {
                                  which_ <- which(model$etamap$canonical == 0)
                                  if (length(which_) > 2) {
                                    for (ii in seq(1, length(which_), by = 2)) {
                                      theta[which_[ii]] <- theta[which_[ii]] * log_fun(network.size(sub_net))
                                    }
                                  } else {
                                    theta[which_[1]] <- theta[which_[1]] * log_fun(network.size(sub_net))
                                  }
                                }
                              }
                              sub_net <- simulate(sub_form,
                                                  nsim = 1,
                                                  coef = theta,
                                                  control = control.simulate(MCMC.burnin = burnin,
                                                                             MCMC.interval = interval))
                              return(sub_net)
                            },
                            parameterization = parameterization,
                            theta = theta, 
                            within_terms = within_terms, 
                            burnin = obj$burnin, 
                            interval = obj$interval, 
                            edge_loc = obj$edge_loc, 
                            mutual_loc = obj$mutual_loc)
  }
  for (k in 1:length(within_nets)) { 
    cur_v <- which(node_memb == unique(node_memb)[k])
    net[cur_v, cur_v] <- within_nets[[k]][ , ] 
  } 
  

  # Simulate between-block edges 
  # Create a list of all between-block edge indices
  if (.Platform$OS.type == "unix") {
    between_edge_indices <- mclapply(1:length(unique(node_memb)), 
                                     function(ind, node_memb, directed) {
                                       u_memb <- unique(node_memb)
                                       K <- length(u_memb)
                                       if (directed) { 
                                         clust_range <- setdiff(1:K, ind)
                                         } else { 
                                         if (ind < K) { 
                                           clust_range <- (ind + 1):K
                                         } else {
                                           return(matrix(0, nrow = 0, ncol = 2))
                                         }
                                       } 
                                       indices <- matrix(0, nrow = 0, ncol = 2) 
                                       clust_nodes <- which(node_memb == u_memb[ind])
                                       for (k in clust_range) { 
                                         cur_nodes <- which(node_memb == u_memb[k])
                                         indices <- rbind(indices, expand.grid(clust_nodes, cur_nodes))
                                       }
                                       return(indices)
                                     },
                                     node_memb = node_memb,
                                     directed = directed, 
                                     mc.cores = obj$par_n_cores)
  } else { 
  between_edge_indices <- parLapply(cl, 
                                    1:length(unique(node_memb)),
                                    function(ind, node_memb, directed) {
                                      u_memb <- unique(node_memb)
                                      K <- length(u_memb)
                                      if (directed) {
                                        clust_range <- setdiff(1:K, ind)
                                      } else {
                                        if (ind < K) {
                                          clust_range <- (ind + 1):K
                                        } else {
                                          return(matrix(0, nrow = 0, ncol = 2))
                                        }
                                      }
                                      indices <- matrix(0, nrow = 0, ncol = 2)
                                      clust_nodes <- which(node_memb == u_memb[ind])
                                      for (k in clust_range) {
                                        cur_nodes <- which(node_memb == u_memb[k])
                                        indices <- rbind(indices, expand.grid(clust_nodes, cur_nodes))
                                      }
                                      return(indices)
                                    },
                                    node_memb = node_memb, 
                                    directed = directed)
    stopCluster(cl)
  }
  between_edge_indices <- as.matrix(Reduce(rbind, between_edge_indices))

  if (!is.null(between_prob)) { 
    net[between_edge_indices] <- rbinom(nrow(between_edge_indices), 1, between_prob)
  } else if (!is.null(between_form)) { 
    hollow_net <- network.initialize(network.size(net), is.directed(net))
    hollow_net[between_edge_indices] <- net[between_edge_indices]
    between_form_sim <- as.formula(paste0("hollow_net ~ ", paste(all.vars(between_form), collapse = " + ")))
    hollow_net <- simulate(between_form_sim, 
                           coef = between_theta, 
                           constraints = ~ fixallbut(as.edgelist(between_edge_indices, n = length(node_memb))), 
                           control = control.simulate(MCMC.burnin = obj$burnin, 
                                                      MCMC.interval = obj$interval),
                           nsim = 1)
    net[between_edge_indices] <- hollow_net[between_edge_indices] 
    rm(hollow_net); clean_mem() 

  } else {
    # No between edge processes, so set all edges to zero 
    net[between_edge_indices] <- 0 
  } 

  return(net) 
}






