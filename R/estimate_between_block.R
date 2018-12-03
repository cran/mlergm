
estimate_between_block <- function(obj)  {
  
  if (obj$verbose > 0) {  
    cat("\n\nEstimating between block model.")
  }
  
  # Make a network for between block edges
  between_network <- obj$net$net 
  for (k in 1:obj$net$num_clust) { 
    v_inds <- which(obj$net$block_memb == k)
    between_network[v_inds, v_inds] <- 0  
  }

  # Get between block edge indices
  if (.Platform$OS.type == "unix") {  
    between_edge_indices <- mclapply(1:obj$net$num_clust,
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
                                     node_memb = obj$net$block_memb,
                                     directed = obj$net$directed_flag,
                                     mc.cores = obj$est$par_n_cores)
  } else { 
    cl <- makeCluster(obj$est$par_n_cores) 
    clusterEvalQ(cl, library(mlergm))
    between_edge_indices <- parLapply(cl, 
                                      1:obj$net$num_clust, 
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
                                      node_memb = obj$net$block_memb,
                                      directed = obj$net$directed_flag)
    stopCluster(cl)
  }
  between_edge_indices <- as.matrix(Reduce(rbind, between_edge_indices))
  max_edges <- nrow(between_edge_indices) 

  if (obj$net$directed_flag) { 
    stat_val <- summary(between_network ~ edges)
    if (stat_val > 0 & stat_val < max_edges) { 
      between_est <- logit(summary(between_network ~ edges) / max_edges)
      obj$est$between_theta <- between_est
      names(obj$est$between_theta) <- "Between edges" 

      obj$est$between_se <- compute_between_se(between_est, NULL, max_edges) 
      names(obj$est$between_se) <- "Between edges"

      obj$est$between_pvalue <- 2 * pnorm(-abs(obj$est$between_theta / obj$est$between_se))
      names(obj$est$between_pvalue) <- "Between edges"

    } else { 
      obj$est$between_theta <- "Between block MLE does not exist." 
    }
  } else { 
    stat_val <- summary(between_network ~ edges) 
    if (stat_val > 0 & stat_val < max_edges) { 
      between_est <- logit(summary(between_network ~ edges) / nrow(between_edge_indices))
      obj$est$between_theta <- between_est
      names(obj$est$between_theta) <- "Between edges"

      obj$est$between_se <- compute_between_se(between_est, NULL, max_edges) 
      names(obj$est$between_se) <- "Between edges"

      obj$est$between_pvalue <- 2 * pnorm(-abs(obj$est$between_theta / obj$est$between_se))
      names(obj$est$between_pvalue) <- "Between edges"

    } else { 
      obj$est$between_theta <- "Between block MLE does not exist." 
    } 
  }
  return(obj)
}


