compute_initial_estimate <- function(obj) {
  
  # Find the initial point for a non-curved ERGM 
  if (!is.curved(obj$net$model)) {
    net <- reorder_block_matrix(obj$net$net_list)
    if (summary(net ~ edges) == 0) { 
      stop("Network provided is empty and has no edges. Maximum likelihood estimator will not exist.", 
            call. = FALSE) 
    }
    form <- as.formula(paste0("net ~ ", obj$net$terms))
    init <- suppressMessages(
              ergm(form, estimate = "MPLE", constraints = ~ blockdiag("node_memb_group"), 
                   verbose = FALSE, eval.loglik = FALSE)
            )
    if (obj$est$parameterization == "offset") { 
      obj$est$theta <- coef(init) 
      if (!is.null(obj$net$edge_loc)) { 
        obj$est$theta[obj$net$edge_loc] <- obj$est$theta[obj$net$edge_loc] + log(median(obj$net$clust_sizes))
      }
      if (!is.null(obj$net$mutual_loc)) { 
        obj$est$theta[obj$net$mutual_loc] <- obj$est$theta[obj$net$mutual_loc] - log(median(obj$net$clust_sizes))
      }
      obj$est$theta_0 <- obj$est$theta 
    } else if (obj$est$parameterization == "size") { 
      cur_theta <- coef(init)
      which_canonical <- which(obj$net$etamap$canonical != 0)
      cur_theta[which_canonical] <- cur_theta[which_canonical] / log_fun(median(obj$net$clust_sizes))
      if (sum(obj$net$etamap$canonical == 0) > 0) {
        which_ <- which(obj$net$etamap$canonical == 0)
        if (length(which_) > 2) {
          for (ii in seq(1, length(which_), by = 2)) {
            cur_theta[which_[ii]] <- cur_theta[which_[ii]] / log_fun(median(obj$net$clust_sizes))
          }
        } else {
          cur_theta[which_[1]] <- cur_theta[which_[1]] / log_fun(median(obj$net$clust_sizes))
        }
      }
      obj$est$theta <- cur_theta 
      obj$est$theta_0 <- obj$est$theta 
    } else {  
      obj$est$theta <- coef(init)
      obj$est$theta_0 <- obj$est$theta 
    }

  # Find the initial point for a curved ERGM
  } else {
    net <- reorder_block_matrix(obj$net$net_list)
    form <- as.formula(paste0("net ~ ", obj$net$terms)) 
    #model <- ergm_model(form, net) 
    #fixed_form   <- fix.curved(form, rep(0.25, length(model$etamap$canonical)))$formula
    init <- suppressMessages(
              ergm(form, estimate = "MPLE", constraints = ~ blockdiag("node_memb_group"),
                   verbose = FALSE, eval.loglik = FALSE)
            )
    cur_theta <- coef(init) 
    if (obj$est$parameterization == "offset") {
      obj$est$theta <- coef(init)
      if (!is.null(obj$net$edge_loc)) {
        obj$est$theta[obj$net$edge_loc] <- obj$est$theta[obj$net$edge_loc] + log(median(obj$net$clust_sizes))
      }
      if (!is.null(obj$net$mutual_loc)) {
        obj$est$theta[obj$net$mutual_loc] <- obj$est$theta[obj$net$mutual_loc] - log(median(obj$net$clust_sizes))
      }
      obj$est$theta_0 <- obj$est$theta 
    } else if (obj$est$parameterization == "size") { 
      which_canonical <- which(obj$net$etamap$canonical != 0)
      cur_theta[which_canonical] <- cur_theta[which_canonical] / log_fun(median(obj$net$clust_sizes))
      if (sum(obj$net$etamap$canonical == 0) > 0) {
        which_ <- which(obj$net$etamap$canonical == 0)
        if (length(which_) > 2) {
          for (ii in seq(1, length(which_), by = 2)) {
            cur_theta[which_[ii]] <- cur_theta[which_[ii]] / log_fun(median(obj$net$clust_sizes))
          }
        } else {
          cur_theta[which_[1]] <- cur_theta[which_[1]] / log_fun(median(obj$net$clust_sizes))
        }
      }
      obj$est$theta <- cur_theta
      obj$est$theta_0 <- obj$est$theta
    } else {
      obj$est$theta <- cur_theta 
      obj$est$theta_0 <- obj$est$theta 
    }
  }

  return(obj) 
}






