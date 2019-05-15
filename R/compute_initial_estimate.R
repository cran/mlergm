compute_initial_estimate <- function(obj) {
  
  # Find the initial point for a non-curved ERGM 
  if (!is.curved(obj$net$model)) {
    net <- reorder_block_matrix(obj$net$net_list)
    if (summary(net ~ edges) == 0) { 
      stop("Network provided is empty and has no edges. Maximum likelihood estimator will not exist.", 
            call. = FALSE) 
    }
    form <- as.formula(paste0("net ~ ", as.character(obj$net$model$formula[3])))
    init <- suppressMessages(
              ergm(form, estimate = "MPLE", constraints = ~ blockdiag("node_memb_group"), 
                   verbose = FALSE, eval.loglik = FALSE)
            )
    if (obj$est$parameterization == "offset") { 
      obj$est$theta <- init$coef 
      if (!is.null(obj$net$edge_loc)) { 
        obj$est$theta[obj$net$edge_loc] <- obj$est$theta[obj$net$edge_loc] + log(median(obj$net$clust_sizes))
      }
      if (!is.null(obj$net$mutual_loc)) { 
        obj$est$theta[obj$net$mutual_loc] <- obj$est$theta[obj$net$mutual_loc] - log(median(obj$net$clust_sizes))
      }
      obj$est$theta_0 <- obj$est$theta 
    } else {  
      obj$est$theta <- init$coef
      obj$est$theta_0 <- obj$est$theta 
    }

  # Find the initial point for a curved ERGM
  } else {
    net <- reorder_block_matrix(obj$net$net_list)
    form <- as.formula(paste0("net ~ ", as.character(obj$net$model$formula[3])))  
    model <- ergm_model(form, net) 
    fixed_form   <- fix.curved(form, rep(0.25, length(model$etamap$canonical)))$formula
    init <- suppressMessages(
              ergm(fixed_form, estimate = "MPLE", constraints = ~ blockdiag("node_memb_group"),
                   verbose = FALSE, eval.loglik = FALSE)
            )
    if (obj$est$parameterization == "offset") {
      obj$est$theta <- init$coef
      if (!is.null(obj$net$edge_loc)) {
        obj$est$theta[obj$net$edge_loc] <- obj$est$theta[obj$net$edge_loc] + log(median(obj$net$clust_sizes))
      }
      if (!is.null(obj$net$mutual_loc)) {
        obj$est$theta[obj$net$mutual_loc] <- obj$est$theta[obj$net$mutual_loc] - log(median(obj$net$clust_sizes))
      }
      obj$est$theta_0 <- obj$est$theta 
    } else {
      obj$est$theta <- init$coef
      obj$est$theta_0 <- obj$est$coef 
    }
    theta <- numeric(0)
    num_curved_terms <- length(obj$net$model$etamap$curved)
    curved_params <- numeric(0)
    for (l in 1:num_curved_terms) {
      add_base_param_loc <- obj$net$model$etamap$curved[[l]]$from[seq(1, length(obj$net$model$etamap$curved[[l]]$from), by = 2)] 
      curved_params <- c(curved_params, add_base_param_loc) 
    }
    skip_flag <- FALSE
    full_param_names <- get_coef_names(obj$net$model, FALSE)
    iter <- 1
    for (m in 1:length(full_param_names)) {
      if ((m %in% curved_params) & !skip_flag) { 
        theta <- c(theta, obj$est$theta[iter], 0.25)
        iter <- iter + 1 
        skip_flag <- TRUE
      } else if (!skip_flag) {
        theta <- c(theta, obj$est$theta[iter])
        iter <- iter + 1 
      } else { 
        skip_flag <- FALSE
      }
    }
    names(theta) <- full_param_names
    obj$est$theta  <- obj$est$theta_0 <- theta
  }

  return(obj) 
}






