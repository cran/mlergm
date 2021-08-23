
initialize_object <- function(net_list,
                              net, 
                              form, 
                              terms, 
                              block_memb,
                              param_list,
                              theta_init,
                              sim_param,
                              est_param,
                              verbose,
                              parameterization) {

  ##### Create object top structure
  obj <- list(net = NULL,
              sim = NULL,
              est = NULL)
  obj$verbose <- verbose; 

  ##### Construct obj$net
  obj$net <- list(num_clust = length(net_list),
                  clust_sizes = numeric(length(net_list)),
                  model = param_list$model,
                  terms = terms, 
                  net_list = net_list,
                  net = net,
                  block_memb = block_memb, 
                  dims = param_list$block_dims,
                  num_terms = param_list$model_dim,
                  directed_flag = is.directed(net_list[[1]]), 
                  na_flag = FALSE,
                  na_clust_flag = vector(length = length(net_list)),
                  obs_stats = NULL,
                  obs_stats_step = NULL,
                  edge_loc = param_list$edge_loc,
                  mutual_loc = param_list$mutual_loc,
                  etamap = param_list$eta_map,
                  theta_names = NULL)
  obj$net$model$formula <- form 

  # Fill in the number of nodes in each cluster
  for (i in 1:obj$net$num_clust) {
    obj$net$clust_sizes[i] <- network.size(obj$net$net_list[[i]])
  }

  # Check if there are missing edges
  for (i in 1:obj$net$num_clust) {
    if (network.naedgecount(obj$net$net_list[[i]]) > 0) {
      obj$net$na_clust_flag[i] <- TRUE
    } else {
      obj$net$na_clust_flag[i] <- FALSE
    }
  }

  if (sum(obj$net$na_clust_flag == TRUE) > 0) {
    obj$net$na_flag <- TRUE
  }

  if (!obj$net$na_flag) {
    obj$net$obs_stats <- numeric(obj$net$num_terms)
  }

  # Set the theta_names for printing 
  theta_names <- get_coef_names(obj$net$model, !is.curved(obj$net$model$formula))
  num_chars <- nchar(theta_names)
  max_char <- max(num_chars)
  for (i in 1:length(theta_names)) { 
    theta_names[i] <- paste0(paste(rep(" ", max_char - num_chars[i]), collapse = ""), theta_names[i]) 
  }
  obj$net$theta_names <- theta_names  


  ##### Construct obj$sim
  obj$sim <- sim_param
  obj$sim$seq <- TRUE
  obj$sim$stats <- rep(list(matrix(0, nrow = obj$sim$num_obs,
                             ncol = obj$net$num_terms)), obj$net$num_clust)
  for (i in 1:length(obj$sim$stats)) { 
    colnames(obj$sim$stats[[i]]) <- param_list$statistic_names; 
  }

  if (obj$net$na_flag) {
    obj$sim$cond_stats <- obj$sim$stats
  }
  


  ##### Construct obj$est 
  obj$est           <- est_param
  obj$est$theta     <- theta_init
  obj$est$theta_0   <- theta_init
  obj$est$sizes     <- numeric(length(obj$net$net_list));
  obj$est$gamma     <- 1
  obj$est$score_val <- NULL
  for (i in 1:obj$net$num_clust) {
    obj$est$sizes[i] <- network.size(obj$net$net_list[[i]])
  }
  obj$est$parameterization <- parameterization;
  if (obj$net$num_clust < obj$est$par_n_cores) { 
    obj$est$par_n_cores <- obj$net$num_clust;
  }

  return(obj)
}







