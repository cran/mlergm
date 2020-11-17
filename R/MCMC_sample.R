MCMC_sample <- function(obj) {
  
  # Get network sizes 
  sizes <- numeric(obj$net$num_clust) 
  for (i in 1:obj$net$num_clust) { 
    sizes[i] <- obj$net$net_list[[i]]$gal$n
  } 

  # Sample sufficient statistics for each block
  obj$sim$stats <- rep(list(matrix(0, nrow = obj$sim$num_obs, ncol = obj$net$num_terms)), obj$net$num_clust)

  if (obj$net$na_flag) {
    obj$sim$cond_stats <- obj$sim$stats
  } else {
    obj$net$obs_stats <- numeric(obj$net$num_terms)
  }

  dims <- obj$net$dims
  for (i in 1:obj$net$num_clust) {
    if (obj$verbose == 2) {
      cat(paste0(" ", i))
    }

    # Setup net and formula for ergm::simulate simulation
    cur_net <- obj$net$net_list[[i]]
    if (obj$est$parameterization == "standard") {
      cur_theta <- obj$est$theta_0

    } else if (obj$est$parameterization == "size") {
      which_canonical <- which(obj$net$etamap$canonical != 0)
      cur_theta <- obj$est$theta_0
      cur_theta[which_canonical] <- cur_theta[which_canonical] * log_fun(obj$est$sizes[i])
      if (sum(obj$net$etamap$canonical == 0) > 0) {
        which_ <- which(obj$net$etamap$canonical == 0)
        if (length(which_) > 2) {
          for (ii in seq(1, length(which_), by = 2)) { 
            cur_theta[which_[ii]] <- cur_theta[which_[ii]] * log_fun(obj$est$sizes[i])
          }
        } else {
          cur_theta[which_[1]] <- cur_theta[which_[1]] * log_fun(obj$est$sizes[i])
        }
      }
    } else if (obj$est$parameterization == "offset") {
      cur_theta <- obj$est$theta_0
      if (is.numeric(obj$net$edge_loc)) { 
        cur_theta[obj$net$edge_loc] <- cur_theta[obj$net$edge_loc] - log_fun(obj$est$sizes[i]) 
      }
      if (is.numeric(obj$net$mutual_loc)) {  
        cur_theta[obj$net$mutual_loc] <- cur_theta[obj$net$mutual_loc] + log_fun(obj$est$sizes[i])
      }
    }
    form <- as.formula(paste0("cur_net ~ ", obj$net$terms))
    
    # Simulate sufficient statistics
    cur_sim_stat <- simulate(form,
                             output = "stats", 
                             sequential = obj$sim$seq,
                             control = control.simulate(MCMC.burnin = obj$sim$burnin,
                                                        MCMC.interval = obj$sim$interval),
                             nsim = obj$sim$num_obs,
                             coef = cur_theta)
    if (obj$est$parameterization == "size") {
      cur_sim_stat <- cur_sim_stat * log_fun(obj$est$sizes[i])
    }

    if (is.curved(obj$net$model)) {
      num_curved <- sum(obj$net$model$etamap$canonical == 0) / 2
      mod_temp <- ergm_model(form, cur_net)
      cur_curved_ind <- numeric(0)
      curved_ind <- numeric(0)
      for (cur_t in 1:num_curved) { 
        curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
        curved_ind <- c(curved_ind, curved_ind_)
        
        cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
        cur_curved_ind <- c(cur_curved_ind, cur_curved_ind_)
        
        cur_len <- length(cur_curved_ind_)
        
        obj$sim$stats[[i]][ , curved_ind_[1:cur_len]] <- cur_sim_stat[ , cur_curved_ind_]
      }
      obj$sim$stats[[i]][ ,-curved_ind] <- cur_sim_stat[ , -cur_curved_ind]
    } else {
      obj$sim$stats[[i]][ , dims[i, ]] <- cur_sim_stat
    }
    
    # Impute observed sufficient statistics if missing data
    if (obj$net$na_flag) {
      # If there are missing edges in the network perform a simulation
      if (obj$net$na_clust_flag[i]) {
        cur_sim_cond <- simulate(form,
                                 output = "stats", 
                                 sequential = obj$sim$seq, 
                                 control = control.simulate(MCMC.burnin = obj$sim$burnin,
                                                            MCMC.interval = obj$sim$interval),
                                 nsim = obj$sim$num_obs,
                                 coef = cur_theta,
                                 constraints = ~ observed)

        if (obj$est$parameterization == "size") {
          cur_sim_cond <- cur_sim_cond * log_fun(obj$est$sizes[i])
        }

        if (is.curved(obj$net$model)) {
          for (cur_t in 1:num_curved) { 
            cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
            curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to

            cur_len <- length(cur_curved_ind_)

            obj$sim$cond_stats[[i]][, curved_ind_[1:cur_len]] <- cur_sim_cond[ , cur_curved_ind_]
          }
          obj$sim$cond_stats[[i]][ ,-curved_ind] <- cur_sim_cond[ , -cur_curved_ind]
        } else {
          obj$sim$cond_stats[[i]][ , dims[i, ]] <- cur_sim_cond
        }


      } else {
        # If there are no missing edges in the network, the network is fixed
        cur_sim_cond <- rep_row(summary(form), obj$sim$num_obs)
        if (obj$est$parameterization == "size") {
          cur_sim_cond <- cur_sim_cond * log_fun(obj$est$sizes[i])
        }

        if (is.curved(obj$net$model)) {
          for (cur_t in 1:num_curved) { 
            curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
            cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to

            cur_len <- length(cur_curved_ind_)

            obj$sim$cond_stats[[i]][ , curved_ind_[1:cur_len]] <- cur_sim_cond[ , cur_curved_ind_]
          }
          obj$sim$cond_stats[[i]][ , -curved_ind] <- cur_sim_cond[ , -cur_curved_ind]
        } else {
          obj$sim$cond_stats[[i]][ , dims[i, ]] <- cur_sim_cond
        }

      }
    } else {
      sum_val <- summary(form)
      if (obj$est$parameterization == "size") {
        sum_val <- sum_val * log_fun(obj$est$sizes[i])
      }

      if (is.curved(obj$net$model)) {

        for (cur_t in 1:num_curved) { 
          curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
          cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
          cur_len <- length(cur_curved_ind_)
          obj$net$obs_stats[curved_ind_[1:cur_len]] <- obj$net$obs_stats[curved_ind_[1:cur_len]] + 
                                                          sum_val[cur_curved_ind_]
        }
        obj$net$obs_stats[-curved_ind] <- obj$net$obs_stats[-curved_ind] + sum_val[-cur_curved_ind]
      } else {
        obj$net$obs_stats[dims[i, ]] <- obj$net$obs_stats[dims[i, ]] + sum_val
      }
    }
  }
  if (obj$net$na_flag) {  
    obj$net$obs_stats <- apply(Reduce("+", obj$sim$cond_stats), 2, mean, na.rm = TRUE)
  }
  return(obj)
}



MCMC_sample_par <- function(obj) {
  
  # Compute block sizes 
  sizes <- numeric(obj$net$num_clust) 
  for (i in 1:obj$net$num_clust) { 
    sizes[i] <- obj$net$net_list[[i]]$gal$n
  }

  # Sample sufficient statistics for each block
  obj$sim$stats <- rep(list(matrix(0, nrow = obj$sim$num_obs, ncol = obj$net$num_terms)), obj$net$num_clust)

  if (obj$net$na_flag) {
    obj$sim$cond_stats <- obj$sim$stats
  } else {
    obj$net$obs_stats <- numeric(obj$net$num_terms)
  }


  dims <- obj$net$dims
  splits_ <- msplit(1:obj$net$num_clust, 1:obj$est$par_n_cores)
  split_sizes <- numeric(obj$est$par_n_cores)
  split_ind <- rep(list(NULL), length(split_sizes))
  for (i in 1:length(split_sizes)) { 
    split_sizes[i] <- length(splits_[[i]])
  }
  for (i in 1:length(split_sizes)) {
    if (i == 1) {  
      split_ind[[i]] <- 1:split_sizes[i]
    } else { 
      split_ind[[i]] <- (1 + sum(split_sizes[1:(i-1)])):(sum(split_sizes[1:i]))
    } 
  }
  
  if (.Platform$OS.type == "windows") { 
    cl <- makeCluster(obj$est$par_n_cores)
    clusterEvalQ(cl, library(mlergm))
    sim_stats <- parLapply(cl, 1:length(split_ind), par_sim_fun, obj, split_ind = split_ind)
    stopCluster(cl)
  } else { 
    sim_stats <- mclapply(1:length(split_ind), par_sim_fun, obj = obj, 
                            split_ind = split_ind,
                            mc.cores = obj$est$par_n_cores)
  }
  sim_stats <- unlist(sim_stats, recursive = FALSE)

  for (i in 1:obj$net$num_clust) {
    obj$sim$stats[[i]] <- sim_stats[[i]]$stat_matrix
    if (obj$net$na_flag) {
      obj$sim$cond_stats[[i]] <- sim_stats[[i]]$cond_matrix
    } else {
      obj$net$obs_stats <- obj$net$obs_stats +
                                        sim_stats[[i]]$obs_stats
    }
  }
  if (obj$net$na_flag) { 
    obj$net$obs_stats <- apply(Reduce("+", obj$sim$cond_stats), 2, mean, na.rm = TRUE)
  }
  return(obj)
}


par_sim_fun <- function(i, obj, split_ind) {


  # Setup net and formula for ergm::simulate simulation
  dims <- obj$net$dims
  inds_ <- split_ind[[i]]
  stat_list <- rep(list(NULL), length(inds_))
  for (ll in 1:length(inds_)) { 
    cur_ind <- inds_[ll]
    cur_net <- obj$net$net_list[[cur_ind]]


    if (obj$est$parameterization == "standard") {
      cur_theta <- obj$est$theta_0
  
    } else if (obj$est$parameterization == "size") {
      cur_theta <- obj$est$theta_0
      which_canonical <- which(obj$net$etamap$canonical != 0)
      cur_theta[which_canonical] <- cur_theta[which_canonical] * log_fun(obj$est$sizes[cur_ind])
      if (sum(obj$net$etamap$canonical == 0) > 0) {
        which_ <- which(obj$net$etamap$canonical == 0)
        if (length(which_) > 2) {
          for (ii in seq(1, length(which_), by = 2)) { 
            cur_theta[which_[ii]] <- cur_theta[which_[ii]] * log_fun(obj$est$sizes[cur_ind])
          }
        } else {
          cur_theta[which_[1]] <- cur_theta[which_[1]] * log_fun(obj$est$sizes[cur_ind])
        }
      }
    } else if (obj$est$parameterization == "offset") {
      cur_theta <- obj$est$theta_0 
      if (is.numeric(obj$net$edge_loc)) {  
        cur_theta[obj$net$edge_loc] <- cur_theta[obj$net$edge_loc] - log_fun(obj$est$sizes[cur_ind])
      }
      if (is.numeric(obj$net$mutual_loc)) { 
        cur_theta[obj$net$mutual_loc] <- cur_theta[obj$net$mutual_loc] + log_fun(obj$est$sizes[cur_ind])
      }
    }

    stat_matrix <- matrix(0, nrow = obj$sim$num_obs,
                             ncol = obj$net$num_terms)
    form <- as.formula(paste("cur_net ~ ", obj$net$terms)) 
  
    if (obj$net$na_flag) {
      cond_matrix <- stat_matrix
    } else {
      obs_stats <- numeric(length(obj$net$obs_stats))
    }

    # Simulate sufficient statistics
    cur_sim_stat <- simulate(form,
                             output = "stats",
                             sequential = obj$sim$seq,
                             control = control.simulate(MCMC.burnin = obj$sim$burnin,
                                                        MCMC.interval = obj$sim$interval),
                             nsim = obj$sim$num_obs,
                             coef = cur_theta)

    if (obj$est$parameterization == "size") {
      cur_sim_stat <- cur_sim_stat * log_fun(obj$est$sizes[cur_ind])
    }

    if (is.curved(obj$net$model)) {
      num_curved <- sum(obj$net$model$etamap$canonical == 0) / 2
      mod_temp <- ergm_model(form, cur_net)
      cur_curved_ind <- numeric(0)
      curved_ind <- numeric(0)
      for (cur_t in 1:num_curved) { 
        curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
        curved_ind <- c(curved_ind, curved_ind_)

        cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to
        cur_curved_ind <- c(cur_curved_ind, cur_curved_ind_)

        cur_len <- length(cur_curved_ind_)

        stat_matrix[ , curved_ind_[1:cur_len]] <- cur_sim_stat[ , cur_curved_ind_] 
      }
      stat_matrix[ , -curved_ind] <- cur_sim_stat[ , -cur_curved_ind]
    } else {
      stat_matrix[ , dims[cur_ind, ]] <- cur_sim_stat
    }

    # Impute observed sufficient statistics if missing data
    if (obj$net$na_flag) {
      # If there are missing edges in the network perform a simulation
      if (obj$net$na_clust_flag[cur_ind]) {
        cur_sim_cond <- simulate(form,
                                 output = "stats",
                                 sequential = obj$sim$seq, 
                                 control = control.simulate(MCMC.burnin = obj$sim$burnin,
                                                            MCMC.interval = obj$sim$interval),
                                 nsim = obj$sim$num_obs,
                                 coef = cur_theta,
                                 constraints = ~ observed)

        if (is.curved(obj$net$model)) {
          for (cur_t in 1:num_curved) { 
            curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
            cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to

            cur_len <- length(cur_curved_ind_)

            cond_matrix[ , curved_ind_[1:cur_len]] <- cur_sim_cond[ , cur_curved_ind_]
          }
          cond_matrix[ ,-curved_ind] <- cur_sim_cond[ , -cur_curved_ind]
        } else {
          cond_matrix[ , dims[cur_ind, ]] <- cur_sim_cond
        }

      } else {
      # If there are no missing edges in the network, the network is fixed
        cur_sim_cond <- rep_row(summary(form), obj$sim$num_obs)

        if (is.curved(obj$net$model)) {
          for (cur_t in 1:num_curved) { 
            curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
            cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to

            cur_len <- length(cur_curved_ind_)

            cond_matrix[ , curved_ind_[1:cur_len]] <- cur_sim_cond[ , cur_curved_ind_]
          }
          cond_matrix[ , -curved_ind] <- cur_sim_cond[ , -cur_curved_ind]
        } else {
          cond_matrix[ , dims[cur_ind, ]] <- cur_sim_cond
        }
      }
      if (obj$est$parameterization == "size") {
        cond_matrix <- cond_matrix * log_fun(obj$est$sizes[cur_ind])
      }
    } else {
      obs_ <- summary(form)
  
      if (is.curved(obj$net$model)) {
        for (cur_t in 1:num_curved) { 
          curved_ind_ <- obj$net$model$etamap$curved[[cur_t]]$to
          cur_curved_ind_ <- mod_temp$etamap$curved[[cur_t]]$to

          cur_len <- length(cur_curved_ind_)

          obs_stats[curved_ind_[1:cur_len]] <- obs_[cur_curved_ind_]
        }
        obs_stats[-curved_ind] <- obs_[-cur_curved_ind]
      } else {
        obs_stats[dims[cur_ind, ]] <-  obs_
      }

      if (obj$est$parameterization == "size") {
        obs_stats <- obs_stats * log_fun(obj$est$sizes[cur_ind])
      }
    }
    if (obj$net$na_flag) { 
      stat_list[[ll]] <- list(stat_matrix = stat_matrix, cond_matrix = cond_matrix)
    } else { 
      stat_list[[ll]] <- list(stat_matrix = stat_matrix, obs_stats = obs_stats)
    }
  }
  return(stat_list)
}
