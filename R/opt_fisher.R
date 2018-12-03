opt_fisher <- function(obj) {

  # Reset iter counters and status
  obj$est$NR_iter <- 1
  obj$est$NR_status <- FALSE
  obj$est$max_iter_flag <- FALSE
  obj$est$step_err <- 1
  obj$est$NR_step_len_ <- obj$est$NR_step_len

  # Compute iterations while termination criterion is not met
  if (obj$verbose > 1) {
    cat("\n")
  }
  
  if (any(is.nan(ergm.eta(obj$est$theta, obj$net$etamap)))) {
    obj$est$ML_status_fail <- TRUE
  }
  
  while (!obj$est$NR_status & (obj$est$NR_iter <= obj$est$NR_max_iter) & !obj$est$ML_status_fail) {
    
    # Compute weights
    if (obj$verbose > 1) {
      cat("\n    Computing step.")
    }
    weights  <- comp_weights(obj)
    if (weights == "broken") { 
      obj$est$ML_status_fail <- TRUE
    }

    # Compute information matrix
    if (!obj$est$par_flag & (obj$est$ML_status_fail == FALSE)) {
      info_mat <- comp_info_mat(obj, weights)
    } else if (obj$est$par_flag & (obj$est$ML_status_fail == FALSE)) {
      info_mat <- comp_info_mat_par(obj, weights)
    }
    if (obj$est$ML_status_fail == FALSE) { 
      obj$est$info_mat <- info_mat
    }

    # Compute the step
    if (obj$est$ML_status_fail == FALSE) { 
      step_list <- comp_step(obj, weights, info_mat)
      step <- step_list$step
      if (!is.null(step_list$ML_status_fail)) {
        obj$est$ML_status_fail <- step_list$ML_status_fail
      }
      if (obj$est$adaptive_step_len) {
        obj$est$NR_step_len_ <- 1 / (1 + sum(abs(step)^2))
      }
      if (any(is.nan(ergm.eta(obj$est$theta + step * obj$est$NR_step_len_, obj$net$etamap)))) {
        obj$est$ML_status_fail <- TRUE
      }
    }

    if (!obj$est$ML_status_fail) {
      if (obj$verbose > 1) {
        val <- sum(abs(step))
        print_1 <- val
        print_1[val < 0.0001] <- "<.0001"
        print_1[val >= 0.0001] <- formatC(val[val >= 0.0001], digits = 4, format = "f")
        cat(paste0("\n      L1-norm of increment of parameter vector: ", print_1))
      }
      obj$est$score_val <- step_list$score_val

      # Check if lowest point
      if (obj$est$NR_iter > 1) { 
        if (norm(obj$est$score_val, type = "2") < norm(obj$est$score_min, type = "2")) { 
          obj$est$theta_min <- obj$est$theta
          obj$est$score_min <- obj$est$score_val
        }
      } else { 
        obj$est$theta_min <- obj$est$theta
        obj$est$score_min <- obj$est$score_val
      }

      # Iterate the next step
      obj$est$theta <- obj$est$theta + step * obj$est$NR_step_len_
      if (any(is.nan(ergm.eta(obj$est$theta, obj$net$etamap)))) {
        obj$est$ML_status_fail <- TRUE
      }

      if (obj$verbose > 1) {
        cat(paste0("\n      Iteration: ", 
                   obj$est$NR_iter, 
                   ",  ", 
                   obj$net$theta_names, 
                   " = ", 
                   round(obj$est$theta, digits = 4)))
        cat("\n")
      }

      # Update status and iterations
      if (!obj$est$ML_status_fail) {
        obj$est$NR_iter <- obj$est$NR_iter + 1
        obj <- check_convergence(obj, step)
      }

    } else {
      if (obj$est$step_err < 3 & !obj$est$adaptive_step_len) {
        if (obj$verbose > 1) {
          cat("\n    Optimization did not converge. Decreasing step length and restarting.\n")
        }
        obj$est$step_err <- obj$est$step_err + 1
        obj$est$NR_step_len_ <- obj$est$NR_step_len_ * obj$est$NR_step_len_multiplier 
        obj$est$theta <- obj$est$theta_0
        obj$est$NR_iter <- 1
        obj$est$ML_status_fail <- FALSE
      } else { 
        cat("\n    Optimization failed to converge.")
        cat("\n    Proposed model may be near-degenerate or the MCMLE may be unstable or not exist.") 
        cat("\n    Decreasing the step length (argument 'NR_step_len' in set_options())") 
        cat(" or increasing the MCMC sample-size may help.")
      }
    }
  }

  if (obj$est$NR_iter >= obj$est$NR_max_iter) {
    if (obj$verbose > 0) { 
      cat("\n    NOTE: Optimization reached maximum number of allowed iterations.")
      cat(paste0("\n\n      - Minimum theta found:"))
      cat(paste0("\n        ",
                 obj$net$theta_names,
                 " = ",
                 round(obj$est$theta_min, digits = 4)))
      
      cat("\n\n")
      cat("      - L2-norm of score at this theta: ")
      cat(paste(round(sum(abs(obj$est$score_min^2)), digits = 4)))
    }
    obj$est$max_iter_flag <- TRUE
  }
  obj$est$theta_0 <- obj$est$theta
  obj$est$NR_conv_thresh <- sqrt(sum(step^2))
  return(obj)
}





#### comp_weights
comp_weights <- function(obj) {

  if (obj$net$na_flag) {
    weights <- list(weight_full = rep(list(numeric(obj$sim$num_obs)),
                                      obj$net$num_clust),
                    weight_cond = rep(list(numeric(obj$sim$num_obs)),
                                      obj$net$num_clust))
  } else {
    weights <- list(weight_full = rep(list(numeric(obj$sim$num_obs)),
                                      obj$net$num_clust))
  }

  theta_diff <- ergm.eta(obj$est$theta, obj$net$etamap) - ergm.eta(obj$est$theta_0, obj$net$etamap)
  for (i in 1:obj$net$num_clust) {
    adjust_flag <- FALSE
    if (any(is.na(obj$sim$stats[[i]] %*% theta_diff))) { 
      return("broken")
    }
    norm_const_full <- sum(exp_fun(obj$sim$stats[[i]] %*% theta_diff))
    if (norm_const_full == Inf) {
      norm_const_full <- .Machine$double.xmax
      adjust_flag <- TRUE
    }
    if (norm_const_full == 0) {
      norm_const_full <- .Machine$double.xmin
      adjust_flag <- TRUE
    }
    weights$weight_full[[i]] <- exp_fun(obj$sim$stats[[i]] %*% theta_diff) /
                                   norm_const_full
    if (adjust_flag) {
      weights$weight_full[[i]] <- weights$weight_full[[i]] / sum(weights$weight_full[[i]])
    }
    if (obj$net$na_flag) {
      if (any(is.na(obj$sim$cond_stats[[i]] %*% theta_diff))) { 
        return("broken")
      }
      norm_const_cond <- sum(exp_fun(obj$sim$cond_stats[[i]] %*% theta_diff))
      weights$weight_cond[[i]] <- exp_fun(obj$sim$cond_stats[[i]] %*% theta_diff) /
                                    norm_const_cond
    }
  }
  return(weights)
}




#### comp_info_mat
comp_info_mat <- function(obj, weights) {
  info_mat <- rep(list(NULL), obj$net$num_clust)
  theta_grad_val <- t(ergm.etagrad(obj$est$theta, obj$net$etamap))
  for (i in 1:obj$net$num_clust) {
    term_1_full <- numeric(obj$net$num_terms)
    term_2_full <- numeric(obj$net$num_terms)

    if (obj$net$na_flag) {
      term_1_cond <- numeric(obj$net$num_terms)
      term_2_cond <- numeric(obj$net$num_terms)
    }

    list_ <- as.list(data.frame(t(obj$sim$stats[[i]])))
    outers_ <- lapply(list_, outer_fun)
    term_1_full <- Reduce("+", Map("*", outers_, weights$weight_full[[i]]))
    term_2_full <-  as.vector(t(obj$sim$stats[[i]]) %*% weights$weight_full[[i]])

    J_full <- t(theta_grad_val) %*%
              (term_1_full - outer(term_2_full, term_2_full)) %*%
              theta_grad_val

    if (obj$net$na_flag) {
      list_ <- as.list(data.frame(t(obj$sim$cond_stats[[i]])))
      outers_ <- lapply(list_, outer_fun)
      term_1_cond <- Reduce("+", Map("*", outers_, weights$weight_cond[[i]]))
      term_2_cond <-  as.vector(t(obj$sim$cond_stats[[i]]) %*%
                                  weights$weight_cond[[i]])

      J_cond <- t(theta_grad_val) %*%
                (term_1_cond - outer(term_2_cond, term_2_cond)) %*%
                theta_grad_val
    }

    if (obj$net$na_flag) {
      info_mat[[i]] <- J_full - J_cond
    } else {
      info_mat[[i]] <- J_full
    }
  }
  return(Reduce("+", info_mat))
}


comp_info_mat_par <- function(obj, weights) {

  theta_grad_val <- t(ergm.etagrad(obj$est$theta, obj$net$etamap))

  par_fun <- function(i, obj, weights, theta_grad_val, job_split) {
    split_ <- job_split[[i]]
    info_mat_list <- rep(list(NULL), length(split_))
    for (ll in 1:length(split_)) {
      cur_ind <- split_[ll]
      term_1_full <- numeric(obj$net$num_terms)
      term_2_full <- numeric(obj$net$num_terms)

      if (obj$net$na_flag) {
        term_1_cond <- numeric(obj$net$num_terms)
        term_2_cond <- numeric(obj$net$num_terms)
      }

      list_ <- as.list(data.frame(t(obj$sim$stats[[cur_ind]])))
      outers_ <- lapply(list_, outer_fun)
      term_1_full <- Reduce("+", Map("*", outers_, weights$weight_full[[cur_ind]]))
      term_2_full <-  as.vector(t(obj$sim$stats[[cur_ind]]) %*% weights$weight_full[[cur_ind]])

      J_full <- t(theta_grad_val) %*%
                (term_1_full - outer(term_2_full, term_2_full)) %*%
                theta_grad_val

      if (obj$net$na_flag) {
        list_ <- as.list(data.frame(t(obj$sim$cond_stats[[cur_ind]])))
        outers_ <- lapply(list_, outer_fun)
        term_1_cond <- Reduce("+", Map("*", outers_, weights$weight_cond[[cur_ind]]))
        term_2_cond <-  as.vector(t(obj$sim$cond_stats[[cur_ind]]) %*%
                                    weights$weight_cond[[cur_ind]])

        J_cond <- t(theta_grad_val) %*%
                  (term_1_cond - outer(term_2_cond, term_2_cond)) %*%
                  theta_grad_val
      }

      if (obj$net$na_flag) {
        info_mat <- J_full - J_cond
      } else {
        info_mat <- J_full
      }
      info_mat_list[[ll]] <- info_mat
    }
    return(info_mat_list)
  }

  # Split the jobs for available cores
  if (obj$est$par_n_cores > obj$net$num_clust) {
    split_num <- obj$net$num_clust
  } else if (obj$net$num_clust >= obj$est$par_n_cores) {
    split_num <- obj$est$par_n_cores
  }
  split_ <- msplit(1:obj$net$num_clust, 1:split_num)
  
  if (.Platform$OS.type == "windows") { 
    cl <- makeCluster(split_num) 
    clusterEvalQ(cl, library(mlergm))
    info_mats <- clusterApply(cl, 1:length(split_), par_fun, obj = obj, weights = weights,
                              theta_grad_val = theta_grad_val, job_split = split_) 
    stopCluster(cl) 
  } else { 
    info_mats <- mclapply(1:length(split_), par_fun,
                          obj, weights, theta_grad_val, split_, mc.cores = split_num)
  }
  return(Reduce("+", unlist(info_mats, recursive = FALSE)))
}




### comp_step
comp_step <- function(obj, weights, info_mat) {
  if (is.nan(norm(info_mat))) {
    ML_status_fail <- TRUE
    cat("\n    Information matrix is near-singular.")
  } else {
    ML_status_fail <- FALSE
  }
  if (!ML_status_fail) {
    info_mat_inv <- tryCatch({ solve(info_mat) },
                             error = function(e) {
                                  return("error")
                             })

    if (is.character(info_mat_inv)) {
      ML_status_fail <- TRUE
      cat("\n    Information matrix is near-singular.")
    } 
  }
  if (!ML_status_fail) {
    theta_grad_val <- t(ergm.etagrad(obj$est$theta, obj$net$etamap))
    exp_approx <- as.vector(Reduce("+", Map("%*%", lapply(weights$weight_full, t),
                                          obj$sim$stats)))
    if (obj$net$na_flag) {
      obs_approx <- as.vector(Reduce("+",
                                  Map("%*%", lapply(weights$weight_cond, t),
                                              obj$sim$cond_stats)))
    } else {
      obs_approx <- obj$net$obs_stats_step
    }

    if ((norm(theta_grad_val) == Inf) | is.nan(norm(theta_grad_val))) {
      ML_status_fail <- TRUE
      cat("\n    Gradient is not finite. Stopping optimization.")
    }

    if (!ML_status_fail) {
      score_val <- t(theta_grad_val) %*% (obs_approx - exp_approx)
      step <- info_mat_inv %*% score_val
    }

  }
  if (ML_status_fail) {
    step <- NA
    score_val <- NA
  }
  return(list(step = step, score_val = score_val, ML_status_fail = ML_status_fail))
}



### check_convergence
check_convergence <- function(obj, step) {

  step_norm <- sqrt(sum(step^2))

  if (!(step_norm == Inf)) {

    conv_check <- sum(abs(step) > obj$est$adj_NR_tol) == 0

    # Check if step is smaller than some tolerance
    if (conv_check) {
      obj$est$NR_status <- TRUE
    }

  } else {
    obj$est$step_err <- obj$est$step_err + 1
    if (obj$est$step_err < 3) {
      obj$est$NR_step_len_ <- obj$est$NR_step_len_ * obj$est$NR_step_len_multiplier
      obj$est$NR_iter <- 1
    } else {
      obj$est$ML_status_fail <- TRUE
    }
  }
  return(obj)
}
