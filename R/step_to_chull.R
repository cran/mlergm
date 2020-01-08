

step_to_chull <- function(obj) { 
  
  obj$net$obs_stats_step <- obj$net$obs_stats
  obj$est$gamma <- 1
  if (obj$verbose > 0) { 
    cat("\n")
  }

  sim_mat <- Reduce("+", obj$sim$stats)
  check_obs <- t(ergm.etagrad(obj$est$theta, obj$net$etamap) %*% obj$net$obs_stats_step)
  check_sim <- t(ergm.etagrad(obj$est$theta, obj$net$etamap) %*% t(sim_mat))
  check_ <- is.inCHv3.9(check_obs, check_sim)
  if (!check_) { 
    step_flag <- TRUE
    iter <- 1
    sim_mean <- apply(Reduce("+", obj$sim$stats), 2, mean)
    if (obj$verbose > 0) {
      cat("    Observed sufficient statistic vector is not in the interior of the simulated convex hull.")
    }
    obj$est$in_chull <- FALSE
    obj$est$adj_NR_tol <- 0.0001
  } else {
    if (obj$verbose > 0) { 
      cat("    Observed sufficient statistic vector is in the interior of the simulated convex hull.")
    } 
    obj$est$in_chull <- TRUE
    obj$est$adj_NR_tol <- obj$est$NR_tol 
    step_flag <- FALSE
  }

  adjust_flag <- FALSE
  max_grid_flag <- FALSE
  zero_flag <- FALSE
  step_grid <- seq(1, 0.000001, length.out = 200)
  while (!check_ & step_flag & !zero_flag) {
    
    if (iter > 200) { 
      obj$est$gamma <- exp(- iter / 4)
      max_grid_flag <- TRUE
    }
    
    if (obj$verbose > 1 & iter == 1) { 
      cat("\n    Stepping towards observed sufficient statistics vector.")
    } 

    if (!max_grid_flag) { 
      obj$est$gamma <- step_grid[iter]
    }
    
    if (obj$est$gamma == 0) {
      if (obj$verbose > 1) {  
        cat("\n\nWARNING: gamma = 0 required to move observation into the simulated convex hull. This could be due to a poor initial estimate, insufficient MCMC parameters indicating the Markov Chain did not mix well, or possible model misspecification.")
      }
      zero_flag <- TRUE
    }
    obj$net$obs_stats_step <- obj$net$obs_stats * 1.05 * obj$est$gamma + (1 - 1.05 * obj$est$gamma) * sim_mean
    check_obs <- t(ergm.etagrad(obj$est$theta, obj$net$etamap) %*% obj$net$obs_stats_step)
    check_ <- suppressMessages(is.inCHv3.9(check_obs, check_sim))
    iter <- iter + 1
  }
  if (step_flag) { 
    obj$net$obs_stats_step <- obj$net$obs_stats * obj$est$gamma + (1 - obj$est$gamma) * sim_mean
  }

  if (obj$verbose > 0) {
    if (step_flag) {
      cat("\n")
      cat("    Using adjusted observed sufficient statistics vector with gamma = ")
      cat(round(obj$est$gamma, digits = 4))
      cat(". L1 norm of difference between observation and adjusted observation: ")
      L1_diff <- sum(abs(obj$net$obs_stats_step - obj$net$obs_stats))
      if (round(L1_diff, digits = 3) < 0.001) { 
        cat("<0.001")
      } else { 
        cat(round(L1_diff, digits = 3))
        cat(".")
      }
    } 
  }
  return(obj)
}


