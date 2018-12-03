
MCMLE <- function(obj) {

  # Run the estimation procedure until termination criterions are met
  conv_flag_1 <- FALSE
  conv_flag_2 <- FALSE
  obj$est$inCH_counter <- 0
  while (!obj$est$MCMLE_status & (obj$est$MCMLE_iter <= obj$est$MCMLE_max_iter) & !obj$est$ML_status_fail) {

    if (obj$verbose > 0) {
      if (obj$est$MCMLE_iter == 1) {  
        cat("\n")
      } else { 
        cat("\n\n")
      }
      cat(paste0("MCMLE Iteration: ", obj$est$MCMLE_iter))
    }

    # Get MCMC sample
    if (obj$verbose > 0) { 
      cat("\n")
      cat("  Sampling networks:")
    } 

    if (obj$est$par_flag) {
      obj <- MCMC_sample_par(obj)
    } else {
      obj <- MCMC_sample(obj)
    }

    # Check that the distribution is not degenerate
    degen_flag <- FALSE 
    sd_ <- lapply(obj$sim$stats, function(x) { apply(x, 2, sd) })
    sd_check <- matrix(0, nrow = obj$net$num_clust, ncol = obj$net$num_terms)
    for (i in 1:obj$net$num_clust) { 
      sd_check[i, ] <- sd_[[i]]
    }
    sd_check_sum <- colSums(sd_check)
    sd_check_sum <- ergm.etagrad(obj$est$theta, obj$net$etamap) %*% sd_check_sum
    if (sum(sd_check_sum == 0) > 0) { 
      warning("Standard deviation of some simulated statistics are zero. Proposed model may be nearly degenerate.")
      degen_flag <- TRUE
    }
    
    if (!degen_flag) { 
      # Use stepping algorithm to put observed vector within the convex hull of the sim sufficient statistics 
      obj <- step_to_chull(obj)

      if (obj$est$gamma == 1) { 
        obj$est$inCH_counter <- obj$est$inCH_counter + 1
      }
  
      # Optimize log-likelihood approximation
      if (obj$verbose > 0) { 
        cat("\n\n  Optimizing likelihood:")
      }
      obj <- opt_fisher(obj)
  
      # If two consecutive samples have been within the convex hull, accept convergence 
      if (obj$est$inCH_counter == 2) { 
        obj$est$MCMLE_status <- TRUE
      }
   
      if (obj$est$inCH_counter < 2) { 
        obj$est$MCMLE_iter <- obj$est$MCMLE_iter + 1
      }
        
      if ((obj$est$MCMLE_iter <= obj$est$MCMLE_max_iter) & !obj$est$ML_status_fail &
          !obj$est$max_iter_flag) { 
        if (obj$verbose > 0) { 
          cat("\n         - Current estimate: ")
          cat(paste0("\n            ",
                     obj$net$theta_names,
                     " = ",
                     round(obj$est$theta, digits = 4)))
          
          score_val <- sum(abs(obj$est$score_val^2))
          print_1 <- ifelse(score_val < 0.000001, "<.000001", as.character(round(score_val, digits = 6)))
          cat("\n\n         - L2-norm of score function at estimate: ", print_1)
        }
      } else if (obj$est$MCMLE_status) {
        # Bypass this 
      } 
    } else {
      obj$est$ML_status_fail <- TRUE  
    }
  }
  return(obj)
}



