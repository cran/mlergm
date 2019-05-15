#' Evaluate the goodness-of-fit of an estimated model. 
#' 
#' Performs a Goodness-of-Fit procedure along the lines of 
#' Hunter, Goodreau, and Handcock (2008). Statistics are simulated 
#' from an fitted \code{mlergm} object, which can then be plotted and visualized.  
#' An example is given in the documentation of \code{\link{mlergm}}.
#'
#' @param object An object of class \code{\link{mlergm}}, likely produced by function \code{\link{mlergm}}.  
#' @param \dots Additional arguments to be passed if necessary. 
#' @param options See \code{\link{set_options}} for details.   
#' @param seed A seed to be provided to ensure reproducibility of results. 
#' @param gof_form A formula object of the form \code{~ term1 + term2 + ...} for statistics to compute for the GOF procedure. 
#' @name gof.mlergm 
#' @return 
#' \code{\link{gof.mlergm}} returns an object of class \code{gof_mlergm} which is a list containing:
#' \item{obs_stats}{The GOF statistic values of the observed network.}
#' \item{gof_stats}{The GOF statistic values simulated from the the estimated \code{mlergm} objeect provided.}
#' 
#' @references
#'
#' Hunter, D. R.,  Goodreau, S. M., and Handcock, M. S. (2008). 
#' Goodness of fit of social network models. 
#' Journal of the American Statistical Association, 103(481), 248-258.
#'
#' @seealso \code{\link{plot.gof_mlergm}}
#' @keywords estimation 
#' @export 
gof.mlergm <- function(object, ..., options = set_options(), seed = NULL, gof_form = NULL) { 
 
  if (!is.mlergm(object)) { 
    stop("Argument is not an 'mlergm' object. See 'help(mlergm') for details.\n", call. = FALSE)
  }

  if (missing(object)) { 
    stop("\nArgument 'object' not provided. An object of class 'mlergm' must be provided.\n", call. = FALSE)
  }

  if (!is.null(seed)) {
    check_integer(seed, "seed")
    set.seed(seed, "L'Ecuyer")
  }


  # Extract network 
  memb_list <- check_and_convert_memb(object$node_memb) 
    memb_internal <- memb_list$memb_internal
    rm(memb_list)
  net_list <- make_net_list(object$net, memb_internal)
    rm(memb_internal); clean_mem()
  max_mem <- max(table(object$node_memb))



  if (!is.null(gof_form)) { 
    if (length(as.character(gof_form)) != 2) { 
      msg <- "\nArgument 'gof_form' must be of the form '~ term1 + term2 ...' to be valid." 
      msg <- paste(msg, "Please see 'ergm.terms' for descriptions of valid term specifications.\n")
      stop(msg, call. = FALSE)
    }
    test_net <- net_list[[1]]
    test <- summary(as.formula(paste0("test_net ~ ", as.character(gof_form)[2])))
    dim_num <- length(test) 
    rm(test_net); rm(test); clean_mem()
  } else { 
    if (is.directed(object$network)) { 
      gof_form <- paste0("~ ", 
                         "esp(0:", max_mem - 2, ") + ",
                         "dsp(0:", max_mem - 2, ") + ", 
                         "odegree(0:", max_mem - 1, ") + ",
                         "idegree(0:", max_mem - 1, ") + ",
                         "triangles")
      gof_form <- as.formula(gof_form)
      test_net <- net_list[[1]] 
      test <- summary(as.formula(paste0("test_net ~ ", as.character(gof_form)[2])))
      dim_num <- length(test) 
      rm(test_net); rm(test); clean_mem()
    } else { 
      gof_form <- paste0("~ ", 
                         "esp(0:", max_mem - 2, ") + ",
                         "dsp(0:", max_mem - 2, ") + ",
                         "degree(0:", max_mem - 1, ") + ",
                         "triangles")
      gof_form <- as.formula(gof_form)
      test_net <- net_list[[1]] 
      test <- summary(as.formula(paste0("test_net ~ ", as.character(gof_form)[2])))
      dim_num <- length(test) 
      rm(test_net); rm(test); clean_mem() 
    }
  }
   if (.Platform$OS.type == "unix") { 
    gof_results <- mclapply(net_list, 
                            gof_slave_fun,
                            form = object$formula,
                            theta = object$theta, 
                            burnin = options$sim_param$burnin, 
                            interval = options$sim_param$interval, 
                            sample_size = options$sim_param$num_obs,
                            gof_form = gof_form, 
                            dim_num = dim_num,
                            max_mem = max_mem, 
                            mc.cores = options$est_param$par_n_cores)
  } else { 
    cl <- makeCluster(options$est_param$par_n_cores) 
    clusterEvalQ(cl, library(mlergm))
    gof_results <- parLapply(cl, 
                             net_list, 
                             gof_slave_fun,
                             form = object$formula, 
                             theta = object$theta, 
                             burnin = options$sim_param$burnin, 
                             interval = options$sim_param$interval, 
                             sample_size = options$sim_param$num_obs, 
                             gof_form = gof_form, 
                             dim_num = dim_num, 
                             max_mem = max_mem)
    stopCluster(cl)
  }


  gof_stats <- rep(list(NULL), length(gof_results))
  obs_stats <- rep(list(NULL), length(gof_results))
  for (i in 1:length(gof_results)) { 
    gof_stats[[i]] <- gof_results[[i]]$gof_stats
    obs_stats[[i]] <- gof_results[[i]]$obs_stats 
  }
  gof_stats <- Reduce("+", gof_stats) 
  obs_stats <- Reduce("+", obs_stats) 

  different_stats <- unique(str_replace_all(colnames(gof_stats), "[0-9]", ""))
  return_list <- list(gof_stats = gof_stats,
                      obs_stats = obs_stats,
                      stat_names = different_stats)
  class(return_list) <- "gof_mlergm"

  return(return_list)
}




# Slave function to compute GOF on individual subnetworks  
gof_slave_fun <- function(net, form, theta, burnin, interval, sample_size, gof_form, dim_num, max_mem) { 
  
  # Simulate networks 
  form <- as.formula(paste0("net ~ ", as.character(form)[3]))
  sim_nets <- simulate(form, 
                       nsim = sample_size, 
                       control = control.simulate(MCMC.burnin = burnin, MCMC.interval = interval), 
                       coef = theta,
                       monitor = gof_form) 

  stat_dim <- ncol(attr(sim_nets, "stats"))
  gof_stats_ind <- (stat_dim - dim_num + 1):stat_dim
  gof_stats <- attr(sim_nets, "stats")[ , gof_stats_ind] 

  geo_stat <- matrix(0, nrow = sample_size, ncol = max_mem - 1)
  geo_net_stats  <- lapply(sim_nets, function(x) { ergm.geodistdist(x) })
  geo_net_stats <- Reduce(rbind, geo_net_stats); rownames(geo_net_stats) <- NULL
  geo_stat[ , 1:(network.size(net) - 1)] <- geo_net_stats[ , -ncol(geo_net_stats)]
  colnames(geo_stat) <- paste0("geo", 1:ncol(geo_stat))
  gof_stats <- cbind(gof_stats, geo_stat)

  obs_form <- as.formula(paste0("net ~ ", as.character(gof_form)[2]))
  obs_stats <- summary(obs_form)
  obs_geo <- numeric(max_mem - 1)
  geo_comp <- ergm.geodistdist(net)
  obs_geo[1:(network.size(net) - 1)] <- geo_comp[-length(geo_comp)]
  names(obs_geo) <- paste0("geo", 1:length(obs_geo))
  obs_stats <- c(obs_stats, obs_geo)

  return_list <- list(obs_stats = obs_stats, gof_stats = gof_stats)
  
  return(return_list)
}



