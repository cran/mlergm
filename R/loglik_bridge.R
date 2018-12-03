
bridge_fun <- function(net, form, theta, offset, burnin, interval, num_bridges, sample_size) { 
  form <- as.formula(paste0("net ~ ", as.character(form)[3]))
  model <- ergm_model(form, net)
  etamap <- model$etamap
  coef_names <- get_coef_names(model, !is.curved(model))
  if (offset == TRUE) { 
    if ("edges" %in% coef_names) {
        edge_loc <- which(coef_names == "edges")
        theta[edge_loc] <- theta[edge_loc] - log(network.size(net))
    }
    if ("mutual" %in% coef_names) {
        mutual_loc <- which(coef_names == "mutual")
        theta[mutual_loc] <- theta[mutual_loc] + log(network.size(net))
    }
  }
  bridge_val <- suppressMessages(
                  ergm.bridge.llr(form, 
                                  to = theta, 
                                  from = rep(0, length(theta)), 
                                  llronly = TRUE,
                                  control = control.ergm.bridge(MCMC.samplesize = sample_size, 
                                                                MCMC.interval = interval,
                                                                MCMC.burnin = burnin, 
                                                                nsteps = num_bridges)))

  return(bridge_val) 
}

lik_fun <- function(form, memb, theta, bridge_num = 10, ncores = 3, offset = FALSE, 
                    burnin = NULL, interval = NULL, sample_size = NULL) {


  # Make net_list + compute obs
  network <- ergm.getnetwork(form)
  obs <- summary(form)
  net_list <- rep(list(NULL), length(unique(memb))) 
  u_memb <- unique(memb)
  if (.Platform$OS.type == "unix") { 
    net_list <- mclapply(u_memb, 
                         function(x, network, memb) {  
                           get.inducedSubgraph(network, v = which(memb == x))
                         }, network = network, memb = memb,
                         mc.cores = ncores)
  } else { 
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library(mlergm))
    net_list <- parLapply(cl, 
                          u_memb, 
                          function(x, network, memb) { 
                            get.inducedSubgraph(network, v = which(memb == x))
                          }, network = network, memb = memb)
  }
  terms <- as.character(form)[3]
 
  # Simulate bridges
  if (.Platform$OS.type == "unix") { 
    bridges <- mclapply(net_list,
                        bridge_fun, 
                        theta = theta, offset = offset, form = form, num_bridges = bridge_num, 
                        burnin = burnin, interval = interval, sample_size = sample_size,
                        mc.cores = ncores) 
  } else { 
    bridges <- parLapply(cl, 
                         net_list,
                         bridge_fun, 
                         theta = theta, offset = offset, form = form, num_bridges = bridge_num, 
                         burnin = burnin, interval = interval, sample_size = sample_size)
    stopCluster(cl)
  }

  null_bridge <- Reduce("+", lapply(net_list, 
                                    function(net) {
                                      if (is.directed(net)) {
                                        2 * choose(network.size(net), 2) * log(2)
                                      } else { 
                                        choose(network.size(net), 2) * log(2)
                                      }
                                    }))

  lik_val <- Reduce("+", bridges) - null_bridge
  
  return(lik_val)
}


get_coef_names <- function(model_obj, is_canonical) { 
  if(is_canonical) {
    model_obj$coef.names
  } else { 
    unlist(lapply(model_obj$terms,
                  function(term) {
                    find_first_non_null(names(term$params),  term$coef.names)
                  }))
  }
}


find_first_non_null <- function(...) { 
  for (x in list(...)) {
    if (!is.null(x)) {
      break
    }
  }
  x
}












