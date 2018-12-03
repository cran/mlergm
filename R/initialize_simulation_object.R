
initialize_simulation_object <- function(options, net_list, edge_loc, mutual_loc) { 
  
  sim_obj <- list(burnin = options$sim_param$burnin, 
                  interval = options$sim_param$interval, 
                  net_list = net_list,
                  edge_loc = edge_loc, 
                  mutual_loc = mutual_loc,
                  par_n_cores = options$est_param$par_n_cores)

  return(sim_obj)

}
