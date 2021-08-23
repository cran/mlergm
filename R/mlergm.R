#' Multilevel Exponential-Family Random Graph Models
#'
#'  This function estimates an exponential-family random graph model for multilevel network data. At present, \code{mlergm} covers network data where the set of nodes is nested within known blocks (see, e.g., Schweinberger and Handcock, 2015). An example is groups of students nested within classrooms, which is covered in the \code{\link{classes}} data set. It is assumed that the node membership, that to which block each node is associated, is known (or has been previously estimated).
#'
#' The estimation procedures performs Monte-Carlo maximum likelihood for the specified ERGM using a version of the Fisher scoring method detailed by Hunter and Handcock (2006). Settings governing the MCMC procedure (such as \code{burnin}, \code{interval}, and \code{sample_size}) as well as more general settings for the estimation procedure can be adjusted through \code{\link{set_options}}. The estimation procedure uses the the stepping algorithm of Hummel, et al., (2012) for added stability. 
#'
#' @param form Formula of the form:  \code{network ~ term1 + term2 + ...}; allowable model terms are a subset of those in R package ergm,
#' see \code{\link{ergm.terms}}.  
#' @param node_memb Vector (length equal to the number of nodes in the network) indicating to which  block or group the nodes belong.  
#' If the network provided in \code{form} is an object of class \code{mlnet}, 
#' then \code{node_memb} can be exctracted directly from the network and need not be provided. 
#' @param parameterization Parameterization options include 'standard', 'offset', or 'size'. 
#' \itemize{
#' \item 'standard' : Does not adjust the individual block parameters for size. 
#' \item 'offset' : The offset parameterization uses edge and mutual offsets along the lines of Krivitsky, Handcock, and Morris (2011) and Krivitsky and Kolaczyk (2015). The edge parameter is offset by \eqn{-log n(k)} and the mutual parameter is offset by \eqn{+log n(k)}, where \eqn{n(k)} is the size of the kth block.  
#' \item 'size' : Multiplies the block parameters by \eqn{log n(k)}, where \eqn{n(k)} is the size of the kth block.
#' }
#' @param options See \code{\link{set_options}} for details. 
#' @param theta_init Parameter vector of initial estimates for theta to be used. 
#' @param verbose Controls the level of output. A value of \code{0} corresponds to no output, except for warnings; a value of \code{1} corresponds to minimal output, and a value of \code{2} corresponds to full output. 
#' @param eval_loglik (Logical \code{TRUE} or \code{FALSE}) If set to \code{TRUE}, the bridge estimation procedure of Hunter and Handcock (2006) is used to estimate the loglikelihood for BIC calculations, otherwise the loglikelihood and therefore the BIC is not estimated.   
#' @param seed For reproducibility, an integer-valued seed may be specified.
#' @return 
#' \code{\link{mlergm}} returns an object of class \code{\link{mlergm}} which is a list containing:
#' \item{theta}{Estimated parameter vector of the exponential-family random graph model.}
#' \item{between_theta}{Estimated parameter vector of the between group model.}
#' \item{se}{Standard error vector for theta.}
#' \item{between_se}{Standard error vector for between_theta.}
#' \item{pvalue}{A vector of p-values for the estimated parameter vector.}
#' \item{between_pvalue}{A vector of p-values for the estimated parameter vector.}
#' \item{logLikval}{The loglikelihood for at the estimated MLE.}
#' \item{bic}{The BIC for the estimated model.} 
#' \item{mcmc_chain}{The MCMC sample used in the final estimation step, which can be used to diagnose non-convergence.}
#' \item{estimation_status}{Indicator of whether the estimation procedure had \code{succcess} or \code{failed}.}
#' \item{parameterization}{The model parameterization (either \code{standard} or \code{offset}).}
#' \item{formula}{The model formula.}
#' \item{network}{The network for which the model is estimated.}
#' \item{node_memb}{Vector indicating to which group or block the nodes belong.}
#' \item{size_quantiles}{The quantiles of the block sizes.}
#'
#' @references
#'
#' Schweinberger, M. and Stewart, J. (2019) 
#' Concentration and consistency results for canonical and curved exponential-family random graphs. 
#' The Annals of Statistics, to appear.
#'
#'
#' Schweinberger, M. and Handcock, M. S. (2015).
#' Local dependence in random graph models: characterization, properties and statistical inference.
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(3), 647-676.
#'
#' Hunter, D. R., and Handcock, M. S. (2006). 
#' Inference in curved exponential family models for networks. 
#' Journal of Computational and Graphical Statistics, 15(3), 565-583.
#'
#' Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012). 
#' Improving simulation-based algorithms for fitting ERGMs. 
#' Journal of Computational and Graphical Statistics, 21(4), 920-939.
#'
#' Krivitsky, P. N., Handcock, M. S., & Morris, M. (2011). 
#' Adjusting for network size and composition effects in exponential-family random graph models. 
#' Statistical methodology, 8(4), 319-339.
#'
#' Krivitsky, P.N, and Kolaczyk, E. D. (2015). 
#' On the question of effective sample size in network modeling: An asymptotic inquiry. 
#' Statistical science: a review journal of the Institute of Mathematical Statistics, 30(2), 184.
#'
#' Hunter D., Handcock M., Butts C., Goodreau S., and Morris M. (2008).
#' ergm: A Package to Fit, Simulate and Diagnose Exponential-Family Models for Networks.
#' Journal of Statistical Software, 24(3), 1-29.
#'
#' Butts, C. (2016).
#' sna: Tools for Social Network Analysis.
#' R package version 2.4. \url{https://CRAN.R-project.org/package=sna}.
#'
#' Butts, C. (2008).
#' network: a Package for Managing Relational Data in R.
#' Journal of Statistical Software, 24(2). 
#'
#' Stewart, J., Schweinberger, M., Bojanowski, M., and M. Morris (2019). 
#' Multilevel network data facilitate statistical inference for curved {ERGM}s with geometrically weighted terms. 
#' Social Networks, 59, 98-119.
#'
#' Schweinberger, M., Krivitsky, P. N., Butts, C.T. and J. Stewart (2018). 
#' Exponential-family models of random graphs: Inference in finite-, super-, and infinite-population scenarios. 
#' https://arxiv.org/abs/1707.04800
#'
#' @seealso \code{\link{gof.mlergm}}, \code{\link{mlnet}}
#' @keywords estimation
#' @export 
#' @importFrom stats median sd as.formula update simulate update.formula pnorm quantile coef  
#' @importFrom parallel stopCluster mclapply makeCluster clusterEvalQ clusterApply parLapply 
#' @importFrom Matrix bdiag
#' @importFrom stringr str_match str_split str_trim str_replace_all 
#' @importFrom cowplot plot_grid 
#' @importFrom reshape2 melt
#' @importFrom plyr is.formula 
#' @importFrom methods is
#' @importFrom graphics plot 
#' @import ergm 
#' @import network
#' @examples 
#' \donttest{
#' ### Load the school classes data-set 
#' data(classes) 
#'
#' # Estimate a curved multilevel ergm model with offset parameter 
#' # Approximate run time (2 cores): 1.2m, Run time (3 cores): 55s 
#' model_est <- mlergm(classes ~ edges + mutual + nodematch("sex") +  gwesp(fixed = FALSE), 
#'                     seed = 123, 
#'                     options = set_options(number_cores = 2))
#'
#' # To access a summary of the fitted model, call the 'summary' function 
#' summary(model_est)
#'
#' # Goodness-of-fit can be run by calling the 'gof.mlergm' method 
#' # Approximate run time (2 cores): 48s, Run time (3 cores): 34s  
#' gof_res <- gof(model_est, options = set_options(number_cores = 2))
#' plot(gof_res, cutoff = 15)
#' } 
mlergm <- function(form,
                   node_memb,
                   parameterization = "standard",
                   options = set_options(),  
                   theta_init = NULL,
                   verbose = 0,
                   eval_loglik = TRUE,
                   seed = NULL) {

 

  # Check that required arguments are provided 
  if (missing(form)) {
    stop("\nArgument 'form' not provided. A formula object must be provided.\n", call. = FALSE)
  } else { 
    check_formula(form)
  }
  if (missing(node_memb)) {
    # Check if network provided is of class 'mlnet' 
    net <- get_network_from_formula(form) 
    if (inherits(net, "mlnet")) { 
      node_memb <- get.vertex.attribute(net, "node_memb")
    } else { 
      stop("\nArgument 'node_memb' not provided. The node memberships must be provided.\n", call. = FALSE)
    }
  }
  # Check that the formula and terms requested are valid
  check_terms(form, K = length(unique(node_memb))) 
  # If a seed is provided, set it
  if (!is.null(seed)) { 
    check_integer(seed, "seed")
    set.seed(seed, "L'Ecuyer")
  } 

  # Adjust formula if necessary 
  form <- adjust_formula(form) 
  

  # Parse formula to get network and model
  net <- get_network_from_formula(form)
  check_net(net)
  terms <- get_terms_from_formula(form, net)
 
  if (verbose > 0) { 
    cat("\nBeginning estimation procedure for:")
    cat(paste0("\n       ", terms))
  }  

  # Since net_list was not provided, we need to make it 
  memb_list <- check_and_convert_memb(node_memb)
    memb_labels <- memb_list$memb_labels
    memb_internal <- memb_list$memb_internal 
    rm(memb_list); clean_mem()
  net_list <- make_net_list(net, memb_internal) 
  

  # Determine the parameterization type and compute necessary quantities
  model_temp <- ergm_model(form, net)
  param_list <- check_parameterization_type(net_list, terms, parameterization, model_temp)
  statistic_names <- param_list$statistic_names
  which_largest <- param_list$which_largest

  # Initialize object
  obj <- initialize_object(net_list = net_list,
                           net = net,
                           form = form, 
                           terms = terms,  
                           block_memb = memb_internal, 
                           theta_init = theta_init,
                           param_list = param_list, 
                           sim_param = options$sim_param,
                           est_param = options$est_param,
                           verbose = verbose,
                           parameterization = parameterization)
  
  # Remove objects that are no longer needed 
  rm(param_list)
  rm(options)
  clean_mem()
  
  # Initialize estimate if an initial estimate is not provided 
  cd_flag <- TRUE
  if (is.null(obj$est$theta)) {
    if (verbose > 0) { 
      cat("\n\nComputing initial estimate.")  
    }
    obj <- compute_initial_estimate(obj)
    if (verbose > 0) {
      cat("\n    Initial estimate:")
      cat(paste("\n     ", obj$net$theta_names, " = ", formatC(obj$est$theta, digits = 4, format = "f")))
    }
  } else { 
    obj$est$theta <- theta_init
  }



  # Call MCMLE to perform estimation
  if (verbose > 0) { 
    cat("\n\n\nBeginning Monte-Carlo maximum likelihood estimation\n")
    cat("===================================================")
    cat("\n")
  }
  obj <- MCMLE(obj)
 
  # Estimate the between block model (if possible)
  obj <- estimate_between_block(obj) 


  # Evalaute results of MCMLE procedure 
  if (verbose > 0 & !obj$est$ML_status_fail) {
    cat(paste0("\n\nComputing approximate loglikelihood at estimate using ", 
               obj$est$bridge_num, " bridges."))
    cat("\n\n")

  } else if (verbose > 0 & obj$est$ML_status_fail) {
    cat("\n\nEstimation procedure stopping. Estimation unsuccesful.\n\n", call. = FALSE)
  }


  # Compute quantiles of block sizes 
  quantiles_of_sizes <- quantile(obj$net$clust_sizes)

  # Make structure to be returned
  if (!obj$est$ML_status_fail) {
    if (obj$est$inCH_counter > 0) {
      obj <- compute_pvalue(obj)
      if (eval_loglik) { 
        obj$likval <- lik_fun(form = form, memb = node_memb, theta = obj$est$theta, 
                              bridge_num = obj$est$bridge_num, ncores = obj$est$par_n_cores,
                              offset = obj$est$parameterization == "offset",
                              burnin = obj$sim$bridge_burnin, 
                              interval = obj$sim$bridge_interval, 
                              sample_size = obj$sim$bridge_sample_size,
                              size = obj$est$parameterization == "size") 
        obj$bic <- compute_bic(obj) 
      } else { 
        obj$likval <- NULL
        obj$bic <- NULL
      }
      mcmc_path <- Reduce("+", obj$sim$stats)
      colnames(mcmc_path) <- statistic_names
      obj$est$theta <- as.numeric(obj$est$theta)
      names(obj$est$theta) <- get_coef_names(obj$net$model, FALSE)
      names(obj$se) <- names(obj$est$theta)
      names(obj$pvalue) <- names(obj$est$theta)
      estimates <- list(theta = obj$est$theta,
                        between_theta = obj$est$between_theta, 
                        between_se = obj$est$between_se,
                        se = obj$se,
                        pvalue = obj$pvalue, 
                        between_pvalue = obj$est$between_pvalue,  
                        bic = obj$bic, 
                        logLikval = obj$likval,  
                        mcmc_chain = mcmc_path,
                        estimation_status = "success",
                        parameterization = obj$est$parameterization,
                        formula = form,
                        network = net,
                        node_memb = node_memb,
                        size_quantiles = quantiles_of_sizes)
      class(estimates) <- "mlergm" 
      rm(mcmc_path); clean_mem()
    } else if (obj$est$inCH_counter == 0) { 
      cat("\n\nWarning: Maximum number of iterations reached without the observation lying in the")
      cat(" interior of the simulated convex hull. Parameters not estimated.\n\n")
      estimates <- list(theta = NA,
                        se = NA,
                        formula = form, 
                        network = net,
                        node_memb = node_memb, 
                        size_quantiles = quantiles_of_sizes,
                        mcmc_chain = NULL,
                        estimation_status = "failed")
      class(estimates) <- "mlergm"
    }
  } else {
    estimates <- list(theta = NA, 
                      se = NA,
                      formula = form, 
                      network = net, 
                      node_memb = node_memb,
                      size_quantiles = quantiles_of_sizes, 
                      mcmc_chain = NULL,
                      estimation_status = "failed")
    class(estimates) <- "mlergm" 
  }

  if (verbose > 0) { 
    summary(estimates) 
  }
  return(estimates)
}
