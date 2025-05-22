#' Set and adjust options and settings. 
#'
#' Function allows for specification of options and settings for simulation and estimation procedures. 
#'
#' The main simulation settings are \code{burnin}, \code{interval}, and \code{sample_size}. For estimation of the loglikelihood value, options include \code{bridge_num} which controls the number of bridges to be used for approximating the loglikelihood (see, e.g., Hunter and Handcock (2006) for a discussion). The main estimation settings and options include \code{NR_tol}, \code{NR_max_iter}, \code{MCMLE_max_iter}, \code{adaptive_step_len}, and \code{step_len}. Parameters \code{NR_tol} and \code{NR_max_iter} control the convergence tolerance and maximum number of iterations for the Newton-Raphson, or Fisher scoring, optimization. When the L2 norm of the incremenet in the Newton-Raphson procedure is under the specified tolerance \code{NR_tol} convergence is reached; and, no more than \code{NR_max_iter} iterations are performed. The MCMLE procedure uses the stepping algorithn of Hummel, et al., (2012) to give stabiity to the estimation procedure. Each MCMLE iteration draws samples from an MCMC chain, and \code{MCMLE_max_iter} controls how many iterations are performed before termination. Most functions support parallel computing for efficiency; by default \code{do_parallel} is \code{TRUE}. The number of computing cores can be adjusted by \code{number_cores}, and the default is one less than the number of cores available. 
#'
#'
#' @param burnin The burnin length for MCMC chains. 
#' @param interval The sampling interval for MCMC chains. 
#' @param sample_size The number of points to sample from MCMC chains for the MCMLE procedure. 
#' @param NR_tol The convergence tolerance for the Newton-Raphson optimization (implemented as Fisher scoring). 
#' @param NR_max_iter The maximum number of Newton-Raphson updates to perform. 
#' @param MCMLE_max_iter The maximum number of MCMLE steps to perform. 
#' @param do_parallel (logical) Whether or not to use parallel processesing (defaults to TRUE). 
#' @param number_cores The number of parallel cores to use for parallel computations. 
#' @param adaptive_step_len (logical) If \code{TRUE}, an adaptive steplength procedure is used 
#' for the Newton-Raphson procedure. Arguments \code{NR_step_len} and \code{NR_step_len_multiplier} 
#' are ignored when \code{adaptive_step_len} is \code{TRUE}.
#' @param step_len_multiplier The step_len adjustment multplier when convergence fails. 
#' @param step_len The step length adjustment default to be used for the Newton-Raphson updates. 
#' @param bridge_num The number of bridges to use for likelihood computations. 
#' @param bridge_burnin The burnin length for the bridge MCMC chain for approximate likelihood computation. 
#' @param bridge_interval The sampling interval for the brdige MCMC chain for approximate likelihood computation. 
#' @param bridge_sample_size The number of points to sample from the bridge MCMC chain for approximate likelihood computation. 
#'
#' @references 
#' Hunter, D. R., and Handcock, M. S. (2006).
#' Inference in curved exponential family models for networks.
#' Journal of Computational and Graphical Statistics, 15(3), 565-583.
#'
#' Hummel, R. M., Hunter, D. R., and Handcock, M. S. (2012).
#' Improving simulation-based algorithms for fitting ERGMs.
#' Journal of Computational and Graphical Statistics, 21(4), 920-939.
#'
#' @export
#' @importFrom parallel detectCores
set_options <- function(burnin = 1e+5,
                        interval = 2000,
                        sample_size = 1000,
                        NR_tol = 1e-4,
                        NR_max_iter = 50,
                        MCMLE_max_iter = 10,
                        do_parallel = TRUE,
                        number_cores = detectCores(all.tests = FALSE, logical = TRUE) - 1,
                        adaptive_step_len = TRUE,
                        step_len_multiplier = 0.5,
                        step_len = 1,
                        bridge_num = 10,
                        bridge_burnin = 1e+4,
                        bridge_interval = 500,
                        bridge_sample_size = 5000) {


  sim_param <- list(burnin = burnin,
                    interval = interval,
                    num_obs = sample_size,
                    stats = NULL,
                    cond_stats = NULL,
                    bridge_burnin = bridge_burnin,
                    bridge_interval = bridge_interval,
                    bridge_sample_size = bridge_sample_size)


  est_param <- list(eta = NULL,
                    eta_0 = NULL,
                    eta_grad = NULL,
                    eta_fun = NULL,
                    score_val = NULL,
                    NR_tol = NR_tol,
                    NR_iter = 1,
                    NR_max_iter = NR_max_iter,
                    NR_status = FALSE,
                    step_err = 0,
                    MCMLE_iter = 1,
                    MCMLE_max_iter = MCMLE_max_iter,
                    MCMLE_status = FALSE,
                    info_mat = NULL,
                    bridge_num = bridge_num,
                    adaptive_step_len = adaptive_step_len,
                    NR_step_len = step_len,
                    NR_step_len_multiplier = step_len_multiplier,
                    NR_conv_thresh = NULL,
                    MCMLE_conv_thresh = NULL,
                    par_flag = do_parallel,
                    par_n_cores = number_cores,
                    ML_status_fail = FALSE)

  return(list(sim_param = sim_param,
              est_param = est_param))
}
