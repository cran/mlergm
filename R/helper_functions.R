
###########################
###
### Basic util functions
###
###########################

clean_mem <- function() {
  invisible(gc());
}

norm2 <- function(vec) {
  val <- norm(as.matrix(vec), type = "2")
  return(val)
}

norm1 <- function(vec) {
  val <- norm(as.matrix(vec), type = "1")
  return(val)
}


rep_row <- function(x, n) {
       matrix(rep(x, each = n), nrow = n)
}


exp_fun <- function(x) {
  val <- exp(x)
  if (length(val) > 1) {
    for (i in 1:length(val)) {
      if (val[i] > .Machine$double.xmax) {
        val[i] <- .Machine$double.xmax
      }
    }
  } else {
   if (val > .Machine$double.xmax) {
     val <- .Machine$double.xmax
   }
 }
 return(val)
}


log_fun <- function(x) {
  val <- log(x)
  if (abs(val) > .Machine$double.xmax) {
    val <- sign(val) * .Machine$double.xmax
  }
  if (abs(val) < .Machine$double.xmin) {
    val <- sign(val) * .Machine$double.xmin
  }
  return(val)
}


outer_fun <- function(v) {
  outer(v, v)
}


# Sort and order the blocks by size
sort_blocks <- function(net_list, mod_names) {

  # Get the sizes of each block 
  n_ <- numeric(length(net_list))
  for (i in 1:length(net_list)) {
    n_[i] <- net_list[[i]]$gal$n
  }

  sort_order <- order(n_)
  net_list_order <- rep(list(NULL), length(net_list))
  mod_names_order <- rep(list(NULL), length(mod_names))
  for (i in 1:length(net_list)) {
    net_list_order[[i]] <- net_list[[sort_order[i]]]
    mod_names_order[[i]] <- mod_names[[sort_order[i]]]
  }
  return(list(net_list = net_list_order, mod_names = mod_names_order, sort_order = sort_order))
}

dim_fun <- function(n_obj, n_groups, eta_len) {
  sizes <- split_blocks(n_obj, n_groups)
  dims <- matrix(0, nrow = n_obj, ncol = eta_len)
  for (i in 1:n_groups) {
    num_ <- sizes[i]
    dim_1 <- 1 + sum(sizes[seq_spec(i, adjust = -1)])
    dim_2 <- sum(sizes[1:i])
    dims[dim_1:dim_2, ] <-
          rep_row(seq(1 + (i - 1) * eta_len, eta_len * i), sizes[i])
  }
  return(list(dims = dims, sizes = sizes))
}

split_blocks <- function(n_obj, n_groups) {
  base_ <- n_obj %/% n_groups
  left_ <- n_obj %% n_groups
  sizes <- rep(base_, n_groups) + c(rep(1, left_), rep(0, n_groups - left_))
  return(sizes)
}

seq_spec <- function(i, adjust = 0) {
  if (i == 1) {
    return(numeric(0))
  } else {
    return(1:(i + adjust))
  }
}

make_eta_fun <- function(num_group, parameterization) {
  if (parameterization == "multi_group") {
    eta_fun <- function(eta) {
      num_ <- 1
      len_ <- length(eta) / num_
      eta_base <- eta[1:len_]
      eta_val <- eta_base
      for (i in 2:num_) {
        dim_1 <- 1 + len_ * (i - 1)
        dim_2 <- len_ * i
        cur_val <- eta_base + eta[dim_1:dim_2]
        eta_val <- c(eta_val, cur_val)
      }
      return(eta_val)
    }
    body(eta_fun)[[2]] <- substitute(num_ <- num_group,
                                     list(num_group = num_group))

  } else if (parameterization == "size") {
    eta_fun <- function(eta) {
      return(eta)
    }
  }
  return(eta_fun)
}


make_eta_grad <- function(num_group, parameterization) {
  if (parameterization == "multi_group") {
    eta_grad <- function(eta) {
      num_ <- 1
      len_ <- length(eta) / num_
      eta_grad_val <- diag(len_)
      for (i in 2:num_) {
        eta_grad_val <- as.matrix(bdiag(eta_grad_val, diag(len_)))
      }
      eta_grad_val[ , 1:len_] <- rbind(t(matrix(rep(diag(len_), num_group),
                                         nrow = len_,
                                         ncol = num_group * len_)))
      return(eta_grad_val)
    }
    body(eta_grad)[[2]] <- substitute(num_ <- num_group,
                                      list(num_group = num_group))

  } else if (parameterization == "size") {
    eta_grad <- function(eta) {
      return(diag(length(eta)))
    }
  }
  return(eta_grad)
}


assign_labels <- function(K, sizes) {
  labels <- numeric(K)
  size_ <- c(0, sizes)
  for (i in 1:K) {
    labels[i] <- max(which(i > cumsum(size_)))
  }
  return(labels)
}


make_return_obj <- function(obj, labels, sort_order) {
  n_ <- length(unique(labels))
  return_list <- rep(list(NULL), n_)
  len_ <- length(obj$est$eta) / n_
  names(return_list) <- sprintf("group%i", 1:n_)
  grad <- obj$est$eta_grad(obj$est$eta)
  info_mat <- t(solve(grad)) %*% obj$est$info_mat %*% solve(grad)
  se_vec <- sqrt(diag(solve(info_mat)))
  for (i in 1:n_) {
    return_list[[i]] <- list(labels = NULL, estimates = NULL, se = NULL)
    return_list[[i]]$labels <- sort(sort_order[labels == i])

    dim_1 <- 1 + len_ * (i - 1)
    dim_2 <- len_ * i
    return_list[[i]]$estimates <- obj$est$eta_fun(obj$est$eta)[dim_1:dim_2]
    return_list[[i]]$se <- se_vec[dim_1:dim_2]
  }
  return(return_list)
}


check_extensions <- function(mod_names) {
  L <- length(mod_names)
  for (i in 1:L) {
    mod_names[[i]] <- strsplit(as.character(mod_names[[i]]), "_ijk")
    mod_names[[i]] <- strsplit(as.character(mod_names[[i]]), "_ij")
  }
  return(mod_names)
}


##################################################################
###
### tryCatch functions and others for error handling / checking
###
##################################################################

get_network_from_formula <- function(form) {
  result <- tryCatch(
              expr = {
                ergm.getnetwork(form);
              },
              error = function(err) {
                cat("\n")
                msg <- paste("The formula object provided to mlergm does not",
                             "contain a 'network' class object.\n",
                             "Formulas are specified:  net ~ term1 + term2 + ...")
                stop(msg)
              },
              warning = function(warn) {
                warning(warn);
              })
  return(result)
}

get_terms_from_formula <- function(form, net) {
  update.formula(form, net ~ .)
  result <- tryCatch(
              expr = {
                terms <- as.character(form)[3];
                sum_test <- summary(form)
                return(terms)
              },
              error = function(err) {
                bad_term <- str_match(as.character(err), "ERGM term (.*?) ")[2];
                if (is.na(bad_term)) {
                  bad_covariate <- str_match(as.character(err), "ergm(.*?): (.*?) is")[3];
                  err$message <- paste0("Covariate ", bad_covariate, " not a valid covariate.",
                                      " Please make sure that ", bad_covariate, " is a covariate of your network.")
                } else {
                  err$message <- paste0("Model term ", bad_term, " not a valid model term.",
                                        " Please reference 'help(ergm.terms)' for a list of",
                                        " valid model terms.");
                }
                cat("\n");
                stop(err);
              },
              warning = function(warn) {
                warning(warn);
              });
  return(terms);
}


check_and_convert_memb <- function(memb) {
  # Check if memb is a vector or can be converted to a vector
  if (!is.vector(memb)) {
    vec_memb <- tryCatch(
                  expr = {
                    as.vector(memb);
                  },
                  error = function(err) {
                    err$message <- paste0("Provided block memberships 'memb' not of class",
                                          " 'vector' and not convertable to class 'vector'.");
                    cat("\n")
                    stop(err);
                  },
                  warning = function(warn) {
                    warning(warn);
                  });
  } else {
    vec_memb <- memb;
  }

  # Now convert membership to numeric integer representation
  converted_memb <- vec_memb;
  unique_labels <- unique(vec_memb);
  iter <- 1;
  for (block_label in unique_labels) {
    which_match <- which(block_label == vec_memb);
    converted_memb[which_match] <- iter;
    iter <- iter + 1;
  }
  return_list <- list(memb_labels = unique_labels,
                      memb_internal = converted_memb)
  return(return_list);
}

check_net <- function(net) {
  if (!is.network(net)) {
    cat("\n")
    stop("Left-hand side of provided formula does not contain a valid object of class 'network'.");
  }
}

make_net_list <- function(net, memb_internal) {

  # Check that the dimensions of memb and net match
  if (network.size(net) != length(memb_internal)) {
    cat("\n")
    stop("Number of nodes in network and length of block membership vector are not equal.")
  }

  list_block_ind <- as.numeric(unique(memb_internal));
  net_list <- rep(list(NULL), length(list_block_ind));
  for (block_ind in list_block_ind) {
    nodes_in_cur_block <- which(block_ind == memb_internal);
    sub_net <- get.inducedSubgraph(net, v = nodes_in_cur_block);
    net_list[[block_ind]] <- sub_net;
  }
  return(net_list);
}


check_parameterization_type <- function(net_list, terms, parameterization, model) {

  # Check sufficient statistic sizes for each block
  block_statistic_dimensions <- numeric(length(net_list));
  for (i in 1:length(net_list)) {
    cur_net <- net_list[[i]];
    form_ <- as.formula(paste("cur_net ~", terms));
    block_statistic_dimensions[i] <- length(summary(form_));
  }
  which_largest <- which.max(block_statistic_dimensions)
  largest_block <- net_list[[which_largest]];
  form_ <- update(form_, largest_block ~ .);
  statistic_names <- names(summary(form_));
  model <- ergm_model(form_, largest_block);
  eta_map <- ergm.etamap(model);
  model_dimension <- max(block_statistic_dimensions);

  if (parameterization %in% c("standard", "offset")) {
    block_dims <- rep_row(rbind(seq(1, model_dimension)), length(net_list));
  }

  if (parameterization %in% c("offset")) {
    param_names <- get_coef_names(model, !is.curved(model))
    edge_ind <- which(param_names == "edges")
    mutual_ind <- which(param_names == "mutual")
    edge_loc <- ifelse(length(edge_ind) > 0, edge_ind, 0)
    mutual_loc <- ifelse(length(mutual_ind) > 0, mutual_ind, 0)
    if (edge_loc == 0) {
      edge_loc <- NULL
    }
    if (mutual_loc == 0) {
      mutual_loc <- NULL
    }

  } else {
    edge_loc <- NULL
    mutual_loc <- NULL
  }

  return_list <- list(model_dim = model_dimension,
                      model = model,
                      block_dims = block_dims,
                      eta_map = eta_map,
                      statistic_names = statistic_names,
                      edge_loc = edge_loc,
                      mutual_loc = mutual_loc,
                      which_largest = which_largest);

  return(return_list);
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


find_first_non_null <-  function(...) {
  for (x in list(...)) {
    if (!is.null(x)) {
      break
    }
  }
  x
}


check_integer <- function(val, name) {
  if (!is.numeric(val)) {
    cat("\n")
    stop(paste(name, "must be numeric."));
  }
  if (length(val) != 1) {
    cat("\n")
    stop(paste(name, "must be a single integer. Cannot supply multiple integers."));
  }
  if (!(val %% 1) == 0) {
    cat("\n")
    stop(paste(name, "must be an integer."));
  }
  if ((abs(val) > .Machine$integer.max)) {
    cat("\n")
    stop(paste(name, "provided is not a valid integer."));
  }
}



msplit <- function(x, y) {
  val <- suppressWarnings(split(x, y))
  return(val)
}



remove_between_block_edges <- function(net, memb) {
  index_mat <- matrix(TRUE, nrow = network.size(net), ncol = network.size(net))
  u_memb <- unique(memb)
  for (k in 1:length(u_memb)) {
    v_ind <- which(memb == u_memb[k])
    index_mat[v_ind, v_ind] <- FALSE
  }
  net[index_mat] <- 0
  return(net)
}


reorder_block_matrix <- function(net_list) {
  memb_vec <- numeric(0)
  attr_names <- list.vertex.attributes(net_list[[1]])
  v_attr <- rep(list(numeric(0)), length(attr_names))
  net_mat <- matrix(0, nrow = 0, ncol = 0)
  for (k in 1:length(net_list)) {
    sub_net <- net_list[[k]]
    for (i in 1:length(attr_names)) {
      v_attr[[i]] <- c(v_attr[[i]], get.vertex.attribute(sub_net, attr_names[i]))
    }
    memb_vec <- c(memb_vec, rep(k, network.size(sub_net)))
    net_mat <- bdiag(net_mat, sub_net[ , ])
  }
  net_mat <- as.matrix(net_mat)
  net <- network(net_mat, directed = is.directed(net_list[[1]]))
  for (i in 1:length(attr_names)) {
    set.vertex.attribute(net, attr_names[i], v_attr[[i]])
  }
  set.vertex.attribute(net, "node_memb_group", memb_vec)
  return(net)
}



adjust_formula <- function(form) {
  all_vars <- str_trim(str_split(as.character(form)[3], "\\+")[[1]])

  # Check if gw* terms are included without modifier
  if (any(all_vars == "gwesp")) {
    location <- which(all_vars == "gwesp")
    all_vars[location] <- "gwesp(fixed = FALSE)"
  }
  if (any(all_vars == "gwodegree")) { 
    location <- which(all_vars == "gwodegree")
    all_vars[location] <- "gwodegree(fixed = FALSE)"
  }
  if (any(all_vars == "gwidegree")) { 
    location <- which(all_vars == "gwidegree")
    all_vars[location] <- "gwidegree(fixed = FALSE)"
  }
  if (any(all_vars == "gwdegree")) { 
    location <- which(all_vars == "gwdegree")
    all_vars[location] <- "gwdegree(fixed = FALSE)"
  }

  # Put all the pieces back together
  right_side_change <- paste("~", paste0(all_vars, collapse = " + "))
  form <- update.formula(form, right_side_change)

  return(form)
}



compute_pvalue <- function(obj) {
  se <- sqrt(diag(solve(obj$est$info_mat)))
  obj$se <- se
  theta_est <- obj$est$theta

  z_val <- theta_est / se
  pvalue <- 2 * pnorm(-abs(z_val))
  pvalue <- as.numeric(pvalue)

  obj$pvalue <- pvalue
  return(obj)
}

format_form_for_cat <- function(form, len = 10) {
  all_vars <- str_trim(str_split(as.character(form)[3], "\\+")[[1]])
  char_lens <- nchar(all_vars)
  print_form <- paste0(as.character(form)[2] , " ~ ")
  base_len <- nchar(print_form)
  cur_len <- base_len
  for (i in 1:length(all_vars)) {
    print_form <- paste0(print_form, all_vars[i])
    cur_len <- cur_len + char_lens[i]
    if ((cur_len > 50) & (i < length(all_vars))) {
      print_form <- paste0(print_form, "\n")
      if (i < length(all_vars)) {
        print_form <- paste0(print_form, paste0(rep(" ", base_len + len), collapse = ""), "+ ")
        cur_len <- base_len
      } else {
        print_form <- paste0(print_form, paste0(rep(" ", base_len + len), collapse = ""))
      }
    } else {
      if (i < length(all_vars)) {
        print_form <- paste0(print_form, " + ")
        cur_len <- cur_len + 3
      }
    }
  }
  print_form <- paste0(print_form, "\n")
  return(print_form)
}


compute_bic <- function(obj) {
  total_edges <- sapply(obj$net$clust_sizes,
                        function(x, dir_flag ) {
                          if (dir_flag) {
                            2 * choose(x, 2)
                          } else {
                            choose(x, 2)
                          }
                        },
                        dir_flag = obj$net$directed_flag)
  total_edges <- sum(total_edges)
  bic_val <- log(total_edges) * length(obj$est$theta) - 2 * obj$likval
  return(bic_val)
}


compute_between_se <- function(eta1, eta2, num_dyads) {
  if (!is.null(eta2)) {
    covar_val <- matrix(0, nrow = 2, ncol = 2)
    covar_val[1, 1] <- (2 * exp(eta1) + 2 * exp(2 * eta1 + eta2) + exp(3 * eta1 + eta2)) /
                       (1 + 2 * exp(eta1) + exp(2 * eta1 + eta2))^2
    covar_val[2, 2] <- (exp(2 * eta1 + eta2) + 2 * exp(3 * eta1 + eta2)) /
                       (1 + 2 * exp(eta1) + exp(2 * eta1 + eta2))^2
    covar_val[1, 2] <- covar_val[2, 1] <- covar_val[2, 2]
  } else {
    covar_val <- matrix(0, nrow = 1, ncol = 1)
    covar_val[1, 1] <- exp(eta1) / (1 + exp(eta1))^2
  }

  covar_tot <- covar_val * num_dyads

  se_val <- as.numeric(sqrt(diag(solve(covar_tot))))
  return(se_val)
}



logit <- function(p) {
  val <- log_fun(p / (1 - p))
  return(val)
}




boxplot_fun <- function(dat_mat, line_dat = NULL, cutoff = NULL,
                        x_labels = NULL, x_angle = 0, 
                        x_axis_label = NULL, y_axis_label = "Count",
                        plot_title = "", title_size = 18, 
                        x_axis_size = NULL, y_axis_size = NULL, 
                        axis_size = 12, axis_label_size = 14,
                        x_axis_label_size = NULL, y_axis_label_size = NULL,
                        line_size = 1, stat_name = NULL, pretty_x = FALSE) {

  if (!is.null(line_dat)) { 
    if (length(line_dat) != ncol(dat_mat)) { 
      msg <- "Dimensions of 'line_dat' and 'dat_mat' must match;" 
      msg <- paste(msg, "'line_dat' must be a vector of length equal") 
      msg <- paste(msg, "to the number of columns of 'dat_mat'.\n") 
      stop(msg)
    }
  }

  if (!is.numeric(x_angle)) { 
    stop("Argument 'x_angle' must be numeric.\n")
  } else if (length(x_angle) != 1) { 
    stop("Argument 'x_angle' must be of length 1.\n")
  }

  if (!is.numeric(line_size)) { 
    stop("Argument 'line_size' must be numeric.\n")
  } else if (length(line_size) != 1) { 
    stop("Argument 'line_size' must be of length 1.\n")
  } else if (line_size < 0) { 
    stop("Argument 'line_size' must be non-negative.\n")
  }


  if (is.null(x_axis_label)) { 
    x_axis_label <- stat_name 
  } 
  if (!(length(x_axis_label) == 1)) { 
    stop("Argument 'x_axis_label' is not a single character string.\n") 
  } else if (!is.character(x_axis_label)) { 
    stop("Argument 'x_axis_label' is not a character string.\n")
  }
  if (!(length(y_axis_label) == 1)) {
    stop("Argument 'y_axis_label' is not a single character string.\n")
  } else if (!is.character(y_axis_label)) { 
    stop("Argument 'y_axis_label' is not a character string.\n")
  }
  if (!(length(plot_title) == 1)) { 
    stop("Argument 'plot_title' is not a single character string.\n")
  } else if (!is.character(plot_title)) { 
    stop("Argument 'plot_title' is not a character string.\n")
  }

  if (!is.numeric(title_size)) { 
    stop("Argument 'title_size' must be numeric.\n")
  } else if (length(title_size) != 1) { 
    stop("Argument 'title_size' must be of length 1.\n")
  } else if (title_size <= 0) { 
    stop("Argument 'title_size' must be a positive number.\n")
  }

  if (!is.numeric(axis_label_size)) {
    msg <- "Argument 'axis_label_size' must be a positive number." 
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_label_size' and 'y_axis_label_size.\n")
    stop(msg)
  }
  if (axis_label_size <= 0) { 
    stop("Argument 'axis_label_size' must be a positive number.\n")
  }
  
  if (!is.numeric(axis_size)) { 
    msg <- "Argument 'axis_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_size' and 'y_axis_size.\n")
    stop(msg)
  }
  if (axis_size <= 0) { 
    stop("Argument 'axis_size' must be a positive number.\n")
  }
  if (!is.numeric(axis_label_size)) {
    msg <- "Argument 'axis_label_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_label_size' and 'y_axis_label_size.\n")
    stop(msg)
  }
  if (axis_label_size <= 0) {
    stop("Argument 'axis_label_size' must be a positive number.\n")
  }

  if (!is.numeric(axis_size)) {
    msg <- "Argument 'axis_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_size' and 'y_axis_size.\n")
    stop(msg)
  }
  if (axis_size <= 0) {
    stop("Argument 'axis_size' must be a positive number.\n")
  }
  if (is.null(x_axis_label_size)) { 
    x_axis_label_size <- axis_label_size  
  } else { 
    if (!is.numeric(x_axis_label_size)) { 
      warning("Argument 'x_axis_label_size' not numeric. Using 'axis_label_size' instead.\n")
      x_axis_label_size <- axis_label_size
    } else if (!(length(x_axis_label_size) == 1)) { 
      warning("Argument 'x_axis_label_size' is not of length 1. Using 'axis_label_size instead.\n")
      x_axis_label_size <- axis_label_size
    } else if (x_axis_label_size <= 0) { 
      warning("Argument 'x_axis_label_size' not a positive number. Using 'axis_label_size' instead.\n")
      x_axis_label_size <- axis_label_size
    }
  }
  if (is.null(y_axis_label_size)) { 
    y_axis_label_size <- axis_label_size 
  } else { 
    if (!is.numeric(y_axis_label_size)) {
      warning("Argument 'y_axis_label_size' not numeric. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    } else if (!(length(y_axis_label_size) == 1)) {
      warning("Argument 'y_axis_label_size' is not of length 1. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    } else if (y_axis_label_size <= 0) {
      warning("Argument 'y_axis_label_size' not a positive number. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    }
  }

  if(is.null(x_axis_size)) { 
    x_axis_size <- axis_size 
  } else { 
    if (!is.numeric(x_axis_size)) { 
      warning("Argument 'x_axis_size' not numeric. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    } else if (!(length(x_axis_size) == 1)) { 
      warning("Argument 'x_axis_size' is not of length 1. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    } else if (x_axis_size <= 0) { 
      warning("Argument 'x_axis_size' not a positive number. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    }
  }
  if (is.null(y_axis_size)) { 
    y_axis_size <- axis_size 
  } else { 
    if (!is.numeric(y_axis_size)) { 
      warning("Argument 'y_axis_size' not numeric. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size 
    } else if (!(length(y_axis_size) == 1)) { 
      warning("Argument 'y_axis_size' is not of length 1. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size 
    } else if (y_axis_size <= 0) { 
      warning("Argument 'y_axis_size' is not a positive number. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size
    }
  } 

  first_colname <- colnames(dat_mat)[1]
  if (!is.null(x_labels) & !is.null(cutoff)) {
    if (cutoff != length(x_labels)) {
      stop("Value of argument 'cutoff' must be equal to length of 'x_labels'.\n")
    }
    if (grepl("0", first_colname)) { 
      dat_mat <- dat_mat[ , 1:(cutoff + 1)]
    } else { 
      dat_mat <- dat_mat[ , 1:cutoff] 
    }
  } else if (!is.null(x_labels)) {
    if (length(x_labels) != ncol(dat_mat)) {
      msg <- "Dimensions of 'x_labels' and 'dat_mat' must match;"
      msg <- paste(msg, "'x_labels' must be a vector of character labels equal")
      msg <- paste(msg, "to the number of columns of 'dat_mat'.\n")
      stop(msg) 
    }
    x_breaks <- 1:ncol(dat_mat)
  } else {
    if (!is.null(cutoff)) { 
      if (grepl("0", first_colname)) { 
        dat_mat <- dat_mat[ , 1:(cutoff + 1)]
      } else { 
        dat_mat <- dat_mat[ , 1:cutoff]
      }
    } 
    x_breaks <- 1:ncol(dat_mat)
    if (grepl("0", first_colname)) { 
      x_labels <- as.character(0:(length(x_breaks - 1)))
    } else { 
      x_labels <- as.character(x_breaks)
    }
    if (pretty_x) { 
      pretty_labels <- as.character(pretty(as.numeric(x_labels), n = 5))
      x_labels[!(x_labels %in% pretty_labels)] <- "" 
    }
  }
  dat_mat_colnames <- colnames(dat_mat) 

  if (is.null(cutoff)) { 
    cutoff <- ncol(dat_mat)
  }
  dat_mat <- melt(dat_mat)[ , 2:3]
  colnames(dat_mat) <- c("group", "values")
  dat_mat$group <- factor(dat_mat$group, levels = dat_mat_colnames)

  if (!is.null(line_dat)) {
    if (length(line_dat) > cutoff) {
      if (grepl("0", first_colname)) {
        line_dat <- line_dat[1:(cutoff + 1)]
       } else {  
        line_dat <- line_dat[1:cutoff]
      }
    }
  } else { 
    line_dat <- matrix(0, nrow = 0, ncol = ncol(dat_mat))
  }
  names(line_dat) <- dat_mat_colnames 
  line_dat <- melt(t(as.matrix(line_dat)))[ , 2:3]
  colnames(line_dat) <- c("group", "values")


  y_breaks <- pretty(dat_mat$values)
  y_labels <- as.character(y_breaks)
  geom_id <- c(rep("box", nrow(dat_mat)), rep("line", nrow(line_dat)))
  box_dat <- as.data.frame(cbind(rbind(dat_mat, line_dat), geom_id))

  # NULL out aes() inputs to appease CRAN check 
  group <- values <- NULL

  plot_ <- ggplot() +
              geom_boxplot(data = subset(box_dat, geom_id == "box"), 
                           aes(x = group, y = values), outlier.color = "NA") +
              geom_line(data = subset(box_dat, geom_id == "line"),
                        aes(x = 1:length(x_breaks), y = values),
                        color = "red", size = line_size) +  
              theme_classic() +
              labs(title = plot_title) +
              xlab(x_axis_label) +
              ylab(y_axis_label) +
              theme(axis.title.x = element_text(family = "Times",
                                                size = x_axis_label_size,
                                                colour = "Black",
                                                vjust = 0.5)) +
              theme(axis.title.y = element_text(family = "Times",
                                                size = y_axis_label_size,
                                                colour = "Black",
                                                margin = margin(r = 10))) +
              theme(plot.title = element_text(family = "Times",
                                              size = title_size,
                                              colour = "Black",
                                              vjust = 1)) +
              theme(axis.text.x = element_text(color = "black",
                                               family = "Times",
                                               size = x_axis_size,
                                               angle = x_angle,
                                               vjust = 0.2,
                                               hjust = 0.8)) +
              theme(axis.text.y = element_text(color = "black", 
                                               size = y_axis_size,
                                               family = "Times")) +
              theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
              theme(legend.position = "none") +
              scale_x_discrete(labels = x_labels) +  
              scale_y_continuous(expand = c(0, 1),
                                 breaks = y_breaks)

  return(plot_)
}



histplot_fun <- function(dat_mat, line_dat = NULL, 
                         x_axis_label = NULL, y_axis_label = "Count",  
                         plot_title = "", title_size = 18,
                         axis_label_size = 16, axis_size = 14, line_size = 1, 
                         x_axis_label_size = NULL, y_axis_label_size = NULL, 
                         x_axis_size = NULL, y_axis_size = NULL, 
                         stat_name = NULL) { 


  if (!is.numeric(dat_mat)) { 
    stop("Argument 'dat_mat' must be numeric.\n")
  } else if (!is.vector(dat_mat)) { 
    stop("Argument 'dat_mat' must be a vector.")
  }

  if (!is.null(line_dat)) { 
    if (!is.numeric(line_dat)) { 
      stop("Argument 'line_dat' must be numeric.\n")
    } else if (!is.vector(line_dat)) { 
      stop("Argument 'line_dat' must be a single number.\n")
    } else if (length(line_dat) != 1) { 
      stop("Argument 'line_dat' must be a single number.\n")
    }
  }

  if (!is.numeric(line_size)) { 
    stop("Argument 'line_size' must be numeric.\n")
  } else if (length(line_size) != 1) { 
    stop("Argument 'line_size' must be of length 1.\n")
  } else if (line_size < 0) { 
    stop("Argument 'line_size' must be non-negative.\n")
  }


  if (is.null(x_axis_label)) { 
    x_axis_label <- stat_name 
  } 
  if (!(length(x_axis_label) == 1)) { 
    stop("Argument 'x_axis_label' is not a single character string.\n") 
  } else if (!is.character(x_axis_label)) { 
    stop("Argument 'x_axis_label' is not a character string.\n")
  }
  if (!(length(y_axis_label) == 1)) {
    stop("Argument 'y_axis_label' is not a single character string.\n")
  } else if (!is.character(y_axis_label)) { 
    stop("Argument 'y_axis_label' is not a character string.\n")
  }
  if (!(length(plot_title) == 1)) { 
    stop("Argument 'plot_title' is not a single character string.\n")
  } else if (!is.character(plot_title)) { 
    stop("Argument 'plot_title' is not a character string.\n")
  }

  if (!is.numeric(title_size)) { 
    stop("Argument 'title_size' must be numeric.\n")
  } else if (length(title_size) != 1) { 
    stop("Argument 'title_size' must be of length 1.\n")
  } else if (title_size <= 0) { 
    stop("Argument 'title_size' must be a positive number.\n")
  }

  if (!is.numeric(axis_label_size)) {
    msg <- "Argument 'axis_label_size' must be a positive number." 
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_label_size' and 'y_axis_label_size.\n")
    stop(msg)
  }
  if (axis_label_size <= 0) { 
    stop("Argument 'axis_label_size' must be a positive number.\n")
  }

  if (!is.numeric(axis_size)) { 
    msg <- "Argument 'axis_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_size' and 'y_axis_size.\n")
    stop(msg)
  }
  if (axis_size <= 0) { 
    stop("Argument 'axis_size' must be a positive number.\n")
  }
  if (!is.numeric(axis_label_size)) {
    msg <- "Argument 'axis_label_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_label_size' and 'y_axis_label_size.\n")
    stop(msg)
  }
  if (axis_label_size <= 0) {
    stop("Argument 'axis_label_size' must be a positive number.\n")
  }

  if (!is.numeric(axis_size)) {
    msg <- "Argument 'axis_size' must be a positive number."
    msg <- paste(msg, "If you want to change the individual axis font sizes")
    msg <- paste(msg, "then you should use specify 'x_axis_size' and 'y_axis_size.\n")
    stop(msg)
  }
  if (axis_size <= 0) {
    stop("Argument 'axis_size' must be a positive number.\n")
  }
  if (is.null(x_axis_label_size)) { 
    x_axis_label_size <- axis_label_size  
  } else { 
    if (!is.numeric(x_axis_label_size)) { 
      warning("Argument 'x_axis_label_size' not numeric. Using 'axis_label_size' instead.\n")
      x_axis_label_size <- axis_label_size
    } else if (!(length(x_axis_label_size) == 1)) { 
      warning("Argument 'x_axis_label_size' is not of length 1. Using 'axis_label_size instead.\n")
      x_axis_label_size <- axis_label_size
    } else if (x_axis_label_size <= 0) { 
      warning("Argument 'x_axis_label_size' not a positive number. Using 'axis_label_size' instead.\n")
      x_axis_label_size <- axis_label_size
    }
  }
  if (is.null(y_axis_label_size)) { 
    y_axis_label_size <- axis_label_size 
  } else { 
    if (!is.numeric(y_axis_label_size)) {
      warning("Argument 'y_axis_label_size' not numeric. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    } else if (!(length(y_axis_label_size) == 1)) {
      warning("Argument 'y_axis_label_size' is not of length 1. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    } else if (y_axis_label_size <= 0) {
      warning("Argument 'y_axis_label_size' not a positive number. Using 'axis_label_size' instead.\n")
      y_axis_label_size <- axis_label_size
    }
  }

  if(is.null(x_axis_size)) { 
    x_axis_size <- axis_size 
  } else { 
    if (!is.numeric(x_axis_size)) { 
      warning("Argument 'x_axis_size' not numeric. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    } else if (!(length(x_axis_size) == 1)) { 
      warning("Argument 'x_axis_size' is not of length 1. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    } else if (x_axis_size <= 0) { 
      warning("Argument 'x_axis_size' not a positive number. Using 'axis_size' instead.\n")
      x_axis_size <- axis_size
    }
  }
  if (is.null(y_axis_size)) { 
    y_axis_size <- axis_size 
  } else { 
    if (!is.numeric(y_axis_size)) { 
      warning("Argument 'y_axis_size' not numeric. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size 
    } else if (!(length(y_axis_size) == 1)) { 
      warning("Argument 'y_axis_size' is not of length 1. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size 
    } else if (y_axis_size <= 0) { 
      warning("Argument 'y_axis_size' is not a positive number. Using 'axis_size' instead.\n")
      y_axis_size <- axis_size
    }
  } 

  # Obtain histogram breaks using David Scott's binwidth rule
  hist_values <- hist(dat_mat, plot = FALSE, breaks = "Scott")
  hist_breaks <- diff(hist_values$breaks)[1]
  if (is.null(line_dat)) { 
    line_dat <- matrix(0, nrow = 0, ncol = ncol(dat_mat))
  }

  y_breaks <- pretty(hist_values$counts)
  y_labels <- as.character(y_breaks)
  x_breaks <- pretty(dat_mat)
  x_labels <- as.character(x_breaks)
  geom_id <- c(rep("hist", length(dat_mat)), rep("line", 1))
  hist_values <- c(dat_mat, line_dat)
  hist_dat <- as.data.frame(cbind(hist_values, geom_id))
  rownames(hist_dat) <- NULL
  colnames(hist_dat) <- c("values", "geom_id")
  hist_dat$values <- as.numeric(levels(hist_dat$values))[hist_dat$values]
 
  
  # NULL out the aes() inputs to appease CRAN check 
  values <- NULL

  plot_ <- ggplot() +
              geom_histogram(data = subset(hist_dat, geom_id == "hist"), 
                             aes(values), binwidth = hist_breaks, 
                             fill = "grey75", color = "grey25") +
              geom_vline(data = subset(hist_dat, geom_id == "line"),
                        aes(xintercept = values),
                        color = "red", size = line_size) +
              theme_classic() +
              labs(title = plot_title) +
              xlab(x_axis_label) +
              ylab(y_axis_label) +
              theme(axis.title.x = element_text(family = "Times",
                                                size = x_axis_label_size,
                                                colour = "Black",
                                                vjust = 0.5)) +
              theme(axis.title.y = element_text(family = "Times",
                                                size = y_axis_label_size,
                                                colour = "Black",
                                                margin = margin(r = 10))) +
              theme(plot.title = element_text(family = "Times",
                                              size = title_size,
                                              colour = "Black",
                                              vjust = 1)) +
              theme(axis.text.x = element_text(color = "black",
                                               family = "Times",
                                               size = x_axis_size,
                                               vjust = 0.2,
                                               hjust = 0.8)) +
              theme(axis.text.y = element_text(color = "black", 
                                               size = y_axis_size,
                                               family = "Times")) +
              theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank()) +
              theme(legend.position = "none") + 
              scale_x_continuous(breaks = x_breaks, 
                                 labels = x_labels) +  
              scale_y_continuous(expand = c(0, 0),
                                 breaks = y_breaks)
  
  return(plot_)
}


check_terms <- function(form, K) { 

  check_formula(form) 
  all_vars <- all.vars(form, functions = TRUE)
  all_vars <- all_vars[!(all_vars %in% c("-", "+", "~"))]
  all_vars <- all_vars[-1]

  allowable_terms <- c("edges",
                       "mutual",
                       "gwesp",
                       "gwdegree",
                       "gwodegree",
                       "gwidegree", 
                       "triangle", 
                       "nodematch",
                       "transitiveties",
                       "cycle",
                       "ddsp",
                       "degree",
                       "desp",
                       "gwdsp",
                       "dsp",
                       "esp",
                       "isolates",
                       "kstar", 
                       "istar", 
                       "nodefactor",
                       "nodeifactor",
                       "nodeofactor",
                       "nodemix", 
                       "idegree",
                       "odegree",
                       "ostar", 
                       "twopath")
  if (K == 1) { 
    allowable_terms <- c(allowable_terms, "sender", "receiver", "sociality")
  }

  
  check_terms <- all_vars %in% allowable_terms 
  
  if (any(check_terms == FALSE)) { 
    location <- which(check_terms == FALSE) 
    msg <- "The following terms are not supported at this time: "
    for (i in 1:length(location)) { 
      cur_loc <- location[i]
      if (i < length(location)) { 
        msg <- paste0(msg, all_vars[cur_loc], ", ")
      } else { 
        msg <- paste0(msg, all_vars[cur_loc], ".\n")
      }
    }
    stop(msg)
  }
}




check_formula <- function(form) { 
  
  if (!is.formula(form)) { 
    stop("Argument 'form' must be a 'formula' class object.\n")
  }
  
  can_get_network <- tryCatch(ergm.getnetwork(form), 
                              error = function(err) { return(err) })

  if (!is.network(can_get_network)) { 
    stop("Cannot extract network from formula provided. Check that a valid formula was specified.")
  }

}





