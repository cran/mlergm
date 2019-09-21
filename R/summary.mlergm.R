#' @describeIn mlergm
#'
#' Prints a summary of the estimated \code{mlergm} model.
#'
#' @param object An object of class \code{mlergm}, probably produced by \code{\link{mlergm}}.
#' @param \dots Additional arguments to be passed if necessary.  
#'
#' @export 
summary.mlergm <- function(object, ...) { 

  if (!is.mlergm(object)) { 
    stop("Argument must be an 'mlergm' object. See 'help(mlergm)' for details.\n", call. = FALSE)
  }

  if (object$estimation_status == "success") { 
    cat("\n")
    cat("============================")
    cat(" Summary of model fit ")
    cat("============================")
    cat("\n\n")
    cat("Formula:  ")
    form_cat <- format_form_for_cat(object$form)
    cat(form_cat)
    if (object$parameterization == "offset") { 
      cat("Parameterization set to 'offset'\n")
    }
    if (object$parameterization == "size") { 
      cat("Parameterization set to 'size'\n")
    }
    cat("\n")
    cat("Number of blocks:  ") 
    cat(paste(length(unique(object$node_memb))))
    cat("\n\n")
    cat("Quantiles of block sizes:")
    cat("\n")
    print(object$size_quantiles) 
    cat("\n\n")
    cat("Monte Carlo MLE Results:\n")
    if ((object$parameterization == "offset") & ("edges" %in% names(object$theta))) {
      cat("    Within-block edge parameter    =  edges  - log(Block size)\n")
    } 
    if ((object$parameterization == "offset") & ("mutual" %in% names(object$theta))) { 
      cat("    Within-block mutual parameter  =  mutual + log(Block size)\n")
    }
    if (object$parameterization == "size") { 
      cat("    Within-lock parameter = parameter * log n(k),   n(k) is the size of block k\n")
    }
    cat("\n")
    theta_names <- names(object$theta)
    max_char <- max(nchar(c(theta_names, "between edges", "between mutual")))
    name_space <- paste(rep(" ", max_char), collapse = "")
    cat(name_space)
    cat("    Estimate   Std. Error    p-value    Sig.\n")
    for (i in 1:length(theta_names)) { 
      cur_name <- theta_names[i] 
      white_space <- paste(rep(" ", max_char - nchar(cur_name)), collapse = "")
      cat(cur_name)
      cat(white_space)
      cat("     ")
      if (object$theta[i] >= 0) { 
        cat(" ")
      }
      theta_val <- sprintf("%.4f", signif(object$theta[i], 4))
      cat(theta_val)
      cat("    ") 
      if (object$se[i] >= 0) { 
        cat(" ")
      }
      se_val <- sprintf("%.4f", signif(object$se[i], 4))
      cat(se_val)
      cat("      ")
      p_value <- object$pvalue[i]
      if (p_value < 0.00001) { 
        p_value <- "<0.00001" 
      } else { 
        p_value <- sprintf("%.5f", round(p_value, digits = 5))
        cat(" ")
      }
      cat(p_value)
      cat("    ")
      if (p_value >= 0.1) { 
        cat("   ")
      } else if (p_value >= 0.05) { 
        cat(".  ")
      } else if (p_value >= 0.01) { 
        cat("*  ")
      } else if (p_value >= 0.001) { 
        cat("** ")
      } else { 
        cat("***")
      }
      cat("\n") 
    }
    cat("-----------------------------------------------------------------") 
    cat("\n")
    if (is.character(object$between_theta)) { 
      cat(object$between_theta)
    } else { 
      for (i in 1:length(object$between_theta)) { 
        cur_name <- names(object$between_theta)[i]
        cat(cur_name) 
        white_space <- paste(rep(" ", max_char - nchar(cur_name)), collapse = "")
        cat("      ")
        if (object$between_theta[i] >= 0) {
          cat("")
        }
        theta_val <- sprintf("%.4f", signif(object$between_theta[i], 4))
        cat(theta_val)
        cat("    ")
        if (object$between_se[i] >= 0) {
          cat(" ")
        }
        se_val <- sprintf("%.4f", signif(object$between_se[i], 4))
        cat(se_val)
        
        cat("      ")
        p_value <- object$between_pvalue[i]
        if (p_value < 0.00001) {
          p_value <- "<0.00001"
        } else {
          p_value <- sprintf("%.5f", round(p_value, digits = 5))
          cat(" ")
        }
        cat(p_value)
        cat("    ")
        if (p_value >= 0.1) {
          cat("   ")
        } else if (p_value >= 0.05) {
          cat(".  ")
        } else if (p_value >= 0.01) {
          cat("*  ")
        } else if (p_value >= 0.001) {
          cat("** ")
        } else {
          cat("***")
        }
        cat("\n")
      }
    }
    cat("\n\n")
    cat("Sig. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
    cat("\n\n")
    cat("BIC:  ")
    if (!is.null(object$bic)) { 
      cat(paste(round(object$bic, digits = 3)))
      cat("\n* Note: BIC is based on the within-block model, and ignores the between-block model.")
    } else { 
      cat("Not estimated.")
    }
  } else { 
    cat("\n")
    cat("============================")
    cat(" Summary of model fit ")
    cat("============================")
    cat("\n\n")
    cat("Formula:  ")
    form_cat <- format_form_for_cat(object$form)
    cat(form_cat)
    cat("\n")
    cat("Model estimation unsuccesful.")
  }
  cat("\n\n")
}


