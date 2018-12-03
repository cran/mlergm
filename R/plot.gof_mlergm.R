#' Plot goodness-of-fit results
#' 
#' Produces goodness-of-fit plots for a \code{gof_mlergm} object in order to visualize and assess the fit of an estimated model produced by \code{\link{mlergm}}. 
#'
#' @param x An object of class \code{gof_mlergm}, produced by \code{\link{gof.mlergm}}. 
#' @param \dots Additional argument to be passed if necessary. 
#' @param individual_plots (Logical \code{TRUE} or \code{FALSE}) If \code{TRUE}, individual gof plots are produced. Defaults to 
#' \code{FALSE}.  
#' @param save_plots (Logical \code{TRUE} or \code{FALSE}) If \code{TRUE}, the individual GOF plots are saved.
#' @param show_plots (Logical \code{TRUE} or \code{FALSE}) If \code{TRUE}, the plots are printed to the screen, and if \code{FALSE} no plots are displayed. This may be helpful when the only desire is to save the individual GOF plots.
#' @param width If \code{save_plots == TRUE}, controls the plot width dimension saved. 
#' @param height If \code{save_plots == TRUE}, controls the plot height dimension saved. 
#' @param cutoff For statistics that are distributions (e.g., degree distributions), specifies a cutoff point. Dimensions past the cutoff are ignored and not plotted. 
#' @param x_labels Character vector specifying the statistic names or labels. 
#' @param x_angle Adjusts the angle of the x axis tick labels (typically the statistic names). 
#' @param x_axis_label Label for the x axis. 
#' @param y_axis_label Label for the y aixs. 
#' @param plot_title Title for the plot. 
#' @param title_size Font size for the plot title. 
#' @param axis_label_size Font size for the axis labels. Individual axes label sizes can be changed using \code{x_axis_label_size} and \code{y_axis_label_size} which are detailed below. 
#' @param axis_size Font size for the axis tick labels. Individual axes tick label sizes can be changed using \code{x_axis_size} and \code{y_axis_size} which are detailed below. 
#' @param line_size (Numeric, non-negative) If \code{line_size} is positive, then a red line will be plotted to indicate the observed network value of the statistic. If \code{line_size} is equal to zero, then the observed data line will not be plotted. 
#' @param x_axis_label_size The font size of the x axis label. When \code{NULL}, \code{axis_label_size} is used. Defaults to \code{NULL}.
#' @param y_axis_label_size The font size of the y axis label. When \code{NULL}, \code{axis_label_size} is used. Defaults to \code{NULL}. 
#' @param x_axis_size The font size of the x axis tick labels. When \code{NULL}, \code{axis_size} is used. Defaults to \code{NULL}. 
#' @param y_axis_size The font size of the y acis tick labels. When \code{NULL}, \code{axis_size} is used. Defaults to \code{NULL}.
#' @param pretty_x (Logical \code{TRUE} or \code{FALSE}) If set to \code{TRUE}, the \code{link{pretty}} function will be called to format the x-axis breaks. This can be useful for when the x-axis range is large.  
#'
#' @export
#' @import ggplot2
#' @importFrom ggplot2 ggsave ggplot geom_boxplot geom_line theme_classic labs xlab ylab theme element_text margin element_line element_blank scale_x_discrete scale_y_continuous geom_histogram geom_vline scale_x_continuous aes
#' @importFrom graphics hist 
plot.gof_mlergm <- function(x, ..., 
                            individual_plots = FALSE, save_plots = FALSE, show_plots = TRUE,
                            width = 8, height = 4.5, 
                            cutoff = NULL, x_labels = NULL, x_angle = 0, 
                            x_axis_label = NULL, y_axis_label = "Count",
                            plot_title = "", title_size = 18, 
                            axis_label_size = 14, axis_size = 10, line_size = 1, 
                            x_axis_label_size = NULL, y_axis_label_size = NULL,
                            x_axis_size = NULL, y_axis_size = NULL, pretty_x = TRUE) {  

  if (!is.gof_mlergm(x)) { 
    stop("Argument must be an 'gof_mlergm' object. See 'help(mlergm)' and 'help(gof.mlergm)' for details.\n")
  }

  # Check main arguments (plot arguments checked in boxplot_fun) 
  if (!is.logical(individual_plots)) { 
    stop("Argument 'individual_plots' must be TRUE or FALSE.\n")
  }
  if (!is.logical(save_plots)) { 
    stop("Argument 'save_plots' must be TRUE or FALSE.\n")
  }
  if (!is.numeric(width)) { 
    stop("Argument 'width' must be numeric.\n")
  } else if (width <= 0) { 
    stop("Argument 'width' must be positive.\n")
  }
  if (!is.numeric(height)) { 
    stop("Argument 'height' must be numeric.\n")
  } else if (height <= 0) { 
    stop("Argument 'height' must be positive.\n")
  }

  # Parse the different model terms
  gof_stats <- x$gof_stats 
  different_stats <- x$stat_names 
  obs_stats <- x$obs_stats
  rm(x); clean_mem() 


  hist_flag <- box_flag <- FALSE
  for (i in 1:length(different_stats)) {
    cur_stat <- different_stats[i]
    which_cur <- grepl(cur_stat, colnames(gof_stats))
    if (sum(which_cur) == 1) { 
      hist_flag <- TRUE
    } else { 
      box_flag <- TRUE
    }
    cur_gof_data <- gof_stats[ , which_cur]
    cur_obs_data <- obs_stats[which_cur]
    
    if (hist_flag) { 
      cur_plot <- histplot_fun(dat_mat = cur_gof_data, 
                               line_dat = cur_obs_data,
                               x_axis_label = x_axis_label, 
                               y_axis_label = y_axis_label,
                               plot_title = plot_title, 
                               title_size = title_size, 
                               axis_label_size = axis_label_size,
                               axis_size = axis_size, 
                               line_size = line_size, 
                               x_axis_label_size = x_axis_label_size,
                               y_axis_label_size = y_axis_label_size,
                               x_axis_size = x_axis_size, 
                               y_axis_size = y_axis_size,
                               stat_name = cur_stat)


    } else if (box_flag) { 
      cur_plot <- boxplot_fun(dat_mat = cur_gof_data, 
                              line_dat = cur_obs_data,
                              cutoff = cutoff, 
                              x_labels = x_labels,
                              x_angle = x_angle, 
                              x_axis_label = x_axis_label,
                              y_axis_label = y_axis_label, 
                              plot_title = plot_title, 
                              title_size = title_size, 
                              axis_label_size = axis_label_size, 
                              axis_size = axis_size, 
                              line_size = line_size, 
                              x_axis_label_size = x_axis_label_size, 
                              y_axis_label_size = y_axis_label_size, 
                              x_axis_size = x_axis_size, 
                              y_axis_size = y_axis_size,
                              stat_name = cur_stat,
                              pretty_x = pretty_x) 
    } 
    assign(paste0("plot", i), cur_plot)

    # Rest plot flags 
    hist_flag <- box_flag <- FALSE
  }

  if (individual_plots) { 
    for (i in 1:length(different_stats)) {
      plot_name <- paste0("plot", i)
      cur_plot <- get(plot_name)
      if (show_plots) { 
        print(cur_plot) 
        readline(prompt="Press [enter] to continue")
      }
      if (save_plots) {
        cur_filename <- paste0("gof_", different_stats[i], ".pdf") 
        ggsave(filename = cur_filename, 
               plot = cur_plot, 
               device = "pdf", 
               width = width, 
               height = height)
      }
    }
  } else { 
    if (length(different_stats) == 1) { 
      cur_plot <- get("plot1") 
      if (show_plots) { 
        print(cur_plot)
      }
      if (save_plots) { 
        cur_filename <- paste0("gof_", different_stats[i], ".pdf")
        ggsave(filename = cur_filename,
               plot = cur_plot, 
               device = "pdf", 
               width = width, 
               height = height)
      }
    } else { 
      plot_command <- "plot_grid("
      for (i in 1:length(different_stats)) { 
        plot_name <- paste0("plot", i)
        cur_plot <- get(plot_name)
        if (i < length(different_stats)) { 
          plot_command <- paste0(plot_command, plot_name, ", ")
        } else { 
          plot_command <- paste0(plot_command, plot_name, ")")
        }
      }
      plot_ <- eval(parse(text = plot_command))
      if (show_plots) { 
        print(plot_)
      }
      if (save_plots) { 
        cur_filename <- paste0("gof.pdf")
        ggsave(filename = cur_filename, 
               plot = plot_, 
               device = "pdf", 
               width = width, 
               height = height)
      }
    }
  }
}






