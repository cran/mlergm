#' @describeIn mlnet
#'
#' Plots network objects of type \code{mlnet}. 
#'
#' @param x An object of class \code{mlnet}, possibly produced by \code{\link{mlnet}} or \code{\link{simulate_mlnet}}. 
#' @param node_size Controls the size of nodes. 
#' @param memb_colors Specifies the named colors to be used for the membership colors.  
#' @param palette If package \code{RColorBrewer} is installed, then the name of an R color brewer pallete can be specified and used for the block colors. See \code{brewer.pal} for details on RColorBrewer palletes. 
#' @param arrow.gap (Directed graphs only) Controls the amount of space between arrowheads and the nodes. 
#' @param arrow.size (Directed graphs only) Controls the size of the arrowhead. 
#' @param legend (\code{TRUE} or \code{FALSE}) Controls whether the block membership legend is printed. 
#' @param legend.position The position of the legend in the plot. Defaults to the "right" position. 
#' @param color_legend_title Name for the node color legend title. 
#' @param layout_type Viable layout options. See \code{gplot.layout} for options. 
#' @param \dots Additional arguments to be passed to \code{\link{ggnet2}}.
#'
#' @export 
#' @importFrom GGally ggnet2
#' @importFrom grDevices rainbow 
#' @importFrom sna gplot.layout.kamadakawai
plot.mlnet <- function(x, node_size = 2.5, palette = NULL, memb_colors = NULL, 
                       arrow.gap = 0.015, arrow.size = 4, 
                       color_legend_title = "", legend = TRUE, legend.position = "right", 
                       layout_type = "kamadakawai", ...) {

  if (!is.mlnet(x)) { 
    stop("Argument must be an 'mlnet' object. See 'help(mlnet)' for details.\n", call. = FALSE)
  }

  # Rename 'x' to a more informative name 
  mlnet_object <- x 
  rm(x); clean_mem() 

  # Find out how many blocks there are in the network 
  node_memb <- get.vertex.attribute(mlnet_object, "node_memb") 
  K <- length(unique(node_memb))

  # Set the color palette 
  if (is.null(palette) & is.null(memb_colors)) { 
    if (K <= 9) { 
      memb_palette <- "Set1"
    } else { 
      memb_palette <- rainbow(K)
      # getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      # memb_palette <- getPalette(K)
      names(memb_palette) <- as.character(unique(node_memb))
    }
  } else if (!is.null(palette) & is.null(memb_colors)) { 
    memb_palette <- palette
  }

  if (!is.null(memb_colors) & (is.null(names(memb_colors))))  { 
    memb_palette <- memb_colors 
    names(memb_palette) <- as.character(unique(node_memb))
  }

  # Check directedness 
  if (!is.directed(mlnet_object)) { 
    arrow.gap <- 0
    arrow.size <- 0
  }

  # Check legend 
  if (legend == FALSE) { 
    legend.position = "n"
  } 


  ggnet2(mlnet_object, 
         size = node_size, 
         color = "node_memb",
         palette = memb_palette,
         arrow.gap = arrow.gap, 
         arrow.size = arrow.size,
         legend.position = legend.position,
         color.legend = color_legend_title, 
         mode = layout_type, 
         ...)
}




