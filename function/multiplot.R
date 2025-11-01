#' Arrange multiple plots on a single page (with an optional title row)
#'
#' Draws several plots on the current graphics device in a grid layout. The
#' layout can be supplied explicitly via a matrix, or it will be constructed
#' from the requested number of columns. If a non-empty \code{title} is given,
#' an extra top row is reserved for the title text.
#'
#' @param ... One or more plot objects (typically \pkg{ggplot2} plots or other
#'   grid-compatible grobs) to draw.
#' @param plotlist list or \code{NULL}. An optional list of plots to include in
#'   addition to those passed via \code{...}. This is useful when the plots are
#'   stored in a list.
#' @param file character. \strong{Deprecated/ignored.} Present for historical
#'   compatibility; the function does not write to a file.
#' @param cols integer(1). Number of columns in the layout when \code{layout}
#'   is not supplied. The number of rows is computed from the number of plots.
#' @param layout integer matrix or \code{NULL}. A matrix specifying the
#'   arrangement of plots. Each positive integer in the matrix identifies the
#'   position(s) for the corresponding plot (e.g., all cells with value \code{2}
#'   are occupied by the second plot). When a non-empty \code{title} is used,
#'   a top row of zeros is prepended internally to reserve space for the title.
#' @param title character(1). An optional page title placed in a dedicated top
#'   row spanning all columns.
#' @param fontsize numeric(1). Title font size (points) passed to \code{grid::gpar()}.
#' @param fontfamily character(1). Title font family passed to \code{grid::gpar()}.
#'
#' @details
#' The function relies on the \pkg{grid} system: it opens a new page
#' (\code{grid::grid.newpage()}), sets up a \code{grid::grid.layout()} and then
#' prints each plot into its assigned viewport. If only a single plot is
#' provided, it is printed directly without creating a new page layout.
#'
#' The \code{layout} matrix governs plot placement and spanning. For example, a
#' matrix such as \code{rbind(c(1, 1), c(2, 3))} places plot 1 across the top
#' row (spanning two columns), and plots 2 and 3 in the bottom-left and
#' bottom-right cells, respectively.
#'
#' @return Invisibly returns \code{NULL}. The function is called for its side
#'   effect of drawing plots on the active graphics device.
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL, title = "",
                      fontsize = 12, fontfamily = "Helvetica") {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
      ncol = cols, nrow = ceiling(numPlots / cols)
    )
  }

  if (nchar(title) > 0) {
    layout <- rbind(rep(0, ncol(layout)), layout)
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout),
      ncol(layout),
      heights = if (nchar(title) > 0) {
        unit(c(0.5, rep(5, nrow(layout) - 1)), "null")
      } else {
        unit(c(rep(5, nrow(layout))), "null")
      }
    )))

    # Make each plot, in the correct location
    if (nchar(title) > 0) {
      grid.text(title,
        vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)),
        gp = gpar(fontsize = fontsize, fontfamily = fontfamily)
      )
    }

    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}
