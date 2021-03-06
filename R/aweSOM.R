################################################################################
################################################################################
## aweSOM funtion: launch shiny interface

#' aweSOM interface
#'
#' Launches the (offline) web-based interface for training and visualizing
#' self-organizing maps.
#' 
#' @references Kohonen T. (2001) \emph{Self-Organizing Maps}, 3rd edition,
#'   Springer Press, Berlin.
#'   
aweSOM <- function() {
  shiny::runApp(system.file('shiny', package='aweSOM'))
}

################################################################################
################################################################################
## S3 methods

print.somQual <- function(x, ...) {
  cat("\n## Quality measures:\n", 
      "* Quantization error     : ", x$err.quant, "\n",
      "* (% explained variance) : ", x$err.varratio, "\n",
      "* Topographic error      : ", x$err.topo, "\n",
      "* Kaski-Lagus error      : ", x$err.kaski, "\n", 
      "\n## Number of obs. per map cell:")
  print(x$cellpop)
  
}

