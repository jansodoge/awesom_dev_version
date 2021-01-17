
#' Initialize SOM prototypes
#'
#' Prototypes are the artificial points in data space that are used to cluster
#' observations: each observation is assigned to the cluster of its closest
#' prototype. In self-organizing maps, each cell of the map has its own
#' prototype, and training is performed by iteratively adjusting the prototypes.
#' This function creates an initial guess for the prototypes of a SOM grid, to
#' be used as the `init` argument to the `kohonen::som` function (see example).
#'
#' @param traindat Matrix of training data, that will also be used to train the
#'   SOM.
#' @param nrows Number of rows on the map.
#' @param ncols Number of columns on the map.
#' @param method Method used, see Details. "pca" or "random"
#'
#' @return A matrix of prototype coordinates.
#' 
#' @details The default method "pca.sample" takes as prototypes the observations
#'   that are closest to the nodes of a 2d grid placed along the first two
#'   components of a PCA. The "pca" method uses the nodes instead of the
#'   observations. The "random" method samples random observations.
#'
#' @examples
#' dat <- iris[, c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")]
#' ### Scale training data
#' dat <- scale(dat)
#' ## Train SOM
#' ### Initialization (PCA grid)
#' init <- somInit(dat, 4, 4)
#' the.som <- kohonen::som(dat, grid = kohonen::somgrid(4, 4, 'hexagonal'), 
#'                         rlen = 100, alpha = c(0.05, 0.01), 
#'                         radius = c(6.08,-6.08), init = init, 
#'                         dist.fcts = 'sumofsquares')
somInit <- function(traindat, nrows, ncols, method= c("pca.sample", "pca", "random")) {
  method <- match.arg(method)
  if (method == "random") {
    init <- traindat[sample(nrow(traindat), nrows * ncols, replace= T), ]
  } else if (method %in% c("pca.sample", "pca")) {
    # the most detailed grid axis is assigned to the first component
    if (nrows >= ncols) {
      x.ev <- 1
      y.ev <- 2
    } else {
      x.ev <- 2
      y.ev <- 1
    }
    # perform PCA
    traindata.pca <- prcomp(traindat, center= F, scale.= F)
    init.x <- seq(from= quantile(traindata.pca$x[,x.ev], .025), 
                  to= quantile(traindata.pca$x[,x.ev], .975),
                  length.out= nrows)
    init.y <- seq(from= quantile(traindata.pca$x[,y.ev], .025), 
                  to= quantile(traindata.pca$x[,y.ev], .975),
                  length.out= ncols)
    init.base <- as.matrix(expand.grid(x= init.x, y= init.y)) # here a hex variant could be created instead if hex topology
    
    if (method == "pca.sample") {
      ## As in SOMbrero, init to observations closest to a 2D PCA grid
      closest.obs <- apply(init.base, 1, function(point) 
        which.min(colSums((t(traindata.pca$x[,c(x.ev,y.ev)])-point)^2)))
      init <- traindat[closest.obs,]
    } else if (method == "pca") {
      ## Pure PCA grid
      init <- tcrossprod(init.base, traindata.pca$rotation[, 1:2])
    }
  } 
  init
}


#' Distance measures of a SOM
#'
#' Several distance measures between cells or prototypes of a trained SOM (in
#' grid space, in data space).
#'
#' @param som ```kohonen``` object, a SOM created by the ```som``` function.
#'
#' @return A list with distance measures
#'
somDist <- function(som){
  if (is.null(som)) return(NULL)
  proto.gridspace.dist <- kohonen::unit.distances(som$grid, F)
  proto.dataspace.dist <- as.matrix(dist(som$codes[[1]]))
  neigh <- round(proto.gridspace.dist, 3) == 1
  proto.dataspace.dist.neigh <- proto.dataspace.dist
  proto.dataspace.dist.neigh[!neigh] <- NA
  list(proto.grid.dist= proto.gridspace.dist, 
       neigh.matrix= neigh, 
       proto.data.dist= proto.dataspace.dist, 
       proto.data.dist.neigh= proto.dataspace.dist.neigh)
}


#' Calculate SOM quality measures
#'
#' @param som `kohonen` object, a SOM created by the `kohonen::som` function.
#' @param traindat matrix containing the training data
#'
#' @references Kohonen T. (2001) \emph{Self-Organizing Maps}, 3rd edition,
#'   Springer Press, Berlin.
#'
#' Kaski, S. and Lagus, K. (1996) Comparing Self-Organizing Maps. In C. von der
#' Malsburg, W. von Seelen, J. C. Vorbruggen, and B. Sendho (Eds.)
#' \emph{Proceedings of ICANN96, International Conference on Articial Neural
#' Networks , Lecture Notes in Computer Science} vol. 1112, pp. 809-814.
#' Springer, Berlin.
#'
#' @return A list with quality measures : quantization error, share of explained
#'   variance, topographic error and Kaski-Lagus error.
#'   
somQuality <- function(som, traindat){
  if(is.null(som)) return(NULL)
  ok.dist <- somDist(som)
  
  ## BMU, Squared distance from obs to BMU
  bmu <- som$unit.classif
  sqdist <- rowSums((traindat - som$codes[[1]][bmu, ])^2)
  
  ## Quantization error
  err.quant <- mean(sqdist)
  
  ## Interclass variance ratio
  totalvar <- sum(apply(traindat, 2, var)) * 
    (nrow(traindat) - 1) / nrow(traindat)
  err.varratio <- 100 - round(100 * err.quant / totalvar, 2)
  
  ## Topographic error
  bmu2 <- apply(traindat, 1, function(row) {
    dist <- colSums((t(som$codes[[1]]) - row)^2)
    order(dist)[2]
  })
  err.topo <- mean(!ok.dist$neigh.matrix[cbind(bmu, bmu2)])
  
  ## Kaski-Lagus error
  err.kaski <- e1071::allShortestPaths(ok.dist$proto.data.dist.neigh)$length[cbind(bmu, bmu2)]
  err.kaski <- mean(err.kaski + sqrt(sqdist))
  
  ## Distribution of individuals in cells
  cellpop <- table(factor(som$unit.classif, levels= 1:nrow(som$grid$pts)))
  
  res <- list(err.quant= err.quant, err.varratio= err.varratio, 
              err.topo= err.topo, err.kaski= err.kaski, cellpop= cellpop)
  class(res) <- "somQual"
  res
}

