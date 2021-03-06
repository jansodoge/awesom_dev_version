% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{aweSOMsilhouette}
\alias{aweSOMsilhouette}
\title{Silhouette plot}
\usage{
aweSOMsilhouette(som, clust)
}
\arguments{
\item{som}{`kohonen` object, a SOM created by the `kohonen::som` function.}

\item{clust}{object containing the result of the super-clustering of the SOM
prototypes (either a `hclust` or a `pam` object).}
}
\description{
Plots a silhouette plot, used to assess the quality of the super-clustering
of SOM prototypes into superclasses. Available for both PAM and
hierarchical clustering.
}
\examples{
## Build training data
dat <- iris[, c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")]
### Scale training data
dat <- scale(dat)
## Train SOM
### RNG Seed (for reproducibility)
### Initialization (PCA grid)
init <- somInit(dat, 4, 4)
ok.som <- kohonen::som(dat, grid = kohonen::somgrid(4, 4, 'hexagonal'),
                       rlen = 100, alpha = c(0.05, 0.01),
                       radius = c(6.08,-6.08), init = init,
                       dist.fcts = 'sumofsquares')
## Group cells into superclasses (PAM clustering)
superclust <- cluster::pam(ok.som$codes[[1]], 2)
superclasses <- superclust$clustering
aweSOMsilhouette(ok.som, superclasses)
}
