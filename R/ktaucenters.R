#' ktaucenters
#'
#' Robust Clustering algorithm.
#' @param x numeric matrix of size n x p, or an object that can be coerced to a matrix (such as a numeric vector or a data frame with all numeric columns).
#' @param centers either the number of clusters, say *k*, or a matrix of initial (distinct) cluster centers. If a number, a random set of (distinct) rows in x is chosen as the initial centers.
#' @param nstart if centers is a number, how many random sets should be chosen?
#' @param use_kmeans use kmeans centers as starting point?
#' @param use_robin use robin algorithm centers as starting point?
#' @param max_iter the maximum number of iterations allowed.
#' @param max_tol maximum tolerance parameter used for the algorithm as stopping rule.
#' @param cutoff quantile of chi-square distribution to be used as a threshold for outliers detection, defaults to 0.999

#' @return A list including the estimated k centers and labels for the observations
##' \itemize{
##'  \item{\code{centers}}{:   matrix of size K x p, with the estimated K centers.}
##'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
##'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
##'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights
##' associated to each observation.}
##'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some
##' iteration there were clusters totally empty}
##'  \item{\code{niter}}{: number of iterations until convergence is achieved
##' or maximun number of iteration is reached}
##'  \item{\code{di}}{: distance of each observation to its assigned cluster-center}
##'  \item{\code{outliers}}{: indices observation that can be considered as outliers}
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
#'
#' @importFrom stats kmeans
#' @export
ktaucenters <- function(x,
                        centers,
                        nstart = 1L,
                        use_kmeans = TRUE,
                        use_robin = TRUE,
                        max_iter = 100L,
                        max_tol = 1e-6,
                        cutoff = 0.999) {
  # Parameters check
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  if (missing(centers)) {
    stop("'centers' must be a number or a matrix")
  }
  
  if (length(centers) == 1L) {
    n_clusters <- centers
    ## To avoid duplicates
    if (nstart >= 2L) {
      unique_x <- unique(x)
      n_unique <- nrow(unique_x)
      if (n_unique < n_clusters)
        stop("more cluster centers than distinct data points.")
    }
  } else {
    centers <- as.matrix(centers)
    
    if (any(duplicated(centers)))
      stop("initial centers are not distinct")
    n_clusters <- nrow(centers)
    if (n < n_clusters)
      stop("more cluster centers than data points")
    if (ncol(x) != ncol(centers)) {
      stop("'x' and 'centers' must have same number of columns")
    }
  }
  
  # Set up center initialization
  if (is.matrix(centers)) {
    init_centers <- list(centers)
  } else {
    init_centers <-
      replicate(nstart, unique_x[sample.int(n_unique, n_clusters),], simplify = FALSE)
  }
  
  if (use_kmeans) {
    kmeans_centers <-
      kmeans(x, n_clusters, nstart = 20)$centers
    init_centers <- append(init_centers, list(kmeans_centers))
  }
  
  if (use_robin) {
    robin_centers_idx <- robinden(distance(x), n_clusters, 10)$centers + 1
    robin_centers <- x[robin_centers_idx,]
    init_centers <- append(init_centers, list(robin_centers))
  }
  
  # Runs
  best_tau <- Inf
  for (iter in seq_along(init_centers)) {
    current_run <-
      ktaucenters_run(x, init_centers[[iter]], max_tol, max_iter)
    if (current_run$tau < best_tau) {
      best_run <- current_run
      best_tau <- best_run$tau
    }
  }
  
  # Outlier detection
  best_run <- flag_outliers(cutoff, 0.5, best_run)
  
  return(best_run)
}