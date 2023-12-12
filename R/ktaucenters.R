#' ktaucenters
#'
#' Robust Clustering algorithm based on centers, a robust and efficient version of KMeans.
#' @param X numeric matrix  of size n x p.
#' @param K the number of cluster.
#' @param centers a matrix of size K x p containing the K initial centers,
#'  one at each matrix-row. If centers is NULL a random set of (distinct) rows in \code{X}
#'  are chosen as the initial centres.
#' @param tolmin a tolerance parameter used for the algorithm stopping rule
#' @param NiterMax a maximum number of iterations used for the algorithm stopping rule
#' @param nstart the number of trials that the base algorithm is run.
#' If it is greater than 1 and center is not set as NULL, a random set of (distinct) rows
#' in \code{X} will be chosen as the initial centres.
#' @param startWithKmeans  TRUE if kmean centers values is included as starting point.
#' @param startWithROBINPD TRUE if ROBINDEN estimator is included as starting point.
#' @param cutoff optional argument for outliers detection - quantiles of chi-square to be used as a threshold for outliers detection, defaults to 0.999
#' @return A list including the estimated K centers and labels for the observations
##' \itemize{
##'  \item{\code{centers}}{: matrix of size K x p, with the estimated K centers.}
##'  \item{\code{cluster}}{: array of size n x 1  integers labels between 1 and K.}
##'  \item{\code{tauPath}}{: sequence of tau scale values at each iterations.}
##'  \item{\code{Wni}}{: numeric array of size n x 1 indicating the weights associated to each observation.}
##'  \item{\code{emptyClusterFlag}}{: a boolean value. True means that in some iteration there were clusters totally empty}
##'  \item{\code{niter}}{: number of iterations until convergence is achieved or maximum number of iteration is reached}
##'  \item{\code{di}}{: distance of each observation to its assigned cluster-center}
##'  \item{\code{outliers}}{: indices observation that can be considered as outliers}
##' }
#'
#' @examples
#' # Generate synthetic data (three clusters well separated)
#' Z <- rnorm(600)
#' mues <- rep(c(-3, 0, 3), 200)
#' X <- matrix(Z + mues, ncol = 2)
#'
#' # Generate 60 synthetic outliers (contamination level 20%)
#' X[sample(1:300,60), ] <- matrix(runif( 40, 3 * min(X), 3 * max(X) ),
#'                                 ncol = 2, nrow = 60)
#'
#'sal <- ktaucenters(
#'      X, K = 3, centers = X[sample(1:300, 3), ],
#'      tolmin = 1e-3, NiterMax = 100)
#'
#' oldpar <- par(mfrow = c(1,2))
#' 
#' plot(X,type = "n", main = "ktaucenters (Robust) \n outliers: solid black dots")
#' points(X[sal$cluster==1,], col = 2)
#' points(X[sal$cluster==2,], col = 3)
#' points(X[sal$cluster==3,], col = 4)
#' points(X[sal$outliers, 1], X[sal$outliers, 2], pch = 19)
#'
#' # Classical (non Robust) algortihm
#' sal <- kmeans(X, centers = 3, nstart = 100)
#'
#' plot(X, type = "n", main = "kmeans (Classical)")
#' points(X[sal$cluster==1,], col = 2)
#' points(X[sal$cluster==2,], col = 3)
#' points(X[sal$cluster==3,], col = 4)
#'
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019). 
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198. 
#'
#' @importFrom stats kmeans dist qchisq
#' @export
ktaucenters <- function(X,
                       K,
                       centers = NULL,
                       tolmin = 1e-06,
                       NiterMax = 100,
                       nstart = 1,
                       startWithKmeans = TRUE,
                       startWithROBINPD = TRUE,
                       cutoff = 0.999) {

  warning("In a next version (second semester 2024), this function will be changed.")
  if (!is.matrix(X)) {
    X <- as.matrix(X)
    
  }
  init_centers <- centers
  taumin <- 1e+20
  n <- nrow(X)
  p <- ncol(X)
  
  centers0 <- matrix(0, nrow = K, ncol = p)
  start = 1 * (!startWithKmeans)
  # the start value is one or zero.
  nstartEnd <- nstart + 1 * (startWithROBINPD)
  
  for (trial in start:nstartEnd) {
    # if startWithKmeans its true, start=0, then trial take the zero value.
    if (trial == 0) {
      sal0 <- kmeans(X, centers = K, nstart = 20)
      sal0$labels <- sal0$cluster
      for (jota in 1:K) {
        # when there is a single observation, it is not possible to use apply function
        if (sum(sal0$labels == jota) == 1) {
          centers0[jota,] <- X[sal0$labels == jota,]
        }
        
        if (sum(sal0$labels == jota) > 1) {
          # as.matrix below is necessary because if p=1 function "apply" will not work.
          centers0[jota,] <-
            apply(as.matrix(X[sal0$labels == jota,]), 2, mean)
        }
      }
    }
    
    if (trial >= 1) {
      centers0 <- X[sample(1:dim(X)[1], K),]
    }
    
    if ((trial == 1) & (!is.null(init_centers))) {
      centers0 <- init_centers
      
    }
    if (trial == nstart + 1) {
      retROB <- robinden(as.matrix(dist(X)),
                         n_clusters = K,
                         10)
      
      centers0 <- X[retROB$centers + 1,]
    }
    
    centers <- centers0
    
    ret_ktau <- ktaucenters_run(X, centers, tolmin, NiterMax)

    tau <- ret_ktau$tau;
    niter = ret_ktau$iter
    
    if (tau < taumin) {
      taumin <- tau
      
      best_tau <- tau
      best_ret_ktau <- ret_ktau
    }
  }
  
  newClusters <- best_ret_ktau$cluster
  squaredi <- (best_ret_ktau$di) ^ 2
  robustScale <- Mscale(u = sqrt(squaredi),
                       c = normal_consistency_constants(p),
                       b = 0.5)
  
  outliers <- c()
  value <- qchisq(cutoff, df = p)
  for (j in 1:K) {
    indices <- which(newClusters == j)
    booleansubindices <- (squaredi[indices] / (robustScale ^ 2)) > value
    outliersk <- indices[booleansubindices]
    outliers <- c(outliersk, outliers)
  }
  best_ret_ktau$outliers <- outliers
  best_ret_ktau
}

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
##' }
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
#'
#' @importFrom stats kmeans
#' @export
ktaucentersfast <- function(x,
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
    unique_x <- unique(x)
    n_unique <- nrow(unique_x)
    if (n_unique < n_clusters)
      stop("more cluster centers than distinct data points.")
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
  
  # Setup center initialization
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