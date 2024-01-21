#' improvedktaucenters
#'
#' Robust Clustering algorithm for non-spherical data. This function estimate
#' clusters taking into account that clusters may have
#' different size, volume or orientation.
#' @param X numeric matrix of size n x p.
#' @param K number of clusters.
#' @param cutoff argument for outliers detection - quantiles of chi-square
#' to be used as a threshold for outliers detection, defaults to 0.999.
#' @param nstart number of trials that the base ktaucenters is run at the first stage.
#' If it is greater than 1 and center is not set as NULL, a random set of (distinct) 
#' rows in x is chosen as the initial centres for each trial.
#' @param INITcenters numeric matrix of size K x p indicating the initial centers for
#' that clusters and robust covariance matrices will be computed, if it is set as NULL the
#' algorithm will compute from ktaucenters routine. Set to NULL by default.
#' @return A list with the following components:
#'  \item{\code{centers}}{: Matrix of size K x p, with the estimated K centers.}
#'  \item{\code{cluster}}{: A vector of integer (from 1:k) indicating the cluster to
#' which each point is allocated.}
#'  \item{\code{sigmas}}{: A list containing the k covariance matrices found by the 
#' procedure at its second step.}
#'  \item{\code{outliers}}{: indices observation that can be considered as outliers.}
#'
#' @importFrom GSE GSE getLocation getScatter
#' @importFrom MASS ginv
#' @importFrom stats mahalanobis qchisq
#'
#' @examples
#'
#' # Generate synthetic data (three normal cluster in two dimensions)
#' # Clusters have different shapes and orientation.
#' # The data is contaminated uniformly (level 20%).
#'
#' # Generates base clusters
#' set.seed(1)
#' Z1 <- c(rnorm(100, 0), rnorm(100, 0), rnorm(100, 0))
#' Z2 <- rnorm(300)
#' X <- matrix(0, ncol = 2, nrow = 300)
#' X[, 1] <- Z1
#' X[, 2] <- Z2
#' true.cluster <- c(rep(1, 100), rep(2, 100), rep(3, 100))
#'
#' # Rotate, expand and translate base clusters
#' theta <- pi/3
#' aux1 <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
#' aux2 <- sqrt(4) * diag(c(1, 1/4))
#' B <- aux1 %*% aux2 %*% t(aux1)
#' X[true.cluster == 3, ] <-
#'   X[true.cluster == 3, ] %*% aux2 %*% aux1 + matrix(c(5, 2),
#'                                                   byrow = TRUE,
#'                                                   nrow = 100,
#'                                                   ncol = 2)
#' X[true.cluster == 2, 2] <- X[true.cluster == 2, 2] * 5
#' X[true.cluster == 1, 2] <- X[true.cluster == 1, 2] * 0.1
#' X[true.cluster == 1, ] <- X[true.cluster == 1, ] + matrix(c(-5, -1),
#'                                                           byrow = TRUE,
#'                                                           nrow = 100,
#'                                                           ncol = 2)
#'
#' # Generate 60 synthetic outliers (contamination level 20%)
#'
#' outliers <- sample(1:300, 60)
#' X[outliers, ] <- matrix(runif( 40, 2 * min(X), 2 * max(X) ),
#'                                 ncol = 2, nrow = 60)
#'
#' # Applying the algorithm
#' robust <- improvedktaucenters(X, K = 3, cutoff = 0.999)
#'
#' # Plotting results
#' oldpar <- par(mfrow = c(2, 1))
#' plot(X, main = "Actual clusters")
#' for (j in 1:3){
#'  points(X[true.cluster == j, ], pch = 19, col = j + 1)
#' }
#' points(X[outliers, ], pch = 19, col = 1)
#' plot(X, main = "Clusters estimation")
#' for (j in 1:3){
#'  points(X[robust$cluster == j,], pch = 19, col = j + 1)
#' }
#' points(X[robust$outliers, ], pch = 19)
#' 
#' par(oldpar)
#' @references Gonzalez, J. D., Yohai, V. J., & Zamar, R. H. (2019).
#' Robust Clustering Using Tau-Scales. arXiv preprint arXiv:1906.08198.
#' @export
improvedktaucenters <- function(X,
                                K,
                                cutoff = 0.999,
                                nstart = 5,
                                INITcenters = NULL) {

  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # Determine the best centers
  if (is.null(INITcenters)){
    ret_ktau <- ktaucentersfast(X, K, nstart = nstart)
    centers <- ret_ktau$centers
    newClusters <-ret_ktau$clusters
  } else {
    sphericalDistanceMatrix <- matrix(0, ncol = K, nrow = n)
    centers <- INITcenters
    for (j in 1:K){
    sphericalDistanceMatrix[, j] <- mahalanobis(x = X,
                                                center = centers[j, ],
                                                cov = diag(p))
      
    }  
    newClusters <- apply(sphericalDistanceMatrix, 1,
                        function(x) which(x == min(x))[1])
  }
  
  sigmas <- replicate(K, diag(p), simplify = FALSE)
  newcentersaux <- 0 * centers

  mahalanobisMatrix <- matrix(0, ncol = K, nrow = n)
  for (j in 1:K){
    Xcluster <- X[newClusters == j, ]
    
    # dim works properly when there are two or more observations
    if (sum(newClusters == j) < 2){
      nk <- sum(newClusters == j)
    } else {
      nk <- dim(Xcluster)[1]
    }
    
    # if there is enough observations, we compute
    # robust covariance matrix 
    if(nk > (3*p)){
      sal1 <- GSE(Xcluster, tol = 1e-4, maxiter = 150)
      newcentersaux[j, ] <- getLocation(sal1)
      sigmas[[j]] <- getScatter(sal1)
    } else {
      newcentersaux[j, ] <- centers[j, ]
    }
    mahalanobisMatrix[, j] <- mahalanobis(x = X,
                                          center = newcentersaux[j, ],
                                          cov = ginv(sigmas[[j]]),
                                          inverted = TRUE)
    
  }

  newClusters_Mah <- apply(mahalanobisMatrix, 1,
                          function(x) which(x == min(x))[1])
  value <- qchisq(cutoff, df = p)

  #Improved outliers determination
  outliers <- c()
  for (j in 1:K){
    indices <- which(newClusters_Mah == j)
    outliersk <- indices[mahalanobisMatrix[indices, j] > value]
    outliers <- c(outliersk, outliers)
  }

  list(centers = newcentersaux,
       cluster = newClusters_Mah,
       outliers = outliers,
       sigmas = sigmas)
}
