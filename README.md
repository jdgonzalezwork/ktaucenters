# ktaucenters package: Robust and efficient Clustering

**Juan D. Gonzalez, Victor J. Yohai, Ruben H. Zamar and Douglas Carmona**

## Introduction

This package implements a clustering algorithm similar to kmeans, it has two main advantages:

- The estimator is resistant to outliers, that means that results of estimator are still correct when there are atypical values in the sample.

- The estimator is efficient, roughly speaking, if there are not outliers in the sample (all data is good), results will be similar than those obtained by a classic algorithm (kmeans)

Clustering procedure is carried out by minimizing the overall robust scale so-called tau scale (see Yohai Gonzalez, Yohai and Zamar 2019 [arxiv:1906.08198](https://arxiv.org/abs/1906.08198)).


# How to use the package ktaucenters

### Example 1: behavior when data is clean

Load the package ktaucenters
``` {.r}
rm(list=ls())
library(ktaucenters)
```

Generate synthetic data (three cluster well separated), and apply a classic algorithm (*kmeans*) and the robust *ktaucenters*.

``` {.r}
# Generate synthetic data (three cluster well separated)
set.seed(1)
Z <- rnorm(600);
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol=2)

# Applying the ROBUST algortihm
ktau_output <- ktaucenters(X, K = 3, nstart = 10)

# Applying the classic algortihm
kmeans_output <- kmeans(X, centers = 3, nstart = 10)

# Plotting the center results 
plot(X, main = "Efficiency")
points(ktau_output$centers, pch = 19, col = 2, cex = 2)
points(kmeans_output$centers, pch = 17, col = 3, cex = 2)
legend(-6, 6, pch = c(19,17), col = c(2,3), cex = 1, legend = c("ktau centers", "kmeans centers"))
```

### Example 2: Behavior in presence of outliers

Contaminate the previous data by replacing 60 observations to outliers located in a bounding box that contains the clean data. Then apply *kmeans* and *ktaucenters* algorithms.

``` {.r}
# Generate 60 synthetic outliers (contamination level 20%)
X[sample(1:300,60), ] <- matrix(runif( 40, 2* min(X), 2 * max(X)),
                                ncol = 2, nrow = 60)

# Applying the ROBUST algortihm
ktau_output <- ktaucenters(X, K = 3, nstart = 10)

# Applying the classic algortihm

kmeans_output <- kmeans(X, centers = 3, nstart = 10)

# Plotting the center results 
plot(X, main = "Robustness")
points(ktau_output$centers, pch = 19, col = 2, cex = 2)
points(kmeans_output$centers, pch = 17, col = 3, cex = 2)
legend(-10, 10, pch = c(19, 17), col = c(2, 3), cex = 1, legend = c("ktau centers", "kmeans centers"))
```

After running the code, it can be observed that kmeans centers were very influenced by outliers, while ktaucenters results are still razonable.


### Example 3: Showing clusters and outliers detection procedure

Continuation from Example 2, for outliers recognition purposes we can see the `ktau_output$outliers` that indicates the indices that may be considered as outliers, on the other hand, the labels of each cluster found by the algorithm are coded with integers between 1 and K (in this case K=3), the variable `ktau_output$clusters` contains that information.

``` {.r}
plot(X, main = "Estimated clusters and outliers detection")

# Plotting clusters 
for (j in 1:3){
  points(X[ktau_output$cluster==j, ], col=j+1)
}

# Plotting outliers 
points(X[ktau_output$outliers, ], pch = 19, col = 1, cex = 1)
legend(7, 15, pch = c(1,1,1,19), col = c(2,3,4,1), cex = 1,
       legend = c("cluster 1", "cluster 2", "cluster 3", "detected \n outliers"), bg = "gray")
```
The final figure contains clusters and outliers detected. 

# Improved-ktaucenters

The algorithm *ktaucenters* works well under noisy data, but fails when clusters have different size, shape and orientation, an algorithm suitable for this sort of data is `improvektaucenters`. To show how this algorithm works we use the data set so-called M5data from package `tclust: tclust: Robust Trimmed Clustering`, [tclust](https://CRAN.R-project.org/package=tclust). M5 data were generated by three normal bivariate distributions with different scales, one of the components is very overlapped with another one. A 10% background noise is added uniformly.

First, load the data and run the `improvedktaucenters` function.

``` {.r}
# Load non spherical datadata 

library("tclust")
data("M5data")
X <- M5data[,1:2]
true.clusters <- M5data[,3]

# Estimate clusters
improved_output <- improvedktaucenters(X, K = 3, cutoff = 0.95)
```

Keep the results in the variable `improved_output`, that is a list that contains the fields `outliers` and `cluster`, among others. For example, we can have access to the cluster labeled as 2 by typing

`X[improved_output$cluster==2, ]`.

If we want to know the values of outliers, type

`X[improved_output$outliers, ]`.

By using these commands, it is easy to estimate the original clusters by means of `improvedktaucenters` routine. 

The preprint [arxiv:1906.08198](https://arxiv.org/abs/1906.08198) contains comparison with other robust clustering procedures as well as technical details and applications.   
