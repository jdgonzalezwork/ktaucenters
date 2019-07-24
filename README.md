---
title: 'ktaucenters package: Robust and efficient Clustering'
author: "Juan D. Gonzalez, Victor J. Yohai, Ruben H. Zamar"
date: '2019-07-23'
generator: pandoc
viewport: width=device-width, initial-scale=1
---


# ktaucenters package: Robust and efficient Clustering 
====================================================

#### Juan D. Gonzalez, Victor J. Yohai, Ruben H. Zamar 


## Introduction
------------

This package implements a clustering algorithm similar to kmeans, it has
two main advantages:

-   The estimator is resistant to outliers, that means that results of
    estimator are still correct when there are atipycal values in the
    sample.

-   The estimator is efficient, roughly speaking, if there are not
    outliers in the sample (all data is good), results will be similar
    than those obtained by a classic algorithm (kmeans)

Clustering procedure is carried out by minimizing the overall robust
scale so-called tau scale (see Yohai and Zamar, 1988,
[doi:10.1080/01621459.1988.10478611](https://www.tandfonline.com/doi/abs/10.1080/01621459.1988.10478611)
and Gonzalez, Yohai and Zamar 2019
[arxiv:1906.08198](https://arxiv.org/abs/1906.08198)).



# How to use the package ktaucenters
----------------------------------

### Example 1: behavior when data are clean

First we load the package ktaucenters
``` {.r}
rm(list=ls())
library(ktaucenters)
```

We generate synthetic data (three cluster well separated), and apply a
classic algorithm (kmeans) and the robust ktaucenters. Results and code
are shown below.

``` {.r}
# Generate synthetic data (three cluster well separated)
set.seed(1)
Z <- rnorm(600);
mues <- rep(c(-4, 0, 4), 200)
X <-  matrix(Z + mues, ncol=2)

### Applying the ROBUST algortihm  ####
ktau_output <- ktaucenters(X, K=3,nstart=10)
### Applying the classic algortihm  ####
kmeans_output <- kmeans(X,centers=3,nstart=10)

### plotting the center results 
plot(X,main=" Efficiency")
points(ktau_output$centers,pch=19,col=2,cex=2)
points(kmeans_output$centers,pch=17,col=3,cex=2)
legend(-6,6,pch=c(19,17),col=c(2,3),cex=1,legend=c("ktau centers" ,"kmeans centers"))
```


![alt text](https://github.com/jdgonzalezwork/ktaucenters/imagesPNG/figure1.png)
{.figure}
Clean data. Estimated centers by K-means and KTAU-centers algorithms
{width="480"}

This figure shows that there are no differeces between kmeans and ktaucenters
in clean data.


### Example 2: behavior in the presence of outliers

We contaminate the previous data by replacing 60 observations to
outliers located in a bounding box that contains the clean data. Then we
apply kmeans and ktaucenters algorithms.

``` {.r}
# Generate 60 sintetic outliers (contamination level 20%)
X[sample(1:300,60), ] <- matrix(runif( 40, 2* min(X), 2 * max(X) ),
                                ncol = 2, nrow = 60)

### Applying the ROBUST algortihm  ####
ktau_output <- ktaucenters(X, K=3,nstart=10)
### Applying the classic algortihm  ####
kmeans_output <- kmeans(X,centers=3,nstart=10)
```

``` {.r}
### plotting the estimated centers 
plot(X,main=" Robustness ")
points(ktau_output$centers,pch=19,col=2,cex=2)
points(kmeans_output$centers,pch=17,col=3,cex=2)
legend(-10,10,pch=c(19,17),col=c(2,3),cex=1,legend=c("ktau centers" ,"kmeans centers"))
```

::: {.figure style="text-align: center"} IT REMAINS TO add figure 2
![\\label{fig:fig2} Contaminated data. Estimated centers by K-means and
KTAU-centers
algorithms.](){width="480"}

Contaminated data. Estimated centers by K-means and KTAU-centers
algorithms.
:::

As it can be observed in Figure kmeans center were very influenced by
outliers, while ktaucenters results are still razonable.
:::

::: {#example-3-showing-clusters-and-outliers-detection-procedure .section .level3}
### Example 3: Showing clusters and outliers detection procedure

Continuation from Example 2, for outliers recognition purposes we can
see the `ktau_output$outliers` that indicates the indices that may be
considered as outliers, on the other hand, the labels of each cluster
found by the algorithm are coded with integers between 1 and K (in this
case K=3), the variable `ktau_output$clusters` contains that
information.

``` {.r}
plot(X,main=" Estimated clusters and outliers detection ")
## plottig clusters 
for (j in 1:3){
  points(X[ktau_output$cluster==j, ], col=j+1)
}

## plottig outliers 
points(X[ktau_output$outliers, ], pch=19, col=1, cex=1)
legend(7,15,pch=c(1,1,1,19),col=c(2,3,4,1),cex=1,
       legend=c("cluster 1" ,"cluster 2","cluster 3","detected \n outliers"),bg = "gray")
```
IT REMAINS TO add figure 3
![](){width="480"}


# Improved-ktaucenters 
--------------------

The algorithm ktaucenter works well under noisy data, but fails when
clusters have different size, shape and orientation, an algorithm
suitable for this sort of data is `improvektaucenters`. To show how this
algorithm works we use the data set so-called M5data from package
`tclust: tclust: Robust Trimmed Clustering`,
[tclust](https://cran.r-project.org/web/packages/tclust/index.html). M5
data were generated by three normal bivariate distributions with
different scales, one of the components is very overlapped with another
one. A 10% background noise is added uniformly

::: {#usage .section .level3}
### usage

First we load the data, then, run the `improvedktaucenters` function.

``` {.r}
## load non spherical datadata 
library("tclust")
data("M5data")
X=M5data[,1:2]
true.clusters=M5data[,3]
### done ###### 

#run the function to estimate clusters
improved_output=improvedktaucenters(X,K=3,cutoff=0.95)
```

We keep the results in the variable `improved_output`, that is a list
that contains the fields `outliers` and `cluster`, among others. For
example, we can have access to the cluster labeled as 2 by typing

`X[improved_output$cluster==2, ]`.

If we want to know the values of outliers, type

`X[improved_output$outliers, ]`.

By using these commands, it is easy to estimate the original clusters by
means of `improvedktaucenters` routine.
IT REMAINS TO add figure 4a  and 4b 
![]

# Real data application: finding a screw in Mars
----------------------------------------------

Package has a dataset called `mars_screw`, containing the Intensity and
Saturation pixels values of a picture from Mars taken from Rover
Curiosity.

Data set can be loaded by typing `mars_screw$SI_matrix`, or
`mars_screw$geographic_matrix`, thats variables means.

-   SI\_matrix: A matrix with 5063 rows and 128 columns. Elements 1 to
    64 of each row indicate the Saturation values of pixels in a square
    cell 8 x 8 whereas elements 65 to 128 of each row indicate the
    cell's Intensity values.

-   geographic\_matrix: An integer matrix of dimension 5063 x 2, each
    row indicates each square cell's locations (x-axis y-axis) at the
    picture.

-   screw\_index: the index corresponding to the screw observation
    (screw\_index=4180)

In order to find the screw, in a first stage, we run the clustering
algorithm over the matrix A, with the aim to separate the differents
materials presented in the picture

``` {.r}
A <- mars_screw$SI_matrix;
B <- mars_screw$geographic_matrix
screw_index <-mars_screw$screw_index
## 
ret1=ktaucenters(X=A,K=3);
## 
```

To find the screw candidate, we run a clustering procedure on
geographic\_matrix, for each cluster determined at the previous step.
Outliers will be those that are the furthest observations regarding to
the gruops they were assigned. The procedure to find them is shown
below.

``` {.r}
screw_candidate_index=c()
for (j in 1:3){
  dj=ktaucenters(B[ret1$cluster==j, ],K=5,nstart=1,startWithROBINPD = FALSE)$di; 
  jcandidate=dj==max(dj);
  jcandidate=which(ret1$cluster==j)[dj==max(dj)]
  screw_candidate_index=c(screw_candidate_index,jcandidate)
}
```

Figure , shows locations of three candidates to be the screw.

``` {.r}
plot(B, type="n" )
for (j in 1:3){
#  col1=which(orderInd==j)
  col1=j
  points(B[ret1$cluster==j,],col=col1+1,pch=19)
}
points(B[screw_candidate_index,1],B[screw_candidate_index,2],col=6,pch=1,cex=3,lwd=3)
```
IT REMAINS TO add figure 5
::: {.figure style="text-align: center"}
![\\label{fig:fig3} Plot of Geographic sub matrices, pink circles are
the outliers candidates in the geographic matix
space.](){width="576"}

The images corresponding to the three chanels are quite similar, this is
because the RGB pixels values are usually highly correlated. On the
other hand, if we want to rebuild the data-set from the original
picture, the following code can be used.

``` {.r}
# To reconstruct the data from source image. 
XXX=cbind(myjpg[,,1],myjpg[,,2],myjpg[,,3])
# define functions Intensity and saturation 
Icolor=function(COLOR){((COLOR[1]+COLOR[2]+COLOR[3]))/3}
Scolor=function(COLOR){
  Iaux=Icolor(COLOR)
  ret=0
  if(!(Iaux==0)){ret=1 - min(COLOR[1],COLOR[2],COLOR[3])/Iaux}
  if(Iaux==0){ret=0}
  ret
}

I=apply(XXX,1,Icolor)
S=apply(XXX,1,Scolor)

myjpgHSI=myjpg
myjpgHSI[,,2]=S
myjpgHSI[,,3]=I
scomp=myjpgHSI[,,2];
icomp=myjpgHSI[,,3];
d=8
p=d^2
nrowIm=dim(myjpg)[1] 
ncolIm=dim(myjpg)[2] 
NporM=ncolIm*nrowIm

# transforming each d-square cell into an array of size 2*dxd
A=matrix(0,ncol=2*p,nrow=NporM/p)
B=matrix(0,ncol=2,nrow=NporM/p)
l=1;
for (i in 1:(nrowIm/d)){
  for (j in 1:(ncolIm/d)){
    posi=(i-1)*d;
    posj=(j-1)*d;
    auxs=scomp[posi+(1:d), posj +(1:d)];
    auxi=icomp[posi+(1:d), posj +(1:d)];
    A[l,]=c(auxs[,],auxi[,])
    B[l,]=c(i,j)
    l=l+1;
  }
}

A=A[1:(l-1),]
```

This is a check whether A and mars\_screw\$SI\_matrix are equal

``` {.r}
mean(abs(A-mars_screw$SI_matrix))==0
#> [1] FALSE
```

