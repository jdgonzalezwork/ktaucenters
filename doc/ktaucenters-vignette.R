## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=5, 
  fig.height=5
)

## ---- results = 'hide'---------------------------------------------------
rm(list=ls())
library(ktaucenters)

## ----fig1,fig.cap = "\\label{fig:fig1} Clean data. Estimated centers by K-means and KTAU-centers algorithms."----
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

## ------------------------------------------------------------------------
# Generate 60 sintetic outliers (contamination level 20%)
X[sample(1:300,60), ] <- matrix(runif( 40, 2* min(X), 2 * max(X) ),
                                ncol = 2, nrow = 60)

### Applying the ROBUST algortihm  ####
ktau_output <- ktaucenters(X, K=3,nstart=10)
### Applying the classic algortihm  ####
kmeans_output <- kmeans(X,centers=3,nstart=10)

## ----fig2,fig.height = 5, fig.width = 5, fig.align = "center", fig.cap = "\\label{fig:fig2} Contaminated data. Estimated centers by K-means and KTAU-centers algorithms."----
### plotting the estimated centers 
plot(X,main=" Robustness ")
points(ktau_output$centers,pch=19,col=2,cex=2)
points(kmeans_output$centers,pch=17,col=3,cex=2)
legend(-10,10,pch=c(19,17),col=c(2,3),cex=1,legend=c("ktau centers" ,"kmeans centers"))

## ------------------------------------------------------------------------
plot(X,main=" Estimated clusters and outliers detection ")
## plottig clusters 
for (j in 1:3){
  points(X[ktau_output$cluster==j, ], col=j+1)
}

## plottig outliers 
points(X[ktau_output$outliers, ], pch=19, col=1, cex=1)
legend(7,15,pch=c(1,1,1,19),col=c(2,3,4,1),cex=1,
       legend=c("cluster 1" ,"cluster 2","cluster 3","detected \n outliers"),bg = "gray")

## ----fig.show='hold',fig.height = 4.5, fig.width = 4.5-------------------

## load non spherical datadata 
library("tclust")
data("M5data")
X=M5data[,1:2]
true.clusters=M5data[,3]
### done ###### 

#run the function to estimate clusters
improved_output=improvedktaucenters(X,K=3,cutoff=0.95)

## ----echo=FALSE,fig.show='hold',fig.height = 4.5, fig.width = 9----------
par(mfrow=c(1,2))
plot(X,type="n",pch=19,main="(a) M5Data - original clusters \n (outliers in black dots)" )
colors1=c(3,4,2)
for (j in 1:3){
  points(X[true.clusters==j, ], col=colors1[j],pch=j-1)
}
points(X[true.clusters==0, ], col=1,pch=19)

plot(X,type="n",pch=19,main="(b) improved ktaucenter Output \n (estimated outliers in black dots)" )

orderInd=order(improved_output$centers[,1]) 
for (j in 1:3){
  col1=which(orderInd==j)
  points(X[improved_output$cluster==j, ], col=col1+1,pch=col1-1)
}
points(X[improved_output$outliers, ], col=1,pch=19)

## ----fig.show='hold',fig.height = 4.5, fig.width = 4.5-------------------
A <- mars_screw$SI_matrix;
B <- mars_screw$geographic_matrix
screw_index <-mars_screw$screw_index
## 
ret1=ktaucenters(X=A,K=3);
## 

## ----fig.show='hold',fig.height = 4.5, fig.width = 4.5-------------------
screw_candidate_index=c()
for (j in 1:3){
  dj=ktaucenters(B[ret1$cluster==j, ],K=5,nstart=1,startWithROBINPD = FALSE)$di; 
  jcandidate=dj==max(dj);
  jcandidate=which(ret1$cluster==j)[dj==max(dj)]
  screw_candidate_index=c(screw_candidate_index,jcandidate)
}

## ----fig4,fig.height = 5, fig.width = 5, fig.align = "center", fig.cap = "\\label{fig:fig3} Plot of Geographic sub matrices, pink circles are the outliers candidates in the geographic matix space."----

plot(B, type="n" )
for (j in 1:3){
#  col1=which(orderInd==j)
  col1=j
  points(B[ret1$cluster==j,],col=col1+1,pch=19)
}
points(B[screw_candidate_index,1],B[screw_candidate_index,2],col=6,pch=1,cex=3,lwd=3)

## ----fig5,echo=TRUE,fig.show='hold',fig.height = 6 , fig.width = 6-------
####### Reconstruction image from package-data ##### 
d=8;
intensityvar=matrix(0,max(B[,1]+1)*d,max(B[,2])*d)
svar=matrix(0,max(B[,1]+1)*d,max(B[,2])*d)

for (l in 1:dim(A)[1]){
  i=mars_screw$geographic_matrix[l,1]
  j=mars_screw$geographic_matrix[l,2]
  posi=(i-1)*d;
  posj=(j-1)*d;
  intensityvar[posi+(1:d), posj +(1:d)]= mars_screw$SI_matrix[l,(d^2+1):(2*d^2)]
  svar[posi+(1:d), posj +(1:d)]= mars_screw$SI_matrix[l,1:(d^2)]
  #A[l,]=c(auxs[,],auxi[,])
  #B[l,]=c(i,j)
}
svar[svar<0]=0

## ------------------------------------------------------------------------
nrowIm=495 #dim(myjpg)[1] #(nrowIm pixels is equivalent to 1 en the  coordinate system)
ncolIm=664 #dim(myjpg)[2] #(ncolIm pixeles pixels is equivalent to 1 en the  coordinate system)
index2ejx=function(indexj){(1/(ncolIm-1))*(indexj-1) +1 }
index2ejy=function(indexi){((2-1)/(1-nrowIm))*(indexi-nrowIm) +1 }
mapIndex2coord=function(filaij){ c(index2ejx(filaij[2]),index2ejy(filaij[1]))}
#transforming pixel position ij to equivalent xy values at the image 
xyscrew=mapIndex2coord(d*mars_screw$geographic_matrix[screw_index, ])

## ----fig6,echo=TRUE,fig.show='hold',fig.height = 6 , fig.width = 6-------
plot(1:2, type="n",yaxs="i",xaxs="i",main="Intensity variable")
rasterImage(intensityvar, 1, 1, 2, 2)
points(xyscrew[1],xyscrew[2],col=6,pch=1,cex=6,lwd=3)

plot(1:2, type="n",yaxs="i",xaxs="i",main="saturation variable")
rasterImage(svar, 1, 1, 2, 2)
points(xyscrew[1],xyscrew[2],col=6,pch=1,cex=6,lwd=3)

## ----echo=TRUE,fig.show='hold',fig.height = 6 , fig.width = 6------------
### Rebuilding data from image 
library("jpeg")
imageFile="screw2.jpeg"
par(mfrow=c(2,2))
myjpg=readJPEG(imageFile)
plot(1:2, type="n",yaxs="i",xaxs="i",main="original picture")
rasterImage(myjpg, 1, 1, 2, 2)

plot(1:2, type="n",yaxs="i",xaxs="i",main="Red Component")
myjpgR=myjpg; 
myjpgR[,,2]=myjpgR[,,3]=0
rasterImage(myjpgR, 1, 1, 2, 2)

plot(1:2, type="n",yaxs="i",xaxs="i",main="Green Component")
myjpgG=myjpg; 
myjpgG[,,1]=myjpgG[,,3]=0
rasterImage(myjpgG, 1, 1, 2, 2)

plot(1:2, type="n",yaxs="i",xaxs="i",main="Blue Component")
myjpgB=myjpg; 
myjpgB[,,1]=myjpgB[,,2]=0
rasterImage(myjpgB, 1, 1, 2, 2)

## ----echo=TRUE,fig.show='hold',fig.height = 6 , fig.width = 6------------
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
nrowIm=dim(myjpg)[1] #(nrowIm pixeles equivale a 1 en mi sist. de coordenadas)
ncolIm=dim(myjpg)[2] #(ncolIm pixeles equivale a 1 en mi sist. de coordenadas)
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

## ------------------------------------------------------------------------
mean(abs(A-mars_screw$SI_matrix))==0

