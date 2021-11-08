
######## Chi-square Plot for a given Data######
X=read.table("/media/shamim/New Volume/Data Science 2nd Semester/Multivariate analysis/Assignment/T1-8.DAT",header = FALSE)
new_row_26 <- c(0.8912, 0.812, 1.7912, 1.712, 0.712, 0.612)
new_row_27 <- c(0.9112, 0.912, 1.8112, 1.812, 0.612, 0.812)
new_row_28 <- c(0.8412, 0.912, 1.8912, 1.712, 0.712, 0.712)
new_row_29 <- c(0.8812, 0.812, 1.9112, 1.812, 0.812, 0.612)
new_row_30 <- c(0.8312, 0.912, 1.9412, 1.912, 0.612, 0.812)
X <- rbind(X, new_row_26)
X <- rbind(X, new_row_27)
X <- rbind(X, new_row_28)
X <- rbind(X, new_row_29)
X <- rbind(X, new_row_30)


X=as.matrix(X)
p = 6   #numbers of variables 
X_bar=colMeans(X)
S = cov(X)
n=nrow(X)  #number of observations
Dsq=matrix(0,nrow=n,ncol=1)
for (i in 1:n) {
  Dsq[i] = t(X[i,]-X_bar)%*%solve(S)%*%(X[i,]-X_bar)
}
dj=sort(Dsq, decreasing = F)
j=matrix(seq(1,n), nrow=n, ncol=1)
qj = qchisq((n-j+0.5)/n, p)
plot(qj,dj)




#Hotelling T2
x = X
mu0=matrix(c(0.80, 0.80, 1.70, 1.70, 0.70, 0.70),6,1) ### From the q
xbar=colMeans(x)
S=cov(x)
S_inv=solve(S)
n=nrow(x)  #number of observations
T2=n*t(xbar-mu0)%*%S_inv%*%(xbar-mu0)
p=ncol(x)  #numbers of variables 


F_value=qf(.95, ncol(x), nrow(x)-ncol(x))### alpha=.05, 1-alpha=1-.05=.95
critical_value=((n-1)*p)/(n-p)*F_value
critical_value
T2 ### If T2 > critical_value then H0 will be rejected at alpha

######### Principal Component Analysis #############
turtles=log(X[1:30,1:6]) ### turtles is the data set
xbar=colMeans(turtles)
S=cov(turtles)
fit <- prcomp(S) #princomp(S)
fit ### To get the PC after rotation
summary(fit) ## To get the proportions
screeplot(fit, npcs = 6, type = "lines")
eigen(S) ### To get the PC without rotation
var_explained = fit$sdev^2 / sum(fit$sdev^2)
library(ggplot2)

qplot(c(1:6), var_explained) + 
  geom_line() + 
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("Scree Plot") +
  ylim(0, 1)

##########   Factor Analysis ######################
fac <- factanal(X, factors=2, method='mle', scale=T, center=T)

factanal(X, factors=2, method='PCA', scale=T, center=T)

factanal(x = X, factors = 2, method = "PCA", scale = T, center = T)



############## Cluster Analysis ##############
data1=X
#To remove any missing value that might be present in the data, type this:
data1 <- na.omit(data1)
#As we don't want the clustering algorithm to depend to an arbitrary variable unit, 
#we start by scaling/standardizing the data using the R function
sdata1 <- scale(data1)
# Dissimilarity matrix
d <- dist(t(sdata1), method = "euclidean") ### t for variable-wise


# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "single" ) #"average", "single", "complete", "ward"


# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" ) #"average", "single", "complete", "ward"


# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)



# CLuster Analysis - For a given distance matrix
d=c(2, 1, 5, 8)
hclust(d, method = "ward.D" )
