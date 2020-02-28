rm(list=ls())
install.packages("ICSNP")
install.packages("devtools")
install.packages("rrcov")
install.packages("scatterplot2d")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(ICSNP)
getwd()
setwd()
data <- read.csv("/Users/christopherton/Desktop/STA135Project/Admission_Predict_Ver1.1.csv",header=TRUE,sep=",")

#************************************************************************************************************

#One sample inference about a mean vector
xbar <-colMeans( data[1:50,-1])
xvar <- var(data[1:50,-1])

#Testing a multivariate mean vector with Hotelling T^2 (page 214-215)
# testing whether mu = mu_0
p <- 8
n <- 50
nullmean <- c(320,100,4.5,4.6,4.2,8.7,0.7,0.6)
d <- xbar-nullmean

#if t2 > cval, reject nullmean
t2 <- n*t(d)%*%solve(xvar)%*%d
cval <- (n-1)*p/(n-p)*qf(0.95,p,n-p)

#p-value based on t2
t2mod <- (n-p)*t2/(p*(n-1))
pval <- 1- pf(t2mod,p,n-p)
cat("Hotelling T-squared statistic:", t2)
cat("p-value",pval)
#an alternative way
HotellingsT2(X = data[1:50,-1], mu = nullmean)

#************************************************************************************************************

# Confidence Region
data[,2] = data[,2]/340
data[,3] = data[,3]/120
datatemp<-data[1:50,c(2,3),-1]
 
conf.reg<-function(xdata,alpha){
  if(ncol(xdata)!=2) stop("Only for bivariate normal")
  n<-nrow(xdata)
  xbar<-colMeans(xdata)
  S<-cov(xdata)
  es<-eigen(S)
  e1<-es$vec %*% diag(sqrt(es$val))
  r1<-sqrt(qf(alpha,2,n-2))*sqrt(2*(n-1)/(n*(n-2)))
  theta<-seq(0,2*pi,len=250)
  v1<-cbind(r1*cos(theta), r1*sin(theta))
  pts<-t(xbar-(e1%*%t(v1)))
  plot(pts,type="l",main="Confidence Region for Bivariate Normal",xlab=colnames(xdata)[1],ylab=colnames(xdata)[2],asp=1)
  segments(0,xbar[2],xbar[1],xbar[2],lty=2) # highlight the center
  segments(xbar[1],0,xbar[1],xbar[2],lty=2)
  
  th2<-c(0,pi/2,pi,3*pi/2,2*pi)   #adding the axis
  v2<-cbind(r1*cos(th2), r1*sin(th2))
  pts2<-t(xbar-(e1%*%t(v2)))
  segments(pts2[3,1],pts2[3,2],pts2[1,1],pts2[1,2],lty=3)  
  segments(pts2[2,1],pts2[2,2],pts2[4,1],pts2[4,2],lty=3)
  
}

conf.reg(datatemp,alpha=0.95)

#************************************************************************************************************
# Compute 95% simultaneous confidence intervals for the two mean values
p<-ncol(datatemp)
n<-nrow(datatemp)
S<-cov(datatemp)
xbar<-colMeans(datatemp)
mu1.L= xbar[1] - sqrt( ((n-1)*p/(n-p))*qf(0.95,p,n-p) ) * sqrt(S[1,1]/n)
mu1.U= xbar[1] + sqrt( ((n-1)*p/(n-p))*qf(0.95,p,n-p) ) * sqrt(S[1,1]/n)
mu2.L= xbar[2] - sqrt( ((n-1)*p/(n-p))*qf(0.95,p,n-p) ) * sqrt(S[2,2]/n)
mu2.U= xbar[2] + sqrt( ((n-1)*p/(n-p))*qf(0.95,p,n-p) ) * sqrt(S[2,2]/n)
cat("95 CI ",c(mu1.L,mu1.U))
cat("95 CI ",c(mu2.L,mu2.U))

lines(c(mu1.L,mu1.L),c(0.53,mu2.U),lty=2,col=2,lwd=2)
lines(c(mu1.U,mu1.U),c(0.53,mu2.U),lty=2,col=2,lwd=2)
lines(c(0.49,mu1.U),c(mu2.L,mu2.L),lty=2,col=2,lwd=2)
lines(c(0.49,mu1.U),c(mu2.U,mu2.U),lty=2,col=2,lwd=2)

#************************************************************************************************************

# Compute 95% Bonferroni confidence intervals for the two mean values
mu1.LB=xbar[1]-qt(0.05/(2*p),n-1,lower.tail=F)*sqrt(S[1,1]/n)
mu1.UB=xbar[1]+qt(0.05/(2*p),n-1,lower.tail=F)*sqrt(S[1,1]/n)
mu2.LB=xbar[2]-qt(0.05/(2*p),n-1,lower.tail=F)*sqrt(S[2,2]/n)
mu2.UB=xbar[2]+qt(0.05/(2*p),n-1,lower.tail=F)*sqrt(S[2,2]/n)
cat("95 Bonf",c(mu1.LB,mu1.UB))
cat("95 Bonf", c(mu2.LB,mu2.UB))


#************************************************************************************************************

# Plot the confidence intervals together with the confidence ellipse:

lines(c(mu1.LB,mu1.LB),c(0.53,mu2.UB),lty=3,col=3,lwd=2)
lines(c(mu1.UB,mu1.UB),c(0.53,mu2.UB),lty=3,col=3,lwd=2)
lines(c(0.49,mu1.UB),c(mu2.LB,mu2.LB),lty=3,col=3,lwd=2)
lines(c(0.49,mu1.UB),c(mu2.UB,mu2.UB),lty=3,col=3,lwd=2)

#************************************************************************************************************
#### two-sample Hotelling's T2 test  -------

noresdata <- data[data$Research == 0,-1]
nores <- noresdata[-7]

yesresdata <- data[data$Research == 1,-1]
yesres <- yesresdata[-7]

# now we perform the two-sample Hotelling T^2-test

HotellingsT2(nores,yesres)

n<-c(50,50)
p<-7
xmean1<-colMeans(nores)
xmean2<-colMeans(yesres)
d<-xmean1-xmean2
S1<-var(nores)
S2<-var(yesres)
Sp<-((n[1]-1)*S1+(n[2]-1)*S2)/(sum(n)-2)
t2 <- t(d)%*%solve(sum(1/n)*Sp)%*%d
t2

alpha<-0.05
cval <- (sum(n)-2)*p/(sum(n)-p-1)*qf(1-alpha,p,sum(n)-p-1)
cval

#****************************************************************
# since we reject the null, we use the simultaneous confidence intervals
# to check the significant components

# simultaneous confidence intervals
alpha<-0.05
wd<-sqrt(((n[1]+n[2]-2)*p/(n[1]+n[2]-p-1))*qf(1-alpha,p,n[1]+n[2]-p-1))*sqrt(diag(Sp)*sum(1/n))
Cis<-cbind(d-wd,d+wd)
cat("95% simultaneous confidence interval","\n")
Cis

#****************************************************************
#Bonferroni simultaneous confidence intervals
wd.b<- qt(1-alpha/(2*p),n[1]+n[2]-2) *sqrt(diag(Sp)*sum(1/n))
Cis.b<-cbind(d-wd.b,d+wd.b)
cat("95% Bonferroni simultaneous confidence interval","\n")
Cis.b

# alternative way
cbind( abs(d/sqrt(diag(Sp)*sum(1/n)) ))
qt(1-alpha/(2*p),n[1]+n[2]-2)

# both component-wise simultaneous confidence intervals do not contain 0, so they have significant differences. 


#************************************************************************************************************
# Creating a data frame out of a matrix:

data.df <- data.frame(data[-1])
data.dfvar <-sqrt(diag(var(data.df)))

data.pc <- princomp(data.df, cor=T)

# Showing the coefficients of the components:
summary(data.pc,loadings=T)

# Showing the eigenvalues of the correlation matrix:
(data.pc$sdev)^2

# A scree plot:

plot(1:(length(data.pc$sdev)),  (data.pc$sdev)^2, type='b', 
     main="Scree Plot", xlab="Number of Components", ylab="Eigenvalue Size")

# Where does the "elbow" occur? 2
# What seems to be a reasonable number of PCs to use?2

par(pty="s")
plot(data.pc$scores[,1], data.pc$scores[,2], 
     xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
# labeling points with state abbreviations:
text(data.pc$scores[,1], data.pc$scores[,2], cex=0.7, lwd=2)

biplot(data.pc)

#****************************************************************

data.pca <- prcomp(data[,-1], center = TRUE, scale = TRUE)
# exclude the two categorical variables
summary(data.pca)
str(data.pca)

#Plotting the PC scores for the sample data in the space of the first two principal components
ggbiplot(data.pca)
ggbiplot(data.pca, labels=rownames(data))

#dim(nores)
#dim(yesres)
data.research <- c(rep("nores", 220), rep("yesres",280 ))
ggbiplot(data.pca,ellipse=TRUE, groups=data.research)

ggbiplot(data.pca,ellipse=TRUE,choices=c(3,4), groups=data.research)

#contribution to first principal component based on covariance matrix
es<-eigen(var(data.df))
loadings <- c(.374,.371,.349,.350,.318,.393,.265,.390) # from summary
corr <- loadings * sqrt(es$values[1]) / data.dfvar

#************************************************************************************************************
#LDA

ldadata <- rbind(nores,yesres)
ldadataa <-rbind(noresdata,yesresdata)

par(mar=c(4,4,2,1))
plot(data[,2],data[,3],xlab="GRE.score",ylab="TOEFL.score",
     pch=rep(c(18,20),times=c(220,280)),col=rep(c(2,4),times=c(220,280)),main="")
legend("bottomright",legend=c("nores","yesres"),pch=c(18,20),col=c(2,4),cex=1)


ldadata1<-ldadata[1:50,1:2] #class1
ldadata2<-ldadata[221:271,1:2] #class2
# compute sample mean vectors:
ldadata1.mean<-colMeans(ldadata1)
ldadata2.mean<-colMeans(ldadata2)
# compute pooled estimate for the covariance matrix:
S.u<-49*(var(ldadata1)+var(ldadata2))/98
w<-solve(S.u)%*%(ldadata1.mean-ldadata2.mean)
w0<-(ldadata1.mean+ldadata2.mean)%*%w/2
lines(ldadata[,2],-(w[1]*ldadata[,2]+w0)/w[2])


library(MASS)
lda.obj<-lda(ldadataa$Research~GRE.Score+TOEFL.Score,data=ldadataa,prior=c(1,1)/2)
plda<-predict(object=lda.obj,newdata=ldadataa)
#determine how well the model fits
table(ldadataa[,4],plda$class)

#plot the decision line
gmean <- lda.obj$prior %*% lda.obj$means
const <- as.numeric(gmean %*%lda.obj$scaling)
slope <- - lda.obj$scaling[1] / lda.obj$scaling[2]
intercept <- const / lda.obj$scaling[2]
#Plot decision boundary
plot(ldadataa[,1:2],pch=rep(c(18,20),times=c(220,280)),col=rep(c(2,4),times=c(220,280)))
abline(intercept, slope)
legend("bottomright",legend=c("nores","yesres"),pch=c(18,20),col=c(2,4))


#************************************************************************************************************

ldadataa.lda <- lda(Research~.,data=ldadataa)
ldadataa.lda
ldadataa.pred <- predict(ldadataa.lda)
#A Stacked Histogram of the LDA Values
ldahist(data = ldadataa.pred$x[,1], g=ldadataa$Research)
# Confusion matrix
table(ldadataa$Research,ldadataa.pred$class)

#****************************************************************

#equal costs and equal priors discriminant function (page588)
#called from ldadataa.lda
coef<- matrix(c(27.40858991,-4.65896494,0.10320122,0.04983113,0.06374563,-0.40450929,4.64907822),1,7)
group1means <- matrix(c(0.9097059 , 0.8665909,2.563636, 2.918182 ,3.095455 ,8.234727,0.6349091),7,1)
group2means <- matrix(c(0.9473739,0.9142262,3.546429, 3.732143, 3.789286 ,8.844929 ,0.7899643),7,1)

ybar1<- coef %*%group1means
ybar2 <- coef %*% group2means

#midpoint between means
m= (1/2)* (ybar1 + ybar2)




