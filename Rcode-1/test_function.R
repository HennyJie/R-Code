source("./transform_function.R")
source("./data_gen.R")
source("./DCOL.r")
library(psych)
library(tsne)

tsne_test <- function(x,is_linear, dem, grps_size)
{
  #x<-data.gen(n.genes=1000, n.samples=300, n.grps=grp, aver.grp.size=grp_size, n.fun.types=6, epsilon=eps, n.depend=0)
  #summary(x)
  d=matrix(data = 0,nr=1000,nc=1000)
  
  if(is_linear==1)
  {
    x$data <- t(x$data)    #transpose the matrix
    d <- cor(x$data)       #calculate the correlative coefficient
    d <- abs(d)            #get the absolute value
    d <- transform1_(d)    #transform it int o distance matrix
  }
  else
  {
    d <- dcol(x$data)      #get the nolinear distance matrix
    sink("intersect.txt")  
    length(d)
    sink()
  }
  
  d<-transform1_(d)      #Get the similarity matrix in higher dimension
  d<-pre_process(d)      #Preprocess
  
  tsne_out<-as.matrix(tsne(d, k=dem, perplexity = grps_size))      #Use tsne to reduce dimension and get the coordinates in lower dimention
  
  Q_dist<-as.matrix(dist(tsne_out))    #Calculate the linear distance in lower dimention
  Q_01dist<-as.matrix(transform01(Q_dist))  #Transform the distances to make them in the range from 0 to 1
  Q<-as.matrix(transform1_(Q_01dist))    #Calculate the similarity to plot ROC curve
  Q_1dem<-matrix1dem(Q)
  image(Q);
  
  return(Q_1dem)
}

pca_test1 <- function(x, is_linear, dem)
{
  # x<-data.gen(n.genes=1000, n.samples=300, n.grps=grp, aver.grp.size=grp_size, n.fun.types=6, epsilon=eps, n.depend=0)
  # summary(x)
  
  d=matrix(data = 0,nr=1000,nc=1000)
  
  if(is_linear==1)
  {
    x$data <- t(x$data)    #transpose the matrix
    d <- cor(x$data)       #calculate the correlative coefficient
    d <- abs(d)            #get the absolute value
    d <- transform1_(d)    #transform it into distance matrix
  }
  else
  {
    d <- dcol(x$data)      #get the nolinear distance matrix
    sink("intersect.txt")  
    length(d)
    sink()
  }
  
  p<-princomp((1-d)^12)
  result<-as.matrix(p$loadings[,1:dem])
  #e<-eigen((1-d)^12)
  #result<-e$vector[,1:dem]
  
  Q_dist<-as.matrix(dist(result))    #Calculate the linear distance in lower dimention
  Q_01dist<-as.matrix(transform01(Q_dist))  #Transform the distances to make them in the range from 0 to 1
  Q<-as.matrix(transform1_(Q_01dist))    #Calculate the similarity to plot ROC curve
  Q_1dem<-matrix1dem(Q)
  
  return(Q_1dem)
}

pca_test2 <- function(x, is_linear, dem)
{
  # x<-data.gen(n.genes=1000, n.samples=300, n.grps=grp, aver.grp.size=grp_size, n.fun.types=6, epsilon=eps, n.depend=0)
  # summary(x)
  
  d=matrix(data = 0,nr=1000,nc=1000)
  
  if(is_linear==1)
  {
    x$data <- t(x$data)    #transpose the matrix
    d <- cor(x$data)       #calculate the correlative coefficient
    d <- abs(d)            #get the absolute value
    d <- transform1_(d)    #transform it into distance matrix
  }
  else
  {
    d <- dcol(x$data)      #get the nolinear distance matrix
    sink("intersect.txt")  
    length(d)
    sink()
  }
  
  #p<-princomp((1-d)^12)
  #result<-as.matrix(p$loadings[,1:dem])
  e<-eigen((1-d)^12)
  result<-e$vector[,1:dem]
  
  Q_dist<-as.matrix(dist(result))    #Calculate the linear distance in lower dimention
  Q_01dist<-as.matrix(transform01(Q_dist))  #Transform the distances to make them in the range from 0 to 1
  Q<-as.matrix(transform1_(Q_01dist))    #Calculate the similarity to plot ROC curve
  Q_1dem<-matrix1dem(Q)
  
  return(Q_1dem)
}



