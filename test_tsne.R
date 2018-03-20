source("./transform_function.R")
source("./data_gen.R")
source("./DCOL.r")

library(tsne)


x<-data.gen(n.genes=1000, n.samples=300, n.grps=2, aver.grp.size=200, n.fun.types=6,
            epsilon=0.4, n.depend=0)
summary(x)             #Generate samples

d<-dcol(x$data)        #Use dcol to calculate the nonlinear distance
d<-transform1_(d)      #Get the similarity matrix in higher dimension
d<-pre_process(d)      #Preprocess

tsne_out<-as.matrix(tsne(d, perplexity = 50))      #Use tsne to reduce dimension and get the coordinates in lower dimention
setwd("/home/cuihejie/TSNE/Rcode-1/tsne_data")
plot(tsne_out[,1],tsne_out[,2], type='n' )
dots=matrix(data = 0,nr=1000000,nc=1)
for(i in 1:1000)
  dots[i,1]=i

text(tsne_out[,1],tsne_out[,2], dots)
Q_dist<-as.matrix(dist(tsne_out))    #Calculate the linear distance in lower dimention
Q_01dist<-as.matrix(transform01(Q_dist))  #Transform the distances to make them in the range from 0 to 1
Q<-as.matrix(transform1_(Q_01dist))    #Calculate the similarity to plot ROC curve
Q_1dem<-matrix1dem(Q)

write.table(Q_1dem,"grp2_0.4_2dem_50.csv",sep=",")
