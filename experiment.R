data.gen<-function(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
{
  set.seed(Sys.time())
  link<-function(x, type)
  {
    x<-(x-mean(x))/sd(x)
    if(type == 1) return(x)
    if(type == 2) return(sin(2*x))
    if(type == 3) return(x^2)
    if(type == 4) return(abs(x))
    if(type == 5) return(x^3)
    if(type == 6) return(atan(4*x))
  }
  
  a<-matrix(rnorm(n.genes*n.samples),ncol=n.samples)
  curr.count<-0
  g<-new("list")
  for(i in 1:n.grps)
  {
    #		this.size<-rpois(1, aver.grp.size)
    this.size<-aver.grp.size
    if(this.size < 2) this.size<-2
    
    this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
    this.mat[1,]<-rnorm(n.samples)
    for(j in 2:this.size)
    {
      if(n.depend==0)
      {
        this.basis<-c(1, rep(0,j-2))
      }else{
        #				this.basis<-sample(c(1,0), j-1, replace=T, prob=c(min(1, n.depend/(j-1)), 1-min(1, n.depend/(j-1))))
        if(j-1 <= n.depend) 
        {
          this.basis<-rep(1, j-1)
        }else{
          this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
        }
        
      }
      if(sum(this.basis) > 0)
      {
        x<-rep(0,n.samples)
        for(k in which(this.basis == 1))
        {
          x<-x+link(this.mat[k,], sample(n.fun.types,1))*runif(1,min=-1,max=1)
        }
        #				x[x>quantile(x, 0.95)]<-quantile(x, 0.95)
        #				x[x<quantile(x, 0.05)]<-quantile(x, 0.05)
        this.mat[j,]<-x
        this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
      }else{
        this.mat[j,]<-rnorm(n.samples)
      }
    }
    if(n.depend == 0)
    {
      this.mat[1,]<-link(this.mat[1,], sample(n.fun.types,1))
      this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
    }
    
    if(curr.count+this.size <= n.genes)
    {
      a[(curr.count+1):(curr.count+this.size),]<-this.mat
      g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
    }
    curr.count<-curr.count+this.size		
  }
  a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
  
  g2<-rep(0, nrow(a))
  for(i in 1:length(g)) g2[g[[i]]]<-i
  
  r<-new("list")
  r$data<-a
  r$grps<-g2
  return(r)
}


normrow<-function(array)
{
  m<-apply(array,1,mean,na.rm=T)
  s<-apply(array,1,sd,na.rm=T)
  array<-(array-m)/s
  return(array)
}

gene.specific.null<-function(array, B=500)
{
  null.mat<-matrix(0, nrow=nrow(array), ncol=B)
  l<-ncol(array)
  d.array<-array[,1:(l-1)]
  for(i in 1:B)
  {
    this.order<-sample(l, l, replace=FALSE)
    for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
    null.mat[,i]<-apply(d.array, 1, sum)
  }
  r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
  return(r)
}

scol.matrix.order<-function(array,x) # x is the vector, a is the matrix, find ordered distance of rows.of.a|x
{
  if(is.null(nrow(array)) | nrow(array) == 1)
  {
    array<-as.vector(array)
    array<-array[order(x)]
    d<-array[2:length(array)]-array[1:(length(array)-1)]
    dd<-sum(abs(d),na.rm=T)
  }else{
    array<-array[,order(x)]
    d<-array[,2:ncol(array)]-array[,1:(ncol(array)-1)]
    dd<-apply(abs(d),1,sum,na.rm=T)
  }
  return(dd)
}

scol.matrix<-function(a, direction=2)  # when direction is 1, scol.matrix[i,j] = SCOL(a[i,], a[j,]), j|i
{
  
  rdmat<-matrix(0, ncol=nrow(a), nrow=nrow(a))
  for(j in 1:nrow(a))
  {
    rdmat[j,]<-scol.matrix.order(a, a[j,])
  }
  
  if(direction == 2)
  {
    rdmat.diff<-rdmat-t(rdmat)
    sel<-which(rdmat.diff > 0)
    rdmat[sel]<-t(rdmat)[sel]
  }
  return(rdmat)
}
gene.specific.p<-function(null.distr, new.d)
{
  for(i in 1:length(new.d))
  {
    new.d[i]<-pnorm(new.d[i], mean=null.distr[i,1], sd=null.distr[i,2], lower.tail=TRUE)
  }
  return(new.d)
}

normscore.row<-function(a)
{
  #library(coin)
  b<-t(apply(a, 1, normal_trafo))
  return(b)
}

dcol<-function(array)
{
  null.distr<-gene.specific.null(array)
  
  sim.mat<-scol.matrix(array,direction=1)  ## similarity matrix by SCOL, asymmetric, column given row
  d.tmp.mat<-sim.mat
  
  for(i in 1:nrow(sim.mat)) d.tmp.mat[,i]<-pnorm(sim.mat[,i], mean=null.distr[i,1], sd=null.distr[i,2],lower.tail=TRUE)
  d.tmp.mat
}

x<-data.gen(n.genes=1000, n.samples=300, n.grps=4, aver.grp.size=200, n.fun.types=6,
            epsilon=0.1, n.depend=0)
summary(x)

#y<-dist(x$data)
require(graphics)
u <- as.matrix(dist(x$data))
print(u[1003])
setwd("/home/cuihejie/TSNE/Rcode-1/tsne_pca")
write.table(u,"dist1.csv",sep=",")
source("/home/cuihejie/TSNE/Rcode-1/transform0-1.R")
u<-transform01(u)
setwd("/home/cuihejie/TSNE/Rcode-1/tsne_pca")
write.table(u,"dist.csv",sep=",")

#d<-dcol(x$data)
#library(tsne)

#tsne_out<-tsne(d)

#plot(tsne_out[,1],tsne_out[,2])


#sink("intersect.txt")  
#length(d)
#sink()


