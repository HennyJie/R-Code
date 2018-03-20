transform01<-function(a){
  a<-as.matrix(a)
  maxv=0.0
  for(i in 1:1000000)
  {
      if(a[i]>maxv)
      {
        maxv=a[i]
      }
  }
  for(i in 1:1000000)
  {
      a[i]=a[i]/maxv
  }
  return(a)
}

transform1_ <-function(a){
  a<-as.matrix(a)
  for(i in 1:1000000)
  {
    a[i]=1.0-a[i]
  }
  return(a)
}

pre_process <-function(a){
  for(i in 1:1000)
  {
    for(j in 1:1000)
    {
      if(a[i,j]<a[j,i])
      {
        a[j,i]=a[i,j]
      }
    }
  }
  for(i in 1:1000)
  {
    for(j in 1:1000)
    {
      a[i,j]=a[i,j]^10;
    }
  }
  return(a)
}

matrix1dem <- function(mtx){
  k=1
  ans1=matrix(data = 0,nr=1000000,nc=1)
  for(i in 1:1000){
    for(j in 1:1000){
      ans1[k,1]=mtx[i,j]
      k=k+1
    }
  }
  return(ans1)
}

gen01 <- function(grp, grp_size)
{
  ans=matrix(data = 0, nr=1000, nc=1000)
  for(i in 1:grp)
  {
    num=(i-1)*grp_size
    for(j in 1:grp_size)
      for(k in 1:grp_size)
      {
        ans[num+j,num+k]=1
      }
  }
  
  ans<-matrix1dem(ans)
  
  return(ans)
}

myROC <- function(Q, T)
{
  x=quantile(Q,seq(0,1,by=0.01));
  xx<-matrix(nrow=1,ncol=100);
  yy<-matrix(nrow=1,ncol=100);
  for(i in 1:100)
  {
    B<- (Q>=x[i]);
    xx[i]=sum((1-T)*B)/sum(1-T);
    yy[i]=sum(T*B)/sum(T);
  }
  
  #plot(xx,yy);
  
  s=0;
  for(i in 1:99)
  {
    s=s+(xx[i+1]-xx[i])*yy[i];
  }
  return(abs(s))
}
