library(methods)
source("./DCOL.r")
source("./test_function.r")
source("./transform_function.R")
source("./data_gen.R")
library(psych)
library(tsne)
library(pROC)
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

n.node=3
circle_time=3
cl<-makeCluster(n.node)
registerDoParallel(cl)

m01_2_200 <- gen01(2, 200)
m01_4_200 <- gen01(4, 200)
m01_8_100 <- gen01(8, 100)

grps_vec <- c(2,4,8)
grps_size_vec <- c(200,200,100)
eps_vec <- c(0.5,1,2)
linear_vec <- c(0, 1)
dem_vec <- c(2,4,6,8)

for(a in 1:3)   #choose the number of groups from grps_vec
{
  grps_num=grps_vec[a]
  grps_size=grps_size_vec[a]
  m01=matrix(data = 0, nr=1000000, nc=1)
  if(grps_num==2)
  {
    m01=m01_2_200
  }
  else if(grps_num==4)
  {
    m01=m01_4_200
  }
  else
    m01=m01_8_100
  
  for(b in 1:3)  #choose epsilon
  {
    eps=eps_vec[b]
    
    for(d in 1:4) #choose the demension
    {
      dem=dem_vec[d]
      
      for(c in 1:2) #linear or nolinear
      {
        is_linear=linear_vec[c]
        cat("No., T-SNE, PCA_eigen, kernal_PCA_degree2, kernal_PCA_degree4, kernal_PCA_degree6, kernal_PCA_rbf","\n",file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"),append=FALSE)
        reslist<-foreach(e = 1:circle_time,.combine="rbind")%dorng% 
        {
          library(tsne)
          library(pROC)
          library(kernlab)
          #generate data
          x<-data.gen(n.genes=1000, n.samples=300, n.grps=grps_num, aver.grp.size=grps_size, n.fun.types=6, epsilon=eps, n.depend=0)
          summary(x)
          Q<-tsne_test(x, is_linear, dem, grps_size)
          Q1<-pca_test(x, is_linear, dem)

          temp<-t(x$data)
          Q2<-kernel_pca_degree(temp, dem, 2)
          Q3<-kernel_pca_degree(temp, dem, 4)
          Q4<-kernel_pca_degree(temp, dem, 6)
          Q5<-kernel_pca_rbf(temp, dem)
          
          #result<-roc(unlist(m01), unlist(Q), plot=FALSE, print.thres=TRUE, print.auc=TRUE)
          #res=auc(result)
          res<-myROC(Q,m01)
          res1<-myROC(Q1,m01);
          res2<-myROC(Q2,m01);
          res3<-myROC(Q3,m01);
          res4<-myROC(Q4,m01);
          res5<-myROC(Q5,m01);
          c(res, res1, res2, res3, res4, res5)
        }
        res=0.0
        res1=0.0
        res2=0.0
        res3=0.0
        res4=0.0
        res5=0.0
        for(ii in 1:circle_time)
        {
          cat(ii, reslist[ii,1], reslist[ii,2], reslist[ii,3],reslist[ii,4], reslist[ii,5],reslist[ii,6],"\n", file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"), append=TRUE)
          res=res+reslist[ii,1];
          res1=res1+reslist[ii,2];
          res2=res2+reslist[ii,3];
          res3=res3+reslist[ii,4];
          res4=res4+reslist[ii,5];
          res5=res5+reslist[ii,6];
        }
        res=res/circle_time
        res1=res1/circle_time
        res2=res2/circle_time
        res3=res3/circle_time
        res4=res4/circle_time
        res5=res5/circle_time
        cat("Average", res, res1, res2, res3, res4, res5, "\n", file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"), append=TRUE)
      }
    }
  }
}
stopCluster(cl)

