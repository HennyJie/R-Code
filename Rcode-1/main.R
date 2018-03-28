library(methods)
source("./DCOL.r")
source("./test_function.R")
source("./transform_function.R")
source("./data_gen.R")
library(psych)
library(tsne)
library(pROC)
library(parallel)
library(foreach)
library(doParallel)

n.node=25
circle_time=25
cl<-makeCluster(n.node)
registerDoParallel(cl)
# registerDoRNG(123)
# cl.cores <- detectCores()-1
# cl <- makeCluster(cl.cores,type = "FORK")

m01_2_200 <- gen01(2,200)
m01_4_200 <- gen01(4,200)
m01_8_100 <- gen01(8,100)

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
        #tsne tests
        #timestart<-Sys.time()
library(doRNG)        
cat("No., T-SNE, PCA_prin, PCA_eigen","\n",file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"),append=TRUE)
        reslist<-foreach(e = 1:circle_time,.combine="rbind")%dorng% 
        {
          library(tsne)
          library(pROC)
          #generate data
          x<-data.gen(n.genes=1000, n.samples=300, n.grps=grps_num, aver.grp.size=grps_size, n.fun.types=6, epsilon=eps, n.depend=0)
          summary(x)
          
          Q<-tsne_test(x, is_linear, dem, grps_size)
          # Q1<-pca_test1(x, is_linear, dem)
          Q2<-pca_test2(x, is_linear, dem)
          #result<-roc(unlist(m01), unlist(Q), plot=FALSE, print.thres=TRUE, print.auc=TRUE)
          #res=auc(result)
          res<-myROC(Q,m01)
          # res1<-myROC(Q1,m01);
          res2<-myROC(Q2,m01);
          c(res, res2)
        }
        res=0.0
        # res1=0.0
        res2=0.0
        for(ii in 1:circle_time)
        {
          cat(ii, reslist[ii,1], reslist[ii,2],  "\n", file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"), append=TRUE)
          res=res+reslist[ii,1];
          # res1=res1+reslist[ii,2];
          res2=res2+reslist[ii,2];
        }
        res=res/circle_time
        # res1=res1/circle_time
        res2=res2/circle_time
        cat("Average", res, res2, "\n", file=paste("./grps_",grps_num,"_grp_size_",grps_size,"_epsilon_",eps,"_demention_",dem,"_islinear_",is_linear,".txt"), append=TRUE)
        # cat("tsne: grps=", grps_num, ", epsilon=", eps, ", is_linear=", is_linear, ", demension=", dem, ", the result is ", res,"\n", file="./test.txt", append=TRUE)
        
        # res2=0.0
        # for(kk in 1:circle_time)
        # {
        #   res2=res2+reslist[[kk]][2];
        # }
        # res2=res2/circle_time
        # cat("pca:  grps=", grps_num, ", epsilon=", eps, ", is_linear=", is_linear, ", demension=", dem, ", the result is ", res2,"\n", file="./test.txt", append=TRUE)
        # 
        #timeend<-Sys.time()
        #runningtime<-timeend-timestart
        #cat("tsne_runtime: ", runningtime,"\n", file="/home/ilab/tsne_pca/Rcode-1/test.txt", append=TRUE)      
        
        #pca tests
        #timestart<-Sys.time()
        # reslist1<-foreach(e = 1:3)%dopar%
        # {
        #   Q1<-pca_test1(grps_num, grps_size, eps, dem, is_linear)
        #   #result<-roc(unlist(m01), unlist(Q), plot=FALSE, print.thres=TRUE, print.auc=TRUE)
        #   #res=auc(result)
        #   res1=myROC(Q1,m01);
        # }
        # reslist2<-foreach(e = 1:3)%dopar%
        # {
        #   Q2<-pca_test2(grps_num, grps_size, eps, dem, is_linear)
        #   #result<-roc(unlist(m01), unlist(Q), plot=FALSE, print.thres=TRUE, print.auc=TRUE)
        #   #res=auc(result)
        #   res2=myROC(Q2,m01);
        # }
        # res1=0.0
       
        # for(jj in 1:3)
        # {
        #   res1=res1+reslist1[[jj]][1];
        # }
        # res1=res1/3
        # cat("pca_princomp:  grps=", grps_num, ", epsilon=", eps, ", is_linear=", is_linear, ", demension=", dem, ", the result is ", res1,"\n", file="./test.txt", append=TRUE)
        
        #timeend<-Sys.time()
        #runningtime<-timeend-timestart
        #cat("pca_runtime: ", runningtime,"\n", file="/home/ilab/tsne_pca/Rcode-1/test.txt", append=TRUE)   
      }
    }
  }
}
stopCluster(cl)

