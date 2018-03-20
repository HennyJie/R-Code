
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
    
