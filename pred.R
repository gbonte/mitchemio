library(rcdk)
library(gbcode)

sumw<-function(x,w){
        return(sum(x*w))
}

spl<-function(s){
as.numeric(strsplit(s,"\\.")[[1]][1])

    }
B=50
ntree=1000
nfs=20
K=20
U=5
test=FALSE
datafiles=c("dmetricsfs1.Rdata","descr.Rdata","mix") #,"dmetrics2.Rdata","data7.Rdata")#,"data4.Rdata","data3.Rdata")
aPts=NULL

for( dat in datafiles[2])
    for (featsel in c(FALSE)){

        if (dat=="mix"){
            load(datafiles[1])
            aX=X
            load(datafiles[2])
            X=cbind(aX,X)
        } else
        load(dat)
        
        print(dat)
        N=length(Y)
        
        library(foreach)
        library(doParallel)
        ncores=40
        cl <- makeForkCluster(ncores)
        registerDoParallel(cl)
        
        
        
        aBER=NULL
        aAUC=NULL
        aAUC2=NULL
        set.seed(0)
                                        #X=X[,1:300]
        sdX=apply(X,2,sd)
        w<-which(is.na(sdX))
        if (length(w)>0){
            X=X[,-w]
            Xts=Xts[,-w]
        }
        
        sdX=apply(X,2,sd)
        
        w<-which(sdX<0.05)
        if (length(w)>0){
            X=X[,-w]
            Xts=Xts[,-w]
        }
        
        set.seed(0)
        sX=scale(X)
        cat(dat,"dim=",dim(X),"\n")
        n=NCOL(X)
        
        aY=NULL
        aYh2=NULL
        aYh=NULL
        aAUC=NULL
        aAUC2=NULL
        afs=numeric(n)

        for (i in 1:K){
            set.seed(i)
     
            Itr=sample(N,round(8*N/10))
            Its=setdiff(1:N,Itr)
            
            print(length(Itr))
            
            ##cat(length(which(Y[Itr]==1)),length(which(Y[Its]==1)),"\n")
            I1=Itr[which(Y[Itr]==1)]
            I0=Itr[which(Y[Itr]==0)]
            I0b=sample(I0,length(I1))
            
            aP=NULL
            wP=NULL
            FF<-foreach (b=1:B) %dopar%{
                set.seed(b)
                I0b=sample(I0,min(length(I0),sample(max(1,U-3):U,1)*length(I1)))
                I0bts=setdiff(I0,I0b)
                fs=1:n
                if (featsel)
                    fs=mrmr(sX[c(I0b,I1),],as.numeric(Y[c(I0b,I1)]),nmax=min(c(nfs,n-10)))
                
                P=pred("rf",sX[c(I0b,I1),fs],factor(Y[c(I0b,I1)]),sX[Its,fs],
                       class=TRUE,ntree=ntree)
               
                list(P=P$prob[,2],fs=fs)
               
            }
            for ( f in 1:length(FF)){
                aP=cbind(aP,FF[[f]]$P)
                afs[FF[[f]]$fs]=afs[FF[[f]]$fs]+seq(length(FF[[f]]$fs),1,by=-1)
            }
            
            ##print(AUC(Y[Its],P$prob[,2]))
            ##   wP=wP+1e-5
            ##   wP=wP/sum(wP)
            ##   Phat=apply(aP,1,sumw,wP)
            Phat=apply(aP,1,mean)
            Phat2=Phat
            aY=c(aY,as.character(Y[Its]))
            aYh=c(aYh,Phat)
            aYh2=c(aYh2,Phat2)
            iAUC=AUC(factor(aY),aYh)
            aAUC<-c(aAUC,iAUC)
            iAUC2=AUC(factor(aY),aYh2)
            aAUC2<-c(aAUC2,iAUC2)
            aBER=BER(factor(aY),factor(round(aYh)))
            cat("i=",i,"BER=",aBER,"AUC=",iAUC,mean(aAUC),"AUC2=",iAUC2,mean(aAUC2),"\n")
            
            
            if (i==K)
                cat("data=",dat,":featsel=",featsel, ":",summary(aAUC),"\n")
        }
        
        
        
    }


       
