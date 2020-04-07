library(rcdk)
library(gbcode)
library(foreach)
library(doParallel)
sumw<-function(x,w){
        return(sum(x*w))
}

spl<-function(s){
as.numeric(strsplit(s,"\\.")[[1]][1])

    }
B=50

ncores=min(10,B)
cl <- makeForkCluster(ncores)
registerDoParallel(cl)
ntree=100
nfs=40
K=20
U=5
test=FALSE
datafiles=c("dmetrics1.Rdata") #,"dmetrics2.Rdata","data7.Rdata")#,"data4.Rdata","data3.Rdata")
aPts=NULL
RQ=read.csv(file="SMILES_COVID19-Sheet1.csv",sep=',')
drugs=RQ[,2]
for( dat in datafiles)
    for (featsel in c(TRUE)){
        load(dat)
        print(dat)
        N=length(Y)    
        
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
        if (!test)
            for (i in 1:K){
                set.seed(i)
                #I=read.table(paste("split_",i-1,".csv",sep=""),sep=",",skip=1)
                Itr=sample(N,round(9*N/10))
                Its=setdiff(1:N,Itr)

                print(length(Itr))
                
                ##cat(length(which(Y[Itr]==1)),length(which(Y[Its]==1)),"\n")
                I1=Itr[which(Y[Itr]==1)]
                I0=Itr[which(Y[Itr]==0)]
                I0b=sample(I0,length(I1))
                
                aP=NULL
                wP=NULL
                FF<-foreach (b=1:B) %dopar%{
                ##  for (b in 1:B){
                set.seed(b)
                    I0b=sample(I0,min(length(I0),sample(max(1,U-3):U,1)*length(I1)))
                    I0bts=setdiff(I0,I0b)
                    fs=1:n
                    if (featsel)
                        fs=mrmr(sX[c(I0b,I1),],as.numeric(Y[c(I0b,I1)]),nmax=min(c(8*nfs,n-10)))
                    
                    X=sX[c(I0b,I1),fs]
                    YY=factor(Y[c(I0b,I1)])
                    d<-data.frame(YY,X)
                    
                    names(d)[1]<-"Y"
                    mod.rf<-randomForest(Y~.,data=d,ntree=ntree,importance=featsel)

                    if (TRUE){
                        fs=fs[sort(importance(mod.rf,type=1),decr=TRUE,index=TRUE)$ix[1:nfs]]
                        print(fs)
                        X=sX[c(I0b,I1),fs]
                        #n=NCOL(X)
                        d<-data.frame(YY,X)
                        names(d)[1]<-"Y"
                        mod.rf<-randomForest(Y~.,data=d,ntree=ntree)
                    }
                    X.ts=sX[Its,fs]
                    d.ts<-data.frame(X.ts)
                    names(d.ts)[1:NCOL(X)]<-names(d)[2:(NCOL(X)+1)]
                    P<-predict(mod.rf,d.ts,type="prob")
                    
                    list(P=P[,2],fs=fs)
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
                L=length(which(afs>0))
                cc=colnames(X)[sort(afs,decr=TRUE,index=TRUE)$ix[1:min(L,50)]]
                
                ndrugs=unique(unlist(lapply(cc,spl)))
                seldrugs=drugs[ndrugs]
                print(seldrugs)
                cat("Lsel=", length(ndrugs),
                    "Lneg=",length(which(ndrugs>130)),
                    " phyper=", phyper(length(which(ndrugs>130)),130,70,length(ndrugs)),
                    "\n")
                print(cc)
                
                if (i==K)
                    cat("data=",dat,":featsel=",featsel, ":",summary(aAUC),"\n")
            }
        
        
        if (test){
            sXts=scale(Xts,attr(sX,"scaled:center"),attr(sX,"scaled:scale"))
            I1=which(Y==1)
            I0=which(Y==0)
            
            wP=NULL
            aP=NULL
            FF<-foreach (b=1:B) %dopar%{
                
                I0b=sample(I0,min(length(I0),U*length(I1)))
                P=pred("rf",X[c(I0b,I1),],factor(Y[c(I0b,I1)]),Xts,class=TRUE,ntree=ntree)
                aP=cbind(aP,P$prob[,2])
                ##I0bts=setdiff(I0,I0b)
                
                ##P0=pred("rf",sX[c(I0b,I1),],factor(Y[c(I0b,I1)]),sX[I0bts,],class=TRUE,ntree=4000)
                ##wP=c(wP,length(which(P0$prob[,1]>0.75))/length(I0bts)) ## the higher the better
                
                list(P=P$prob[,2])
            }
            for ( f in 1:length(FF)){
                aP=cbind(aP,FF[[f]]$P)
                aPts=cbind(aPts,FF[[f]]$P)
            }
            ## wP=wP+1e-5
            ##wP=wP/sum(wP)
            #R3<-read.table("test_predictions_sample.csv",sep=",",comment.char="%",skip=1)
            #R3[,2]=apply(aP,1,mean)
            #colnames(R3)=c("smiles","activity")
            #write.table(R3,file=paste("subm",dat,featsel,"csv",sep="."),
             #           row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")
                                        #print("saved")
            print(dim(aPts))
        }
    }
    

       
if (test){        
    R3<-read.table("test_predictions_sample.csv",sep=",",comment.char="%",skip=1)
    
    print(dim(aPts))
    R3[,2]=apply(aPts,1,mean)
    colnames(R3)=c("smiles","activity")
    write.table(R3,file=paste("newsubm","csv",sep="."),
                row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")
}
