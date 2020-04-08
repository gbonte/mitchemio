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
ntree=500
nfs=40
K=20
U=5

datafiles=c("dmetrics1.Rdata") 
aPts=NULL
RQ=read.csv(file="SMILES_COVID19-Sheet1.csv",sep=',')
drugs=RQ[,2]
dat=datafiles

load(dat)
print(dat)
N=length(Y)    

aBER=NULL
aAUC=NULL
aAUC2=NULL
set.seed(0)
           
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

fX=X  ## X without the constant vars
fXts=Xts ## Xts without the constant vars

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
        
        fs=mrmr(sX[c(I0b,I1),],as.numeric(Y[c(I0b,I1)]),nmax=min(c(8*nfs,n-10)))
        
        X=sX[c(I0b,I1),fs]
        YY=factor(Y[c(I0b,I1)])
        d<-data.frame(YY,X)
        
        names(d)[1]<-"Y"
        mod.rf<-randomForest(Y~.,data=d,ntree=ntree,importance=TRUE)
        
        if (TRUE){
            fs=fs[sort(importance(mod.rf,type=1),decr=TRUE,index=TRUE)$ix[1:nfs]]
            print(fs)
            X=sX[c(I0b,I1),fs]
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
    print(unique(seldrugs))
    cat("Lsel=", length(ndrugs),
        "Lneg=",length(which(ndrugs>130)),
        " phyper=", phyper(length(which(ndrugs>130)),130,70,length(ndrugs)),
        "\n")
    print(cc)
    
    if (i==K)
        cat("data=",dat, ":",summary(aAUC),"\n")
}

FS=which(afs>0)
X=fX[,FS]
Xts=fXts[,FS]
save(file='demetricsfs1.Rdata',list=c("X","Y","Xts"))

                
 

    

       
