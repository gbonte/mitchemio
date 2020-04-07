library(rcdk)
library(gbcode)
R<-read.table("train.csv",sep=",",comment.char="%",skip=1)
R2<-read.table("test.csv",sep=",",comment.char="%",skip=1)
Y=factor(R[,3])
library(foreach)
library(doParallel)
ncores=50
cl <- makeForkCluster(ncores)
registerDoParallel(cl)

## tanimoto euclidean hamming meanHamming
## soergel patternDifference, variance size shape jaccard dice mt simple russelrao rodgerstanimoto
## cosine achiai carbo baroniurbanibuser kulczynski2 robust hamann yule pearson mcconnaughey stiles simpson
## petke tversky


descriptor<-function(smiles){
   
   mols <- parse.smiles(as.character(smiles))
   dc <- get.desc.categories()
   troubles=c("org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor",
              "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor",
              "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor",
              "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor")
   
   dc <- get.desc.categories()
   descs=NULL
   for (dci in 1:length(dc)) {
      descNames <- get.desc.names(dc[dci])
      
      for (nn in descNames){
         if (!is.element(nn,troubles)){
            cat("dci=",dci,":",nn,"\n")
            descsi <- eval.desc(mols, nn)
            if (!is.null(descs))
               descs<-cbind(descs,descsi)
            else
               descs<-descsi
            print(dim(descs))
         }
      }
   }
   
   print(dim(descs))
   return(descs)
}  
dmetrics<-function(smiles,queries){   
   descs=NULL
   
   metrics=c('tanimoto','euclidean','hamming','meanHamming','patternDifference','variance','size','shape',
             'jaccard','dice','mt','simple','russelrao','rodgerstanimoto','achiai','carbo','baroniurbanibuser',
             'hamann','yule','pearson','stiles','simpson','petke','kulczynski2')
   # problems with  'robust' 'tversky'  'mcconnaughey',
  
  
   print(length(queries))
   types= c('circular','extended','graph','maccs','shortestpath','hybridization','estate')
   cnt=1
   for ( q in queries){
      print(q)
      query.mol <- parse.smiles(q)[[1]]
      target.mols <- parse.smiles(as.character(smiles))
      fd=fingerprint::distance
      for (ty in types){
         print(ty)
         
         query.fp <- get.fingerprint(query.mol, type=ty)
         target.fps <- lapply(target.mols, get.fingerprint, type=ty)
         
         
         FF<-foreach(me=metrics) %dopar%{
            library(fingerprint)
            sims <- data.frame(sim=do.call(rbind, lapply(target.fps,
                                                         fingerprint::distance,
                                                         fp2=query.fp, method=me)))
            #print(sims[1])
            list(sm=sims[1]) 
         }
         for (f in 1:length(FF))
            if (is.null(descs)){
               descs=FF[[f]]$sm
               colnames(descs)=paste(cnt,ty,metrics[f],sep=".")
            }else{
               qdescs=FF[[f]]$sm
               colnames(qdescs)=paste(cnt,ty,metrics[f],sep=".")
               descs=cbind(descs,qdescs) #sims[,1])
               
            }
         
         print(dim(descs))
        
      }
      cat(".")
      
      cnt=cnt+1
   } ## for q
   
   
   
   return(descs)
}



RQ=read.csv(file="SMILES_COVID19-Sheet1.csv",sep=',')
queries=as.character(RQ[,3])



I1=which(Y==1)
smiles1=as.character(R[I1,2])

N1=NROW(R)
mX=dmetrics(c(as.character(R[,2]),as.character(R2[,1])),queries)
X=mX[1:N1,]
Xts=mX[(N1+1):NROW(mX),]
save(file="dmetrics1.Rdata", list=c("N1","X","Y","Xts"))


if (FALSE){
    dX=descriptor(c(as.character(R[,2]),as.character(R2[,1])))
    X=dX[1:N1,]
    Xts=dX[(N1+1):NROW(dX),]
    save(file="descr.Rdata", list=c("N1","X","Y","Xts"))
    
    
    
    
    mX2=dmetrics(c(as.character(R[,2]),as.character(R2[,1])),smiles1)
    X=mX2[1:N1,]
    Xts=mX2[(N1+1):NROW(mX2),]
    
    save(file="dmetrics2.Rdata", list=c("N1","X","Y","Xts"))
}
