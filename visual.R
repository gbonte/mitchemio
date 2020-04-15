load("dmetricsfs1.Rdata")
I0=which(Y==0)
I1=which(Y==1)
RQ=read.csv(file="SMILES_COVID19-Sheet1.csv",sep=',')
drugs=RQ[,2]
snames=c("51.extended.petke" , "52.extended.petke","51.shortestpath.petke",
         "51.shortestpath.tanimoto" , "51.extended.tanimoto","51.shortestpath.jaccard" ,
         "52.shortestpath.petke" ,"51.shortestpath.dice" ,
         "51.shortestpath.baroniurbanibuser", "51.extended.dice" ,"51.shortestpath.achiai",
         "12.circular.baroniurbanibuser",     "147.extended.hamming",
         "51.extended.jaccard"    ,           "12.extended.petke" , "179.maccs.simpson",
         "12.circular.achiai" ,"147.extended.hamann" , "50.circular.baroniurbanibuser",
         "51.extended.achiai", "52.shortestpath.tanimoto",          "51.shortestpath.mt" ,
         "12.circular.petke", "52.shortestpath.mt", "83.maccs.russelrao",
         "51.extended.mt", "12.extended.tanimoto", "12.extended.achiai",
         "100.shortestpath.petke","51.shortestpath.kulczynski2",
         "52.extended.mt", "96.shortestpath.hamming","52.shortestpath.kulczynski2",
         "75.extended.hamming","100.shortestpath.tanimoto",
         "52.extended.tanimoto","69.extended.hamming",
         "134.extended.hamming", "100.extended.tanimoto",
         "52.circular.achiai", "110.maccs.kulczynski2",
         "110.maccs.tanimoto", "52.circular.petke",
         "52.circular.dice", "115.maccs.simpson"   ,
         "91.maccs.simpson","75.hybridization.tanimoto",
         "163.extended.tanimoto" )
spl<-function(s){
   as.numeric(strsplit(s,"\\.")[[1]][1])
   
}
ndrugs=unlist(lapply(snames,spl))
seldrugs=drugs[ndrugs]
print(unique(seldrugs))



for (i in 1:length(ndrugs)){
   s=snames[i]
   boxplot(X[I0, s],X[I1, s],main=paste(seldrugs[i],s))
   browser()
}

if (FALSE){
   snames=colnames(X)
   
   ndrugs=unlist(lapply(snames,spl))
   seldrugs=drugs[ndrugs]
   print(unique(seldrugs))
   
   
   
   for (i in 1:length(ndrugs)){
      s=snames[i]
      boxplot(X[I0, s],X[I1, s],main=paste(seldrugs[i],ndrugs[i],s))
      browser()
   }
}
