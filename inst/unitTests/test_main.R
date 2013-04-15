
library(eiR)
library(snow)

options(warn=2)
test_dir="test_workspace"
r<- 50
d<- 40
N<- 100
j=1
runDir<-file.path(test_dir,paste("run",r,d,sep="-"))
fpDir=file.path(test_dir,"fp_test")
descType="ap"


test_aa.eiInit <- function() {
#	DEACTIVATED("slow")

   cleanup()

	checkData <- function(cids,dir=test_dir){
		checkTrue(file.exists(file.path(dir,"data","chem.db")))
		checkTrue(file.exists(file.path(dir,"data","main.iddb")))
		i <- readLines(file.path(dir,"data","main.iddb"))
		checkEquals(length(i),N)
		checkEquals(length(cids),N)
		sdfFromDb = getCompounds(initDb(file.path(dir,"data","chem.db")),cids)
		checkEquals(length(sdfFromDb),N)
	}

	setDefaultDistance("dummy",function(x,y) x-y)
	checkTrue(!is.null(eiR:::getDefaultDist("dummy")))

	addTransform("dummy","d2",toObject=function(x)x)
	checkTrue(!is.null(eiR:::getTransform("dummy","d2")))


   data(sdfsample)
   compoundIds = eiInit(sdfsample,descriptorType=descType,dir=test_dir)
	checkData(compoundIds)

	dir.create(fpDir)
	fpCids = eiInit(sdfsample,dir=fpDir,descriptorType="fp")
	checkData(fpCids,fpDir)
}

testRefs <- function(){
	200+c(1,2,5,8,9,10,11,17,18,19,20,23,24,25,26,29,31,33,34,36,38,43,45,46,47,48,49,51,53,66,67,70,71,72,73,74,75,77,78,79,80,81,82,83,87,88,89,91,99,100)
}
test_ba.eiMakeDb <- function() {

	#DEACTIVATED("slow")
   runChecks = function(){
      checkMatrix(".cdb$",r,1)
      checkMatrix(".cdb.distmat$",r,r)
      checkMatrix(".cdb.distmat.coord$",r,d)
      checkMatrix(".cdb.distances$",N,r)
      checkMatrix(sprintf("coord.%d-%d",r,d),N,d)
      checkMatrix(sprintf("coord.query.%d-%d",r,d),20,d)
      checkTrue(file.info(file.path(runDir,sprintf("matrix.%d-%d",r,d)))$size>0)
      checkTrue(file.info(file.path(runDir,sprintf("matrix.query.%d-%d",r,d)))$size>0)
      Map(function(x)
         checkTrue(!file.exists(file.path(runDir,paste(r,d,x,sep="-")))),1:j)
      Map(function(x)
         checkTrue(!file.exists(file.path(runDir,paste("q",r,d,x,sep="-")))),1:j)
   }

	print("by file name")
   refFile = file.path(test_dir,"reference_file.cdb")
	eiR:::writeIddb((1:r)+200,refFile)
   eiMakeDb(refFile,d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""),descriptorType=descType,dir=test_dir)
   runChecks()
	unlink(runDir,recursive=TRUE)

	print("by number")
   eiMakeDb(r,d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""),descriptorType=descType,dir=test_dir)
   runChecks()
	unlink(runDir,recursive=TRUE)

	print("by vector")
   eiMakeDb(testRefs(),d,numSamples=20,cl=makeCluster(j,type="SOCK",outfile=""),descriptorType=descType,dir=test_dir)
   runChecks()

	
}
test_ca.eiQuery <- function(){

	#DEACTIVATED("slow")
   data(sdfsample)
   refIddb = findRefIddb(runDir)
   results = eiQuery(r,d,refIddb,sdfsample[1:2],K=15,descriptorType=descType,dir=test_dir)
   checkTrue(length(results$distance) != 0)
   checkTrue(all(results$distance <= 1))
   checkEquals(results$distance[16],0)


	results=eiQuery(r,d,refIddb,203:204,format="compound_id",K=15,descriptorType=descType,dir=test_dir)
   checkEquals(results$distance[1],0)

	results=eiQuery(r,d,refIddb,c("650002","650003"), format="name",K=15,descriptorType=descType,dir=test_dir)
   checkEquals(results$distance[1],0)
   #checkEquals(results$distance[9],0) # not reliable

}

test_da.eiPerformanceTest <- function() {
	#DEACTIVATED("slow")
   eiPerformanceTest(r,d,K=22,descriptorType=descType,dir=test_dir)
   checkMatrix("chemical-search.results$",20, N,file.path(test_dir,"data"))
   checkMatrix(sprintf("eucsearch.%d-%d",r,d),20,N)
   checkMatrix("^indexed$",20,22)
   checkMatrix("indexed.performance",20,1)
}
test_ea.eiAdd<- function(){

	#DEACTIVATED("slow")
   data(example_compounds)
   cat(paste(paste(example_compounds,collapse="\n"),"\n",sep=""),file=file.path(test_dir,"example_compounds.sdf"))
   options(warn=-1)
   examples=read.SDFset(file.path(test_dir,"example_compounds.sdf"))
   options(warn=2)
   eiAdd(r,d,findRefIddb(runDir),examples[1:2],descriptorType=descType,dir=test_dir)

   results = eiQuery(r,d,findRefIddb(runDir),examples[1:2],descriptorType=descType,dir=test_dir)
   print(results)
   checkEquals(results$distance[1],0)

   eiAdd(r,d,findRefIddb(runDir),examples[4:8],descriptorType=descType,dir=test_dir)
   results = eiQuery(r,d,findRefIddb(runDir),examples[4],descriptorType=descType,dir=test_dir)
   checkEquals(results$distance[1],0)
   print(results)
}
test_fa.eiCluster <- function(){
#	DEACTIVATED("off")
	numNbrs=5
	minNbrs=2
	cutoff=0.5


	clustering=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs,cutoff=1-cutoff,descriptorType=descType,dir=test_dir)
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff

	conn = initDb(file.path(test_dir,"data","chem.db"))
	compoundIds=names(clustering)
	compoundNames=getCompoundNames(conn,compoundIds)
	names(clustering)=compoundNames
	sizes= clusterSizes(clustering)
	print(sizes)
	checkTrue(nrow(sizes) %in% c(9,10)) # 10 if eiAdd has run
	checkTrue(all(sizes[,2]==2))


	nnm=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs,type="matrix",cutoff=1-cutoff,descriptorType=descType,dir=test_dir)
   clustering = jarvisPatrick(nnm,k=numNbrs)
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff

	#print(sort(clustering))
}
test_fn.cluster_comparison <- function(){

	#DEACTIVATED("off")
	numNbrs=10
	minNbrs=2
	fast=TRUE
	cutoff=0.8

	dir=test_dir
	#dir="/home/khoran/runs/drug_bank_1000"
	#r=300
	#d=100


	clustering=eiCluster(r,d,K=numNbrs,minNbrs=minNbrs,dir=dir,cutoff=1-cutoff,descriptorType=descType)
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff

	conn = initDb(file.path(dir,"data","chem.db"))
	compoundIds=names(clustering)
	compoundNames=getCompoundNames(conn,compoundIds)
	names(clustering)=compoundNames
	#print(sort(clustering))


	#### non lsh clustering
	preProcess = eiR:::getTransform(descType)$toObject
	desc = preProcess(eiR:::getDescriptors(conn,descType,compoundIds))
	aps=if(descType=="ap") as(desc,"APset")
		 else if(descType=="fp"){
			 x=as(sapply(1:length(desc),function(i) as(desc[[i]],"character")),"FPset")
			 cid(x) = compoundNames
			 x
		 }

	cl2nnm = nearestNeighbors(aps,numNbrs=numNbrs,cutoff=cutoff)

	#print(tail(cl2nnm))


	cl2 = jarvisPatrick_c(cl2nnm$indexes,minNbrs,fast=fast)
	names(cl2)=compoundNames
	#print(cl2)


	#print(sort(clustering))
	#print(sort(cl2))


	source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/clusterIndex.R")
	#print(clusterSizes(clustering))
	#print(clusterSizes(cl2))
	ci <- cindex(clV1=clustering, clV2=cl2, minSZ=0, method="jaccard")$Jaccard_Index

	rand <- cindex(clV1=clustering, clV2=cl2, minSZ=0,
						method="rand")
	rand=rand$Rand_Index
	arand <- cindex(clV1=clustering, clV2=cl2, minSZ=0,
						 method="arand")
	arand=arand$Adjusted_Rand_Index
	print(paste("cluster similarity:",ci,rand,arand))
	checkEquals(ci,1)

}

test_fo.nnm_test  <- function(){
	DEACTIVATED("off")
	numNbrs=10
	minNbrs=2
	fast=TRUE
	

	dir="/home/khoran/runs/drug_bank_1000"
	r=300
	d=100
	#dir="."
	cutoff = 0.5

	compoundIds=eiR:::readIddb(file.path(dir,eiR:::Main))

	#compoundIds = compoundIds[1:100]

	print("computeing lsh results")
	lshnnm=lshNnm(r=r,d=d,K=numNbrs,minNbrs=minNbrs,dir=dir,cutoff=cutoff)

#	lshnnm = queriedNnm(compoundIds,r,d,numNbrs,dir)
	#print(lshnnm)

	print("lsh NNM:")
	print(tail(lshnnm,n=30))
	#print(lshnnm)

	print("computing true results")
	truennm = trueNnm(compoundIds,numNbrs,minNbrs,dir=dir,cutoff=cutoff)

	print("true NNM: ")
	print(tail(truennm,n=30))
	#print(truennm)

	#compute precsion/recall
	results=sapply(1:nrow(lshnnm),function(i){
			intersection =  intersect(truennm[i,],lshnnm[i,])
			t=length( intersection[!is.na(intersection)])
			#numFetched = sum(lshnnm[i,] != -1)
			numFetched = sum(!is.na(lshnnm[i,])) 
			numTrue = sum(!is.na(truennm[i,])) 
			#print(paste(t,numFetched,numTrue))
			p=t/numFetched
			r=t/numTrue
			#print(paste(t,numFetched,numTrue,p,r))
			c(p,r,  2* (p*r)/(p+r)  )
				 })
	print(results)
	print(paste(sum(results[1,])/length(results[1,]), sum(results[2,])/length(results[2,])
					,  sum(results[3,])/length(results[3,])))

}
lshNnm <- function(...){

	nnm=eiCluster(...,type="matrix",descriptorType=descType)
	nnm
}
queriedNnm <- function(compoundIds,r,d,numNbrs,dir){

		refIddb=findRefIddb(file.path(dir,paste("run",r,d,sep="-")))
		results = eiQuery(r,d,refIddb,compoundIds,format="compound_id",K=numNbrs,dir=dir,descriptorType=descType)
		cidToPosition = 1:length(compoundIds)
		names(cidToPosition) = as.character(compoundIds)
		print(results)

		t(sapply(seq(along=compoundIds),function(i){
					#print(paste("i:",i,"cid: ",compoundIds[i],paste(which(results$query==compoundIds[i]),collapse=",")))
					cidToPosition[as.character(results[results$query==compoundIds[i],"target_ids"])]
				 }))


}
trueNnm <- function(compoundIds,numNbrs,minNbrs,dir,cutoff=NA){

	conn = initDb(file.path(test_dir,"data","chem.db"))
	preProcess = eiR:::getTransform(descType)$toObject
	aps=as(preProcess(eiR:::getDescriptors(conn,descType,compoundIds)),"APset")
	#cid(aps)=compoundNames
	#cid(aps)=as.character(compoundIds)
	cid(aps)=as.character(1:length(compoundIds))


	nnm = nearestNeighbors(aps,cutoff=cutoff,numNbrs=numNbrs)
	d=dim(nnm)
	nnm=as.numeric(nnm)
	dim(nnm)=d
	rownames(nnm)=cid(aps)

	nnm$indexes
}

clusterSizes <- function(clustering) {
	sizes=Reduce(rbind,lapply(unique(clustering),function(cid)
							  cbind(cid=cid,size=sum(clustering==cid))))
	sizes[sizes[,2]>1,]
}



cleanup<- function(){
   unlink(test_dir,recursive=T)
   dir.create(test_dir)
#   setwd(test_dir) # this breaks check
   #junk <- c("data","example_compounds.sdf","example_queries.sdf",paste("run",r,d,sep="-"),fpDir)
   #junk <- c("example_compounds.sdf","example_queries.sdf",paste("run",r,d,sep="-"))
   #unlink(junk,recursive=T)
}
findRefIddb <- function(runDir){
   matches<-dir(runDir,pattern=".cdb$",full.names=T)
   checkEquals(length(matches),1)
   matches[1]
}
checkMatrix <- function(pattern,x,y,dir=runDir){
	#print(paste("searching for ",pattern))
   matches<-dir(dir,pattern=pattern,full.names=T)
	#print(matches)
   checkEquals(length(matches),1)
   file <- matches[1]
   checkTrue(file.info(file)$size>0)
   checkEquals(dim(read.table(file)),c(x,y))
}

