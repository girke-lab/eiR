
library(eiR)
library(snow)
library(DBI)

options(warn=2)
options(error=dump.frames)
test_dir="test_workspace"
r<- 50
d<- 40
N<- 100
j=1  #can only use 1 when using sqlite
runDir<-file.path(test_dir,paste("run",r,d,sep="-"))
fpDir=file.path(test_dir,"fp_test")
descType="ap"

sqliteSource=function(){
	require(eiR)
	require(RSQLite)
	path = file.path(test_dir,"data")
	if(!file.exists(path))
		dir.create(path,recursive=TRUE)
	conn=initDb(file.path(path,"chem.db"))
	#checkTrue(file.exists(file.path(test_dir,"data","chem.db")))
	conn
}

resetDb <- function(conn){
   tables = dbListTables(conn)
   #print(tables)
   for(table in tables)
      dbGetQuery(conn, paste("DROP TABLE ",table," CASCADE"))
   initDb(conn)
}

pgSource=function(){
	require(eiR)
	require(RPostgreSQL)
	conn=dbConnect(dbDriver("PostgreSQL"),user="chemminer_tests",password="40ersfdv90erijgfk",
				 dbname="chemminer_tests",host="girke-db-1.bioinfo.ucr.edu")
	initDb(conn)
	eiR:::setDefaultConn(conn)
	conn
}

connSource = sqliteSource
#connSource = pgSource

lastRunId=0

debug=TRUE


test_aa.eiInit <- function() {
#	DEACTIVATED("slow")

   cleanup()
	connSource() #sets default connection

	checkData <- function(cids,dir=test_dir){
		checkTrue(file.exists(file.path(dir,"data","main.iddb")))
		i <- readLines(file.path(dir,"data","main.iddb"))
		checkEquals(length(i),N)
		checkEquals(length(cids),N)
		sdfFromDb = getCompounds(connSource(),cids)
		checkEquals(length(sdfFromDb),N)
	}

	message("eiInit tests")

	setDefaultDistance("dummy",function(x,y) x-y)
	checkTrue(!is.null(eiR:::getDefaultDist("dummy")))

	addTransform("dummy","d2",toObject=function(x,conn=NA,dir=".")x)
	checkTrue(!is.null(eiR:::getTransform("dummy","d2")))


	message("default descriptor")
   data(sdfsample)
   compoundIds = eiInit(sdfsample,descriptorType=descType,dir=test_dir)
	#checkData(compoundIds)

	message("fp descriptor")
	dir.create(fpDir)
	write.SDF(sdfsample[1:30],file.path(fpDir,"f1"))
	write.SDF(sdfsample[31:60],file.path(fpDir,"f2"))
	write.SDF(sdfsample[61:100],file.path(fpDir,"f3"))
	# we just test with one node as SQLite does not support parallel writes
	cl=makeCluster(1,type="SOCK",outfile=file.path(test_dir,"eiInit.snow"))
	fpCids = eiInit(file.path(fpDir,c("f1","f2","f3")),dir=fpDir,descriptorType="fp",cl=cl,
						 connSource=connSource )
	stopCluster(cl)
	#checkData(fpCids)
}

testRefs <- function(){
	200+c(1,2,5,8,9,10,11,17,18,19,20,23,24,25,26,29,31,33,34,36,38,43,45,46,47,48,49,51,53,66,67,70,71,72,73,74,75,77,78,79,80,81,82,83,87,88,89,91,99,100)
}

test_ba.parDist <- function(){

	#DEACTIVATED("slow")
	conn = connSource()
	distance = eiR:::getDefaultDist("ap") 
	require(snow)
	cl = makeSOCKcluster(3,outfile="")

	set1=eiR:::readIddb(conn,file.path(test_dir,"data","main.iddb"))
	set2 = testRefs()

	if(debug){
		print("set1")
		print(set1)

		print("set2")
		print(set2)
	}
   
	eiR:::IddbVsIddbDist(conn,set1,set2, distance,"ap",file=file.path(test_dir,"parDist.final"),
								cl=cl, connSource= connSource )
	stopCluster(cl)

	checkTrue(file.exists(file.path(test_dir,"parDist.final")))


}
test_bb.eiMakeDb <- function() {

	#DEACTIVATED("slow")

	message("  eiMakeDb   ")

	conn = connSource()
	runDbChecks = function(rid){

		parameters = eiR:::runQuery(conn,paste("SELECT e.name,e.embedding_id,dimension,num_references FROM runs as r JOIN embeddings as e USING(embedding_id) 
										WHERE r.run_id = ",rid))
		#print(parameters)
		checkEquals(d,parameters$dimension)
		checkEquals(r,parameters$num_references)

		numRefs = eiR:::runQuery(conn,paste("SELECT count(*) FROM runs as r JOIN embeddings as e USING(embedding_id) 
											JOIN compound_group_members as refs ON(e.references_group_id=refs.compound_group_id)
										WHERE r.run_id = ",rid))[[1]]
		message("found ",numRefs," refs")
		checkEquals(r,numRefs)
		numCompounds = eiR:::runQuery(conn,paste("SELECT count(*) FROM runs as r  
											JOIN compound_group_members as cgm ON(r.compound_group_id=cgm.compound_group_id)
										WHERE r.run_id = ",rid))[[1]]

		checkEquals(N,numCompounds)
		numSamples = eiR:::runQuery(conn,paste("SELECT count(*) FROM runs as r  
											JOIN compound_group_members as cgm ON(r.sample_group_id=cgm.compound_group_id)
										WHERE r.run_id = ",rid))[[1]]
		checkEquals(20,numSamples)


		numDescriptors = eiR:::runQuery(conn,paste("SELECT count(distinct descriptor_id) FROM runs as r JOIN embedded_descriptors as ed USING(embedding_id) 
										WHERE r.run_id = ",rid))[[1]]
		#checkEquals(N,numDescriptors)
		#this should be 96 now because some descriptors are identical
		# and we now only store a unique set
		checkEquals(96,numDescriptors)

		for(i in 1:numDescriptors)
			checkDescriptor(conn,rid,descriptorIndex=i)

      checkTrue(file.info(file.path(runDir,sprintf("matrix.%d-%d",r,d)))$size>0)
      checkTrue(file.info(file.path(runDir,sprintf("matrix.query.%d-%d",r,d)))$size>0)


	}

	cl=makeCluster(j,type="SOCK",outfile="")

	print("by file name")
   refFile = file.path(test_dir,"reference_file.cdb")
	eiR:::writeIddbFile((1:r)+200,refFile)
   rid=eiMakeDb(refFile,d,numSamples=20,cl=cl,descriptorType=descType,dir=test_dir,
		connSource=connSource	)
   runDbChecks(rid)
	unlink(runDir,recursive=TRUE)

	print("by number")
   rid=eiMakeDb(r,d,numSamples=20,cl=cl,descriptorType=descType,dir=test_dir,
		connSource=connSource	)
   runDbChecks(rid)
	unlink(runDir,recursive=TRUE)

	print("by vector")
   rid=eiMakeDb(testRefs(),d,numSamples=20,cl=cl,descriptorType=descType, dir=test_dir,
		connSource=connSource	)
	stopCluster(cl)
   runDbChecks(rid)

	lastRunId<<-rid
	
}

test_ca.eiQuery <- function(){

	#DEACTIVATED("slow")
	message("eiQuery")
	conn=connSource()
   data(sdfsample)

	runId = lastRunId
	message("using runId ",runId)


	message("eiQuery test 1")
	results = eiQuery(runId,sdfsample[1:2],K=15,asSimilarity=TRUE,dir=test_dir)
   checkTrue(length(results$similarity) != 0)
   checkTrue(all(results$similarity>= 0))
   checkEquals(results$similarity[16],1)

	message("eiQuery test 2")
	results=eiQuery(runId,203:204,format="compound_id",K=15,asSimilarity=TRUE,dir=test_dir)
   checkEquals(results$similarity[1],1)

	message("eiQuery test 3")
	results=eiQuery(runId,c("650002","650003"), format="name",K=15,dir=test_dir)
   
   checkEquals(results$distance[1],0)
   #checkEquals(results$distance[9],0) # not reliable

	message("eiQuery test 4")
	lshData = loadLSHData(r,d,dir=test_dir)
	results=eiQuery(runId,c("650002","650003"), format="name",K=15,lshData=lshData,dir=test_dir)
	freeLSHData(lshData)


}

test_da.eiPerformanceTest <- function() {
	#DEACTIVATED("slow")
	runId = lastRunId
   eiPerformanceTest(runId,K=22,dir=test_dir)
   checkMatrix("chemical-search.results$",20, N,file.path(test_dir,"data"))
   checkMatrix(sprintf("eucsearch.%d-%d",r,d),20,N)
   #checkMatrix("^indexed$",20,22)
   checkMatrix("indexed.performance",20,1)
}
test_ea.eiAdd<- function(){

	#DEACTIVATED("slow")
	conn = connSource()

   data(example_compounds)
   cat(paste(paste(example_compounds,collapse="\n"),"\n",sep=""),
		 file=file.path(test_dir,"example_compounds.sdf"))
   options(warn=-1)
   examples=read.SDFset(file.path(test_dir,"example_compounds.sdf"))
   options(warn=2)

	runId = lastRunId
	message("using runId ",runId)

	message("eiAdd test1")
   newCompoundIds = eiAdd(runId,examples[1:2],dir=test_dir)
   results = eiQuery(runId,examples[1:2],dir=test_dir)
   print(results)
   checkEquals(results$distance[1],0)

	for(i in newCompoundIds)
		checkDescriptor(conn,runId,compoundId=i,printAll=F)

	message("eiAdd test2")
   newCompoundIds2 = eiAdd(runId,examples[4:8],dir=test_dir)
   results = eiQuery(runId,examples[4],dir=test_dir)
   checkEquals(results$distance[1],0)
   print(results)

	for(i in newCompoundIds2)
		checkDescriptor(conn,runId,compoundId=i,printAll=F)



}
test_fa.eiCluster <- function(){
	#DEACTIVATED("off")
	numNbrs=5
	minNbrs=2
	cutoff=0.5

	runId = lastRunId

	clustering=eiCluster(runId,K=numNbrs,minNbrs=minNbrs,cutoff=1-cutoff,dir=test_dir)
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff
	print(byCluster(clustering))

	conn=connSource()
	compoundIds=names(clustering)
	compoundNames=getCompoundNames(conn,compoundIds)
	names(clustering)=compoundNames
	sizes= clusterSizes(clustering)
	print(sizes)
	checkTrue(nrow(sizes) %in% c(8,9,10)) # 10 if eiAdd has run
	checkTrue(all(sizes[,2]==2))


	nnm=eiCluster(runId,K=numNbrs,minNbrs=minNbrs,type="matrix",cutoff=1-cutoff,dir=test_dir)

   clustering = jarvisPatrick(nnm,k=minNbrs,mode="a1b")
	sizes= clusterSizes(clustering)

	#print(clusterSizes(clustering))
	checkTrue(length(clustering) >= N) #eiAdd will add some stuff
	checkTrue(nrow(sizes) %in% c(8,9,10)) # 10 if eiAdd has run
	checkTrue(all(sizes[,2]==2))
}
test_fn.cluster_comparison <- function(){

	DEACTIVATED("off, out of date")
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

	conn = connSource()
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
	DEACTIVATED("off, local only,  out of date")
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

	conn = connSource()
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
	conn=connSource()
	if(inherits(conn,"PostgreSQLConnection"))
		resetDb(conn)

   #unlink(test_dir,recursive=T) #doesn't work
	system(paste("rm -rf ",test_dir))
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
	if(debug) print(paste("searching for ",pattern," expected dims: ",x,y))
   matches<-dir(dir,pattern=pattern,full.names=T)
	if(debug) print(matches)
   checkEquals(length(matches),1)
   file <- matches[1]
   checkTrue(file.info(file)$size>0)
	d = dim(read.table(file))
	if(debug) print(paste("found dims: ",d))
   checkEquals(d,c(x,y))
}

checkDescriptor = function(conn,rid,descriptorIndex=NULL,
									compoundId = eiR:::readIddb(conn,
																		 file.path(test_dir,eiR:::Main))[descriptorIndex],
									printAll=FALSE) {

		parameters = eiR:::runQuery(conn,paste("SELECT e.name,e.embedding_id,dimension,num_references FROM runs as r JOIN embeddings as e USING(embedding_id) 
										WHERE r.run_id = ",rid))
		#check that embedded descriptor values are stored in correct order
		message("checking descriptor of compound id ",compoundId)


		descId = eiR:::runQuery(conn,paste("SELECT descriptor_id FROM compound_descriptors where
												 compound_id=",compoundId))[[1]]

		desc = eiR:::getDescriptors(conn,descType,compoundId)[[1]]
		embeddedDesc = eiR:::embedDescriptor(conn,r,d,parameters$name,desc,
														 descriptorType=descType,dir=test_dir)



		dbEmbeddedDesc = eiR:::runQuery(conn,paste("SELECT value FROM embedded_descriptors WHERE
									 embedding_id=",parameters$embedding_id," and descriptor_id=",descId  ,
									 " ORDER BY ordering"))[[1]]

		#if(!all(as.vector(embeddedDesc) == dbEmbeddedDesc)){
		if(printAll){
			print(desc)
			print(as.vector(embeddedDesc))
			print(dbEmbeddedDesc)
		}
		checkEquals(as.vector(embeddedDesc),dbEmbeddedDesc)


	}

