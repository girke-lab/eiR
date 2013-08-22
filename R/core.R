
library(snow)

DataDir = "data"
TestQueries = file.path(DataDir,"test_query.iddb")
TestQueryResults=file.path(DataDir,"chemical-search.results")
ChemPrefix="chem"
ChemDb = file.path(DataDir,paste(ChemPrefix,".db",sep=""))
ChemIndex = file.path(DataDir,paste(ChemPrefix,".index",sep=""))
Main = file.path(DataDir,"main.iddb")


#debug=TRUE
debug=FALSE

# Notes
#  Need function to produce descriptors from sdf or smile
#  Need function to compute distances between descriptors

cdbSize <- function(dir=".") {
	#TODO: make this more efficient
	length(readIddb(file.path(dir,Main)))
}
embedCoord <- function(s,len,coords) 
	.Call("embedCoord",s,as.integer(len),as.double(coords))

embedCoordTest <- function(r,d,refCoords,coords) 
	.Call("embedCoordTest",as.integer(r),as.integer(d),as.double(refCoords),as.double(coords))


lshPrep <- function(matrixFile,
	W=NA,H=NA,M=NA,L=NA,K=NA,T=NA,R=NA,dir=".") {

	if(!file.exists(matrixFile)) stop(paste("could not find matrix file:",matrixFile))

	#TODO: see about removeing old .lshindex files
	matrixHash = digest(file = matrixFile,serialize=FALSE)
	matrixDir=dirname(matrixFile)
	indexName = file.path(matrixDir,paste(paste(matrixHash,W,H,M,L,T,R,sep="-"),"lshindex",sep="."))
	if(debug) print(paste("index name: ",indexName))
	indexName
}

# requires one query per column, not per row
lshsearch <- function(queries,matrixFile,
	W=NA,H=NA,M=NA,L=NA,K=NA,T=NA,R=NA) 
{
	indexFile = lshPrep(matrixFile,W,H,M,L,K,T,R)
	.Call("lshsearch",queries,as.character(matrixFile),indexFile,
		as.double(W),as.integer(H),as.integer(M),as.integer(L),
		as.integer(K),as.integer(T), as.double(R))
}
lshsearchAll <- function(matrixFile,
	W=NA,H=NA,M=NA,L=NA,K=NA,T=NA,R=NA) 
{

	indexFile = lshPrep(matrixFile,W,H,M,L,K,T,R)

	.Call("lshsearchAll",as.character(matrixFile),indexFile,
		as.double(W),as.integer(H),as.integer(M),as.integer(L),
		as.integer(K),as.integer(T), as.double(R))
}



# functions needed for sql backend:
# distance

# descriptorStr  raw format (sdf,smile) -> descriptor object -> string
# str2Descriptor  string -> descriptor object
# also need descrptor type, ie, "ap" "fpap", etc.

# and compound format, ei, "SDF", "SMILE", etc.
# X optionally: compound -> string and string -> compound

eiInit <- function(compoundDb,dir=".",format="sdf",descriptorType="ap",append=FALSE,
						 conn=defaultConn(dir))
{
	if(!file.exists(file.path(dir,DataDir)))
		if(!dir.create(file.path(dir,DataDir)))
			stop("failed to create data directory ",file.path(dir,DataDir))

	descriptorFunction = function(set)
		data.frame(descriptor=getTransform(descriptorType,"sdf")$toString(set),
					  descriptor_type=descriptorType)
	

	print("CONNECTION:")
	print(str(conn))
	if(is.null(conn))
		conn = initDb(file.path(dir,ChemDb))
	print(str(conn))
	
	if(tolower(format) == "sdf"){
		compoundIds = loadSdf(conn,compoundDb, descriptors=descriptorFunction)
	}else if(tolower(format) == "smiles" || tolower(format)=="smi"){
		stop("smiles are not yet supported")
		compoundIds = loadSmiles(conn,compoundDb,descriptors=descriptorFunction)
	}else{
		stop(paste("unknown input format:",format," supported formats: SDF, SMILE"))
	}
	print(paste(length(compoundIds)," loaded by eiInit"))

	writeIddb(compoundIds,file.path(dir,Main),append=append)
	compoundIds
}
eiMakeDb <- function(refs,d,descriptorType="ap",distance=getDefaultDist(descriptorType), 
				dir=".",numSamples=cdbSize(dir)*0.1,conn=defaultConn(dir),
				cl=makeCluster(1,type="SOCK"))
{
	workDir=NA
	createWorkDir <- function(r){
		workDir<<-file.path(dir,paste("run",r,d,sep="-"))
		if(!file.exists(workDir))
			if(!dir.create(workDir))
				stop("Could not create run directory ",workDir)
	}

	if(is.null(conn))
		stop("no database connection given")

	if(is.character(refs)){ #assume its a filename
		refIds=readIddb(refs)
		r=length(refIds)
		createWorkDir(r)
		refIddb=file.path(workDir,basename(refs))
		file.copy(refs,workDir,overwrite=TRUE)
	}else if(is.numeric(refs)){
		if(length(refs)==0){ #assume its the number of refs to use
			stop(paste("variable refs must be posative, found ",refs))
		}else if(length(refs)==1){ #assume its the number of refs to use
			r=refs
			createWorkDir(r)
			refIddb=genRefName(workDir)
			refIds=genRefs(r,refIddb,dir)
		}else{ #refs is a vector of compound indexes to use a referances
			refIds=refs
			r=length(refIds)
			createWorkDir(r)
			refIddb=genRefName(workDir)
			writeIddb(refIds,refIddb)
		}
	}else{
		stop(paste("don't know how to handle refs:",str(refs)))
	}

	matrixFile = file.path(workDir,sprintf("matrix.%d-%d",r,d))

	if(d >= length(refIds))
		stop("d must be less than the number of reference compounds")
	if(file.exists(matrixFile))
		stop(paste("found existing",matrixFile),"stopping")
	if(!file.exists(file.path(dir,Main)))
		stop(file.path(dir,Main)," not found. Did you run eiInit first?")
	

	queryIds=genTestQueryIds(numSamples,dir,refIds)
	mainIds <- readIddb(file.path(dir,Main))
	#print("queryids")
	#print(queryIds)


	selfDistFile <- paste(refIddb,"distmat",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	ref2AllDistFile <- paste(refIddb,"distances",sep=".")
	embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))
	embeddedQueryFile <- file.path(workDir,sprintf("coord.query.%d-%d",r,d))

	#embed references in d dimensional space 
	coords <- if(file.exists(coordFile)){
		message(paste("re-using coordfile",coordFile))
		as.matrix(read.table(coordFile))
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			message("generating selfDistFile")
			IddbVsIddbDist(conn,refIds,refIds,distance,descriptorType,file=selfDistFile)
		}
		selfDist<-read.table(selfDistFile)
		#print(head(selfDist))
		#print(tail(selfDist))

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}
	#compute dist between refs and all compounds
	if(!file.exists(ref2AllDistFile))
		IddbVsIddbDist(conn,readIddb(file.path(dir,Main)),
							refIds,distance,descriptorType,file=ref2AllDistFile)
	
	#each job needs: R, D, coords, a chunk of distance data
	solver <- getSolver(r,d,coords)	
	distConn <- file(ref2AllDistFile,"r")


	numJobs=length(cl)
	jobSize = as.integer(cdbSize(dir) / numJobs + 1) #make last job short

	dataBlocks = Map(function(x)
		strsplit(readLines(distConn,jobSize),"\\s+"),1:numJobs)

	clusterApply(cl,1:numJobs, 
		function(i) { # job i has indicies [(i-1)*jobSize+1, i*jobSize]
			solver <- getSolver(r,d,coords)	
			data = sapply(dataBlocks[[i]],function(x) 
								embedCoord(solver,d,as.numeric(x)))

			write.table(t(data),
				file=file.path(workDir,paste(r,d,i,sep="-")),
				row.names=F,col.names=F)

			#list indexes for this job, see which of them are queries, 
			#then shift indexes back to this jobs range before selecting 
			#from data.
			#print(mainIds[((i-1)*jobSize+1):(i*jobSize)] %in% queryIds)
			selected = which( mainIds[((i-1)*jobSize+1):(i*jobSize)]  
										%in% queryIds ) - ((i-1)*jobSize)
			#selected = queryIds[queryIds %in% 
									  #mainIds[((i-1)*jobSize+1):(i*jobSize)]] - 
								#((i-1)*jobSize)
			#print(selected)
			# R magically changes the data type depending on the size, yay!
			qd = if(length(selected)==1) 
						t(data[,selected]) else t(data)[selected, ]
			write.table(qd ,
				file=file.path(workDir,paste("q",r,d,i,sep="-")),
				row.names=F,col.names=F)
		})
	close(distConn)

	unlink(c(embeddedFile,embeddedQueryFile))
	for(x in 1:numJobs){
		cat(scan(file.path(workDir,paste(r,d,x,sep="-")),what="raw",sep="\n"),
			 sep="\n",file=embeddedFile, append=TRUE)
		cat(scan(file.path(workDir,paste("q",r,d,x,sep="-")),what="raw",sep="\n"),
			 sep="\n",file=embeddedQueryFile, append=TRUE)
	}

	Map(function(x) unlink(file.path(workDir,paste(r,d,x,sep="-"))),1:numJobs)
	Map(function(x) unlink(file.path(workDir,paste("q",r,d,x,sep="-"))),1:numJobs)

	binaryCoord(embeddedFile,matrixFile,d)
	binaryCoord(embeddedQueryFile,
		file.path(workDir,sprintf("matrix.query.%d-%d",r,d)),d)

	#file.path(workDir,sprintf("matrix.%d-%d",r,d))
	refIddb
}
eiQuery <- function(r,d,refIddb,queries,format="sdf",
		dir=".",descriptorType="ap",distance=getDefaultDist(descriptorType),
		asSimilarity=FALSE, K=200, W = 1.39564, M=19,L=10,T=30)
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		refIds = readIddb(refIddb)


		descriptorInfo = getTransform(descriptorType,format)$toObject(queries,dir)
		queryDescriptors = descriptorInfo$descriptors
		numQueries = length(queryDescriptors)
		queryNames = descriptorInfo$names
		#print("queryNames")
		#print(queryNames)
		stopifnot(length(queryNames)==length(queryDescriptors))

		#embed queries in search space
		embeddedQueries = embedFromRefs(r,d,refIddb, 
								  t(IddbVsGivenDist(conn,refIds,queryDescriptors,distance,descriptorType)))

		#search for nearby compounds
		#if(debug) print(embeddedQueries)
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		hits = search(embeddedQueries,matrixFile,
							queryDescriptors,distance,dir,descriptorType=descriptorType,K=K,W=W,M=M,L=L,T=T)
		#if(debug) print("hits")
		#if(debug) print(hits)

		targetIds=unlist(lapply(1:length(hits),function(x) hits[[x]][,1]))
		#targetIds=targetIds[targetIds!=-1]
		targetIds=targetIds[!is.na(targetIds)]
		targetNames=as.matrix(getNames(targetIds,dir))
		rownames(targetNames)=targetIds
		#print(paste(targetIds,targetNames))


		#numHits=sum(sapply(hits,function(x) sum(x[,1]!=-1)))
		numHits=sum(sapply(hits,function(x) !is.na(sum(x[,1]))))
		#print(paste("numHits:",numHits))
		#fetch names for queries and hits and put in a data frame
		results = data.frame(query=rep(NA,numHits),
								  target = rep(NA,numHits),
								  distance=rep(NA,numHits),
								  target_ids = rep(NA,numHits))
		i=1
		lapply(1:numQueries,function(queryIndex)
			lapply(1:(length(hits[[queryIndex]][,1])),function(hitIndex){
					results[i,"query"]<<-queryNames[queryIndex]
					results[i,"target"]<<-targetNames[ as.character(hits[[queryIndex]][hitIndex,1]),1]
					results[i,"target_ids"]<<-hits[[queryIndex]][hitIndex,1]
					d=hits[[queryIndex]][hitIndex,2]
					results[i,"distance"]<<- if(asSimilarity) 1-d else d
					i<<-i+1
			}))

		if(asSimilarity)
			names(results)=c("query","target","similarity","target_ids")

		if(debug) print("results:")
		if(debug) print(results)
		results
}

eiAdd <- function(r,d,refIddb,additions,dir=".",format="sdf",
						conn=defaultConn(dir), descriptorType="ap",
						distance=getDefaultDist(descriptorType))
{
		tmpDir=tempdir()
		workDir=file.path(dir,paste("run",r,d,sep="-"))

		# add additions to database
		compoundIds = eiInit(additions,dir,format,descriptorType,append=TRUE)
		additionDescriptors=getDescriptors(conn,descriptorType,compoundIds)
		numAdditions = length(compoundIds)
		refIds = readIddb(refIddb)

		#embed queries in search space
		embeddedAdditions= embedFromRefs(r,d,refIddb,
									t(IddbVsGivenDist(conn,refIds,additionDescriptors,distance,descriptorType)))
		#if(debug) print(dim(embeddedAdditions))
		#if(debug) print(embeddedAdditions)
		embeddedFile <- file.path(workDir,sprintf("coord.%d-%d",r,d))

		#add additions to existing coord and names files
		write.table(t(embeddedAdditions),
				file=embeddedFile, append=TRUE,row.names=F,col.names=F)
		binaryCoord(embeddedFile,file.path(workDir,sprintf("matrix.%d-%d",r,d)),d)
}

eiCluster <- function(r,d,K,minNbrs, dir=".",cutoff=NULL,
							 descriptorType="ap",distance=getDefaultDist(descriptorType),
							 conn=defaultConn(dir),
							  W = 1.39564, M=19,L=10,T=30,type="cluster",linkage="single"){

		workDir=file.path(dir,paste("run",r,d,sep="-"))
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		mainIndex = readIddb(file.path(dir,Main))
		neighbors = lshsearchAll(matrixFile,K=2*K,W=W,M=M,L=L,T=T)


		ml=length(mainIndex)
#		neighbors = array(matrix(1:ml,nrow=ml,
#								 ncol=ml,byrow=TRUE),dim=c(ml,ml,2))
#		neighbors[,,2]=-2.0
		#print(neighbors)

		refinedNeighbors=array(NA,dim=c(length(mainIndex),K))
		if(type=="matrix")
			similarities = array(NA,dim=c(length(mainIndex),K))
	
		#print("refining")
		batchByIndex(mainIndex,function(indexSet){

			#print("indexset:"); print(indexSet)
			descriptors = getDescriptors(conn,descriptorType,indexSet)

			lapply(1:length(indexSet),function(i){
			#	print(neighbors[i,,])
				#nonNegs=neighbors[i,,1]!=-1
				nonNegs=!is.na(neighbors[i,,1])
			#	print(nonNegs)
			   n=neighbors[i,nonNegs,]
			#	print(dim(n))
				#must set this manually because if only one row is
				#left after filtering, R decides its not really an
				#array anymore
				dim(n)=c(sum(nonNegs) ,2)
		      #translate from matrixFile index space to database index space
				reverseIndex=n[,1]
				names(reverseIndex)=mainIndex[n[,1]]
			   n[,1] = mainIndex[n[,1]]
				#print(reverseIndex)
				#print(n)
				#print(paste("refining",i))
				refined = refine(n,descriptors[i],K,distance,dir,descriptorType=descriptorType,cutoff=cutoff,conn=conn)
				dim(refined)=c(min(sum(nonNegs),K) ,2)
				#print("refined: ")
				#print(refined)
				#print(paste(mainIndex[i],paste(refined[,1],collapse=",")))
				refinedNeighbors[i,1:nrow(refined)]<<-
							reverseIndex[as.character(refined[,1])]
				if(type=="matrix")
					similarities[i,1:nrow(refined)] <<- 1 - refined[,2]
				#print(refinedNeighbors[i,1:nrow(refined)])
			})
		 })

		

		rownames(refinedNeighbors)=1:ml  ##
		#print("refined:")
		#print((refinedNeighbors))

		if(type=="matrix")
			return(list(indexes=refinedNeighbors,
							names=rownames(refinedNeighbors),
							similarities=similarities))

		#print("clustering")
		rawClustering = jarvisPatrick_c(refinedNeighbors,minNbrs,fast=TRUE)
		clustering = mainIndex[rawClustering]
		names(clustering) = mainIndex
		clustering
}

#expects one query per column
search <- function(embeddedQueries,matrixFile,queryDescriptors,distance,K,dir,descriptorType,...)
{
		neighbors = lshsearch(embeddedQueries,matrixFile,K=2*K,...)
		mainIds <- readIddb(file.path(dir,Main))
		#print(paste("got ",paste(dim(neighbors),callapse=","),"neighbors back from lshsearch"))
		#print("neighbors:")
		#print(neighbors)

		#compute distance between each query and its candidates	
		Map(function(i) {
			 #nonNegs=neighbors[i,,1]!=-1
			 nonNegs = ! is.na(neighbors[i,,1])
			 #print(nonNegs)
			 n=neighbors[i,nonNegs,]
			 dim(n)=c(sum(nonNegs) ,2)
			 #print(sum(nonNegs))
			 #print(n)
			 n[,1] = mainIds[n[,1]]
		#	 print("neighbors:")
		#	 print(n)
			 refine(n,queryDescriptors[i],K,distance,dir,descriptorType=descriptorType)
		  }, 1:(dim(embeddedQueries)[2]))
}
#fetch coords from refIddb.distmat.coord and call embed
embedFromRefs <- function(r,d,refIddb,query2RefDists)
{
		coordFile=paste(refIddb,"distmat","coord",sep=".")
		coords = as.matrix(read.table(coordFile))
		embed(r,d,coords,query2RefDists)
}
#take referenace coords, and distance matrix
#return the emedding of each row in the matrix
embed <- function(r,d,coords, query2RefDists)
{
		solver = getSolver(r,d,coords)
		embeddedQueries = apply(query2RefDists,c(1),
			function(x) embedCoord(solver,d,x))
}
refine <- function(lshNeighbors,queryDescriptors,limit,distance,dir,descriptorType,cutoff=NULL,
						 conn=defaultConn(dir))
{

	d = t(IddbVsGivenDist(conn,lshNeighbors[,1],
								 queryDescriptors,distance,descriptorType))

	#if(debug) print("result distance: ")
	#if(debug) print(str(d))
	lshNeighbors[,2]=d 

	if(!is.null(cutoff))
		lshNeighbors[lshNeighbors[,2] > cutoff,]=NA

	limit = min(limit,length(lshNeighbors[,2]))
	#print(paste("num dists:",length(lshNeighbors[,2]), "limit:",limit,"cutoff: ",cutoff))
	lshNeighbors[order(lshNeighbors[,2])[1:limit],]
}
getNames <- function(indexes,dir,conn=defaultConn(dir))
	getCompoundNames(conn,indexes)

writeIddb <- function(data, file,append=FALSE)
		write.table(data,file,quote=FALSE,append=append,col.names=FALSE,row.names=FALSE)
readIddb <- function(file) as.numeric(readLines(file))
readNames <- function(file) as.numeric(readLines(file))


# randomly select n reference compounds. Also sample and stash away
# numSamples query compounds that are not references for later
# testing
genTestQueryIds <- function(numSamples,dir,refIds=c())
{
	testQueryFile <-file.path(dir,TestQueries)
	mainIds <- readIddb(file.path(dir,Main))
	set=setdiff(mainIds,refIds)
	if(numSamples < 0 || numSamples > length(set)) 
		stop(paste("trying to take more samples than there are compounds available",numSamples,length(set)))
	queryIds <- sort(sample(set,numSamples))
	writeIddb(queryIds,testQueryFile)
	queryIds
}
genRefs <- function(n,refFile,dir,queryIds=c())
{
	mainIds <- readIddb(file.path(dir,Main))
	set=setdiff(mainIds,queryIds)
	if(n < 0 || n > length(set)) stop(paste("found more refereneces than compound candidates",n,length(set)))
	refIds = sort(sample(set,n))
	writeIddb(refIds,refFile)
	refIds
}
genRefName <- function(workDir)
	file.path(workDir,
				 paste(paste(sample(c(0:9,letters),32,replace=TRUE),
								 collapse=""),
						 "cdb",sep="."))
genTestQueryResults <- function(distance,dir,descriptorType,conn=defaultConn(dir))
{
	if(file.exists(file.path(dir,TestQueryResults)))
		return()

	out=file(file.path(dir,TestQueryResults),"w")
	d=IddbVsIddbDist(conn,
		readIddb(file.path(dir,TestQueries)),
		readIddb(file.path(dir,Main)),distance,descriptorType)
	if(debug) print(paste("dim(d): ",dim(d)))
	maxLength=min(dim(d)[2],50000)
	for(i in 1:(dim(d)[1]))
		cat(paste(
				paste(1:dim(d)[2],d[i,],sep=":")[order(d[i,])[1:maxLength]],
				collapse=" "),"\n",file=out)	
	close(out)
}
eiPerformanceTest <- function(r,d,distance=getDefaultDist(descriptorType),descriptorType="ap",
										conn=defaultConn(dir),
										dir=".",K=200, W = 1.39564, M=19,L=10,T=30)
{
	workDir=file.path(dir,paste("run",r,d,sep="-"))
	eucsearch=file.path(workDir,sprintf("eucsearch.%s-%s",r,d))
	genTestQueryResults(distance,dir,descriptorType)
	eucsearch2file(file.path(workDir,sprintf("matrix.%s-%s",r,d)),
				 file.path(workDir,sprintf("matrix.query.%s-%s",r,d)),
				 50000,eucsearch)

	#evaluator TestQueryResuts eucsearch-r-d recall
	evaluator(file.path(dir,TestQueryResults),eucsearch,
		file.path(workDir,"recall"))

	matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
	coordQueryFile =file.path(workDir,sprintf("coord.query.%d-%d",r,d))

	testQueryDescriptors=getDescriptors(conn,descriptorType,readIddb(file.path(dir,TestQueries)) )
	embeddedTestQueries = t(as.matrix(read.table(coordQueryFile)))
	hits = search(embeddedTestQueries,matrixFile,
						testQueryDescriptors,distance,descriptorType=descriptorType,dir,K=K,W=W,M=M,L=L,T=T)
	indexed=file.path(workDir,"indexed")
	out=file(indexed,"w")
	#if(debug) print(hits)
	for(x in hits)
		cat(paste(x[,1],x[,2],sep=":",collapse=" "),"\n",file=out)
	close(out)

	#indexed_evalutator TestQueryResults indexed indexed.performance
	write.table(compareSearch(file.path(dir,TestQueryResults),indexed),
			file=file.path(workDir,"indexed.performance"),
			row.names=F,col.names=F,quote=F)
}


#distance functions
# subset vs descriptors  db,iddb,descriptors
# subset vs subset       db,iddb1,iddb2

desc2descDist <- function(desc1,desc2,dist)
	as.matrix(sapply(desc2,function(x) sapply(desc1,function(y) dist(x,y))))


IddbVsGivenDist<- function(conn,iddb,descriptors,dist,descriptorType,file=NA){

	preProcess = getTransform(descriptorType)$toObject
	descriptors=preProcess(descriptors)

	process = function(record){
		batchByIndex(iddb,function(ids){
			outerDesc=preProcess(getDescriptors(conn,descriptorType,ids))
			record(desc2descDist(outerDesc,descriptors,dist))
		})
	}
	output(file,length(iddb),length(descriptors),process)
}
IddbVsIddbDist<- function(conn,iddb1,iddb2,dist,descriptorType,file=NA){

	#print(paste("iddb2:",paste(iddb2,collapse=",")))
	preProcess = getTransform(descriptorType)$toObject
	descriptors = preProcess(getDescriptors(conn,descriptorType,iddb2))
	#print("descriptors")
	#print(str(descriptors))
	process = function(record){
		batchByIndex(iddb1,function(ids){
			outerDesc = preProcess(getDescriptors(conn,descriptorType,ids))
			record(desc2descDist(outerDesc,descriptors,dist))
		})
	}
	output(file,length(iddb1),length(iddb2),process)
}

#choose whether to output a file or a matrix
output <- function(filename,nrows,ncols,process)
	if(!is.na(filename)){ #return result as matrix
		toFile(filename,process)	
	}else{ #write result to file
		toMatrix(nrows,ncols,process)
	}

#send data produced by body to a file
toFile <- function(filename,body){
	f = file(filename,"w")
	body(function(data) write.table(data,file=f,quote=F,sep="\t",row.names=F,col.names=F))
	close(f)
}

#send data produced by body to a matrix
toMatrix <- function(nrows,ncols,body){

	allDists = matrix(NA,nrows,ncols)
	rowCount = 1

	body(function(data){
#		if(debug) print(paste("recording data. rowCount: ",rowCount,", dim:",paste(dim(data),collapse=",")))
#		if(debug) print(dim(allDists))
		allDists[rowCount:(rowCount+dim(data)[1]-1),] <<- data
		rowCount <<- rowCount + dim(data)[1]
	})

	allDists
}

selectDescriptors <- function(type,ids){
	q=paste("SELECT compound_id, descriptor FROM descriptors JOIN descriptor_types USING(descriptor_type_id) WHERE ",
				" descriptor_type='",type,"' AND compound_id IN (", paste(ids,collapse=","),") ORDER
				BY compound_id",sep="")
	#print(q)
	q
}
getDescriptors <- function(conn,type,idList){
	data = selectInBatches(conn,idList,function(ids) selectDescriptors(type,ids))
	n=data$descriptor
	if(length(n) != length(idList)){
		if(debug) print(idList)
		stop(paste("missing some descriptors! Found only",
					  length(n),"out of",length(idList),"given ids"))
	}
	names(n)=data$compound_id
	ordered=n[as.character(idList)]
	#write.table(n,file="descriptors.out")
	ordered
}
