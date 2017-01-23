
library(snow)

DataDir = "data"
TestQueries = file.path(DataDir,"test_query.iddb")
DescriptorSpaceDists=file.path(DataDir,"chemical-search.results")
ChemPrefix="chem"
ChemDb = file.path(DataDir,paste(ChemPrefix,".db",sep=""))
ChemIndex = file.path(DataDir,paste(ChemPrefix,".index",sep=""))
Main = file.path(DataDir,"main.iddb")

tmessage = function(...) message(Sys.time(),": ",...)

debug=TRUE
#debug=FALSE

embedCoord <- function(s,len,coords){
	.Call("embedCoord",s,as.integer(len),as.double(coords))
}

embedCoordTest <- function(r,d,refCoords,coords) 
	.Call("embedCoordTest",as.integer(r),as.integer(d),as.double(refCoords),as.double(coords))


annoySearch <- function(queries,matrixFile,dimension,numNeighbors,searchK=-1)
{
	message("creating annoy index")
	index <-new(AnnoyEuclidean,dimension)
	message("loading matrix file")
	index$load(matrixFile)
	message("done. num items in index: ",index$getNItems())

	numQueries = dim(queries)[2]
	message("found ",numQueries," queries")
	result = array(0,c(numQueries,numNeighbors,2))
	for(i in seq(numQueries)){
		neighbors = index$getNNsByVectorList(queries[,i],numNeighbors,searchK,TRUE)
		numResults = length(neighbors$item)
		#message("got ",numResults," neighbors for query ",i)
		if(numResults < numNeighbors){ #pad any missing neighbors with NA
			neighbors$item[(numResults+1):numNeighbors]=NA
			neighbors$distance[(numResults+1):numNeighbors]=NA
		}
		result[i,,1] = neighbors$item + 1 #shift to 1 indexed
		result[i,,2] = neighbors$distance
	}
	result
}
#find neighbors for everything in database
annoySearchAll <- function(matrixFile,dimension,numNeighbors,searchK=-1)
{
	message("creating annoy index
			  ",matrixFile,dimension,numNeighbors,searchK)
	index <-new(AnnoyEuclidean,dimension)
	message("loading matrix file")
	index$load(matrixFile)
	numItems = index$getNItems()
	message("done. num items in index: ",numItems)

	result = array(0,c(numItems,numNeighbors,2))
	for(i in seq(numItems)){
		neighbors = index$getNNsByItemList(i-1,numNeighbors,searchK,TRUE)
		numResults = length(neighbors$item)
		#message("got ",numResults," neighbors for query ",i)
		if(numResults < numNeighbors){ #pad any missing neighbors with NA
			neighbors$item[(numResults+1):numNeighbors]=NA
			neighbors$distance[(numResults+1):numNeighbors]=NA
		}
		result[i,,1] = neighbors$item + 1 #shift to 1 indexed
		result[i,,2] = neighbors$distance
	}
	result
}

eiInit <- function(inputs,dir=".",format="sdf",descriptorType="ap",append=FALSE,
						 conn=defaultConn(dir,create=TRUE),updateByName=FALSE,
						 cl=NULL,connSource=NULL)
{

	if(!file.exists(file.path(dir,DataDir)))
		if(!dir.create(file.path(dir,DataDir)))
			stop("failed to create data directory ",file.path(dir,DataDir))

	descriptorFunction = function(set)
		data.frame(descriptor=getTransform(descriptorType,"sdf")$toString(set,conn,dir),
					  descriptor_type=descriptorType)
	
	if(debug) message("input type: ",class(inputs))

	if(is.null(conn))
		stop("no database connection found")

	ensureSchema(conn)
	
	if(tolower(format) == "sdf"){
		loadFormat=loadSdf
	}else if(tolower(format) == "smiles" || tolower(format)=="smi"){
		loadFormat=loadSmiles
	}else{
		stop(paste("unknown input format:",format," supported formats: SDF, SMILE"))
	}

	connSource
	loadInput = function(input){
		withConnection(connSource,function(conn){
			if(is.character(input)) message("loading ",input)
			ids = loadFormat(conn,input, descriptors=descriptorFunction,updateByName=updateByName)
			if(is.character(input))
				message("loaded ",length(ids)," compounds from ",input)
			else
				message("loaded ",length(ids)," compounds")
			ids
		})
	}

		

	if(!is.null(cl) && is.character(inputs)){ #if its a list of filenames, use the cluster
		if(is.null(connSource))
			stop("a connSource must be provided when using a cluster")
		message("using cluster")
		withConnection(connSource,function(conn)
							addDescriptorType(conn,descriptorType))
		compoundIds = unlist(clusterApplyLB(cl,inputs, loadInput))
	}else{
		if(debug) message("loading locally")
		connSource= conn
		if(is.character(inputs) && length(inputs) > 1) # list of filenames
			compoundIds=unlist(lapply(inputs,loadInput))
		else
			compoundIds = loadInput(inputs)
	}


	print(paste(length(compoundIds)," loaded by eiInit"))

	writeIddb(conn,compoundIds,file.path(dir,Main),append=append)
	if(length(compoundIds)!=0 ){
		descIds = getDescriptorIds(conn,compoundIds,descriptorType)
		setPriorities(conn,forestSizePriorities,descIds,cl=cl,connSource=connSource)
	}
	compoundIds
}
eiMakeDb <- function(refs,d,descriptorType="ap",distance=getDefaultDist(descriptorType), 
				dir=".",numSamples=getGroupSize(conn,name=file.path(dir,Main))*0.1,conn=defaultConn(dir),
				cl=makeCluster(1,type="SOCK",outfile=""),connSource=NULL,numTrees=100)
{
	conn
	workDir=NA
	createWorkDir <- function(r){
		workDir<<-file.path(dir,paste("run",r,d,sep="-"))
		if(!file.exists(workDir))
			if(!dir.create(workDir))
				stop("Could not create run directory ",workDir)
		workDir <<- normalizePath(workDir)
	}
	if(debug) message("createWorkDir envir: ",ls(environment(createWorkDir)))

	if(is.null(conn))
		stop("no database connection given")

	if(debug) message("reading ",file.path(dir,Main))
	mainIds <- readIddb(conn,file.path(dir,Main))

	if(is.character(refs)){ #assume its a filename
		refIds=readIddbFile(refs)
		r=length(refIds)
		createWorkDir(r)
	}else if(is.numeric(refs)){
		if(length(refs)==0){ #assume its the number of refs to use
			stop(paste("variable refs must be positive, found ",refs))
		}else if(length(refs)==1){ #assume its the number of refs to use
			r=refs
			createWorkDir(r)
			refIds=genRefs(r,mainIds)
		}else{ #refs is a vector of compound indexes to use a referances
			refIds=sort(refs)
			r=length(refIds)
			createWorkDir(r)
		}
	}else{
		stop(paste("don't know how to handle refs:",str(refs)))
	}
	refGroupName = genGroupName(refIds)
	refGroupId = writeIddb(conn,refIds,refGroupName)
	matrixFile = file.path(workDir,sprintf("matrix.%d-%d",r,d))

	if(d >= length(refIds))
		stop("d must be less than the number of reference compounds")
	if(file.exists(matrixFile))
		stop(paste("found existing",matrixFile),"stopping")
	

	#create run and embedding context
	#create embedding
	# need: name, d, r,, descriptor type,refIds
	embeddingId = getEmbeddingId(conn,refGroupName,r,d,descriptorType,refGroupId,create=TRUE)

	if(debug) message("generating test query ids")
	queryIds=genTestQueryIds(numSamples,dir,mainIds,refIds)
	queryGroupId = writeIddb(conn,queryIds,file.path(dir,TestQueries))
	#print("queryids")
	#print(queryIds)
	if(debug) message("done generating test query ids")

	#create run
	# need: name,embedding id, compound group id, sample group id
	mainGroupId = getCompoundGroupId(conn,file.path(dir,Main))
	runId = getRunId(conn,workDir,embeddingId,mainGroupId,queryGroupId,create=TRUE)


	selfDistFile <- paste(file.path(workDir,refGroupName),"distmat",sep=".")
	selfDistFileTemp <- paste(selfDistFile,"temp",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	ref2AllDistFile <- paste(file.path(workDir,refGroupName),"distances",sep=".")
	ref2AllDistFileTemp	 <- paste(ref2AllDistFile,"temp",sep=".")

	#embed references in d dimensional space 
	coords <- if(file.exists(coordFile)){
		message(paste("re-using coordfile",coordFile))
		as.matrix(read.table(coordFile))
	}else{
		#compute pairwise distances for all references
		if(! file.exists(selfDistFile)){
			if(debug) message("generating selfDistFile")
			IddbVsIddbDist(conn,refIds,refIds,distance,descriptorType,file=selfDistFileTemp)
			file.rename(selfDistFileTemp,selfDistFile)
		}
		selfDist<-read.table(selfDistFile)
		#print(head(selfDist))
		#print(tail(selfDist))

		#use MDS and pairwise distances to embed references
		coords <- cmdscale(selfDist,d)
		write.table(coords,file=coordFile,row.names=F,col.names=F)
		coords
	}
	if(is.null(connSource))
		cl=NULL
	embedAll(conn,runId,refIds=refIds,
				coords=coords,distance=distance,
				cl=cl,connSource=connSource)

	writeMatrixFile(conn,runId,dir=dir,cl=cl,connSource=connSource,numTrees=numTrees)
	writeMatrixFile(conn,runId,dir=dir,samples=TRUE,numTrees=numTrees)

	runId
}

#lshData deprecated
eiQuery <- function(runId,queries,format="sdf",
		dir=".",distance=getDefaultDist(descriptorType),
		conn=defaultConn(dir),
		asSimilarity=FALSE,K=200, searchK=-1, lshData=NULL,
		mainIds =readIddb(conn,file.path(dir,Main),sorted=TRUE))
{
		conn
		if(debug) print("eiQuery")


		if(debug) print("getting run info")
		runInfo = getExtendedRunInfo(conn,runId) 
		if(debug) print(runInfo)
		if(nrow(runInfo)==0)
			stop("no information found for ",runId)
		if(debug) print("got run info")
		r=runInfo$num_references
		d=runInfo$dimension
		refGroupName = runInfo$references_group_name
		descriptorType=getDescriptorType(conn,info =runInfo)

		workDir=file.path(dir,paste("run",r,d,sep="-"))
		refIds = readIddb(conn,groupId=runInfo$references_group_id,sorted=TRUE)

		#print("refids: "); print(refIds)

		descriptorInfo = getTransform(descriptorType,format)$toObject(queries,conn,dir)
		queryDescriptors = descriptorInfo$descriptors
		numQueries = length(queryDescriptors)
		queryNames = descriptorInfo$names
		#print("queryNames")
		#print(queryNames)
		stopifnot(length(queryNames)==length(queryDescriptors))

		#embed queries in search space
		embeddedQueries = embedFromRefs(r,d,file.path(workDir,refGroupName), 
								  t(IddbVsGivenDist(conn,refIds,queryDescriptors,distance,descriptorType)))

		#search for nearby compounds
		#if(debug) print(embeddedQueries)
		hits = search(embeddedQueries,runId,
							queryDescriptors,distance,dir=dir,conn=conn,
							lshData=lshData,K=K,
							searchK=searchK)
		#if(debug) print("hits")
		#if(debug) print(hits)

		targetIds=unlist(lapply(1:length(hits),function(x) hits[[x]][,1]))
		#targetIds=targetIds[targetIds!=-1]
		targetIds=targetIds[!is.na(targetIds)]
		targetNames=as.matrix(getNames(targetIds,dir,conn=conn))
		rownames(targetNames)=targetIds
		#print(paste(targetIds,targetNames))


		numHits=sum(sapply(hits,function(x) !is.na(sum(x[,1]))))
		if(debug) print(paste("numHits:",numHits))
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

eiAdd <- function(runId,additions,dir=".",format="sdf",
						conn=defaultConn(dir), 
						distance=getDefaultDist(descriptorType),updateByName=FALSE)
{
		conn

		runInfo = getExtendedRunInfo(conn,runId) 
		if(nrow(runInfo)==0)
			stop("no information found for run ",runId)
		if(debug) print(runInfo)
		if(debug) print("got run info")
		r=runInfo$num_references
		d=runInfo$dimension
		refGroupName = runInfo$references_group_name
		embeddingId = runInfo$embedding_id
		descriptorType=getDescriptorType(conn,info=runInfo)
		if(debug) message("initial compound group size: ", getGroupSize(conn,groupId=runInfo$compound_group_id))

		workDir=file.path(dir,paste("run",r,d,sep="-"))

		#TODO make this work for modified descriptors

		# add additions to database
		compoundIds = eiInit(additions,dir,format,descriptorType,append=TRUE,updateByName=updateByName,conn=conn)
		#print("new compound ids: "); print(compoundIds)
		#message("new compound group size: ", getGroupSize(conn,groupId=runInfo$compound_group_id))
		
		invalidateCache()

		if(length(compoundIds) != 0){
			embedAll(conn,runId,distance,dir=dir)
			writeMatrixFile(conn,runId,dir=dir)
		}
		compoundIds
}

eiCluster <- function(runId,K,minNbrs, compoundIds=c(), dir=".",cutoff=NULL,
							 distance=getDefaultDist(descriptorType),
							 conn=defaultConn(dir),
							  searchK=-1,type="cluster",linkage="single"){

		if(debug) print("starting eiCluster")

		conn
		runInfo = getExtendedRunInfo(conn,runId) 
		if(nrow(runInfo)==0)
			stop("no information found for ",runId)
		if(debug) print(runInfo)
		if(debug) print("got run info")
		r=runInfo$num_references
		d=runInfo$dimension
		descriptorType=getDescriptorType(conn,info=runInfo)
	
		if(length(compoundIds) > 0){
			if(debug) message("using custom set to cluster")
			matrixFile = writeMatrixFile(conn,runId,compoundIds)
		}else{
			workDir=file.path(dir,paste("run",r,d,sep="-"))
			matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))
		}
		mainDescriptorIds = readMatrixIndex(matrixFile)

	
		#neighbors is in matrix space
		neighbors = annoySearchAll(matrixFile,d,numNeighbors=force(2*K),searchK=searchK)

		compIds = descriptorsToCompounds(conn,mainDescriptorIds)

		ml = length(mainDescriptorIds)
		compId2Sequential = 1:ml
		names(compId2Sequential)=compIds

		message("neighbors:")
		print(neighbors)

		refinedNeighbors=array(NA,dim=c(ml,K))
		if(type=="matrix")
			similarities = array(NA,dim=c(ml,K))
	
		i=1
		batchByIndex(mainDescriptorIds,function(indexSet){

			descriptors = getDescriptorsByDescriptorId(conn,indexSet)

			lapply(1:length(indexSet),function(x){ #for each descriptor id, i
				#message("-------------- ",i," -----------------")
				#print(neighbors[i,,])
				nonNAs=!is.na(neighbors[i,,1]) #matrix space
			   n=neighbors[i,nonNAs,,drop=FALSE]   #matrix space
				#must set this manually because if only one row is
				#left after filtering, R decides its not really an
				#array anymore
			   dim(n)=c(sum(nonNAs), 2) #martrix space


				#matrix space -> descriptor space 
			   n[,1] = mainDescriptorIds[n[,1]]
				#print("neighbors: ")
				#print(n[,1])

				#print(paste("refining",i))
				refined = refine(n,descriptors[x],K,distance,dir,descriptorType=descriptorType,cutoff=cutoff,conn=conn)
				#print("refined: ")
				#print(refined)

				#descriptor space -> compound space
				refinedCompounds = compIds[as.character(refined[,1])]
				#print("refined compounds:")
				#print(refinedCompounds)

				# compound space -> sequential id space
				refinedNeighbors[i,1:nrow(refined)]<<-
							compId2Sequential[as.character(refinedCompounds)]
				if(type=="matrix")
					similarities[i,1:nrow(refined)] <<- 1 - refined[,2]
				#print("refined neighbors:")
				#print(refinedNeighbors[i,1:nrow(refined)])
				i<<-i+1
			})
		 },batchSize=1000)

		

		rownames(refinedNeighbors)=1:ml
		print("final refined:")
		print((refinedNeighbors))

		if(type=="matrix")
			return(list(indexes=refinedNeighbors,
							names=compIds,
							similarities=similarities))

		print("clustering")
		rawClustering = jarvisPatrick_c(refinedNeighbors,minNbrs,fast=TRUE)
		# sequential space -> descriptor space -> compound space
		clustering = compIds[as.character(mainDescriptorIds[rawClustering])]
		names(clustering) = compIds 
		clustering
}

searchCache = new.env()
searchCache$descriptorIds=NULL
searchCache$runId=NULL
loadSearchCache <- function(conn,runId,dir) {
		if(debug) message("loading search cache with runId ",runId)
		runInfo = getExtendedRunInfo(conn,runId) 
		if(nrow(runInfo)==0)
			stop("no information found for ",runId)
		r=runInfo$num_references
		d=runInfo$dimension
		descriptorType=getDescriptorType(conn,info=runInfo)
	
		workDir=file.path(dir,paste("run",r,d,sep="-"))
		matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))


		searchCache$runId=runId
		searchCache$dimension=runInfo$dimension
		searchCache$matrixFile = matrixFile
		searchCache$descriptorType = descriptorType
		searchCache$descriptorIds = getRunDescriptorIds(conn,runId)

		if(debug) message("done loading search cache")

}
invalidateCache <- function(){

	searchCache$descriptorIds=NULL
	searchCache$runId=NULL
	searchCache$dimension=NULL
	searchCache$matrixFile = NULL
	searchCache$descriptorType = NULL
}
#expects one query per column
#lshData is deprecated
search <- function(embeddedQueries,runId,queryDescriptors,distance,K,dir,
						 conn=defaultConn(dir),lshData=NULL,searchK=-1)
{
		if(is.null(searchCache$descriptorIds) || 
			(!is.null(searchCache$runId) && searchCache$runId != runId)){
			loadSearchCache(conn,runId,dir)
		}
		message("quries:")
		print(embeddedQueries)
		
		neighbors = annoySearch(embeddedQueries,searchCache$matrixFile,searchCache$dimension,
										numNeighbors=force(2*K),searchK=searchK)

		
		if(debug) print(paste("got ",paste(dim(neighbors),collapse=","),"neighbors back from lshsearch"))

		message("K=",K," neighbors:")
		print(neighbors)

		#select only those positions actually used so we don't have to query
		# all descriptors every time.
		neighborDescIds  = searchCache$descriptorIds[as.vector(neighbors[,,1])]
		neighborDescIds = neighborDescIds[!is.na(neighborDescIds)]
		#map each descriptor back to one compound it could have come from
		compIds = descriptorsToCompounds(conn,neighborDescIds,all=FALSE)


		#compute distance between each query and its candidates	
		Map(function(i) { # for each query 
			 nonNAs = ! is.na(neighbors[i,,1])
		#	 print(nonNAs)
			 n=neighbors[i,nonNAs,]
			 dim(n)=c(sum(nonNAs) ,2)
		#	 print(sum(nonNAs))
		#	 print(n)


		#	 message("length of n[,1]: ",length(n[,1]))
		#	 message("length of descriptorIds: ",length(searchCache$descriptorIds))
		#	 message("length of compIds: ",length(compIds))

		#	 message("mapping to descriptor ids: ")
		#	 print(searchCache$descriptorIds[n[,1]])

		#	 message("mapping to comp ids: ")
		#	 print(compIds[as.character(searchCache$descriptorIds[n[,1]])])

		#	 message("length of final mapping: ",length( compIds[as.character(searchCache$descriptorIds[n[,1]])] ))
		#	 message("descriptor Ids: ")
		#	 print(searchCache$descriptorIds)

		#	 n[,1] = compIds[as.character(searchCache$descriptorIds[n[,1]])]
			 n[,1] = searchCache$descriptorIds[n[,1]]
			 #print("comp id neighbors:")
			 #print(n)
			 refined = refine(n,queryDescriptors[i],K,distance,dir,descriptorType=searchCache$descriptorType,conn=conn)
			 refined[,1] = compIds[as.character(refined[,1])]
			 #print("refined: ")
			 #print(refined)
			 refined
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
		#solver = getSolver(r,d,t(coords))
		embeddedQueries = apply(query2RefDists,c(1),
			#function(x) embedCoord(solver,d,x))
			function(x) newEmbed(x,coords,r,d))
}
refine <- function(lshNeighbors,queryDescriptors,limit,distance,dir,descriptorType,cutoff=NULL,
						 conn=defaultConn(dir))
{

	d = t(IddbVsGivenDistByDescId(conn,lshNeighbors[,1],
								 queryDescriptors,distance,descriptorType))

	#if(debug) print("result distance: ")
	#if(debug) print(str(d))
	lshNeighbors[,2]=d 

	if(!is.null(cutoff))
		lshNeighbors[lshNeighbors[,2] > cutoff,]=NA

	limit = min(limit,length(lshNeighbors[,2]))
	#print(paste("num dists:",length(lshNeighbors[,2]), "limit:",limit,"cutoff: ",cutoff))
	# make sure it stays a matrix, R likes to change things just for the heck of it
	lshNeighbors[order(lshNeighbors[,2])[1:limit],,drop=FALSE]
}
getNames <- function(indexes,dir,conn=defaultConn(dir))
	getCompoundNames(conn,indexes,keepOrder=TRUE,allowMissing=TRUE)

writeIddbFile <- function(data, file,append=FALSE)
		write.table(data,file,quote=FALSE,append=append,col.names=FALSE,row.names=FALSE)

readIddbFile <- function(file){
	binFile=paste(file,".Rdata",sep="")
	if(file.exists(binFile) && file.info(file)$mtime < file.info(binFile)$mtime){
		if(debug) message("reading from binary iddb: ",binFile)
		f=file(binFile,"r")
		x=unserialize(f)
		close(f)
		x
	}else{
		if(debug) message("no binary iddb found, ",binFile)
		x=as.numeric(readLines(file))
		if(length(x) > 1000000){
			message("large iddb found (",length(x),"), generating binary version")
			f=file(binFile,"w")
			serialize(x,f)
			close(f)
		}
		x
	}
}
readNames <- function(file) as.numeric(readLines(file))


# randomly select n reference compounds. Also sample and stash away
# numSamples query compounds that are not references for later
# testing
genTestQueryIds <- function(numSamples,dir,mainIds,refIds=c())
{
	#testQueryFile <-file.path(dir,TestQueries)
	set=setdiff(mainIds,refIds)
	if(numSamples < 0 || numSamples > length(set)) 
		stop(paste("trying to take more samples than there are compounds available",numSamples,length(set)))
	queryIds <- sort(sample(set,numSamples))
#	writeIddb(conn,queryIds,testQueryFile)

	#REMOVE ME!
	queryIds = c( 256,265,269,274,279,280,281,284,290,297)

	queryIds
}
genRefs <- function(n,mainIds,queryIds=c())
{
	set=setdiff(mainIds,queryIds)
	if(n < 0 || n > length(set)) stop(paste("found more refereneces than compound candidates",n,length(set)))
	refIds = sort(sample(set,n))
	#writeIddb(conn,refIds,refFile)
	refIds
}
genGroupName <- function(members)
	digest(paste(members,collapse=""),serialize=FALSE)
genRefName <- function(workDir)
	file.path(workDir,
				 paste(paste(sample(c(0:9,letters),32,replace=TRUE),
								 collapse=""),
						 "cdb",sep="."))
testDescriptorSpaceDists<- function(distance,dir,descriptorType,conn=defaultConn(dir))
{
	if(file.exists(file.path(dir,DescriptorSpaceDists)))
		return()

#	out=file(file.path(dir,DescriptorSpaceDists),"w")
#	print(out)
	print("getting test query ids")
	allIds = readIddb(conn,file.path(dir,Main))
	d=IddbVsIddbDist(conn,
		readIddb(conn,file.path(dir,TestQueries)),
		allIds,distance,descriptorType)

	writeTopN(50000,d,allIds,file.path(dir,DescriptorSpaceDists))

#	if(debug) print(paste("dim(d): ",dim(d)))
#	maxLength=min(dim(d)[2],50000)
#	for(i in 1:(dim(d)[1]))
#		cat(paste(
#				#paste(1:dim(d)[2],d[i,],sep=":")[order(d[i,])[1:maxLength]],
#				paste(allIds[1:dim(d)[2]],d[i,],sep=":")[order(d[i,])[1:maxLength]],
#				collapse=" "),"\n",file=out)	
#	close(out)
}
writeTopN = function(N,d,ids,outputFile){
	message("writing top ",N," to file ",outputFile)
	out=file(outputFile,"w")
	print(out)
	if(debug) print(paste("dim(d): ",dim(d)))
	maxLength=min(dim(d)[2],50000)
	for(i in 1:(dim(d)[1])) # each row
		cat(paste(
				#paste(ids[1:dim(d)[2]],d[i,],sep=":")[order(d[i,])[1:maxLength]],
				paste(1:dim(d)[2],d[i,],sep=":")[order(d[i,])[1:maxLength]],
				collapse=" "),"\n",file=out)	
	close(out)

}
testEmbeddedSpaceDists <- function(outputFile,runInfo,embeddedTestDescriptors,dir, conn=defaultConn(dir)){
	message("computing embedded distances on test quries")

	allIds = readIddb(conn,groupId=runInfo$compound_group_id)
	testIds = readIddb(conn,file.path(dir,TestQueries))
	message(length(allIds)," total descriptors, ",length(testIds)," test quries")
	print(allIds)
	allEmbeddedDescriptors = getEmbeddedDescriptors(conn,runInfo$embedding_id,allIds)
	message("test embbedded descriptors:")
	print(str(embeddedTestDescriptors))
	message("all embbedded descriptors:")
	print(str(allEmbeddedDescriptors))
	
	dists=sapply(1:dim(allEmbeddedDescriptors)[1],function(i){
			 sapply(1:dim(embeddedTestDescriptors)[1],function(j){
				 message("query ",j,": ",paste(format(embeddedTestDescriptors[j,],digits=6),collapse=" "))
				 message("db ",i,": ",paste(format(allEmbeddedDescriptors[i,],digits=6),collapse=" "))
				 dist=sqrt(sum( (allEmbeddedDescriptors[i,]-embeddedTestDescriptors[j,])^2))
				 message("dist=",dist)
				 dist
		})})
	
	message("dists:")
	print(str(dists))
	write.table(format(dists,digits=6),"raw_embedded_dists",quote=FALSE,row.names=FALSE,col.names=FALSE)

	writeTopN(50000,dists,allIds,outputFile)
}
eiPerformanceTest <- function(runId,distance=getDefaultDist(descriptorType),
										conn=defaultConn(dir),
										dir=".", K=200, searchK=-1)
{
	conn
	message("in eiPerformanceTest ========================")

	runInfo = getExtendedRunInfo(conn,runId) 
	if(nrow(runInfo)==0)
		stop("no information found for ",runId)
	r=runInfo$num_references
	d=runInfo$dimension
	descriptorType=getDescriptorType(conn,info=runInfo)
	embeddingId = runInfo$embedding_id
	sampleGroupId = runInfo$sample_group_id

	workDir=file.path(dir,paste("run",r,d,sep="-"))
	eucsearch=file.path(workDir,sprintf("eucsearch.%s-%s",r,d))

	sampleCompoundIds = readIddb(conn,groupId=sampleGroupId)
	embeddedTestQueries = getEmbeddedDescriptors(conn,embeddingId,sampleCompoundIds)
	write.table(embeddedTestQueries,file="embeddedTestQueries",row.names=FALSE,col.names=FALSE,quote=FALSE)
	testEmbeddedSpaceDists(eucsearch,runInfo,embeddedTestQueries,dir,conn=conn)


	testDescriptorSpaceDists(distance,dir,descriptorType,conn=conn)
#	testEmbeddedSpaceDists(eucsearch,embeddingId,dir,conn=conn)

	#computes distance from each query to every other descriptor in db
	#keeps only top 50000 (or db size) nearest results
	#saves them in eucsearch file as "Id:distance" pairs.
#	eucsearch2file(file.path(workDir,sprintf("matrix.%s-%s",r,d)),
#				 file.path(workDir,sprintf("matrix.query.%s-%s",r,d)),
#				 50000,eucsearch)

	#evaluator TestQueryResuts eucsearch-r-d recall
	message("running evaluator ========================")
	evaluator(file.path(dir,DescriptorSpaceDists),eucsearch,
		file.path(workDir,"recall"))

	matrixFile =file.path(workDir,sprintf("matrix.%d-%d",r,d))

	testQueryDescriptors=getDescriptors(conn,descriptorType,sampleCompoundIds )


	hits = search(t(embeddedTestQueries),runId,
						testQueryDescriptors,distance,dir=dir,conn=conn,
						K=K, searchK=searchK)
	indexed=file.path(workDir,"indexed")
	out=file(indexed,"w")
	#if(debug) print(hits)
	
	for(x in hits)
		cat(paste(x[,1],x[,2],sep=":",collapse=" "),"\n",file=out)
	close(out)

	#indexed_evalutator DescriptorSpaceDists indexed indexed.performance
	results = compareSearch(file.path(dir,DescriptorSpaceDists),indexed)
	write.table(results,
			file=file.path(workDir,"indexed.performance"),
			row.names=F,col.names=F,quote=F)
	return(results)
}


#distance functions
# subset vs descriptors  db,iddb,descriptors
# subset vs subset       db,iddb1,iddb2

desc2descDist <- function(desc1,desc2,dist)
	as.matrix(sapply(desc2,function(x) sapply(desc1,function(y) dist(x,y))))


IddbVsGivenDistByDescId <- function(conn,descIds,descriptors,dist,descriptorType,file=NA){
	#message("descIds: ")
	#print(descIds)
	IddbVsGivenDistBase(function(ids) getDescriptorsByDescriptorId(conn,ids) ,
							  conn,descIds,descriptors,dist,descriptorType,file)
}
IddbVsGivenDist<- function(conn,iddb,descriptors,dist,descriptorType,file=NA){
	IddbVsGivenDistBase(function(ids) getDescriptors(conn,descriptorType,ids),
							  conn,iddb,descriptors,dist,descriptorType,file)
}
IddbVsGivenDistBase<- function(descSource,conn,iddb,descriptors,dist,descriptorType,file=NA){

	preProcess = getTransform(descriptorType)$toObject
	descriptors=preProcess(descriptors)

	process = function(record){
		batchByIndex(iddb,function(ids){
			outerDesc=preProcess(descSource(ids))
			record(desc2descDist(outerDesc,descriptors,dist))
		},batchSize=1000)
	}
	output(file,length(iddb),length(descriptors),process)
}


IddbVsIddbDist<- function(conn,iddb1,iddb2,dist,descriptorType,file=NA,cl=NULL,connSource=NULL){

	#print(paste("iddb2:",paste(iddb2,collapse=",")))
	preProcess = getTransform(descriptorType)$toObject
	descriptors = preProcess(getDescriptors(conn,descriptorType,iddb2))
	#print("descriptors")
	#print(str(descriptors))

	if(is.null(cl)){
		process = function(record,recordPart=NULL){
			batchByIndex(iddb1,function(ids){
				outerDesc = preProcess(getDescriptors(conn,descriptorType,ids))
				record(desc2descDist(outerDesc,descriptors,dist))
			},batchSize=1000)
		}
		output(file,length(iddb1),length(iddb2),process)
	}else{
		if(is.null(connSource))
			stop("the connSource parameter is required when using a cluster")


		process = function(record,recordPart=NULL){
				ip=function(ids,jobId){
					f = file(paste("job-",jobId,".out",sep=""),"a")
					message("in indexProcessor: ",jobId)
					cat("in indexProcessor: ",jobId,"\n",file=f); flush(f)
					#print(ids)
					ret=recordPart({
						message("starting")
						cat("recording","\n",file=f); flush(f)

						# this must be done here to ensure connSource() is evaluated
						# before getDescriptors starts to run

						outerDesc = withConnection(connSource, function(conn){
											preProcess(getDescriptors(conn,descriptorType,ids))
										})

				#		conn=NULL
				#		tryCatch({
				#				conn=connSource()
				#				outerDesc = preProcess(getDescriptors(conn,descriptorType,ids))
				#			},
				#			error=function(e) stop(e),
				#			finally= if(!is.null(conn)) dbDisconnect(conn)
				#		)
						cat("got descriptors","\n",file=f); flush(f)
						d2d=desc2descDist(outerDesc,descriptors,dist)
						cat("got distances","\n",file=f); flush(f)
						d2d
					},jobId)
					cat("done with ",jobId,"\n",file=f); flush(f)
					close(f)
					ret
				}

			clusterExport(cl,c("descriptors","getDescriptors","desc2descDist","dist","descriptorType"),envir=environment())

			ipEnv = new.env(parent=globalenv())

			ipEnv$recordPart=recordPart
			ipEnv$connSource=connSource
			ipEnv$preProcess=preProcess
			ipEnv$withConnection=withConnection

			environment(ip)=ipEnv

			parBatchByIndex(iddb1,cl=cl,batchSize=10000,
				indexProcessor=ip ,
				reduce = function(results){
					message("evaluationg results")
					results # force evaluation here
					message("in reduce")
					sapply(results,function(result) {
						if(is.character(result))
							record(read.table(result))
						else
							record(result)
						unlink(result)
						})
				})
		}
		#absPath=file.path(getwd(),file)
		#absPath=normalizePath(file)
		absPath=file #given file is now absolute already
		#print(absPath)

		output(absPath,length(iddb1),length(iddb2),process,mapReduce=TRUE)
	}
}

#choose whether to output a file or a matrix
output <- function(filename,nrows,ncols,process,mapReduce=FALSE)
	if(!is.na(filename)){ #write result to file
		toFile(filename,process,mapReduce)	
	}else{ #return result as matrix
		toMatrix(nrows,ncols,process,mapReduce)
	}


#send data produced by body to a file
toFile <- function(filename,body,mapReduce){
	f = file(filename,"w")

	write = function(data) 
		write.table(data,file=f,quote=F,sep="\t",row.names=F,col.names=F)
	writePart = function(data,jobId) {
		partFilename = paste(filename,".part-",jobId,sep="")
		print(paste("partFilename: ",partFilename))

		if(!file.exists(partFilename) || file.info(partFilename)$size == 0){
			write.table(data,file=partFilename, quote=F,sep="\t",row.names=F,col.names=F)
		}else{
			message("skipping ",partFilename)
			flush(stderr())
		}

		partFilename
	}
	env = new.env(parent=globalenv())
	env$filename=filename
	environment(writePart)=env
	
	if(mapReduce)
		body(write,writePart)
	else
		body(write)
	close(f)
}
testFun <- function() 8
#send data produced by body to a matrix
toMatrix <- function(nrows,ncols,body,mapReduce){

	allDists = matrix(NA,nrows,ncols)
	rowCount = 1

	write = function(data){
		#if(debug) print(paste("recording data. rowCount: ",rowCount,", dim:",paste(dim(data),collapse=",")))
		#if(debug) print(dim(allDists))
		allDists[rowCount:(rowCount+dim(data)[1]-1),] <<- data
		rowCount <<- rowCount + dim(data)[1]
	}
	writePart = function(data)  data

	if(mapReduce)
		body(write,writePart)
	else
		body(write)

	allDists
}
#this is just a utility function used in the unit tests
embedDescriptor <- function(conn,r,d,refName,descriptor,descriptorType="ap",
									 distance=getDefaultDist(descriptorType), dir="."){
	refIddb = file.path(dir,paste("run",r,d,sep="-"),refName)
	refIds = readIddb(conn,refName,sorted=TRUE)
	embedFromRefs(r,d,refIddb,
			t(IddbVsGivenDist(conn,refIds,descriptor,distance,descriptorType)))
}

getCoords <- function(conn,runId, dir="."){ # looks for coord file in current directory
	runInfo = getExtendedRunInfo(conn,runId) 
	refGroupName = runInfo$references_group_name
	r=runInfo$num_references
	d=runInfo$dimension
	workDir = file.path(dir,paste("run",r,d,sep="-"))
	selfDistFile <- paste(file.path(workDir,refGroupName),"distmat",sep=".")
	coordFile <- paste(selfDistFile,"coord",sep=".")
	as.matrix(read.table(coordFile))
}

embedAll <- function(conn,runId, distance,dir=".", 
							refIds=readIddb(conn,groupId=runInfo$references_group_id,sorted=TRUE),
							coords = getCoords(conn,runId,dir),
							cl=NULL,
							connSource=NULL){
	runInfo = getExtendedRunInfo(conn,runId) 
	r=runInfo$num_references
	d=runInfo$dimension
	embeddingId = runInfo$embedding_id
	descriptorType=getDescriptorType(conn,info =runInfo)
	unembeddedDescriptorIds = getUnEmbeddedDescriptorIds(conn,runId)
	if(debug) message("embedding ",length(unembeddedDescriptorIds)," unembedded descriptors")

	embedJob = function(ids,jobId){
		#solver <- getSolver(r,d,t(coords))

		withConnection(connSource,function(conn){
		
			descriptors = getDescriptorsByDescriptorId(conn,ids)
			rawDists = t(IddbVsGivenDist(conn,refIds,descriptors,distance,descriptorType))
			#embeddedDescriptors = apply(rawDists,c(1), function(x) embedCoord(solver,d,x))
			embeddedDescriptors = apply(rawDists,c(1), function(x){
							 desc = newEmbed(x,coords,r,d)
							 #message("comparison| descriptor: ", paste(format(x,digits=5),collapse=","))
							 #message("comparison| embed desc: ", paste(format(desc,digits=5),collapse=","))
							 desc
						})

			insertEmbeddedDescriptors(conn,embeddingId,ids,t(embeddedDescriptors))
		})
	}
	
	if(is.null(cl)){ #don't use cluster
		if(debug) message("embedding locally")
		connSource = conn
		batchByIndex(unembeddedDescriptorIds,embedJob)
	}else{
		if(debug) message("embedding on cluster")
		#ensure we have at least as many jobs as cluster nodes, but if
		# we have a large number of compounds, batch them by no more than 10,000
		num = length(unembeddedDescriptorIds)
		numJobs = max(length(cl), as.integer(num/ 10000))
		jobSize = as.integer(num/ numJobs + 1) #make last job short

		#copy large items to nodes once, then remove them from the closures scope
		# so that they don't get copied to ndoes each time
		clusterExport(cl,c("coords","distance"),envir=environment())
		#rm(coords,unembeddedDescriptorIds,distance,envir=environment(embedJob))
		rm(coords,distance,envir=environment(embedJob))

		x=parBatchByIndex(unembeddedDescriptorIds,embedJob,
							 reduce=identity,cl=cl,batchSize=jobSize)
		if(debug) message("done with cluster embedding")
		x
	}
}
newEmbed = function(rawDists,coords,r,d){
	f =function(x) .Call("calcF",as.double(x),as.double(rawDists),as.double((coords)))
	g =function(x) .Call("calcG",as.double(x),as.double(rawDists),as.double((coords)))

	initValue = rep(0,d)
	result = optim(initValue,f,g,method="L-BFGS-B",
			lower= rep(-0.5,r) ,
			upper=  rep(0.5,r),
			control = list(pgtol=1e-2,factr = 1e12, lmm=3),
			hessian = FALSE)
	#print(str(result))
	result$par

}

checkEmbedding <- function(conn,descriptorIds,runId,distance,dir=".",
							refIds=readIddb(conn,groupId=runInfo$references_group_id,sorted=TRUE),
							coords = getCoords(conn,runId,dir),callback=NA) {
	runInfo = getExtendedRunInfo(conn,runId) 
	message("run info: ")
	print(runInfo)

	descriptorType=getDescriptorType(conn,info =runInfo)
	message("descriptor type:")
	print(descriptorType)

	r= runInfo$num_referances
	d = runInfo$dimension

	#solver <- getSolver(runInfo$num_references,runInfo$dimension,t(coords))

 	embedJob = function(ids){
		print(ids)
		descriptors = getDescriptorsByDescriptorId(conn,ids)
		#print(head(descriptors))
		#print(distance)
		#print(descriptorType)
		rawDists = t(IddbVsGivenDist(conn,refIds,descriptors,distance,descriptorType))
		#print(head(t(rawDists)))
		embeddedDesc = apply(rawDists,c(1), 
									function(x) {
										message("x dims: ",length(x))
										#embedCoord(solver, runInfo$dimension,x)
										newEmbed(x,coords,r,d)
									})
		#print(descriptors)
		#print(embeddedDesc)

		embeddedDesc
 	}

	for(descriptorId in descriptorIds){
		#fetch what is in the db
		dbEmbeddedDescriptor = t(getEmbeddedDescriptors(conn,runInfo$embedding_id,descriptorIds=descriptorId))

		#generate the embedding ourselves
		localEmbeddedDescriptor = embedJob(descriptorId)

		if(!is.na(callback)) 
			callback(descriptorId,dbEmbeddedDescriptor,localEmbeddedDescriptor)

		if( ! isTRUE(all.equal(dbEmbeddedDescriptor, localEmbeddedDescriptor))){
			warning("found mismatched embedding!")
			print("db version: ")
			print(head(dbEmbeddedDescriptor))
			print("local version: ")
			print(head(localEmbeddedDescriptor))
			#print(data.frame(dbEmbeddedDescriptor,localEmbeddedDescriptor))

		}

	}

}
