
withConnection <- function(connSource,connUser){

	#given a connection object directly
	if(inherits(connSource,"DBIConnection")){
		connUser(connSource)
	}else{
		conn=connSource()
		tryCatch({
			connUser(conn)
		  },error = function(e) stop(e),
		  finally = dbDisconnect(conn)
		)
	}
}



ensureSchema <- function(conn) {

	tableList=dbListTables(conn)

#	print(tableList)
	if( ! all(c("compound_groups","compound_group_members","embeddings","runs","embedded_descriptors") 
				 %in% tableList)) {
		if(debug) print("ammending db")

		sqlFile = file.path("schema",if(inherits(conn,"SQLiteConnection")) "data.SQLite" 
								  else if(inherits(conn,"PostgreSQLConnection")) "data.RPostgreSQL")
																	
		nocomments = function(line) !grepl("^\\s*--",line)
		noblank = function(line) !grepl("^\\s*$",line)
		statements = Filter(noblank,unlist(strsplit(paste(
							  Filter(nocomments,readLines(system.file(sqlFile,package="eiR",mustWork=TRUE))),
							  collapse=""),";",fixed=TRUE)))
#		print(statements)

		Map(function(sql) runQuery(conn,sql),statements)
	}

}

getEmbeddingId <- function(conn,name,r,d,descriptorType,refGroupId,create=FALSE){

	embeddingId = getOrCreate(conn,
									  paste("SELECT embedding_id FROM embeddings WHERE name='", name,"'",sep=""),
									  paste("INSERT INTO
											  embeddings(name,dimension,num_references,descriptor_type_id,references_group_id)",
											  "VALUES('",name,"',",d,",",r,
													",(SELECT descriptor_type_id FROM descriptor_types WHERE
															  descriptor_type='",descriptorType,"'),",
													refGroupId,")",sep=""),
									  errorTag=paste("embedding",name),create=create )
	embeddingId
}
getEmbedding <-function(conn,embeddingId){
	runQuery(conn,paste("SELECT embedding_id,name,dimension,num_references,descriptor_type_id,references_group_id ",
								 "FROM embeddings WHERE embedding_id = ",embeddingId))
}
getRunId <- function(conn,name,embeddingId,mainGroupId,queryGroupId,create=FALSE) {
	runId = getOrCreate(conn,
							  paste("SELECT run_id FROM runs WHERE embedding_id
									  =",embeddingId," AND compound_group_id =
									  ",mainGroupId,sep=""),
							  paste("INSERT INTO runs(name,embedding_id,compound_group_id,sample_group_id)", 
									  "VALUES('",name,"',",embeddingId,",",mainGroupId,",",queryGroupId,")",sep=""),
							  errorTag = paste("run "+name),create=create)
	runId
}
getRun <- function(conn,runId){
	runQuery(conn,paste("SELECT run_id,name, embedding_id,compound_group_id,sample_group_id ",
								 "FROM runs WHERE run_id = ",runId))
}
getExtendedRunInfo <-function(conn,runId){
	if( is.null(runId) || length(runId)==0) stop("no runId given in getExtendedRunInfo")
	runQuery(conn,paste("SELECT r.run_id,r.name, r.embedding_id,r.compound_group_id, ccg.name as compound_group_name,",
								 "		r.sample_group_id, rcg.name as sample_group_name, ",
								 "		e.name as embedding_name,e.dimension,e.num_references,e.descriptor_type_id,e.references_group_id,rcg.name as references_group_name ",			
								 "FROM runs as r JOIN embeddings as e USING(embedding_id)  ",
								 "		JOIN compound_groups as ccg USING(compound_group_id) ",
								 "		JOIN compound_groups as rcg ON(rcg.compound_group_id = e.references_group_id) ",
								 "		JOIN compound_groups as scg ON(scg.compound_group_id = r.sample_group_id) ",
								 "WHERE run_id = ",runId))
}
getDescriptorType <- function(conn,runId=NULL,embeddingId=NULL,info = if(!is.null(runId)) getExtendedRunInfo(conn,runId) else NULL ){
	if(debug) print(info)
	result = if(!is.null(info))
		 runQuery(conn,paste("SELECT descriptor_type FROM descriptor_types WHERE descriptor_type_id = ",info$descriptor_type_id))[[1]]
	else if(!is.null(embeddingId))
		 runQuery(conn,paste("SELECT descriptor_type FROM embeddings JOIN descriptor_types USING(descriptor_type_id) WHERE embedding_id = ",embeddingId))[[1]]
	else
		stop("either runId, embeddingId, or info  must be specified")

	if(length(result) ==0)
		stop("could not find the descriptor type for descriptor_type_id ",info$descriptor_type_id)
	result
}
	

writeIddb <- function(conn,ids,name,append=FALSE) {

	dbTransaction(conn,{
		if(debug) message("writeiddb name: ",name)
		groupId = getCompoundGroupId(conn,name,create=TRUE)
		if(!append) # delete existing group
			runQuery(conn,paste("DELETE FROM compound_group_members WHERE compound_group_id = ",groupId))

		if(debug) message("groupid: ",groupId)
		if(debug) message("inserting members")
		#insert ids
		if(length(ids) !=0)
			insertGroupMembers(conn,data.frame(compound_group_id=groupId,compound_id=ids))
		groupId
   })

}
readIddb <- function(conn,name=NULL,groupId=getCompoundGroupId(conn,name),sorted=FALSE) {
	handle = if(is.null(name)) groupId else name
	#if(debug) message("readiddb name: ",handle)
	#groupId = getCompoundGroupId(conn,name)
	if(is.na(groupId))
		stop("compound group ",handle," was not found in the database")
	runQuery(conn,paste("SELECT compound_id FROM compound_group_members WHERE compound_group_id=
								 ",groupId, (if(sorted) " ORDER BY compound_id " else "")     ))[[1]]
}
getGroupSize <- function(conn,name=NULL,groupId=NULL) {
	handle = groupId

	if(is.null(groupId) && !is.null(name)){
		groupId = getCompoundGroupId(conn,name)
		handle=name
		if(length(groupId) == 0)
			stop("could not find compound group ",handle," in 'getGroupSize'")
	}else if(is.null(groupId) && is.null(name))
		stop("either 'groupId' or 'name' must be specified to 'getGroupSize'")

	size = runQuery(conn,paste("SELECT count(*) AS count FROM compound_group_members
										  WHERE compound_group_id = ",groupId,sep=""))$count
	if(length(size) == 0)
		stop("could not find size of compound group ",handle)
	if(debug) message("size of ",handle," is: ",size)
	size
}
getCompoundGroupId<- function(conn,name,create=FALSE) {
	#message("name: ",name)
	getOrCreate(conn, 
					paste("SELECT compound_group_id FROM compound_groups WHERE name = '",name,"'",sep=""), 
					paste("INSERT INTO compound_groups(name) VALUES('",name,"')",sep=""),
					errorTag=paste("compound group",name),create=create)

}

getEmbeddedDescriptors <- function(conn,embeddingId, compoundIds,descriptorIds=NULL){

	if(is.null(descriptorIds)){
		 queryIds = compoundIds
		 queryFn = function(ids) 
			  paste("SELECT cd.compound_id AS ids,value FROM descriptors as d
							JOIN embedded_descriptors ed USING(descriptor_id)
							JOIN embeddings as em USING(embedding_id)
							JOIN compound_descriptors as cd USING(descriptor_id)
						WHERE d.descriptor_type_id = em.descriptor_type_id AND  em.embedding_id = ",embeddingId,
						"  AND cd.compound_id IN (",paste(ids,collapse=","),")
						ORDER BY cd.compound_id,ordering",sep="")
	}
	else {
		 queryIds = descriptorIds
		 queryFn = function(ids) 
				paste("SELECT descriptor_id AS ids,value FROM embedded_descriptors WHERE embedding_id = ",embeddingId,
					"  AND descriptor_id IN (",paste(ids,collapse=","),")
					ORDER BY descriptor_id,ordering",sep="")
	}


	data = selectInBatches(conn,queryIds,queryFn)


	embeddedDescriptors = aggregate(data$value,list(ids=data$ids),identity)$x
	if(!is.matrix(embeddedDescriptors))
		stop("Could not create a matrix from emebdded descriptors. Perhaps they are not all the correct (same) length? ")


	if(nrow(embeddedDescriptors) != length(queryIds)){
		if(debug) print(compoundIds)
		stop(paste("missing some embedded descriptors! Found only",
					  nrow(embeddedDescriptors),"out of",length(queryIds),"given ids"))
	}

	embeddedDescriptors
}
getUnEmbeddedDescriptorIds <- function(conn,runId){
	descriptorIds=runQuery(conn,paste("SELECT descriptor_id FROM unembedded_descriptors WHERE run_id = ",runId,sep=""))
	message("found ",length(descriptorIds$descriptor_id)," un-embedded descriptors for runId ",runId)
	descriptorIds$descriptor_id
}

insertEmbeddedDescriptors <-function(conn,embeddingId,descriptorIds,data){

	numDescriptors = nrow(data)
	if(numDescriptors == 0)
		return()

	descriptorLength = ncol(data)
	assert(numDescriptors == length(descriptorIds))
	data=as.vector(data)

	toInsert = data.frame(embedding_id=embeddingId,descriptor_id=descriptorIds,
				  #ordering = rep(1:descriptorLength,numDescriptors),
				  ordering = as.vector(sapply(1:descriptorLength,function(i) rep(i,numDescriptors))),
				  value = data)

	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, 
			 paste("INSERT INTO embedded_descriptors(embedding_id,descriptor_id,ordering,value) ",
				"VALUES (:embedding_id,:descriptor_id,:ordering,:value)"),bind.data=toInsert)
	}else if(inherits(conn,"PostgreSQLConnection")){

		fields = c("embedding_id","descriptor_id","ordering","value")

		postgresqlWriteTable(conn,"embedded_descriptors",toInsert[,fields],append=TRUE,row.names=FALSE)

	}else{
		stop("database ",class(conn)," unsupported")
	}


}


#deprecated
insertEmbeddedDescriptorsByCompoundId <-function(conn,embeddingId,compoundIds,data){

	descriptorType = getDescriptorType(conn,embeddingId=embeddingId)
	descriptorIds = getDescriptorIds(conn,compoundIds,descriptorType,keepOrder=TRUE)
	numDescriptors = nrow(data)
	descriptorLength = ncol(data)
	assert(numDescriptors == length(descriptorIds))
	data=as.vector(data)
	toInsert = data.frame(embedding_id=embeddingId,descriptor_id=descriptorIds,
				  #ordering = rep(1:descriptorLength,numDescriptors),
				  ordering = as.vector(sapply(1:descriptorLength,function(i) rep(i,numDescriptors))),
				  value = data)

	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, 
			 paste("INSERT INTO embedded_descriptors(embedding_id,descriptor_id,ordering,value) ",
				"VALUES (:embedding_id,:descriptor_id,:ordering,:value)"),bind.data=toInsert)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("embedding_id","descriptor_id","ordering","value")
		postgresqlWriteTable(conn,"embedded_descriptors",toInsert[,fields],append=TRUE,row.names=FALSE)

	}else{
		stop("database ",class(conn)," unsupported")
	}


}

selectDescriptors <- function(type,ids){
	# paste(formatC(c(1,4,10000000,123400000056),format="fg"),collapse=",")
	q=paste("SELECT compound_id, descriptor FROM compound_descriptors 
					JOIN descriptors USING(descriptor_id)
					JOIN descriptor_types USING(descriptor_type_id) WHERE ",
				" descriptor_type='",type,"' AND compound_id IN (", paste(ids,collapse=","),") ORDER
				BY compound_id",sep="")
#	if(debug) message("select descriptors: ",q)
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
getDescriptorsByDescriptorId <- function(conn,descriptorIds){

	data = selectInBatches(conn,descriptorIds, function(ids) 
								  paste("SELECT descriptor_id,descriptor FROM descriptors WHERE descriptor_id IN ("
										  ,paste(ids,collapse=","),")",sep=""))
	if(nrow(data) != length(descriptorIds)){
		if(debug) print(descriptorIds)
		if(debug && nrow(data) < 100) {print("data: "); print(data$descriptor_id)}
		stop(paste("missing some descriptors. Found only ",nrow(data)," out of ",length(descriptorIds), "given descriptor ids"))
	}

	temp=data$descriptor
	names(temp) = data$descriptor_id
	temp[as.character(descriptorIds)]
}

getDescriptorIds <- function(conn,compoundIds,descriptorType=NULL,descriptorTypeId=NULL,keepOrder=FALSE){
	if(length(compoundIds) == 0)
		return(c())

	if(is.null(descriptorType) && is.null(descriptorTypeId))
		stop("either descriptorType or descriptorTypeId must be specified")

	selectClause = if(keepOrder) "descriptor_id, compound_id" else "DISTINCT descriptor_id"
	queryBase = if(!is.null(descriptorType))
					paste("SELECT ",selectClause," FROM compound_descriptors as cd
													  JOIN descriptors USING(descriptor_id)
													  JOIN descriptor_types USING(descriptor_type_id) 
													  WHERE descriptor_type = '",descriptorType,"'",sep="")
			  else #use descriptorTypeId
					paste("SELECT ",selectClause," FROM compound_descriptors as cd
													  JOIN descriptors USING(descriptor_id)
													  WHERE descriptor_type_id = '",descriptorTypeId,"'",sep="")
													  
	data = selectInBatches(conn,compoundIds, function(ids) 
								  paste(queryBase," AND cd.compound_id IN (",paste(ids,collapse=","),")",sep=""))
	descriptorIds =data$descriptor_id
	if(keepOrder){
		names(descriptorIds) = data$compound_id
		as.vector(descriptorIds[as.character(compoundIds)])
	}else
		descriptorIds
}
getRunDescriptorIds <- function(conn,runId){

	data = runQuery(conn,paste("SELECT DISTINCT cd.descriptor_id
										 FROM  runs AS r
												 JOIN compound_groups AS cg USING(compound_group_id)
												 JOIN compound_group_members AS cgm USING(compound_group_id)
												 JOIN compound_descriptors as cd USING(compound_id)
												 JOIN embeddings AS e ON(e.embedding_id = r.embedding_id)
												 JOIN descriptor_types AS dt ON(dt.descriptor_type_id = e.descriptor_type_id)
										WHERE r.run_id = ",runId,
										"ORDER BY cd.descriptor_id "))

	data$descriptor_id
}
getGroupDescriptorIds <- function(conn,groupId,descriptorTypeId){
	runQuery(conn,paste(" SELECT DISTINCT d.descriptor_id
											 FROM  compound_group_members AS cgm 
													 JOIN compound_descriptors USING(compound_id)
													 JOIN descriptors AS d USING(descriptor_id)
											WHERE d.descriptor_type_id = ",descriptorTypeId,
												" AND cgm.compound_group_id = ",groupId))[[1]]
}

writeMatrixFile<- function(conn,runId,compoundIds=c(),dir=".",samples=FALSE,cl=NULL,connSource=NULL){

	message("Regenerating matrix file...")

	runInfo = getExtendedRunInfo(conn,runId)
	descriptorIds=c()
	if(length(compoundIds) == 0){
		matrixFile = file.path(dir,paste("run",runInfo$num_references,runInfo$dimension,sep="-"),
								  paste(if(samples) "matrix.query" else "matrix",".",runInfo$num_references,"-",runInfo$dimension,sep=""))

		viewName = if(samples) "run_sample_embedded_descriptors" else "run_embedded_descriptors"
		descriptorIds= getGroupDescriptorIds(conn,if(samples) runInfo$sample_group_id else runInfo$compound_group_id,
												 runInfo$descriptor_type_id)
		numRows = length(descriptorIds)
	}
	else{
		matrixFile = file.path(if(debug) "." else tempdir(),"matrix")
		descriptorIds = getDescriptorIds(conn,compoundIds,descriptorTypeId=runInfo$descriptor_type_id)
		numRows = length(descriptorIds)
	}
	matrixFileTemp = paste(matrixFile,".temp",sep="")
	matrixFileIndexTemp = paste(matrixFile,".index.temp",sep="")
	if(debug) message("filename: ",matrixFile)

	f = file(matrixFileTemp,"wb")
	floatSize = 4

	numCols = runInfo$dimension
	if(debug)  message("numRows: ",numRows," numCols: ",numCols)

	writeBin(as.integer(floatSize),f,floatSize)
	writeBin(as.integer(numRows),f,floatSize)
	writeBin(as.integer(numCols),f,floatSize)
	

	indexF = file(matrixFileIndexTemp,"w")
	count=0

	writeChunk = function(df){
			if(debug) message("writing chunk. count= ",count)
			# for testing
			#df  = rbind(df,data.frame(descriptor_id=887,value=8.7))
			#df  = rbind(df,data.frame(descriptor_id=887,value=8.7))
			#df  = rbind(df,data.frame(descriptor_id=888,value=8.8))

			if(nrow(df) %% numCols != 0) {#this chunk is not a mulitple of descriptor dimention
				#find out which descriptor is lacking and then fail
				numIncomplete = 0
				for( id in unique(df$descriptor_id)){
					numValues = sum(df$descriptor_id == id)
					if( numValues != numCols){
						numIncomplete = numIncomplete + 1
						message("descriptor ",id," has only ",numValues," values out of ",numCols)
					}
				}
				stop(numIncomplete," descriptors missing values")
			}
			
			for( i in 1:nrow(df)){
				if(count %% numCols == 0)
					cat(paste(df$descriptor_id[i]),file=indexF,sep="\n")
				count <<- count + 1
			}
			writeBin(as.vector(df$value),f,floatSize)
   }


	batchByIndex(descriptorIds,function(ids){
		c = if(is.null(connSource)) conn else connSource()
		writeChunk(dbGetQuery(c,paste("SELECT descriptor_id,value
												  FROM embedded_descriptors 
												  WHERE descriptor_id IN (",paste(ids,collapse=","),")
														  AND embedding_id = ",runInfo$embedding_id,
												  "ORDER BY descriptor_id, ordering")))
		if(!is.null(connSource)) dbDisconnect(c)
	},10)
	if(count/numCols != numRows)
		stop("expected to find ",numRows," but wrote ",count/numCols)
	
	close(f)
	close(indexF)
	file.rename(matrixFileTemp,matrixFile)
	file.rename(matrixFileIndexTemp,paste(matrixFile,".index",sep=""))
	matrixFile
}
readMatrixIndex <- function(matrixFile){
	readIddbFile( paste(matrixFile,".index",sep=""))
}

descriptorsToCompounds <- function(conn,descriptorIds, all=FALSE){


	if(all){
		stop("Returning all compounds for a set of desciptors is not yet supported")
	}else{
		df = selectInBatches(conn,descriptorIds, function(ids)
					paste("SELECT cd.descriptor_id, min(cd.compound_id) AS compound_id 
							FROM compound_descriptors AS cd 
								  JOIN (SELECT descriptor_id, min(priority) AS priority
										  FROM compound_descriptors ",
										  "WHERE descriptor_id IN (",paste(ids,collapse=","),")",
										  "GROUP BY descriptor_id) AS t 
										ON(cd.descriptor_id=t.descriptor_id AND cd.priority=t.priority)",
							#"WHERE cd.descriptor_id IN (",paste(ids,collapse=","),")",
							"GROUP BY cd.descriptor_id",sep=""))
		compIds = df$compound_id
		names(compIds) = df$descriptor_id
		compIds[as.character(descriptorIds)]
	}



#	#message("descriptorIds: ")
#	#print(descriptorIds)
#	df = selectInBatches(conn,descriptorIds, function(ids)
#									paste("SELECT compound_id, descriptor_id FROM compound_descriptors",
#										  " WHERE descriptor_id IN (",paste(ids,collapse=","),")",sep=""))
#	if(all) {
#		compIds = list()
#		for(i in seq(along=df$compound_id)){
#			key = as.character(df$descriptor_id[i])
#			value = df$compound_id[i]
#			compIds[[key]] = if(is.na(compIds[key])) value else  c(compIds[key],value)
#		}
#		compIds
#	}else{
#		compIds = df$compound_id
#		names(compIds) = df$descriptor_id
#		compIds[as.character(descriptorIds)]
#	}
}



getOrCreate <- function(conn,getQuery,createQuery,create=FALSE,errorTag=getQuery){

	#print(getQuery)
	#print(createQuery)

	id = runQuery(conn,getQuery)[[1]]
	if(length(id)==0 || is.na(id)){
		if(!create)
			stop("could not find an entry for ",errorTag)
		runQuery(conn,createQuery)
		id = getOrCreate(conn,getQuery,createQuery,create=FALSE)
		if(length(id)==0 || is.na(id))
			stop("could not find or create an entry for ",errorTag)
	}
	if(length(id) > 1)
		stop("found more than one matches for ",errorTag)
	id
}

insertGroupMembers <- function(conn,data){
	#print("member data: ")
	#print(data)

	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
				"VALUES (:compound_group_id,:compound_id)"),bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("compound_group_id","compound_id")
		postgresqlWriteTable(conn,"compound_group_members",data[,fields],append=TRUE,row.names=FALSE)

		#apply(data[,fields],1,function(row) 
			#runQuery(conn,paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
					#"VALUES( $1,$2)"),row))
	}else{
		stop("database ",class(conn)," unsupported")
	}

}

runQuery <- function(conn,query,...){
	#if(debug) message(query)
	if(is.character(conn)){
		if(debug) print(class(conn))
		print(sys.calls())
	}
	df = dbGetQuery(conn,query,...)
	if(is.null(df))
		return(NULL)

	#the postgres driver insists on returning a data frame
	#with no columns when an empty result is returned.
	#this breaks everything that might use an index
	# so we make up some columns in that case...
	#print(class(df))

	if(ncol(df)==0){	
	   as.data.frame(rep(list(dummy=numeric(0)), 20))
	}else{
		df
	}
}
