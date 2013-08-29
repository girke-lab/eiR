

ensureSchema <- function(conn) {

	tableList=dbListTables(conn)

	if( ! all(c("compound_groups","compound_group_members","embedding","runs","embedded_descriptors") 
				 %in% tableList)) {
		print("ammending db")

		sqlFile = file.path("schema",if(inherits(conn,"SQLiteConnection")) "data.SQLite" 
								  else if(inherits(conn,"PostgreSQLConnection")) "data.RPostgreSQL")
																	
		nocomments = function(line) !grepl("^\\s*--",line)
		noblank = function(line) !grepl("^\\s*$",line)
		statements = Filter(noblank,unlist(strsplit(paste(
							  Filter(nocomments,readLines(system.file(sqlFile,package="eiR",mustWork=TRUE))),
							  collapse=""),";",fixed=TRUE)))
#		print(statements)

		Map(function(sql) dbGetQuery(conn,sql),statements)
	}

}

getEmbeddingId <- function(conn,name,r,d,descriptorType,refGroupId,create=TRUE){

	embeddingId = getOrCreate(conn,
									  paste("SELECT embedding_id FROM embeddings WHERE name='", name,"'",sep=""),
									  paste("INSERT INTO
											  embeddings(name,dimension,num_references,descriptor_type_id,references_group_id)",
											  "VALUES('",name,"',",d,",",r,
													",(SELECT descriptor_type_id FROM descriptor_types WHERE
															  descriptor_type='",descriptorType,"'),",
													refGroupId,")",sep=""),
									  errorTag=paste("embedding",name) )
	embeddingId
}
getEmbedding <-function(conn,embeddingId){
	dbGetQuery(conn,paste("SELECT embedding_id,name,dimension,num_references,descriptor_type_id,references_group_id ",
								 "FROM embeddings WHERE embedding_id = ",embeddingId))
}
getRunId <- function(conn,name,embeddingId,mainGroupId,queryGroupId) {
	runId = getOrCreate(conn,
							  paste("SELECT run_id FROM runs WHERE embedding_id
									  =",embeddingId," AND compound_group_id =
									  ",mainGroupId,sep=""),
							  paste("INSERT INTO runs(name,embedding_id,compound_group_id,sample_group_id)", 
									  "VALUES('",name,"',",embeddingId,",",mainGroupId,",",queryGroupId,")",sep=""),
							  errorTag = paste("run "+name))
	runId
}
getRun <- function(conn,runId){
	dbGetQuery(conn,paste("SELECT run_id,name, embedding_id,compound_group_id,sample_group_id ",
								 "FROM runs WHERE run_id = ",runId))
}
getExtendedRunInfo <-function(conn,runId){
	dbGetQuery(conn,paste("SELECT r.run_id,r.name, r.embedding_id,r.compound_group_id,r.sample_group_id, ",
								 "		e.name as embedding_name,e.dimension,e.num_references,e.descriptor_type_id,e.references_group_id ",			
								 "FROM runs as r JOIN embeddings as e USING(embdding_id)  ",
								 "WHERE run_id = ",runId))
}
	

writeIddb <- function(conn,ids,name,append=FALSE) {

	print("conn:")
	print(conn)
	dbTransaction(conn,{
		message("writeiddb name: ",name)
		groupId = getCompoundGroupId(conn,name)
		if(!append) # delete existing group
			dbGetQuery(conn,paste("DELETE FROM compound_group_members WHERE compound_group_id = ",groupId))

		message("groupid: ",groupId)
		message("inserting members")
		#insert ids
		insertGroupMembers(conn,data.frame(compound_group_id=groupId,compound_id=ids))
		groupId
   })

}
readIddb <- function(conn,name) {
	message("readiddb name: ",name)
	groupId = getCompoundGroupId(conn,name)
	if(is.na(groupId))
		stop("compound group ",name," was not found in the database")
	dbGetQuery(conn,paste("SELECT compound_id FROM compound_group_members WHERE compound_group_id=
								 ",groupId))[[1]]
}
getGroupSize <- function(conn,groupId=NULL,name=NULL) {
	handle = groupId

	if(is.null(groupId) && !is.null(name)){
		groupId = getCompoundGroupId(conn,name,create=FALSE)
		handle=name
		if(length(groupId) == 0)
			stop("could not find compound group ",handle," in 'getGroupSize'")
	}else if(is.null(groupId) && is.null(name))
		stop("either 'groupId' or 'name' must be specified to 'getGroupSize'")

	size = dbGetQuery(conn,paste("SELECT count(*) FROM compound_group_members
										  WHERE compound_group_id = ",groupId,sep=""))
	if(length(size) == 0)
		stop("could not find size of compound group ",handle)
	message("size of ",handle," is: ",size)
	size
}
getCompoundGroupId<- function(conn,name,create=TRUE) { #TODO: change create to default to FALSE
	message("name: ",name)
	getOrCreate(conn, 
					paste("SELECT compound_group_id FROM compound_groups WHERE name = '",name,"'",sep=""), 
					paste("INSERT INTO compound_groups(name) VALUES('",name,"')",sep=""),
					errorTag=paste("compound group",name))

}

insertEmbeddedDescriptors <-function(conn,embeddingId,compoundIds,descriptorType,data){

	descriptorIds = getDescriptorIds(conn,compoundIds,descriptorType)
	numDescriptors = nrow(data)
	descriptorLength = ncol(data)
	assert(numDescriptors == length(descriptorIds))
	#save(embeddingId,descriptorIds,descriptorLength,numDescriptors,data,file="debug.RData")
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
		fields = c("compound_group_id","compound_id")
		apply(toInsert,1,function(row) 
			dbGetQuery(conn,
				paste("INSERT INTO embedded_descriptors(embedding_id,descriptor_id,ordering,value) ",
					"VALUES( $1,$2,$3,$4)"),row))
	}else{
		stop("database ",class(conn)," unsupported")
	}


}
getDescriptorIds <- function(conn,compoundIds,descriptorType){
	data = dbGetQuery(conn,
										paste("SELECT descriptor_id FROM descriptors
													  JOIN descriptor_types USING(descriptor_type_id) 
													  WHERE descriptor_type = '",descriptorType,"'
													  AND compound_id IN (",paste(compoundIds,collapse=","),")",sep=""))
	descriptorIds =data[[1]]
	descriptorIds
}

writeMatrixFile<- function(conn,runId,dir="."){
#writeMatrixFile <- function(matrixFile,data,append=FALSE){

	print("writing matrix file")

	runInfo = getExtendedRunInfo(conn,runId)
	matrixFile = file.path(dir,paste("run",runInfo$dimension,runInfo$num_references,sep="-"),
								  paste("matrix",runInfo$dimension,runInfo$num_references,sep="."))

	f = file(matrixFile,"wb")
	floatSize = 4
	numRows= getGroupSize(conn,groupId=runInfo$compound_group_id)
	numCols = runInfo$dimension

	writeBin(as.integer(floatSize),f,floatSize)
	writeBin(numRows,f,floatSize)
	writeBin(numCols,f,floatSize)




	rs=dbSendQuery(conn,paste("SELECT * FROM run_embedded_descriptors WHERE run_id=",runId))
	bufferResultSet(rs,function(df){
			writeBin(as.vector(df$value),f,floatSize)
   },batchSize = 10000,closeRS=TRUE)
	

	#for(i in seq(1,numRows,length.out=numRows)){
		#for(j in seq(1,numCols,length.out=numCols)){
			#writeBin(data[i,j],f,floatSize)
		#}
	#}
	close(f)
}




getOrCreate <- function(conn,getQuery,createQuery,create=TRUE,errorTag=getQuery){

	print(getQuery)
	print(createQuery)

	id = dbGetQuery(conn,getQuery)[[1]]
	if(length(id)==0 || is.na(id)){
		dbGetQuery(conn,createQuery)
		id = getOrCreate(conn,getQuery,createQuery,create=FALSE)
		if(length(id)==0 || is.na(id))
			stop("could not find or create an entry for ",errorTag)
	}
	if(length(id) > 1)
		stop("found more than one matches for ",errorTag)
	id
}

insertGroupMembers <- function(conn,data){
	#message("member data: ",data)

	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
				"VALUES (:compound_group_id,:compound_id)"),bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("compound_group_id","compound_id")
		apply(data[,fields],1,function(row) 
			dbGetQuery(conn,paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
					"VALUES( $1,$2)"),row))
	}else{
		stop("database ",class(conn)," unsupported")
	}

}
