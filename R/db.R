

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
getGroupSize <- function(conn,name) {
	groupId = getCompoundGroupId(conn,name,create=FALSE)
	if(length(groupId) == 0)
		stop("could not find compound group ",name," while trying to find size")
	size = dbGetQuery(conn,paste("SELECT count(*) FROM compound_group_members
										  WHERE compound_group_id = ",groupId,sep=""))
	if(length(size) == 0)
		stop("could not find size of compound group ",name)
	message("size of ",name," is: ",size)
	size
}
getCompoundGroupId<- function(conn,name,create=TRUE) {
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
	data=as.vector(data)
	toInsert = data.frame(embedding_id=embeddingId,descriptor_id=descriptorIds,
				  order = rep(1:descriptorLength,numDescriptors),
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
										paste("SELECT descriptor_ids FROM descriptors
													  JOIN descriptor_types USING(descriptor_type_id) 
													  WHERE descriptor_type = '",descriptorType,"'
													  AND compound_id IN (",paste(compoundIds,collapse=","),")",sep=""))
	descriptorIds =data[[1]]
	descriptorIds
}

writeMatrixFile<- function(conn,runId){

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
