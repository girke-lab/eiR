

ensureSchema <- function(conn) {

	tableList=dbListTables(conn)

	if( ! all(c("compound_groups","compound_group_members","embedding","runs","embedded_descriptors") 
				 %in% tableList)) {
		print("ammending db")

		sqlFile = file.path("schema",if(inherits(conn,"SQLiteConnection")) "data.SQLite" 
								  else if(inherits(conn,"PostgreSQLConnection")) "data.RPostgreSQL")
																	
		statements = unlist(strsplit(paste(
							  readLines(system.file(sqlFile,package="eiR",mustWork=TRUE)),
							  collapse=""),";",fixed=TRUE))
		#print(statements)

		Map(function(sql) dbGetQuery(conn,sql),statements)
	}

}

writeIddb <- function(conn,ids,name,appned=FALSE) {

	dbTransaction(conn,{
		groupId = compoundGroupId(name)
		#search for existing group
		if(is.na(groupId)){ #create new group if not found
			dbGetQuery(conn,paste("INSERT INTO compound_groups(name) VALUES('",name,"'",sep=""))
			groupId = compoundGroupId(name)
			if(is.na(groupId))
				stop("could not find or create an entry for compound group ",name)
		}
		if(!append){ # delete existing group
			dbGetQuery(conn,paste("DELETE FROM compound_groups WHERE compound_group_id = ",groupId))
		}

		#insert ids
		insertGroupMembers(conn,data.frame(compound_group_id=groupId,compound_id=ids))
   })

}
readIddb <- function(conn,name) {
	groupId = compoundGroupId(name)
	if(is.na(groupId))
		stop("compound group ",name," was not found in the database")
	dbGetQuery(conn,paste("SELECT compound_id FROM compound_groups WHERE compound_group_id=
								 ",groupId))[[1]]
}
compoundGroupId<- function(conn,name) {
	dbGetQuery(conn,paste("SELECT compound_group_id FROM compound_groups WHERE name =
								 '",name,"'",sep=""))[[1]][1]
}

insertGroupMembers <- function(conn,data){
	if(inherits(conn,"SQLiteConnection")){
		dbGetPreparedQuery(conn, paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
				"VALUES (:compound_group_id,:compound_id)"),bind.data=data)
	}else if(inherits(conn,"PostgreSQLConnection")){
		fields = c("compound_group_id","compound_id")
		apply(data[,fields],1,function(row) 
			dbGetQuery(conn,paste("INSERT INTO compound_group_members(compound_group_id,compound_id) ",
					"VALUES( $1,$2)"),row)
	}else{
		stop("database ",class(conn)," unsupported")
	}

}
