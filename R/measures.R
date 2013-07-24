

eiOptions = new.env()
eiOptions$defaultDistances=list()
eiOptions$Translations=list()
eiOptions$defaultConn=NULL



setDefaultDistance <- function(descriptorType,distance)
	eiOptions$defaultDistances[[descriptorType]] =  distance

getDefaultDist <- function(descriptorType){
	d=eiOptions$defaultDistances[[descriptorType]]
	if(is.null(d))
		message("no default distance found for desciptor type ",descriptorType)
	d
}

setDefaultDistance("ap", function(d1,d2) 1-cmp.similarity(d1,d2) )
setDefaultDistance("fp", function(d1,d2) 1-fpSim(d1,d2) )

defaultConn <- function(dir=".",create=FALSE){

	conn=eiOptions$defaultConn
	if(is.null(conn) && (create || file.exists(file.path(dir,ChemDb)))){
		conn = initDb(file.path(dir,ChemDb))
		# if we set this, the examples fail because the DBI packages are unloaded
		# but not re-loaded because we won't re-call initDb if we already have a connection
		#eiOptions$defaultConn=conn
	}
	if(is.null(conn))
		stop("no default connection found, looked for SQLite db in ",file.path(dir,ChemDb))
	
	conn
}
setDefaultConn <- function(conn){
	eiOptions$defaultConn = conn
}


addTransform <- function(descriptorType,compoundFormat=NULL,toString=NULL,toObject){

	name = buildType(descriptorType,compoundFormat)

	if(!is.null(compoundFormat) && is.null(toString))
		toString = function(input,dir=".",conn=defaultConn()) 
			getTransform(descriptorType)$toString(getTransform(
																	descriptorType,compoundFormat)$toObject(input,conn,dir)$descriptors,
															  dir=dir)
	else if(is.null(toString))
		stop("toString function must be specified if compoundFormat is NULL")		


	eiOptions$Translations[[name]] = list(toString=toString,toObject=toObject)

	if( is.null(compoundFormat)){
		# add extra handlers for compound_id and name types
		addTransform(name,"compound_id",
			# compound_id -> ap string
			toString = function(ids,conn,dir="."){
				getDescriptors(conn,name,ids)
			},
			# compound_id -> AP list object
			toObject = function(ids,conn,dir="."){
				descInfo = getTransform(name,"compound_id")$toString(ids,conn,dir)
				list(names=names(descInfo),
					  descriptors=getTransform(name)$toObject(descInfo,conn,dir))
			}
		)

		addTransform(name,"name",
			# name -> ap string
			toString = function(names,conn,dir="."){
				getDescriptors(conn,name,findCompoundsByName(conn,names,keepOrder=TRUE))
			},
			# name -> AP list object
			toObject = function(names,conn,dir="."){
				descInfo = getTransform(name,"name")$toString(names,conn,dir)
				list(names=names,
					  descriptors=getTransform(name)$toObject(descInfo,conn,dir))
			}
		)
	}	



}
getTransform <- function(descriptorType,compoundFormat=NULL){

	name = buildType(descriptorType,compoundFormat)

	t=eiOptions$Translations[[name]]
	if(is.null(t))
		stop("transform ",name," not defined")
	t
}

buildType <- function(descriptorType,compoundFormat) {
	tolower(if(is.null(compoundFormat)) descriptorType else paste(compoundFormat,descriptorType,sep="-"))
}



addTransform("ap","sdf",
	# Any sdf source -> APset
	toObject = function(input,conn=NULL,dir="."){
		sdfset=if(is.character(input) && file.exists(input)){
			read.SDFset(input)
		}else if(inherits(input,"SDFset")){
			input
		}else{
			stop(paste("unknown type for 'input', or filename does not exist. type found:",class(input)))
		}
		list(names=sdfid(sdfset),descriptors=sdf2ap(sdfset))
	}
)

addTransform("fp","sdf",
	# Any sdf source -> FPset
	toObject = function(input,conn=NULL,dir="."){
		apList = getTransform("ap","sdf")$toObject(input,dir=dir)
		apList$descriptors = desc2fp(apList$descriptors)
		apList
	}
)
addTransform("ap",  
   # APset -> string,
	toString = function(apset,conn=NULL,dir="."){
		unlist(lapply(ap(apset), function(x) paste(x,collapse=", ")))
	},
   # string or list -> AP set list
	toObject= function(v,conn=NULL,dir="."){ 
		if(inherits(v,"list") || length(v)==0)
			return(v)

		as( if(!inherits(v,"APset")){
				names(v)=as.character(1:length(v));  
				read.AP(v,type="ap",isFile=FALSE)
			} else v,
			"list")  
	}
)
addTransform("fp",  
   # FPset -> string,
	toString = function(fpset,conn=NULL,dir="."){
		sapply(1:length(fpset), function(i) as(fpset[i],"character") )
	},
   # string or list -> FP set list
	toObject= function(v,conn=NULL,dir="."){ 
		if(inherits(v,"list") || length(v)==0)
			return(v)

		as( if(!inherits(v,"FPset")){
				#names(v)=as.character(1:length(v));  
				read.AP(v,type="fp",isFile=FALSE)
			} else v,
			"FP")  
	}
)

