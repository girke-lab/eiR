

eiOptions = new.env()
eiOptions$defaultDistances=list()
eiOptions$Translations=list()


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



addTransform <- function(descriptorType,compoundFormat=NULL,toString=NULL,toObject){

	name = buildType(descriptorType,compoundFormat)

	if(!is.null(compoundFormat) && is.null(toString))
		toString = function(input,dir=".") 
			getTransform(descriptorType)$toString(getTransform(descriptorType,compoundFormat)$toObject(input)$descriptors)
	else if(is.null(toString))
		stop("toString function must be specified if compoundFormat is NULL")		


	eiOptions$Translations[[name]] = list(toString=toString,toObject=toObject)

	if( is.null(compoundFormat)){
		# add extra handlers for compound_id and name types
		addTransform(name,"compound_id",
			# compound_id -> ap string
			toString = function(ids,dir="."){
				getDescriptors(initDb(file.path(dir,ChemDb)),name,ids)
			},
			# compound_id -> AP list object
			toObject = function(ids,dir="."){
				descInfo = getTransform(name,"compound_id")$toString(ids,dir)
				list(names=names(descInfo),
					  descriptors=getTransform(name)$toObject(descInfo))
			}
		)

		addTransform(name,"name",
			# name -> ap string
			toString = function(names,dir="."){
				conn=initDb(file.path(dir,ChemDb))
				getDescriptors(conn,name,findCompoundsByName(conn,names,keepOrder=TRUE))
			},
			# name -> AP list object
			toObject = function(names,dir="."){
				descInfo = getTransform(name,"name")$toString(names,dir)
				list(names=names,
					  descriptors=getTransform(name)$toObject(descInfo))
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
	toObject = function(input,dir="."){
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
	toObject = function(input,dir="."){
		apList = getTransform("ap","sdf")$toObject(input,dir)
		apList$descriptors = desc2fp(apList$descriptors)
		apList
	}
)
addTransform("ap",  
   # APset -> string,
	toString = function(apset,dir="."){
		unlist(lapply(ap(apset), function(x) paste(x,collapse=", ")))
	},
   # string or list -> AP set list
	toObject= function(v,dir="."){ 
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
	toString = function(fpset,dir="."){
		sapply(1:length(fpset), function(i) as(fpset[i],"character") )
	},
   # string or list -> FP set list
	toObject= function(v,dir="."){ 
		if(inherits(v,"list") || length(v)==0)
			return(v)

		as( if(!inherits(v,"FPset")){
				#names(v)=as.character(1:length(v));  
				read.AP(v,type="fp",isFile=FALSE)
			} else v,
			"FP")  
	}
)

