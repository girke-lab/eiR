\name{addTransform}
\alias{addTransform}
\title{
	Add Transform
}
\description{
	New descriptor types can be added using the \code{addTransform}
	function. These transforms are basically just ways to read descriptors from 
	compound definitions, and to convert descriptors between string and object form. This conversion
	is required because descriptors are stored as strings in the SQL database, but are used
	by the rest of the program as objects.

	There are two main components that need to be added. The
	\code{addTransform} function takes the name of the transform and two functions,
	\code{toString}, and \code{toObject}. These have slightly different
	meanings depending on the component you are adding.
	The first component to add 
	is a transform from a chemical compound format, such as SDF, to a descriptor format,
	such as atom pair (AP), in either string or object form. 
	The toString function should take
	any kind of chemical compound source, such an SDF file, an SDF object or an SDFset, 
	and output a string representation of the descriptors. Since this function
	can be written in terms of other functions that will be defined, you
	can usually accept the default value of this function.  The toObject function
	should take the same kind of input, but output the descriptors as an object.
	The actual return value is a list containing the names of the compounds (in the names field), and
	the actual descriptor objects ( in the descriptors field).

	The second component to add is a transform that converts between string
	and object representations of descriptors. In this case the toString function
	takes descriptors in object form and returns a string representation for each.
	The toObject function performs the inverse operation. It takes descriptors in string
	form and returns them as objects. The objects returned by this function
	will be exactly what is handed to the distance function, so you need to make sure that
	the two match each other.


}
\usage{
addTransform(descriptorType, compoundFormat = NULL, toString = NULL, toObject)
}
\arguments{
  \item{descriptorType}{
		The name of the type of the descritor being added.
	}
  \item{compoundFormat}{
		The format of the compound data the descriptor will be extracted from.
}
  \item{toString}{
		A function with three arguments, the data, an SQL connection object, and a directory name.
		The last two are optional and can be set to a default value of NULL if not used
		in the body of the function. If this parameter is NULL and compoundFormat is not NULL, then a default
		function will be used for this value.
	}
  \item{toObject}{
		A function with three arguments, the data, an SQL connection object, and a directory name.
		The last two are optional and can be set to a default value of NULL if not used
		in the body of the function. 
		If compoundFormat is not NULL, then the return value of this function should
		be a list with the fields "names" and "descriptors", containing the compound names and 
		descriptor objects, respectivly. If compoundFormat is NULL, then the return value should be
		a collection of descriptor objects, in whatever format the distance function for this descrptor type	
		requires.
	}
}
\value{
	No value returned.
}
\author{
	Kevin Horan
}

\seealso{
 \code{\link{setDefaultDistance}}
}
\examples{
	# adding support for atompair (ap) descriptors extracted from
	# sdf formmatted data.

    #first component
    addTransform("ap-example","sdf-example",
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

    #second component
    addTransform("ap-example",  
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
}
