\name{eiQuery}
\alias{eiQuery}
\title{
   Perform a query on an embedded database
}
\description{
   Finds similar compounds for each query.
}
\usage{
	eiQuery(runId,queries,format="sdf",
		dir=".",distance=getDefaultDist(descriptorType),conn=defaultConn(dir),
		asSimilarity=FALSE, K=200, searchK=-1,lshData=NULL,
		mainIds = readIddb(conn,file.path(dir, Main)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{runId}{
		The id number identifying a particular set of settings for a database. This is generally
		the number returned by \code{\link{eiMakeDb}}. If your coming from an older version of eiR, you should
		not use this value instead of
		specifying \code{r}, \code{d},\code{refIddb}, and \code{descriptorType}.
	}

  \item{queries}{
      This can be either an SDFset, or a file containg 1 or more
      query compounds.
   }
	\item{format}{
		The format in which the queries are given. Valid values are: "sdf" when
		\code{queries} is either a filename of an sdf file, or and SDFset object;
		"compound_id" when \code{queries} is a list of id numbers; and "name", when \code{queries}
		is a list of compound names, as returned by \code{cid(apset)}.
	}
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors. The Tanimoto function is used by default.
	}

	\item{conn}{
		Database connection to use.
	}
	\item{asSimilarity}{
		If true, return similarity values instead of distance values. This only works in the given
		distance function returns values between 0 and 1. This is true for the default atom pair and 
		finger print distance functions.
	}
  \item{K}{
     The number of results to return.
   }
  \item{searchK}{
     Tunable Annoy LSH parameter. A larger value will give more accurate results, but will take longer time to return.
	  The default value of -1 will allow the value to chosen automatically, which will set a value of numTrees * (approximate number of nearest neighbors).
	  See Annoy page for details.  \url{https://github.com/spotify/annoy}
   }
	\item{lshData}{
		DEPRECATED. This is no longer used.
	}
	\item{mainIds}{
		A vector of all id numbers in the current database. This is mainly provided as an option
		here to avoid having to re-read the id list multiple times when executing several queries.
		If not supplied it will read it in itself.
	}
}
\details{
	This function identifies the database by the \code{r}, \code{d}, and 
	\code{refIddb} parameters. The queries can be given in a few
	different formats, see the \code{queries} parameter for details. 
	The LSH algorithm is used to quickly identify compounds similar to the
	queries. 
	This function must use a distance function rather than a similarity function.
	However, if the distance function given returns values between 0 and 1, then
	the \code{asSimilarity} parameter may be used to return similarity values rather
	than distance values.
}
\value{
   Returns a data frame with columns 'query', 'target', 'target_ids', and
   'distance'. 'query' and 'target' are the compound names and
   distance is the distance between them, as computed by
   the given distance function.'target_ids' is the compound id of the target.
	Query namess are repeated for each matching target found.
	If \code{asSimilarity} is true then instead of a "distance"
	column there will be a "similarity" column.
}
\author{
   Kevin Horan
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
   \code{\link{eiInit}}
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
}
\examples{

   library(snow)
   r<- 50
   d<- 40

   #initialize
   data(sdfsample)
   dir=file.path(tempdir(),"query")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   runId=eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(1,type="SOCK",outfile=""))

   #find compounds similar two each query
   results = eiQuery(runId,sdfsample[1:2],K=15,dir=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
