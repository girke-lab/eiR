\name{eiQuery}
\alias{eiQuery}
\title{
   Perform a query on an embedded database
}
\description{
   Finds similar compounds for each query.
}
\usage{
	eiQuery(r,d,refIddb,queries,format="sdf",
		dir=".",descriptorType="ap",distance=getDefaultDist(descriptorType),conn=defaultConn(dir),
		K=200, W = 1.39564, M=19,L=10,T=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
      The number of references used to build the database you wish to
      query against.
   }
  \item{d}{
      The number of dimensions used to build the database you wish to
      query against.
   }
  \item{refIddb}{
      An Iddb formatted file containing the reference IDs of the database
      you wish to query against. This should almost always be the file:
      run-\code{r}-\code{d}/<long random string>.cdb.
		The \code{refIddb} value should also be returned by eiMakeDb.
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
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors. The Tanimoto function is used by default.
	}

	\item{conn}{
		Database connection to use.
	}
  \item{K}{
     The number of results to return.
   }
  \item{W}{
     Tunable LSH parameter. See LSHKIT page for details.
     \url{http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{M}{
     Tunable LSH parameter. See LSHKIT page for details.
     \url{http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html}
   }
  \item{L}{
     Number of hash tables
   }
  \item{T}{
     Number of probes
   }
}
\details{
	This function identifies the database by the \code{r}, \code{d}, and 
	\code{refIddb} parameters. The queries can be given in a few
	different formats, see the \code{queries} parameter for details. 
	The LSH algorithm is used to quickly identify compounds similar to the
	queries. 
}
\value{
   Returns a data frame with columns 'query', 'target', 'target_ids', and
   'distance'. 'query' and 'target' are the compound names and
   distance is the distance between them, as computed by
   the given distance function.'target_ids' is the compound id of the target.
	Query namess are repeated for each matching target found.
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
   refIddb=eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(1,type="SOCK",outfile=""))

   #find compounds similar two each query
   results = eiQuery(r,d,refIddb,sdfsample[1:2],K=15,dir=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
