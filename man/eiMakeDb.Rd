\name{eiMakeDb}
\alias{eiMakeDb}
\title{
   Create an embedded database
}
\description{
   Uses the initalized compound data to create an embedded compound
   databbase with \code{r} reference compounds in \code{d} dimensions.
}
\usage{
	eiMakeDb(refs,d,descriptorType="ap",distance=getDefaultDist(descriptorType), 
				dir=".",numSamples=cdbSize(dir)*0.1,conn=defaultConn(dir),
				cl=makeCluster(1,type="SOCK"))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{refs}{
		The reference compounds to use to build the database you wish to query against.
		\code{Refs} can be one of three things. It can be a filename of an iddb file
		giving the index values of the reference compounds to use, it can be vector of 
		index values, or it can be a scalar value giving the number of randomly selected
		references to use.
   }
  \item{d}{
      The number of dimensions used to build the database you wish to
      query against.
   }
  \item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
   }
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}

  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{numSamples}{
      The number of non-reference samples to be chosen now to be used
      later by the eiPerformanceTest function.
   }
	\item{conn}{
		Database connection to use.
	}
  \item{cl}{
     A SNOW cluster can be given here to run this function in
     parrallel.
   }
}
\details{
   This function will embedd compounds from the data
   directory in another space which allows for more
   efficient searching. The main two parameters are r and
   d. r is the number of reference compounds to use and
   d is the dimension of the embedding space.  We have
   found in practice that setting d to around 100 works
   well.  r should be large enough to ``represent'' the
   full compound database. Note that an r by r matrix will be constructed
	during the course of execution, so r should be less than
	about 46,000 to avoid overflowing an integer.
   Since this is the longest running step, a SNOW cluster can be
   provided to parallelize the task.
    
   To help tune these values, \code{eiMakeDb} will pick
   \code{numSamples} non-reference samples which can later be used by the
   \code{eiPerformanceTest} function.

   \code{eiMakdDb} does its job in a job folder, named after the number of reference
   compounds and the number of embedding dimensions. For example, using 300
   reference compounds to generate a 100-dimensional embedding (r=300,
   d=100) will result in a job folder called run-300-100. 
   The embedding result is the file matrix.<r>.<d>. In the above example,
   the output would be run-300-100/matrix.300.100.


}
\value{
   Creates files in \code{dir} ("run-r-d" by default).
	The return value is the name of the refIddb file, which needs to be
	given to other functions such as eiQuery or eiAdd.
}
\author{
   Kevin Horan
}


\seealso{
   \code{\link{eiInit}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}
\examples{
   library(snow)

   r<- 50
   d<- 40

   #initialize 
   data(sdfsample)
   dir=file.path(tempdir(),"makedb")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   refIddb=eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(1,type="SOCK",outfile=""))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
