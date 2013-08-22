\name{eiAdd}
\alias{eiAdd}
\title{
   Add new compounds
}
\description{
   Add additional compounds to and existing database
}
\usage{
	eiAdd(r,d,refIddb,additions,dir=".",format="sdf",conn=defaultConn(dir),
		descriptorType="ap",distance=getDefaultDist(descriptorType),
		updateByName = FALSE)

}
\arguments{
  \item{r}{
      The number of references used to build the database you wish to
      query agains.
   }
  \item{d}{
      The number of dimensions used to build the database you wish to
      query agains.
   }
  \item{refIddb}{
      An Iddb formatted file containing the reference IDs of the database
      you wish to append to. This should almost always be the file:
      run-\code{r}-\code{d}/<long random string>.cdb.
		The \code{refIddb} value should also be returned by eiMakeDb.
   }
  \item{additions}{
      The compounds to add. This can be either a file in sdf format, or 
		and SDFset object.
  }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
  \item{format}{
	  The format of the data given in  \code{additions}. Currently only "sdf" is
	  supported.
   }
	\item{conn}{
		Database connection to use.
	}
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}
	\item{updateByName}{
		If true we make the assumption that all compounds, both in the existing database and the
		given dataset, have unique names. This function will then avoid re-adding existing,
		identical compounds, and will update existing compounds with a new definition if a new
		compound definition with an existing name is given. 

		If false, we allow duplicate compound names to exist in the database, though not
		duplicate definitions. So identical compounds will not be re-added, but if a new version of
		an existing compound is added it will not update the existing one, it will add the modified one 
		as a completely new compound with a new compound id. 
	}

}
\details{
   New Compounds can be added to an existing database, however, the 
   reference compounds cannot be changed. 
	This will also update the matrix file in the run/job directory with the 
	new compounds.
}
\examples{

   library(snow)
   r<- 50
   d<- 40

   #initialize
   data(sdfsample)
   dir=file.path(tempdir(),"add")
   dir.create(dir)
   eiInit(sdfsample[1:99],dir=dir)

   #create compound db
   refIddb=eiMakeDb(r,d,numSamples=20,dir=dir)
      

   #find compounds similar two each query
   eiAdd(r,d,refIddb,sdfsample[100],dir=dir)


}
\author{ Kevin Horan }
\seealso{
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}

