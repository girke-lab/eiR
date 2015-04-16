\name{loadLSHData}
\alias{loadLSHData}
\title{
	Load LSH Data
}
\description{
	Load the LSH index and data. If many queries are going to be performed it is advantageous 
	to load this object first and then hand it to \code{\link{eiQuery}} via the
	\code{lshData} parameter for each query. 
	
	If the data needs to be freed you can call the \code{\link{freeLSHData}} function.
}
\usage{
loadLSHData(r, d, W = NA, M = NA, L = NA, K = NA, T = NA, dir = ".", matrixFile = NULL)
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
  \item{W}{
		See \code{\link{eiQuery}}
	}
  \item{M}{
		See \code{\link{eiQuery}}
	}
  \item{L}{
		See \code{\link{eiQuery}}
	}
  \item{K}{
		See \code{\link{eiQuery}}
	}
  \item{T}{
		See \code{\link{eiQuery}}
	}
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
	}
  \item{matrixFile}{
		The path to the matrix file. If not specified it will
		look for it in the default spot.
	}
}
\value{
	Returns a pointer to the allocated data. This should only be passed to 
	other functions with an \code{lshData} parameter, such as \code{\link{eiQuery}}
}
\author{
Kevin Horan
}


\seealso{
	\code{\link{freeLSHData}}
	\code{\link{eiQuery}}
}
\examples{
	\dontrun{
		lshData = loadLSHData(r,d)
		eiQuery(r,d,refIddb,c("650002","650003"), format="name",K=15,lshData=lshData)
		eiQuery(r,d,refIddb,c("650004","650005"), format="name",K=15,lshData=lshData)
		freeLSHData(lshData)
	}
}
