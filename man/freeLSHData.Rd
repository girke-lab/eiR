\name{freeLSHData}
\alias{freeLSHData}
\title{
	Free LSH Data
}
\description{
	Free the memory allocated by \code{\link{loadLSHData}}.
}
\usage{
freeLSHData(lshData)
}
\arguments{
  \item{lshData}{
		A pointer returned by \code{\link{loadLSHData}}.
	}
}
\value{
	No return value.
}
\author{
	Kevin Horan
}
\seealso{
	\code{\link{loadLSHData}}
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
