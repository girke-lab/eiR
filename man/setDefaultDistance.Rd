\name{setDefaultDistance}
\alias{setDefaultDistance}
\title{
	Set the default distance function for a descriptor type
}
\description{
	Set the default distance function for a descriptor type. This is the distance function
that will be used if none is given for a particular function call.
}
\usage{
setDefaultDistance(descriptorType, distance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{descriptorType}{
		The type of the descriptor to set a distance function for. Built-in values
		are "ap" and "fp". Additional values can be set as well.
  }
  \item{distance}{
		A distance function taking two descriptor objects (as returned by toObject in a descriptor transform, 
		see \code{\ ink{addTransform}} for details), and returning a distance value.
	}
}
\value{
	No return value.
}
\author{
	Kevin Horan
}


\seealso{
  \code{\link{addTransform}}
}
\examples{
	setDefaultDistance("ap", function(d1,d2) 1-cmp.similarity(d1,d2) )
}
