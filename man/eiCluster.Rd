\name{eiCluster}
\alias{eiCluster}
\title{
	Cluster compounds
}
\description{
	Uses Jarvis-Patrick clustering to cluster the compound database using 
	the LSH algorithm to quickly find nearest neighbors.
}
\usage{
	eiCluster(runId,K,minNbrs,compoundIds=c(), dir=".",cutoff=NULL,
							 distance=getDefaultDist(descriptorType),
							 conn=defaultConn(dir), searchK=-1,type="cluster",linkage="single")
}
\arguments{
	\item{runId}{
		The id number identifying a particular set of settings for a database. This is generally
		the number returned by \code{\link{eiMakeDb}}. If your coming from an older version of eiR, you should
		not use this value instead of
		specifying \code{r}, \code{d}, and \code{descriptorType}.
	}

  \item{K}{
	  The number of neighbors to consider for each compound.
	}
  \item{minNbrs}{
	  The minimum number of neighbors that two comopunds must have in common in order to be joined.
	}
	\item{compoundIds}{
		If this variable is set to a vector of compound ids, then clustering will be done
		with just those compounds. If left unset or empty, clustering will apply to all
		compounds in the given run.
	}
	\item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}

  \item{cutoff}{
	  Distance cutoff value. Compounds having a distance larger this this value will not
	  be included in the nearest neighbor table. Note that this is a distance value, not a similarity
	  value, as is often used in other ChemmineR functions.
	}

	\item{conn}{
		Database connection to use.
	}
  \item{searchK}{
     Tunable Annoy LSH parameter. A larger value will give more accurate results, but will take longer time to return.
	  The default value of -1 will allow the value to chosen automatically, which will set a value of numTrees * (approximate number of nearest neighbors).
	  See Annoy page for details.  \url{https://github.com/spotify/annoy}
   }

   \item{type}{
		If "cluster", returns a clustering, else, if "matrix", returns a list in the format 
		expected by the \code{jarvisPatrick} function in ChemmineR. This list contains the 
		nearest neighbor matrix along with the similarity matrix. This allows one to quickly try
		different cutoff values without having to re-compute the whole similarity matrix each time.
		Note that since we are returning similarity values here instead of distance values,
		this will only work if the given distance function returns a value between 0 and 1. This is
		true of the default funtions.
	}
	\item{linkage}{
		Can be one of "single", "average", or "complete", for single linkage, average linkage and complete linkage
		merge requirements, respectively. In the context of Jarvis-Patrick, average linkage means that at least half 
		of the pairs between the clusters under consideration must pass the merge requirement. Similarly, for complete
		linkage, all pairs must pass the merge requirement. Single linkage is the normal case for Jarvis-Patrick and 
		just means that at least one pair must meet the requirement.
	}

}
\details{
	The jarvis patrick clustering algorithm takes a set of items, a distance function, and two
	parameters, \code{K}, and \code{minNbrs}. For each item, it find the \code{K} nearest neighbors
	of that item. Normally this requires computing the distance between every pair of items.
	However, using Locality Sensative Hashing (LSH), the set of nearst neighbors can be found in
	near constant time. Once the nearest neighbor matrix is computed, the algorithm makes one pass
	through the items and merges all pairs that have at least \code{minNbrs} neighbors in common.

	Although not required, it is avisable to specify a \code{cutoff} value. This is the maximum
	distance two items can have from each other and still be considered to be neighbors. 	It is
	thus possible for an item to end up with less than \code{K} neighbors if less than \code{K}
	items are close enough to it. If a cutoff is not specified, it is possible for highly
	un-related items to be listed as neighbors of another item simply because nothing else was
	nearby. This can lead to items being joined into clusters with which they have no true
	connection.

	The \code{type} parameter can be used to return a list which can be used to call the
	\code{jarvisPatrick} function in ChemmineR directly. The advantage of this is that it will
	contain the similarity matrix which can then be used to quickly set different cutoff values
	(using \code{trimNeighbors}) whithout having to re-compute the similarity matrix. Note that this 
	requires that the given distance function return a value between 0 and 1 so it can be converted
	to a similarity function.

}
\value{
		If \code{type} is "cluster", returns a clustering.
		This will be a vector in which the names are the compound names, and the values are the cluster
		labels.
		Otherwise, if \code{type} is "matrix", returns a list with the following components:
      \item{indexes}{index values of nearest neighbors, for each item. }
      \item{names}{The database compound id of each item in the set.}
      \item{similarities}{The similarity values of each neighbor to the item for that row. 
          Each similarity values corresponds to the id number in the 
         same position in the indexes entry}

		If there are not \code{K} neibhbors for a compound, that row will be padded with NAs.
}
\author{
	Kevin Horan
}

\examples{

	library(snow)
   r<- 50
   d<- 40

   #initialize
   data(sdfsample)
   dir=file.path(tempdir(),"cluster")
   dir.create(dir)
   eiInit(sdfsample,dir=dir,skipPriorities=TRUE)

   #create compound db
   runId=eiMakeDb(r,d,numSamples=20,dir=dir, cl=makeCluster(1,type="SOCK",outfile=""))

	eiCluster(runId,K=5,minNbrs=2,cutoff=0.5,dir=dir)

}
