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
	eiCluster(r,d,K,minNbrs, dir=".",cutoff=NULL,
							 descriptorType="ap",distance=getDefaultDist(descriptorType),
							 conn=defaultConn(dir), W = 1.39564, M=19,L=10,T=30,type="cluster",linkage="single")
}
\arguments{
  \item{r}{
      The number of references used to build the database you wish to
      query against.
	}
  \item{d}{
      The number of dimensions used to build the database you wish to
      query against.
	}
  \item{K}{
	  The number of neighbors to consider for each compound.
	}
  \item{minNbrs}{
	  The minimum number of neighbors that two comopunds must have in common in order to be joined.
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
   \item{type}{
		If "cluster", returns a clustering, else, if "matrix", returns a nearest neighbor matrix.
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
	The jarvis patrick clustering algorithm takes a set of items, a distance function, and two parameters, \code{K}, and 
	\code{minNbrs}. For each item, it find the \code{K} nearest neighbors of that item. Normally this requires computing the distance
	between every pair of items. However, using Locality Sensative Hashing (LSH), the set of nearst neighbors can be
	found in near constant time. Once the nearest neighbor matrix is computed, the algorithm makes one pass through
	the items and merges all pairs that have at least \code{minNbrs} neighbors in common.

	Although not required, it is avisable to specify a \code{cutoff} value. This is the maximum distance two items can
	have from each other and still be considered to be neighbors. 	It is thus possible for an item to end up with less
	than \code{K} neighbors if less than \code{K} items are close enough to it. If a cutoff is not specified, it is
	possible for highly un-related items to be listed as neighbors of another item simply because nothing else was
	nearby. This can lead to items being joined into clusters with which they have no true connection.

}
\value{
		If \code{type} is "cluster", returns a clustering.
		This will be a vector in which the names are the compound names, and the values are the cluster
		labels.
		Otherwise, if \code{type} is "matrix", returns a nearest neighbor matrix.
		This will be a matrix with a row for each compound. Each row will contain the index value of
		the neibhboring compounds. If there are not \code{K} neibhbors for a 
		compound, that row will be padded with NAs.

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
   eiInit(sdfsample,dir=dir)

   #create compound db
   eiMakeDb(r,d,numSamples=20,dir=dir, cl=makeCluster(1,type="SOCK",outfile=""))

	print(dir())
	print(dir(dir))
	print(dir(file.path(dir,"data")))
	print(dir(file.path(dir,"run-50-40")))
	eiCluster(r,d,K=5,minNbrs=2,cutoff=0.5,dir=dir)

}
