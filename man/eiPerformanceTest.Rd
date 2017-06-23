\name{eiPerformanceTest}
\alias{eiPerformanceTest}
\title{
   Test the performance of LSH search
}
\description{
   Tests the performance of embedding and LSH.
}
\usage{
	eiPerformanceTest(runId,distance=getDefaultDist(descriptorType),
		conn=defaultConn(dir),dir=".",K=200, searchK=-1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{runId}{
		The id number identifying a particular set of settings for a database. This is generally
		the number returned by \code{\link{eiMakeDb}}. If your coming from an older version of eiR, you should
		not use this value instead of
		specifying \code{r}, \code{d}, and \code{descriptorType}.
	}

	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}
	\item{conn}{
		Database connection to use.
	}
 \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }

  \item{K}{
     Number of search results to use for LSH performance test.
   }
  \item{searchK}{
     Tunable Annoy LSH parameter. A larger value will give more accurate results, but will take longer time to return.
	  See Annoy page for details.
     \url{https://github.com/spotify/annoy}
   }
}
\details{
	This function can be used to tune the two Annoy LSH parameters, \code{numTrees}, and \code{searchK}. 

	NumTrees is provided to the  eiMakeDb function and affects the build time and the index size. A larger value will produce
	more accurate results, but use more disk space. 

	SearchK is given to the eiQuery function, or to this function. A larger value will give more accurate results, but will require more time to run. 


   This function will perform two different tests.
	The first test is how well the embedding is working. When the 
	\link{eiMakeDb} function is run, you can specify the number of test samples
	to use for this test. If not specified, it will default to 10\% of the data set
	size. During this test, we take each sample and compute its distance to 
	every other compund in the dataset using both the given descriptor distance function (e.g., "AP" or "fingerprint"), 
	as well as the euclidean distance computed on the embedded version. We
	then measure how similar the resulting ranks of these lists are using 
	Rank Based Overlap (Webber,2010) (\url{http://www.williamwebber.com/research/papers/wmz10_tois.pdf}).
	The similarity for each sample is output in a file called 'embedding.performance' in the work directory. Each
	line corresponds to one sample. 

	The second test compares the rankings produced using the descriptor distance function, to the
	rankings produced by the final output of the LSH search, for each sample query. Again, rank 
	based overlap (RBO) is used to compare the rankings. The results are output in the same format
	as for the fist test, in a file called 'indexed.performance'.

	RBO is a similarity measure that produces a value in the range of [0,1]. Values closer to 0 are very
	dissimilar, while values closer to 1 are more similar. 

}
\value{
   Returns the results of the indexing test. Each element of the resulting
	vector is the RBO similarity for the coresponding query. 
	Creates files in \code{dir}/run-r-d.
}
\author{
   Kevin Horan
}

\seealso{
   \code{\link{eiInit}}
   \code{\link{eiMakeDb}}
   \code{\link{eiQuery}}
}
\examples{
   library(snow)

   r<- 50
   d<- 40

   #initialize 
   data(sdfsample)
   dir=file.path(tempdir(),"perf")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)

   #create compound db
   runId = eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(1,type="SOCK",outfile=""))

   eiPerformanceTest(runId,dir=dir,K=22)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
