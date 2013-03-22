\name{eiPerformanceTest}
\alias{eiPerformanceTest}
\title{
   Test the performance of LSH search
}
\description{
   Tests the performance of embedding and LSH.
}
\usage{
	eiPerformanceTest(r,d,distance=getDefaultDist(descriptorType),descriptorType="ap",
		dir=".",K=200, W = 1.39564, M=19,L=10,T=30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{r}{
      The number of references used to build the database you wish to
      query agains.
   }
  \item{d}{
      The number of dimensions used to build the database you wish to
      query agains.
   }


	\item{distance}{
		The distance function to be used to compute the distance between two descriptors. A default function is
		provided for "ap" and "fp" descriptors.
	}
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
 \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }

  \item{K}{
     Number of search results to use for LSH performance test.
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
   This will perform two different tests.
   The first tests the embedding results in similarity search. The way this works is by
   approximating 1,000 random similarity searches (determined by
   data/test_queries.iddb) by nearest neighbor search using the coordinates
   from the embedding results. The search results are then compared to the
   reference search results (chemical-search.results.gz). 

   The comparison results are summarized in two types of files. The first
   type lists the recall for different k values, k being the number of
   numbers to retrieve. These files are named as ``recall-ratio-k''.
   For example, if the recall is 70% for top-100
   compound search - 70 of the 100 results are among the real top-100
   compounds - then the value at line 100 is 0.7. Several relaxation ration
   are used, each generating a file in this form. For instance,
   recall.ratio-10 is the file listing the recalls when relaxation
   ratio is 10. The other file, recall.csv, lists recalls of different
   relaxation ratios in one file by limiting to selected k value. In this
   CSV file, the rows correspond to different relaxation ratios, and the
   columns are different k values. You will be able to pick an appropriate
   relaxation ratio for the k values you are interested in.


   The second test measures the performance of the Locality Sensitive Hash (LSH).
   The results for lsh-assisted search will be in
   run-r-d/indexed.performance. It's a 1,000-line files of recall values. Each
   line corresponds to one test query.  LSH search performance is
   highly sensitive to your LSH parameters (K, W, M, L, T). The
   default parameters are listed in the man page for
   \code{eiPerformanceTest}. When you have your embedding result in
   a matrix file, you should follow instruction on
   \url{http://lshkit.sourceforge.net/dd/d2a/mplsh-tune_8cpp.html} to
   find the best values for these parameters.
}
\value{
   No value is returned. Creates files in \code{dir}/run-r-d.
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
   eiMakeDb(r,d,numSamples=20,dir=dir,
      cl=makeCluster(1,type="SOCK",outfile=""))

   eiPerformanceTest(r,d,dir=dir,K=22)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
