
\name{NEWS}
\title{EiR News}
\section{CHANGES IN VERSION 1.6.0 2017-4-14}{
   \subsection{NEW FEATURES}{
      \itemize{
			\item Removed almost all C/C++ code, so compiling on different platforms is easier
			\item New LSH algorithm based on Annoy library produced by Spotify (https://github.com/spotify/annoy).
			\item Performance tests are now based on Rank Based Overlap measure (http://www.williamwebber.com/research/papers/wmz10_tois.pdf).
		}
	}
	\subsection{UPGRADING}{
		\itemize{
			\item LSH parameters have changed. See man pages for eiQuery and eiCluster for details
			\item loadLSHData and freeLSHData are deprecated and no longer do anything. See man pages for details.
		}
	}
}


\section{CHANGES IN VERSION 1.4.0 2014-4-14}{
   \subsection{NEW FEATURES}{
      \itemize{
			\item eiCluster can now cluster subsets of database
			\item use new features in ChemmineR to store duplicate
			descriptors only once. 
			\item embedded descriptors now stored in database and the
			matrix file is written out only as needed by LSH to create an
			index.
		}
	}
	\subsection{UPGRADING}{
		\itemize{
			\item Database schema changes make this version incompatible with
					version 1.2 an earlier. Existing databases will need to be
					re-loaded.
		}
	}
}
\section{CHANGES IN VERSION 1.2.0 2013-9-10}{
   \subsection{NEW FEATURES}{
      \itemize{
			\item speed improvements
			\item  ieInit now accepts a SNOW cluster for parallel inserts
			\item allow compounds to be updated in-place by name
			\item eiQuery can now return similarity values instead of distances
      }
   }
}
\section{CHANGES IN VERSION 1.0.0 2013-5-26}{
   \subsection{NEW FEATURES}{
      \itemize{

         \item The eiR packages introduces efficient methods for accelerating structure
         similarity searches and clustering of very large compound datasets. The
         acceleration is achieved by applying embedding and indexing techniques to
         represent chemical compounds in a high-dimensional Euclidean space and to
         employ ultra-fast pre-screening of the compound dataset using the LSH-assisted
         nearest neighbor search in the embedding space. This method can drastically
         reduce the search time of large databases, by a factor of 40–200 fold when
         searching for the 100 closest compounds to a query. 
      }
   }
}
