\name{eiInit}
\alias{eiInit}
\title{
   Initialize a compound database
}
\description{
   Takes the raw compound database in whatever format the given
   measure supports and creates a "data" directory.
}
\usage{
	eiInit(inputs,dir=".",format="sdf",descriptorType="ap",append=FALSE,
	conn=defaultConn(dir,create=TRUE), updateByName = FALSE, cl = NULL, connSource = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputs}{
	  Either a filename of a file in \code{format} format, or an SDFset. This can
	  also be a vector of  filenames and if \code{cl} is also specified and if you database
	  supports it (SQLite does not), it will load these file in parallel on the cluster.
   }
  \item{dir}{
      The directory where the "data" directory lives. Defaults to the
      current directory.
   }
	\item{format}{
		The format of the data in \code{inputs}. Currenly only "sdf" and "smiles" is 
		supported.
	}
	\item{descriptorType}{
		The format of the descriptor. Currently supported values are "ap" for atom pair, and 
		"fp" for fingerprint.
	}
	\item{append}{
		If true the given compounds will be added to an existing database
		and the <data-dir>/Main.iddb file will be udpated with the new
		compound id numbers. This should not normally be used directly, use
		\code{\link{eiAdd}} instead to add new compounds to a database.
	}
	\item{conn}{
		Database connection to use. If a connection is given, you must ensure that it has been initialized using 
		the \code{\link{initDb}} function from ChemmineR before calling \code{\link{eiInit}}.
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
  \item{cl}{
     A SNOW cluster can be given here to run this function in
     parrallel.
   }
	\item{connSource}{
		A function returning a new database connection. Note that it is not suffient to return a
		reference to an existing connection, it must be a distinct, new connection.  
		This is needed for cluster operations
		that make use of the database as they will each need to craete a new connection.
		If not given, certain parts of this function will not be parrallelized.

		This function can also be used to setup the envrionment on the cluster worker nodes. For
		example, you might need to re-load libraries like RSQLite and such.
	}

}
\details{


   EiInit can take either an SDFset, or a filename.  SDF and SMILES is supported
   by default.
   It might complain if your SDF file does not
   follow the SDF specification. If this happens, you can create an
   SDFset with the \code{read.SDFset} command and then use that
   instead of the filename.  
	
	EiInit will create  a folder called
   'data'. Commands should always be executed in the folder containing
   this directory (ie, the parent directory of "data"), or else
   specify the location of that directory with the \code{dir} option.

}
\value{
   A directory called "data" will have been created in the current working directory.
	The generated compound ids of the given compounds will be returned. These can be used to 
	reference a compound or set of compounds in other functions, such as \code{\link{eiQuery}}.
}
\author{
   Kevin Horan
}


\seealso{
   \code{\link{eiMakeDb}}
   \code{\link{eiPerformanceTest}}
   \code{\link{eiQuery}}
}
\examples{
   data(sdfsample)
   dir=file.path(tempdir(),"init")
   dir.create(dir)
   eiInit(sdfsample,dir=dir)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
