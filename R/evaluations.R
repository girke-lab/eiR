compareListsRBO<- function(file1,file2)
{
	if(debug) print("comparing search results")
   getRanking<- function(line) 
      as.numeric(sapply(line,
                     function(y) strsplit(y,":",fixed=TRUE)[[1]][1]))

   in1=file(file1,"r")
   in2=file(file2,"r")

   p=NA

   results=c()
   while(TRUE){
		suppressWarnings({
			line1=readLines(in1,n=1)
			line2=readLines(in2,n=1)
		})

      if(length(line1)==0 || length(line2)==0)
         break;

      line1=unlist(strsplit(line1,"\\s+"))
      line2=unlist(strsplit(line2,"\\s+"))


      ind1=1:length(line1)
		names(ind1)=as.character(getRanking(line1))

      ind2=1:length(line2)
		names(ind2) = as.character(getRanking(line2))

		#print(ind1)
		#print(ind2)
      ##results=c(results,length(intersect(ind1,ind2))/p)
		rboDiff = rbo(ind1,ind2,p=0.9,side="bottom")
		#message("------------ diff:",rboDiff)
      results=c(results, rboDiff)
   }
   close(in2)
   close(in1)
	#message("average rank difference using RBO: ",mean(results))
   #print(results)
   return(results)
}



## The following three functions were taken from the gespeR package


#' Rank biased overlap (Webber et al., 2010)
#'
#' Evaluates the rank biased overlap (rbo) of two ranked lists based on formula based on (32) from
#' "A Similarity Measure for Indefinite Rankings" (Webber et al.). Two ranked lists with high rbo are
#' very similar, wheras low rbo indicates dissimilar lists. rbo ranges between 0 and 1. In this method
#' the extrapolated version of rbo is implemented.
#'
#' @author Fabian Schmich
#' @export
#'
#' @param s List 1
#' @param t List 2
#' @param p Weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k Evaluation depth for extrapolation
#' @param side Evaluate similarity between the top or the bottom of the ranked lists
#' @param mid Set the mid point to for example only consider positive or negative scores
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return rank biased overlap (rbo)
#'
#' @seealso \code{\link{concordance}}
#'
#' @examples
#' a <- rnorm(26)
#' b <- rnorm(26)
#' names(a) <- names(b) <- LETTERS
#' rbo(a, b, p = 0.95)
rbo <- function(s, t, p, k=floor(max(length(s), length(t))/2), side=c("top", "bottom"), mid=NULL, uneven.lengths = TRUE) {
  side <- match.arg(side)
  if (!is.numeric(s) | !is.numeric(t))
    stop("Input vectors are not numeric.")
  if (is.null(names(s)) | is.null(names(t)))
    stop("Input vectors are not named.")
  ids <- switch(side,
                "top"=list(s=.select.ids(s, "top", mid), t=.select.ids(t, "top", mid)),
                "bottom"=list(s=.select.ids(s, "bottom", mid), t=.select.ids(t, "bottom", mid))
  )
  min(1, .rbo.ext(ids$s, ids$t, p, k, uneven.lengths = uneven.lengths))
}

#' Select top or bottom names of ranked vector
#'
#' @author Fabian Schmich
#' @noRd
#'
#' @param x The ranked list
#' @param side The side to be evaluated ("top" or "bottom" of ranked list)
#' @param mid The mid point to split a list, e.g. to split between positive and negative values choose mid=0
#' @return A vector of selected identifiers
.select.ids <- function(x, side=c("top", "bottom"), mid=NULL) {
  side <- match.arg(side)
  if (side == "top")  {
    x <- sort(x, decreasing=TRUE)
    if (is.null(mid))
      return(names(x))
    else
      return(names(x)[which(x > mid)])
  } else if (side == "bottom") {
    x <- sort(x, decreasing=FALSE)
    if (is.null(mid))
      return(names(x))
    else
      return(names(x)[which(x < mid)])
  }
}

#' Rank biased overlap formula based on (32) from "A Similarity Measure for Indefinite Rankings" (Webber et al.)
#'
#' @author Fabian Schmich
#' @noRd
#'
#' @param x List 1
#' @param y List 2
#' @param p The weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k The evaluation depth
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return The rank biased overlap between x and y
.rbo.ext <- function(x, y, p, k, uneven.lengths = TRUE) {
  if (length(x) <= length(y)) {
    S <- x
    L <- y
  } else {
    S <- y
    L <- x
  }
  l <- min(k, length(L))
  s <- min(k, length(S))

  if (uneven.lengths) {
    Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
    ((1-p) / p) *
      ((sum(Xd[seq(1, l)] / seq(1, l) * p^seq(1, l))) +
         (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
      ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l
  } else {
    #stopifnot(l == s)
    k <- min(s, k)
    Xd <- sapply(1:k, function(i) length(intersect(x[1:i], y[1:i])))
    Xk <- Xd[k]
    (Xk / k) * p^k + (((1-p)/p) * sum((Xd / seq(1,k)) * p^seq(1,k)))
  }
}
