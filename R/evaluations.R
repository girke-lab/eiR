compareListsRBO<- function(file1,file2)
{
	if(debug) print("comparing search results")
   getRanking<- function(line) 
      as.numeric(sapply(line,
                     function(y) strsplit(y,":",fixed=TRUE)[[1]][1]))

   in1=gzfile(file1,"r")
   in2=gzfile(file2,"r")

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
      #results=c(results,length(intersect(ind1,ind2))/p)
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
