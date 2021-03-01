ParseVCF <- function(file, maxres, output){

# your vcf file
con <- file(file)

# this opens a connection to your file
open(con)

# just where you start
current.line <- 1

# this is what we will hold our useful data in
results <- as.data.frame(matrix(,maxres,7))
colnames(results) <- c("chrom", "pos","ref","alt","qual","dp4","stuff")

res.line <- 1
# now we will read one line at a time till we are done
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  #just hold first char to skip comments
  z <- strsplit(line,"")[[1]][1]
  #this if statement skips over commenting lines
  if(z != "#"){
    #i.check = 1 when there is an indel, 0 when there is not
    i.check <- length(which(startsWith(strsplit(strsplit(line,"\t")[[1]][8], ";")[[1]], "INDEL")))
    if(i.check == 0){
      #h.check = 1 when homozygous, less than 1 when heterozygous
      h.check <- as.numeric(substr(strsplit(line,"\t")[[1]][10], 1, 1))/
        as.numeric(substr(strsplit(line,"\t")[[1]][10], 3, 3))
      if(h.check < 1){
        results[res.line, 1:7] <- strsplit(line,"\t")[[1]][c(1,2,4,5,6,8,10)]
        res.line = res.line +1 
        if(res.line %% 2500 == 0){
          print(paste("identified", res.line, "heterozygous sites"))
        }
      }
    }
  }
}

close(con)

write.csv(results, output)
}