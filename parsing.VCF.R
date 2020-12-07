# your vcf file
con <- file('SRR5229975.vcf')

# this opens a connection to your file
open(con)

# just where you start
current.line <- 1

# this is what we will hold our useful data in
results <- as.data.frame(matrix(,1,7))
colnames(results) <- c("chrom", "pos","ref","alt","qual","dp4","stuff")

res.line <- 1
# now we will read one line at a time till we are done
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  #just hold first char to skip comments
  z <- strsplit(line,"")[[1]][1]
  if(z != "#"){

    results[res.line, 1:6] <- strsplit(line,"\t")[[1]][c(1,2,4,5,6,8)]
    info <- results[res.line, 6]
    info <- strsplit(info, ";")[[1]]
    dp4 <- info[which(startsWith(info, "DP4"))]
    ind <- info[which(startsWith(info, "INDEL"))]
    if(length(ind) == 0){
    # pull out any other pieces you want from info
    results[res.line, 6] <- dp4
    # results[res.line, 7] <- blah
    res.line <- res.line + 1
    print(res.line)
    }
  }
  current.line <- current.line + 1
}
close(con)


