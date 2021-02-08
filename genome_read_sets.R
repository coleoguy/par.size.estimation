#read in csv file with genome size and read data set
genome.sets <- data.frame(read.csv("genome_read_sets.csv"), stringsAsFactors = F)

#gathers the read data that we need for each of the ten sets
for(i in 1:nrow(genome.sets)){
  genome.sets[[5]][i] <- genome.sets[[3]][i]/
    genome.sets[[4]][i]
  genome.sets[i,6] <- genome.sets[i,5]*0.1
  genome.sets[i,7] <- genome.sets[i,5]*0.2
  genome.sets[i,8] <- genome.sets[i,5]*0.3
  genome.sets[i,9] <- genome.sets[i,5]*0.4
  genome.sets[i,10] <- genome.sets[i,5]*0.5
  genome.sets[i,11] <- genome.sets[i,5]*0.6
  genome.sets[i,12] <- genome.sets[i,5]*0.7
  genome.sets[i,13] <- genome.sets[i,5]*0.8
  genome.sets[i,14] <- genome.sets[i,5]*0.9
}

for(i in 1:length(genome.sets)){
  #takes the SRA number
  genome <- as.character((genome.sets[i,2]))
  #stores the forward SRA file name
  genome.current.forward <- paste(genome,"_1.fastq.gz", sep = "")
  #stores the reverse SRA file name
  genome.current.reverse <- paste(genome, "_2.fastq.gz", sep = "")
  #count lines in the SRA files
  lines.count <- system(paste("wc -l",genome.current.forward), intern = T)
  #string split the output to isolate the line number
  count <- strsplit(lines.count, split=" ")
  number.lines <- as.numeric(count[[1]][2])
  
  #making the dataset with 10% of the data
  #new forward SRA file name for 10%
  genome.new.forward <- paste(genome, ".10_1.fastq.gz", sep = "")
  #new reverse SRA file name 
  genome.new.reverse <- paste(genome, ".10_2.fastq.gz", sep = "")
  #count the number of lines for 10% of the read depth
  numberlines.10 <- round((number.lines*0.1))
  system(awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < ERR3465834_1.fastq.gz | awk 'NR>1{ printf("%s",$0); n++; if(n%2==0) { printf("\n");} else { printf("\t");} }' | awk -v k=10000 'BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x- 1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' | awk -F"\t" '{print $1"\n"$2 > ERR3465834.10_1.fastq.gz}')
  
  #sample fastq file
  fqsample <- readFastq(dirPath = "/Volumes/seagate_external/sra", pattern = "fastq")
  #making the dataset with 20% of the data
  numberlines.20
  #making the dataset with 30% of the data
  numberlines.30
  #making the dataset with 40% of the data
  numberlines.40
  #making the dataset with 50% of the data
  numberlines.50
  #making the dataset with 60% of the data
  numberlines.60
  #making the dataset with 70% of the data
  numberlines.70
  #making the dataset with 80% of the data
  numberlines.80
  #making the dataset with 90% of the data
  numberlines.90
  
  
}



#uses a reservoire sampling method; k = the number of reads that will be taken 
system("paste forward.fastq reverse.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) 
{ printf("\n");} else { printf("\t");} }' |awk -v k=10000 'BEGIN{srand(systime() + 
PROCINFO["pid"]);}{s=x++<k?x- 1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "forward_sub.fastq";print $2"\n"$4"\n"$6"\n"$8 > 
"reverse_sub.fastq"}'")