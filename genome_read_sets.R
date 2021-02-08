library(ShortRead)
#read in csv file with genome size and read data set
genome.sets <- read.csv("genome_read_sets.csv")
#gathers the read data that we need for each of the ten sets
genome.sets[, 5] <- genome.sets[, 3] / genome.sets[ , 4]
for(i in 1:nrow(genome.sets)){
  genome.sets[i, 6:14] <- genome.sets[i, 5] *
    seq(from = .1, to = .9, by = .1)
}
for(i in 1:nrow(genome.sets)){
  #takes the SRA number
  genome <- genome.sets[i, 2]
  #stores the forward SRA file name
  genome.current.for <- paste(genome,"_1.fastq.gz", sep = "")
  #stores the reverse SRA file name
  genome.current.rev <- paste(genome, "_2.fastq.gz", sep = "")
  freads.av <- as.numeric(countFastq(genome.current.for)[1])
  rreads.av <- as.numeric(countFastq(genome.current.rev)[1])
  system(paste("mkdir", genome))
  for(frac in seq(from = .1, to = .9, by = .1)){
    print(paste("working on", genome.current.for, frac, "sampling"))
    f <- FastqSampler(genome.current.for, n = round(frac * freads.av))
    writeFastq(yield(f), paste(genome,"/", genome, ".", frac,
                               "_1.fastq.gz", sep = ""))
    f <- FastqSampler(genome.current.rev, n = round(frac * rreads.av))
    writeFastq(yield(f), paste(genome,"/", genome, ".", frac,
                               "_2.fastq.gz", sep = ""))
  }
}







#uses a reservoire sampling method; k = the number of reads that will be taken
system("paste forward.fastq reverse.fastq | awk '{ printf("%s",$0); n++; if(n%4==0)
{ printf("\n");} else { printf("\t");} }' |awk -v k=10000 'BEGIN{srand(systime() +
PROCINFO["pid"]);}{s=x++<k?x- 1:int(rand()*x);if(s<k)R[s]=$0}END{for(i in R)print R[i]}' |
awk -F"\t" '{print $1"\n"$3"\n"$5"\n"$7 > "forward_sub.fastq";print $2"\n"$4"\n"$6"\n"$8 >
"reverse_sub.fastq"}'")
