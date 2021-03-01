#***SRA READ DOWNLOAD****
#allows you to change where SRA reads are saved to 
system("./vdb-config-i")
#pulls the forward and reverse reads from the SRA database using the above commands
system("sratoolkit.2.9.6-1-mac64/bin/./fastq-dump --gzip --skip-technical 
       --dumpbase --split-files --clip *SRA_Number*")



#***TRIMMOMATIC (TRIM READS)***
#set working directory 
setwd("~/Documents/GitHub/par")
#read in csv with SRA names
genome.sets <- read.csv("genome_read_sets.csv")
#make a loop that can automatically pull and name files for trimming and 
#runs trimmomatic
for(i in 1:nrow(genome.sets)){
        #takes the SRA number
        genome <- genome.sets[i, 2]
        #store the trimlog file name
#trimlog <- paste("/Volumes/seagate_external/trimmomatic/", 
        #genome,".trimlog", sep = "")
        trimlog <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,".trimlog", sep = "")
        #stores the forward SRA file name
#genome.current.for <- paste("/Volumes/seagate_external/sra/", 
        #genome,"_1.fastq.gz", sep = "")
        genome.current.for <- paste("/Users/michelle/Documents/GitHub/par/sra/", genome,"_1.fastq.gz", sep = "")
        #stores the reverse SRA file name
#genome.current.rev <- paste("/Volumes/seagate_external/sra/", 
        #genome, "_2.fastq.gz", sep = "")
        genome.current.rev <- paste("/Users/michelle/Documents/GitHub/par/sra/", genome, "_2.fastq.gz", sep = "")
        #store the output file names for paired and unpaired
#trimunpair.for <- paste("/Volumes/seagate_external/trimmomatic/", 
        #genome,"_1_unpair.fq", sep = "")
        trimunpair.for <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_1_unpair.fq", sep = "")
#trimunpair.rev <- paste("/Volumes/seagate_external/trimmomatic/", 
        #genome,"_2_unpair.fq", sep = "")
        trimunpair.rev <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_2_unpair.fq", sep = "")
#trimpair.for <- paste("/Volumes/seagate_external/trimmomatic/", 
        #genome,"_1_pair.fq", sep = "")
        trimpair.for <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_1_pair.fq", sep = "")
#trimpair.rev <- paste("/Volumes/seagate_external/trimmomatic/", 
        #genome,"_2_pair.fq", sep = "")        
        trimpair.rev <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_2_pair.fq", sep = "")
        #trims the paired end reads for remaining adapters and assesses quality of reads
        system(paste("java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -trimlog", 
               trimlog,
               genome.current.for,
               genome.current.rev,
               trimpair.for,
               trimunpair.for,
               trimpair.rev,
               trimunpair.rev,
               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 MINLEN:40 MAXINFO:36:0.8 TRAILING:3"))
}

#***BWA (ALIGN READS)***
setwd("/Users/michelle/Documents/GitHub/par/bwa-0.7.17")
#loop to go through BWA
for(i in 1:nrow(genome.sets)){
        #takes the SRA number
        genome <- genome.sets[i, 2]
        #takes the common name
        name <- genome.sets[i,1]
        #create input variables for file names
#forward <- paste("/Volumes/seagate_external/trimmomatic/", genome,
        #"_1_unpair.fq", sep = "")
        forward <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_1_unpair.fq", sep = "")
#reverse <- paste("/Volumes/seagate_external/trimmomatic/", genome,
        #"_2_unpair.fq", sep = "")
        reverse <- paste("/Users/michelle/Documents/GitHub/par/trimmomatic/", genome,"_2_unpair.fq", sep = "")
        #create output variables for file names
#align <- paste("/Volumes/seagate_external/ref/", name, 
        #".fna.gz", sep = "")
        align <- shQuote(paste("/Users/michelle/Documents/GitHub/par/ref/", name, ".fna.gz", sep = ""))
#bwamem <- paste("/Volumes/seagate_external/bwa_mem/", genome,
        #"_aligned.sam", sep = "")
        bwamem <- paste("/Users/michelle/Documents/GitHub/par/bwa_mem/", genome, "_aligned.sam", sep = "")
        #first index the alignment file to align the reads to
        system(paste("./bwa index -p", name, "-a bwtsw", align))
        #unzip the alignement file
        system(paste("gunzip", align))
        #system command to run bwa mem for mapping
        system(paste("./bwa mem", name, forward, reverse, ">", bwamem))
}


#***SAMTOOLS***

setwd("/Users/michelle/Documents/GitHub/par/samtools-1.11")
#loop to go through samtools functions
for(i in 1:nrow(genome.sets)){
        #takes the SRA number
        genome <- genome.sets[i, 2]
        #create output variables for file names
#bwabam <- paste("/Volumes/seagate_external/samtools/", genome,
        #"_aligned.bam", sep = "")
        bwabam <- paste("/Users/michelle/Documents/GitHub/par/samtools/", genome,
                        "_aligned.bam", sep = "")
#bamsort <- paste("/Volumes/seagate_external/samtools/", genome,
        #"_alignsort.bam", sep = "")
        bamsort <- paste("/Users/michelle/Documents/GitHub/par/samtools/", genome,
                        "_alignsort.bam", sep = "")
#depth <- paste("/Volumes/seagate_external/samtools/", genome,
        #".depth", sep = "")        
        depth <- paste("/Users/michelle/Documents/GitHub/par/samtools/", genome,
                       ".depth", sep = "")
        #system command to transform sam into bam
        system(paste("./samtools view -Sb", bwamem, ">", bwabam))
        #system command to sort the bam file
        system(paste("./samtools sort -o", bamsort, bwamem))
        #system command to find depth at each region within the bam file
        system(paste("./samtools depth -a", bwamem, ">", depth))
}




#***BCFTOOLS***
setwd("/Users/michelle/Documents/GitHub/par/bcftools-1.11")
for(i in 1:nrow(genome.sets)){
        #takes the SRA number
        genome <- genome.sets[i, 2]
        #create input variables for creating bcf file
#align_unzip <- paste("/Volumes/seagate_external/ref/", name, 
        #".fna", sep = "")
        align_unzip <- shQuote(paste("/Users/michelle/Documents/GitHub/par/ref/", name, ".fna", sep = ""))
        #create output variables for creating bcf file
#bcf <- paste("/Volumes/seagate_external/bcf/", genome,
        #"_raw.bcf", sep = "")
        bcf <- paste("/Users/michelle/Documents/GitHub/par/bcf/", genome, "_raw.bcf", sep = "")
        #system command using bcftools to creaete bcf file for variant calling
        system(paste("./bcftools mpileup -O b -o", bcf, "-f", align_unzip, bamsort))
       #create output variables for creating vcf file
#bwabam <- paste("/Volumes/seagate_external/bcf/", genome,
        #"_raw.bcf", sep = "")
        vcf <- paste("/Users/michelle/Documents/GitHub/par/bcf/", genome, ".vcf.gz", sep = "")
        #system command using bcftools to call nucleotides and detect SNPs
        system(paste("./bcftools call -vmO z -m", bcf, ">", vcf))
}
        
        
        
        