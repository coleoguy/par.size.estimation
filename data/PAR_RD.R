#Michelle Jonika
#Sliding across the genome to get average read depth in sliding windows


#read in the libraries needed for the analysis
library(data.table)
library(evobiR)
#read in the depth file
#fread is necessary because of the size of the file being read in 
depth <- fread("ERR5101147_1.depth")

#idnetifies the column names for the depth file
colnames(depth) <- c("chr", "locus", "depth")

#code that identifies the two chromosomes
#each species should have an autosome (chromosome 3)
#and the x chromosome
unique(depth$chr)

#NC_000003.12 (autosome, chromosome 3)
#NC_000023.11 (x chromosome)

#separate out chromosome 3 from the x chromosome
chrom.auto <- which(depth$chr == unique(depth$chr)[1])
chrom.auto.depth <- depth[chrom.auto, ]

chrom.x <- which(depth$chr == unique(depth$chr)[2])
chrom.x.depth <- depth[chrom.x, ]

#find the average read depth across the autosome (chr3)
chrom.auto.depth.avg <- mean(chrom.auto.depth$depth)

#standardize the x chromosome read depth data to the average read depth
#across the autosome
standardized <- (chrom.x.depth[,3])/chrom.auto.depth.avg

rm(chrom.auto.depth, chrom.x.depth, depth, chrom.auto, chrom.auto.depth.avg, chrom.x)
standardized[standardized<0.1] <- NA
standardized[standardized>1.9] <- NA

#calculate 3 SDs for the standardized read data to adjust for those sites 
#that have inappropriately high read depth values along the length of the 
#chromosome
sd <- (sd(standardized$depth)) * 3

#perform a sliding window approach to calculate the average
#read depth along chromosome 3

#FUN is the function to be applied across the data
#data is the data that you want the function to be applied to
#window is the size of the window 
#step is the size of the step between windows

FUN <- function(standardized){
  mean_slide <- mean(standardized, na.rm = TRUE)
  return(mean_slide)
}
test <- SlidingWindow(FUN, standardized$depth, window = 50000, step = 2500)
#plot the test results from the sliding window
plot(test,type="l", ylim = c(0,2))

#calculate the average RD across the x chromosome
chrom.x.depth.avg <- mean(standardized$depth, na.rm = TRUE)
#calculate the average RD across the x chromosome in the non-PAR region
chrom.x.depth.avg.NP <- mean(standardized$depth[2600001:156040895], na.rm = TRUE)
#calculate the average RD across the x chromosome in the PAR region
chrom.x.depth.avg.P <- mean(standardized$depth[1:2600000], na.rm = TRUE)

#set a threshold for where to visualize the cutoff 
abline(h = chrom.x.depth.avg.P, col = "green")

abline(h = chrom.x.depth.avg.NP, col = "red")

#calculate the 80% CI
error_80 <- qnorm(0.90)*(sd(standardized$depth, na.rm = TRUE))/sqrt(length(standardized$depth))
left_80 <- (mean(standardized$depth, na.rm = TRUE)) - error_80
right_80 <- (mean(standardized$depth, na.rm = TRUE)) + error_80


#calculate the 90% CI
error_90 <- qnorm(0.95)*(sd(standardized$depth, na.rm = TRUE))/sqrt(length(standardized$depth))
left_90 <- (mean(standardized$depth, na.rm = TRUE)) - error_90
right_90 <- (mean(standardized$depth, na.rm = TRUE)) + error_90

#calculate the 95% CI
error_95 <- qnorm(0.975)*(sd(standardized$depth, na.rm = TRUE))/sqrt(length(standardized$depth))
left_95 <- (mean(standardized$depth, na.rm = TRUE)) - error_95
right_95 <- (mean(standardized$depth, na.rm = TRUE)) + error_95
