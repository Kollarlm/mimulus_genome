# January 23, 2023
# Author: Kollar L.M.
# SNP density plot


## Loading in packages
library(ggplot2)
library(tidyr)
library(tidyverse)

## Loading in snp data in the format of a bed file

snps_L1<-read.table("data/L1_snps_filtered_nohead_chrONLY.bed",header=F,blank.lines.skip=TRUE)
colnames(snps_L1)<-c("chr","start","stop")
snps_L1$geno <- "L1"


snps_S1<-read.table("data/S1_snps_filtered_nohead_chrONLY.bed",header=F,blank.lines.skip=TRUE)
colnames(snps_S1)<-c("chr","start","stop")
snps_S1$geno <- "S1"

snps_all <- rbind(snps_L1,snps_S1)

# Plot the densities of snps in the bed file for each chr seperately

snphistorgram<-snps_all %>% 
  #filter(chr=="chr8") %>% 
  ggplot() + 
  geom_histogram(aes(x=start, group=geno, color=geno),binwidth=1e4) + # pick a binwidth that is not too small 
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Histogram of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density") + 
  theme_bw() # I prefer the black and white theme

snpDensity<-ggplot(snps_all, aes(x=start, group=geno, color=geno)) + 
  geom_density() + # pick a binwidth that is not too small 
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density") +
  theme_bw()

snpDensity_chr8<- snps_all %>% 
  filter(chr=="chr8") %>% 
  ggplot(aes(x=start, group=geno, color=geno)) + 
  geom_density(size = 1) + # pick a binwidth that is not too small 
  geom_vline(xintercept = 6260288, color="pink")+
  geom_vline(xintercept = 5998986, color="blue")+
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density") +
  theme_bw()


snp_chr8<- snps_all %>% 
  filter(chr=="chr8")

snpDensity_chr8<- snp_chr8 %>% 
  filter(start %in% (746467:6612646)) %>% 
  ggplot(aes(x=start, group=geno, color=geno)) + 
  geom_density(size = 1) + # pick a binwidth that is not too small 
  geom_vline(xintercept = 6465310, color="pink")+ #Start large
  geom_vline(xintercept = 858736, color="blue")+ #Start large
  geom_vline(xintercept = 1032334, color="pink",linetype="dotdash" )+ #Start small
  geom_vline(xintercept = 5998986, color="blue", linetype="dotdash")+ #Start small
  geom_vline(xintercept = 1246126, color="pink", linetype="dotdash")+ #End small
  geom_vline(xintercept = 6260288, color="blue", linetype="dotdash")+ #End small
  geom_vline(xintercept = 7604769, color="pink")+ #End large
  geom_vline(xintercept = 6465310, color="blue")+ #End large
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density") +
  theme_bw()

snpDensity_chr8<- snp_chr8 %>% 
  filter(start %in% (846467:7612646)) %>% 
  ggplot(aes(x=start, group=geno, color=geno)) + 
  geom_density(size = 1) + # pick a binwidth that is not too small 
  geom_vline(xintercept = 6465310, color="pink")+
  geom_vline(xintercept = 858736, color="blue")+
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density with large inversion starts") +
  theme_bw()

snpDensity_chr8<- snp_chr8 %>% 
  filter(start %in% (740000:7510000)) %>% 
  ggplot(aes(x=start, group=geno, color=geno)) + 
  geom_density(size = 1) + # pick a binwidth that is not too small 
  geom_vline(xintercept = 7604769, color="pink")+
  geom_vline(xintercept = 6465310, color="blue")+
  facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
  ggtitle("Density of SNPs") +
  xlab("Position in the genome") + 
  ylab("SNP density with large inversion end") +
  theme_bw()


# save the plot to .pdf file
pdf("snp_density.pdf",800,1000)
print(snpDensity)
#dev.off()
