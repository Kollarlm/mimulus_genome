#Set some variables
tandem_cutoff=10 #Must be within this number of genes of a gene with same arrayID to be classified as "tandem"

#Load libraries
library(GENESPACE)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)

runwd <- file.path("/Users/lesliekollar/Desktop/mimulus_genome_R/mimulus_genome/")
setwd("/Users/lesliekollar/Desktop/mimulus_genome_R/mimulus_genome/")

#Read in the pangenome database
pgdb <- fread("results/S1_pangenomeDB.txt.gz",na.strings = c("", "NA"))
pgff <- fread("results/gffWithOgs.txt.gz",na.strings = c("", "NA"))

#Create two new tables for each of the genomes, removing genes not in the pangenome
L1 <- subset(pgdb, genome=="L1" & !is.na(pgChr)) # !is.na = is not blank/empty ... subsetting only the genes that do not have missing values. 
S1 <- subset(pgdb, genome=="S1" & !is.na(pgChr)) 

#Create a separate list of unique genes to each genome, not in pangenome
L1u <- subset(pgdb, genome=="L1" & is.na(pgChr)) # They are not in the pangenome but unique to the specific genotype
S1u <- subset(pgdb, genome=="S1" & is.na(pgChr))

#Combine the files based on the pangenome ID
LS <- full_join(L1,S1,by=("pgID")) # Making an master list of genes for each genotype that are in the pangenome.

#Calculate putative CNV & PAV from pgID
cnv <- dcast(subset(pgdb, !is.na(pgChr)), pgChr + pgOrd + pgID ~ genome, value.var="id", fun.aggregate=length) # Only using rows with complete cases (genes in pangenome), turn the genome variables (L1 and S1) into columns
                                                                                                              # and counts based on the 3 columns to the left of ~. For example it counts the number of occurences of each ID.
#If isDirectSyn & isArrayRep is TRUE for both paired genes in both genomes, it is a syntenic pair
#Subset syntenic genes and write out lists of these genes
synt <- subset(LS, isArrayRep.x & isDirectSyn.x & isArrayRep.y & isDirectSyn.y) #To be syntenic isArrayRep and isDirectSyn have to be true 
write.table(unique(synt$id.x),file="results/L1_syntelogs.txt", 
            quote=FALSE,col.names=FALSE,row.names=FALSE) #Write table of syntenic genes for L1 and S1
write.table(unique(synt$id.y),file="results/S1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get High-confidence 1x1 syntelogs and write out lists of these genes
synt1x1 <- subset(synt, pgID %in% subset(cnv, L1==1 & S1==1)$pgID) ## the list of syntenic genes between L1 and S1 based on the matching pgID in the cnv file where there is only 1 copy number in L1 and S1
write.table(synt1x1$id.x,file="results/L1_1x1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(synt1x1$id.y,file="results/S1_1x1_syntelogs.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
#Get sets of non-syntenic genes and write out lists of these genes
nonsynt <- setdiff(LS, synt) #Removing the synetic genes from the pangenome so that what is left is the non syntenic genes.
L1nonsynt <- subset(nonsynt, !(id.x %in% synt$id.x)) 
S1nonsynt <- subset(nonsynt, !(id.y %in% synt$id.y))
write.table(c(unique(L1nonsynt$id.x),L1u$id),file="results/L1_nonsyntenic.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(c(unique(S1nonsynt$id.y),S1u$id),file="results/S1_nonsyntenic.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

#Get genes in arraypairs arrays
#Add the arrayID to pgdb
pge <- setorder(merge(pgdb,pgff[,c(2,11)],by="id"), pgChr, pgOrd, na.last = T)
#Count count arrays with more than one gene in them
arraycount <- as.data.frame(table(unique(pge[,c(1,15)])$arrayID))
#subset out those genes with an array count of more than 1
arraypairs <- setorder(subset(pge, arrayID %in% subset(arraycount, Freq > 1)$Var1), arrayID)
#Create output dataframe
columns=c("arrayID","id","isTandem","geneCount","arrayGenes")
tandem <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(tandem) <- columns
#Now for a complicated loop to identify tandem arrays
for(i in unique(arraypairs$arrayID)){
  #some arraypairs arrays may belong to multiple orthogroups, collapse these so each gene is represented once
  x <- setorder(unique(subset(arraypairs, arrayID==i)[,c(1,9,10,15)]), ord)
  #Check to make sure all of these are on the same chromosome!
  if(uniqueN(x$chr)==1){
    #Most array pairs are simple arraypairs duplicates and faster to process
    #If the number of genes is 2 and the distance is <= our cutoff classify as arraypairs duplicates
    if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) <= tandem_cutoff){
      for(gene in x$id){
        tandem <- rbind(tandem,data.frame(arrayID=i,id=gene,isTandem=TRUE,geneCount=2,
                                          arrayGenes=paste(x$id,collapse=", ")))
      } 
      #If the number of genes is 2 and distance > cutoff, classify as a dispersed
    } else if(nrow(x)==2 & abs(x[1]$ord - x[2]$ord) > tandem_cutoff) {
      tandem <- rbind(tandem,data.frame(arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
                                        arrayGenes=NA))
      #If the arraypairs have more than one gene, it requires more work
    } else {
      #Create an empty dataframe to calculate pairwise distances between all genes
      x2 <- data.frame(matrix(ncol=nrow(x),nrow=nrow(x)))
      #Set row & column names for the table
      row.names(x2) <- colnames(x2) <- x$id
      #Calculate absolute distances between all genes
      for(gene in 1:nrow(x)){
        x2[,gene] <- abs(x$ord - x[gene]$ord)
      }
      for(gene in x$id){
        #tandem if a gene's closest array pair is greater than the cutoff
        #If so, classify as FALSE
        if(min(x2[x2[,gene] != 0,gene]) > tandem_cutoff){
          tandem <- rbind(tandem,data.frame(arrayID=i,id=gene,isTandem=FALSE,geneCount=NA,
                                            arrayGenes=NA))
          #Otherwise, the gene is part of an arraypairs array
        } else {
          #Get a list of all those genes within the cutoff distance of the arraypairs array
          x3 <- row.names(x2[x2[gene,] <= tandem_cutoff,])
          #Now we are going to loop back over those genes in that initial list and search
          #for other genes which may be within range of the arraypairs array, but not the initial gene
          #these will get added to the list
          newgenes <- 1 #set this to initialize the for loop
          #Keep the loop going until no new genes are added
          while(newgenes != 0){
            #Create x4 list same as x3
            x4 <- x3
            #Loop over genes in x3
            for(gene2 in x3){
              #Extend out the array and add the new genes to x4
              x4 <- unique(c(x4,row.names(x2[x2[gene2,] <= tandem_cutoff,])))
            }
            #Set newgenes to the difference between x4 & x3
            #The loop will terminate when no new genes are added
            newgenes <- length(x4) - length(x3)
            #Set x3 to now be equal to x4
            x3 <- x4
          }
          #Now we wll report it
          tandem <- rbind(tandem,data.frame(arrayID=i,id=gene,isTandem=TRUE,
                                            geneCount=length(x3),arrayGenes=paste(x3,collapse=", ")))
        }
      }
    }	
    #If the genes are not on the same chromosome, they can't be a arraypairs array. Throw a warning.
  } else {
    print(paste("Warning! Array ID",unique(x$arrayID),"genes on different chromosomes. Skipping",sep=" "))
  }
}

table(tandem$geneCount)
table(tandem$isTandem)
table(unique(tandem[,c(1,4)])$geneCount)

# Adding gene id to CNV
dictionary <- pgdb[,1:3] 
dictionary <- cbind(dictionary,pgdb$id) 
dictionary$id <- dictionary$V2 # renaming column


dictionary2 <- dictionary %>% 
  group_by(pgChr, pgOrd, pgID ) %>% #Gather by the columns that are repetative
  summarise_at(vars(-group_cols()), ~toString(.[!is.na(.)], collapse="_")) # Collapse the id column into a single cell seperated by ","

dictionary2$V2 <- NULL #removing column

cnv_withID <- cnv %>% left_join(dictionary2,by=c("pgChr", "pgOrd" ,"pgID"))
write.table(cnv_withID, file="supplemental_tables/cnv.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)


### Sub-setting the cnv dataframe 

# Genes absent in one ecotype but not the other... specific to L1 or S1
absent_s1_present_L1 <- cnv_withID[sapply(lapply(cnv$S1, `%in%`, "0"), any),]
write.table(absent_s1_present_L1, file="supplemental_tables/absent_S1_present_L1.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

absent_L1_present_S1 <- cnv_withID[sapply(lapply(cnv$L1, `%in%`, "0"), any),]
write.table(absent_L1_present_S1, file="supplemental_tables/absent_L1_present_S1.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)

# Copy number of genes in IM62, S1, and L1
cnv_withID$cn_diff <- cnv_withID$L1 - cnv_withID$S1 


# Same number of copies and not a CNV
same_cn <- cnv_withID %>% 
  subset(cn_diff == 0)

# Different number of copies and thus CNV.
# This includes PAVs as well
diff_cn <- cnv_withID %>% 
  subset(cn_diff != 0)
write.table(diff_cn, file = "supplemental_tables/IM62_L1_S1_different_copy_numbers.txt",
            quote=FALSE, col.names = FALSE, row.names = FALSE)

# Tandem arrays but this includes IM62...
tandem_arrays <- as.data.frame(table(unique(tandem[,c(1,4)])$geneCount))

#Plotting tandem arrays
ggplot(data=tandem_arrays, aes(y=tandem_arrays$Freq, x=tandem_arrays$Var1)) +
  geom_bar(stat="identity", fill= "chartreuse4")+
  #geom_text(aes(label=tandem_arrays$Freq), vjust=1.6, color="black", size=3) +
  labs(y= "Frequency", x= "Size of Tandem Array")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
    

