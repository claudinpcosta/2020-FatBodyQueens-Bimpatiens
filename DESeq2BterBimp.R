#################################Bimp <-> Bter Analysis by DESeq2##################################

#downloaded file from Ensembl with orthologs between Apis and Bter 
BterApis <- read.csv("BterApis.csv", header = T, sep = ",")
head(BterApis)
tail(BterApis)
dim(BterApis)

#downloaded file from Amselem et al 2015, readcounts from Bter.
Amselem <- read.csv("readcounts_AmselemBter.csv", header = T, sep = ",")
head(Amselem)
tail(Amselem)
dim(Amselem)

#match the ID between Bter and dataset from Amselem et al. 2015
BterAmselem <- Amselem[which(Amselem$Bter_ID %in% BterApis$Bter_ID),]
head(BterAmselem)
tail(BterAmselem)
dim(BterAmselem)
BterAmselem <- BterAmselem[order(BterAmselem$Bter_ID),]
head(BterAmselem)
dim(BterAmselem)


Amselem_BterApis <- merge(BterAmselem,BterApis, by = "Bter_ID")
dim(Amselem_BterApis)
head(Amselem_BterApis)
A_BterApis <- unique(Amselem_BterApis)
dim(A_BterApis)
head(A_BterApis)

    
#file from our data ortho Apis <-> Bimp 
BimpApis <- read.csv("BimpApis_mod.csv", header = T, sep = ",")
head(BimpApis)
tail(BimpApis)
dim(BimpApis)


BterBimp <- A_BterApis[which(A_BterApis$Amel.GB %in% BimpApis$Amel.GB),]
head(BterBimp)
tail(BterBimp)
dim(BterBimp)


A_BterBimpApis <- merge(BterBimp,BimpApis, by = "Amel.GB")
dim(A_BterBimpApis)
head(A_BterBimpApis)

ReadCounts_Amselem <- A_BterBimpApis[order(A_BterBimpApis$Amel.GB),]
ReadCounts_Amselem$Amel.GB <- NULL 
ReadCounts_Amselem$Bter_ID <- NULL
      #switch the order of column 
      library(dplyr) 
      ReadCounts_Amselem <- ReadCounts_Amselem %>% select(Bimp, everything())
      #Convert the values in first column into row names
      ReadCounts <- ReadCounts_Amselem[-1]
      rownames(ReadCounts) <- make.names(ReadCounts_Amselem[,1], unique = TRUE)
row.names(ReadCounts)
head(ReadCounts)
dim(ReadCounts)

###############################DESeq2 analysis#####################################################

#Get gene lists for Amselem dataset
      
#Loading packages (Whenever you start R)
library("DESeq2") 

readcounts <- ReadCounts
head(readcounts)
dim(readcounts)

coldata <- read.csv("coldata_Amselem.csv", header = T, sep = "," , row.names = 1)
head(coldata)
dim(coldata)

#verify if the rownames in coldata is the same in readcounts
all(rownames(coldata)%in%colnames(readcounts))
all(rownames(coldata)==colnames(readcounts))
      
dds <- DESeqDataSetFromMatrix(countData = readcounts,
                              colData = coldata,
                              design= ~ group)

dds

#pre-filtering, we can use when we would like to filter: keep only rows that have at least 10 reads total.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

dds$group <- relevel(dds$group, ref = "mated")

ddsSeq <- DESeq(dds, betaPrior = TRUE)
resultsNames(ddsSeq)

diapause <- results(ddsSeq, name =  "groupdiapause")
summary(diapause)
diapause <- subset(diapause, padj < 0.05)
summary(diapause)
sum(diapause$padj < 0.05, na.rm=TRUE)
              
        diapauseUP <- subset(diapause, log2FoldChange > 0)
        dim(diapauseUP)
        
        diapauseDOWN <- subset(diapause, log2FoldChange < 0)
        dim(diapauseDOWN)      

mated <- results(ddsSeq, name =  "groupmated")
summary(mated)
mated <- subset(mated, padj < 0.05)
summary(mated)
sum(mated$padj < 0.05, na.rm=TRUE)

        matedUP <- subset(mated, log2FoldChange > 0)
        dim(matedUP)
  
        matedDOWN <- subset(mated, log2FoldChange < 0)
        dim(matedDOWN)      
        
founder <- results(ddsSeq, name =  "groupfounder_diapause")
summary(founder)
founder <- subset(founder, padj < 0.05)
summary(founder)
sum(founder$padj < 0.05, na.rm=TRUE)

        founderUP <- subset(founder, log2FoldChange > 0)
        dim(founderUP)
 
        founderDOWN <- subset(founder, log2FoldChange < 0)
        dim(founderDOWN)      
        

#See overlap among the gene lists from Amsalem et al. 2015
DxM <- diapause[which(row.names(diapause) %in% row.names(mated)),]
head(DxM)
dim(DxM)

      DxMup <- diapauseUP[which(row.names(diapauseUP) %in% row.names(matedUP)),]
      head(DxMup)
      dim(DxMup)
      
      DxMdown <- diapauseDOWN[which(row.names(diapauseDOWN) %in% row.names(matedDOWN)),]
      head(DxMdown)
      dim(DxMdown)

DxF <- diapause[which(row.names(diapause) %in% row.names(founder)),]
head(DxF)
dim(DxF)
              
      DxFup <- diapauseUP[which(row.names(diapauseUP) %in% row.names(founderUP)),]
      head(DxFup)
      dim(DxFup)
                            
      DxFdown <- diapauseDOWN[which(row.names(diapauseDOWN) %in% row.names(founderDOWN)),]
      head(DxFdown)
      dim(DxFdown)

MxF <- diapause[which(row.names(mated) %in% row.names(founder)),]
head(MxF)
dim(MxF)

      MxFup <- matedUP[which(row.names(matedUP) %in% row.names(founderUP)),]
      head(MxFup)
      dim(MxFup)
      
      MxFdown <- matedDOWN[which(row.names(matedDOWN) %in% row.names(founderDOWN)),]
      head(MxFdown)
      dim(MxFdown)
      
      
#See overlap between Amsalem's gene lists and our study

#used results from DESeq2 for our results:
      #go to DESEq2code.R, Wald Test, get the results:
      #e.g.: age12Xdiet75 <- results(ddsSeq, name = "age12d.diet75.")
      
#12days & 75%
Bimp12d75<- read.csv("age12Xdiet75_genelist.csv", header = F) 
head(Bimp12d75)
dim(Bimp12d75)

dBimpBter <- Bimp12d75[which(Bimp12d75$V1 %in% row.names(diapause)),]
head(dBimpBter)
dBimpBter <- as.data.frame(dBimpBter)
dim(dBimpBter)
head(dBimpBter)
          
fBimpBter <- Bimp12d75[which(Bimp12d75$V1 %in% row.names(founder)),]
head(fBimpBter)
fBimpBter <- as.data.frame(fBimpBter)
dim(fBimpBter)
head(fBimpBter)

mBimpBter <- Bimp12d75[which(Bimp12d75$V1 %in% row.names(mated)),]
head(mBimpBter)
mBimpBter <- as.data.frame(mBimpBter)
dim(mBimpBter)
head(mBimpBter)
          

#12 days & 50%
Bimp12d50<- read.csv("age12Xdiet50_genelist.csv", header = T) 
head(Bimp12d50)
dim(Bimp12d50)
          
dBimpBter <- Bimp12d50[which(Bimp12d50$rownames.age12Xdiet50. %in% row.names(diapause)),]
head(dBimpBter)
dBimpBter <- as.data.frame(dBimpBter)
dim(dBimpBter)
head(dBimpBter)


fBimpBter <- Bimp12d50[which(Bimp12d50$rownames.age12Xdiet50. %in% row.names(founder)),]
head(fBimpBter)
fBimpBter <- as.data.frame(fBimpBter)
dim(fBimpBter)
head(fBimpBter)

mBimpBter <- Bimp12d50[which(Bimp12d50$rownames.age12Xdiet50. %in% row.names(mated)),]
head(mBimpBter)
mBimpBter <- as.data.frame(mBimpBter)
dim(mBimpBter)
head(mBimpBter)
          
#9 days and 75%
Bimp9d75<- read.csv("age9Xdiet75_genelist.csv", header = T) 
head(Bimp9d75)
dim(Bimp9d75)

dBimpBter <- Bimp9d75[which(Bimp9d75$rownames.age9Xdiet75. %in% row.names(diapause)),]
head(dBimpBter)
dBimpBter <- as.data.frame(dBimpBter)
dim(dBimpBter)
head(dBimpBter)

fBimpBter <- Bimp9d75[which(Bimp9d75$rownames.age9Xdiet75. %in% row.names(founder)),]
head(fBimpBter)
fBimpBter <- as.data.frame(fBimpBter)
dim(fBimpBter)
head(fBimpBter)

mBimpBter <- Bimp9d75[which(Bimp9d75$rownames.age9Xdiet75. %in% row.names(mated)),]
head(mBimpBter)
mBimpBter <- as.data.frame(mBimpBter)
dim(mBimpBter)
head(mBimpBter)

#######Lists Bter X Bimp foldchange

####UP

#12days & 75%
up12d75<- read.csv("age12d75UP_genelist.csv", header = T) 
head(up12d75)
dim(up12d75)

upOverD <- up12d75[which(up12d75$rownames.age12d75UP. %in% row.names(diapauseUP)),]
head(upOverD)
upOverD <- as.data.frame(upOverD)
dim(upOverD)
head(upOverD)

upOverF <- up12d75[which(up12d75$rownames.age12d75UP. %in% row.names(founderUP)),]
head(upOverF)
upOverF <- as.data.frame(upOverF)
dim(upOverF)
head(upOverF)

upOverM <- up12d75[which(up12d75$rownames.age12d75UP. %in% row.names(matedUP)),]
head(upOverM)
upOverM <- as.data.frame(upOverM)
dim(upOverM)
head(upOverM)

#12days & 50%
up12d50 <- read.csv("age12d50UP_genelist.csv", header = T)
head(up12d50)
dim(up12d50)

upOverD <- up12d50[which(up12d50$V1 %in% row.names(diapauseUP)),]
head(upOverD)
upOverD <- as.data.frame(upOverD)
dim(upOverD)
head(upOverD)

upOverF <- up12d50[which(up12d50$V1 %in% row.names(founderUP)),]
head(upOverF)
upOverF <- as.data.frame(upOverF)
dim(upOverF)
head(upOverF)

upOverM <- up12d50[which(up12d50$V1 %in% row.names(matedUP)),]
head(upOverM)
upOverM <- as.data.frame(upOverM)
dim(upOverM)
head(upOverM)

#9days & 75%
up9d75<- read.csv("age9d75UP_genelist.csv", header = T) 
head(up9d75)
dim(up9d75)

upOverD <- up9d75[which(up9d75$rownames.age9d75UP. %in% row.names(diapauseUP)),]
head(upOverD)
upOverD <- as.data.frame(upOverD)
dim(upOverD)
head(upOverD)

upOverF <- up9d75[which(up9d75$rownames.age9d75UP. %in% row.names(founderUP)),]
head(upOverF)
upOverF <- as.data.frame(upOverF)
dim(upOverF)
head(upOverF)

upOverM <- up9d75[which(up9d75$rownames.age9d75UP. %in% row.names(matedUP)),]
head(upOverM)
upOverM <- as.data.frame(upOverM)
dim(upOverM)
head(upOverM)


####DOWN

#12days & 75%
down12d75<- read.csv("age12d75DOWN_genelist.csv", header = T) 
head(down12d75)
dim(down12d75)

downOverD <- down12d75[which(down12d75$rownames.age12d75DOWN. %in% row.names(diapauseDOWN)),]
head(downOverD)
downOverD <- as.data.frame(downOverD)
dim(downOverD)
head(downOverD)

downOverF <- down12d75[which(down12d75$rownames.age12d75DOWN. %in% row.names(founderDOWN)),]
head(downOverF)
downOverF <- as.data.frame(downOverF)
dim(downOverF)
head(downOverF)


downOverM <- down12d75[which(down12d75$rownames.age12d75DOWN. %in% row.names(matedDOWN)),]
head(downOverM)
downOverM <- as.data.frame(downOverM)
dim(downOverM)
head(downOverM)

#12days & 50% and DOWN
down12d50<- read.csv("age12d50DOWN_genelist.csv", header = T) 
head(down12d50)
dim(down12d50)

downOverD <- down12d50[which(down12d50$rownames.age12d50DOWN. %in% row.names(diapauseDOWN)),]
head(downOverD)
downOverD <- as.data.frame(downOverD)
dim(downOverD)
head(downOverD)

downOverF <- down12d50[which(down12d50$rownames.age12d50DOWN. %in% row.names(founderDOWN)),]
head(downOverF)
downOverF <- as.data.frame(downOverF)
dim(downOverF)
head(downOverF)


downOverM <- down12d50[which(down12d50$rownames.age12d50DOWN. %in% row.names(matedDOWN)),]
head(downOverM)
downOverM <- as.data.frame(downOverM)
dim(downOverM)
head(downOverM)


#9days & 75% 
down9d75<- read.csv("age9d75DOWN_genelist.csv", header = T) 
head(down9d75)
dim(down9d75)

downOverD <- down9d75[which(down9d75$rownames.age9d75DOWN. %in% row.names(diapauseDOWN)),]
head(downOverD)
downOverD <- as.data.frame(downOverD)
dim(downOverD)
head(downOverD)

downOverF <- down9d75[which(down9d75$rownames.age9d75DOWN. %in% row.names(founderDOWN)),]
head(downOverF)
downOverF <- as.data.frame(downOverF)
dim(downOverF)
head(downOverF)


downOverM <- down9d75[which(down9d75$rownames.age9d75DOWN. %in% row.names(matedDOWN)),]
head(downOverM)
downOverM <- as.data.frame(downOverM)
dim(downOverM)
head(downOverM)
