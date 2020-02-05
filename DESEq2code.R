########################FAT BODY TRANSCRIPTOME - DESeq2###############################

#reference to use
#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input

#install the all packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("vsn")

#Loading packages (Whenever you start R)
library("DESeq2")
library("ggplot2")
library("tibble")
library("plotly")
library("dplyr")
library("ggpubr")
library("vsn")
theme_set(theme_pubr())


############import the data.frames#########
    
  #obs. To use DESeqDataSetFromMatrix, the user should provide the counts matrix, the information about the samples (the columns of the count matrix) as a DataFrame or data.frame, and the design formula.
readcounts <- read.csv("ReadCounts.csv", header = T, row.names = 1)
head(readcounts)
dim(readcounts)

coldata <- read.csv("coldata.csv", header = T, row.names = 1, sep = ",")
head(coldata)
dim(coldata)
#for this analysis, I removed the column group   
coldata$group <- NULL 
head(coldata)
dim(coldata)


  #verify if the rownames in coldata is the same in readcounts
all(rownames(coldata)%in%colnames(readcounts))
all(rownames(coldata)==colnames(readcounts))




############Construct a DESeqDataSet############


  ####Wald Test###


#we can put a simple design, and after you can test if the log2 fold change attributable to given a condition is different based on another factor.

###design: ~ age + ~ diet + ~ Colony + ~ age:diet
dds <- DESeqDataSetFromMatrix(countData = readcounts,
                                      colData = coldata,
                                      design= ~ age + ~ diet + ~ Colony + ~ age:diet)

dds

#pre-filtering, we can use when we would like to filter: keep only rows that have at least 10 reads total.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

dds$age <- relevel(dds$age, ref = "03d")
dds$diet <- relevel(dds$diet, ref = "0%")

ddsSeq <- DESeq(dds)
resultsNames(ddsSeq)

#here you can get the paired results, for exemple queens 12 days old that received 75% sucrose diet
age12Xdiet75 <- results(ddsSeq, name = "age12d.diet75.")
summary(age12Xdiet75)
age12Xdiet75 <- subset(age12Xdiet75, padj < 0.05)
summary(age12Xdiet75)
dim(age12Xdiet75)

        #if you need to save the results
          #only ID names
        write.csv(as.data.frame(rownames(age12Xdiet75)),
                  file="age12Xdiet75_genelist.csv")
        
          #all results from table
        write.csv(as.data.frame(age12Xdiet75),
                  file="age12Xdiet75.csv")
        
  #if you want only the DEGs up-regulated
  age12d75UP <- subset(age12Xdiet75, log2FoldChange > 0)
  dim(age12d75UP)

    
          write.csv(as.data.frame(rownames(age12d75UP)),
                    file="age12d75UP_genelist.csv")
          
          write.csv(as.data.frame(age12d75UP),
                    file="age12d75UP.csv")

    #if you want only the DEGs down-regulated
    age12d75DOWN <- subset(age12Xdiet75, log2FoldChange < 0)
    dim(age12d75DOWN)      
    
          write.csv(as.data.frame(rownames(age12d75DOWN)),
                    file="age12d75DOWN_genelist.csv")
          
          write.csv(as.data.frame(age12d75DOWN),
                    file="age12d75DOWN.csv")
        
          ###NOTES: You can get any results, you need only switch "ageTwelve.dietTwenty_five" for what you have interest, e.g., "age12d.diet25."
      
          
              
          
  ###LRT###

          ###NOTES: We ran for each factor: interaction (age vs diet); age, diet and colony
          
ddsData <- DESeqDataSetFromMatrix(countData = readcounts,
                                    colData = coldata,
                                    design= ~ age + ~ diet + ~ Colony + age:diet)

ddsData

keep <- rowSums(counts(ddsData)) >= 10
ddsData <- ddsData[keep,]
ddsData
          
#interaction AGExDIET#


ddsDataAxD <- DESeq(ddsData, test = "LRT", reduced = ~ age:diet)

res05AxD <- results(ddsDataAxD, alpha=0.05); res05AxD
summary(res05AxD)
sum(res05AxD$padj < 0.05, na.rm=TRUE)
res05AxDOrdered <- res05AxD[order(res05AxD$pvalue),]; res05AxDOrdered

resSigAxD <- subset(res05AxDOrdered, padj < 0.05)
summary(resSigAxD)
dim(resSigAxD)

            # #if you need to save the results
            # write.csv(as.data.frame(rownames(resSigAxD)),
            #         file="DegAGExDIET.csv")

#AGE#

ddsDataA <- DESeq(ddsData, test = "LRT", reduced = ~ age)

res05A <- results(ddsDataA, alpha=0.05); res05A
summary(res05A) 
sum(res05A$padj < 0.05, na.rm=TRUE)
res05AOrdered <- res05A[order(res05A$pvalue),]; res05AOrdered

resSigA <- subset(res05AOrdered, padj < 0.05)
summary(resSigA)
dim(resSigA)

            # if you need to save the results
            # write.csv(as.data.frame(rownames(resSigA)),
            #               file="AGE.csv")

#DIET#

ddsDataD <- DESeq(ddsData, test = "LRT", reduced = ~ diet)

res05D <- results(ddsDataD, alpha=0.05); res05D
summary(res05D) 
sum(res05D$padj < 0.05, na.rm=TRUE)
res05DOrdered <- res05D[order(res05D$pvalue),]; res05DOrdered

resSigD <- subset(res05DOrdered, padj < 0.05)
summary(resSigD)
dim(resSigD)

          # if you need to save the results
          # write.csv(as.data.frame(rownames(resSigD)),
          #           file="DIET.csv")

#COLONY#

ddsDataC <- DESeq(ddsData, test = "LRT", reduced = ~ Colony)

res05C <- results(ddsDataC, alpha=0.05); res05C
summary(res05C) 
sum(res05C$padj < 0.05, na.rm=TRUE)
res05COrdered <- res05C[order(res05C$pvalue),]; res05COrdered

resSigC <- subset(res05COrdered, padj < 0.05)
summary(resSigC)
dim(resSigC)

      # if you need to save the results
      # write.csv(as.data.frame(rownames(resSigC)),
      #           file="COLONY.csv")



############Filtering############

        ###NOTES: The DEGs related to natal colony and the interaction of diet and age were subsequently excluded from the Age and Diet DEG lists for enrichment analysis in order to identify genes whose expression is most impacted by these two factors, independently

#Filter by colony

  #AGExDIET#
filterAGExDIET <- resSigAxD[-which(row.names(resSigAxD) %in% row.names(resSigC)),]
head(filterAGExDIET)
dim(filterAGExDIET)

            # if you need to save the results
            # write.csv(as.data.frame(filterAGExDIET),
            #           file="filterAGExDIET.csv")
            # write.csv(as.data.frame(rownames(filterAGExDIET)),
            #           file="filterAxD_genelist.csv")


        #give the number of "genes" with LFC < 0 (down)
        sum(filterAGExDIET$log2FoldChange < 0, na.rm=TRUE)
        AxDdown <- subset(filterAGExDIET, log2FoldChange < 0)
        dim(AxDdown)
        head(AxDdown)
        
                    # if you need to save the results
                    # write.csv(as.data.frame(AxDdown),
                    #           file="AxDdown.csv")
                    # 
                    # write.csv(as.data.frame(rownames(AxDdown)),
                    #           file="AxDdown_genelist.csv")
        
        
        #give the number of "genes" with LFC > 0 (up)
        sum(filterAGExDIET$log2FoldChange > 0, na.rm=TRUE)
        AxDup <- subset(filterAGExDIET, log2FoldChange > 0)
        dim(AxDup)
        head(AxDup)
        
                      # if you need to save the results
                      # write.csv(as.data.frame(AxDup),
                      #           file="AxDup.csv")
                      # 
                      # write.csv(as.data.frame(rownames(AxDup)),
                      #           file="AxDup_genelist.csv")


  #AGE#
filterAGE <- resSigA[-which(row.names(resSigA) %in% row.names(resSigC)),]
head(filterAGE)
dim(filterAGE)

            # if you need to save the results
            # write.csv(as.data.frame(filterAGE),
            #           file="filterAGE.csv")
            # 
            # write.csv(as.data.frame(rownames(filterAGE)),
            #           file="filterA_genelist.csv")


  #DIET#
filterDIET <- resSigD[-which(row.names(resSigD) %in% row.names(resSigC)),]
head(filterDIET)
dim(filterDIET)

            # if you need to save the results
            # write.csv(as.data.frame(filterDIET),
            #           file="filterDIET.csv")
            # 
            # write.csv(as.data.frame(rownames(filterDIET)),
            #           file="filterD_genelist.csv")


#Filter by INTERACTION

  #AGE#
filterAGE1 <- filterAGE[-which(row.names(filterAGE) %in% row.names(filterAGExDIET)),]
head(filterAGE1)
dim(filterAGE1)

            # if you need to save the results
            # write.csv(as.data.frame(filterAGE1),
            #           file="filterAGE1.csv")
            # 
            # write.csv(as.data.frame(rownames(filterAGE1)),
            #           file="filterA1_genelist.csv")


      #give the number of "genes" with LFC < 0 (down)
      sum(filterAGE1$log2FoldChange < 0, na.rm=TRUE)
      Adown <- subset(filterAGE1, log2FoldChange < 0)
      dim(Adown)
      head(Adown)
      
      # if you need to save the results
      # write.csv(as.data.frame(Adown),
      #           file="Adown.csv")
      # 
      # write.csv(as.data.frame(rownames(AxDdown)),
      #           file="Adown_genelist.csv")
      
      
      #give the number of "genes" with LFC > 0 (up)
      sum(filterAGE1$log2FoldChange > 0, na.rm=TRUE)
      Aup <- subset(filterAGE1, log2FoldChange > 0)
      dim(Aup)
      head(Aup)
      
      # if you need to save the results
      # write.csv(as.data.frame(Aup),
      #           file="Aup.csv")
      # 
      # write.csv(as.data.frame(rownames(Aup)),
      #           file="Aup_genelist.csv")



  #DIET#
filterDIET1 <- filterDIET[-which(row.names(filterDIET) %in% row.names(filterAGExDIET)),]
head(filterDIET1)
dim(filterDIET1)

            # if you need to save the results
            # write.csv(as.data.frame(filterDIET1),
            #           file="filterDIET1.csv")
            # 
            # write.csv(as.data.frame(rownames(filterDIET1)),
            #           file="filterD1_genelist.csv")


        #give the number of "genes" with LFC < 0 (down)
        sum(filterDIET1$log2FoldChange < 0, na.rm=TRUE)
        Ddown <- subset(filterDIET1, log2FoldChange < 0)
        dim(Ddown)
        head(Ddown)
        
        # if you need to save the results
        # write.csv(as.data.frame(Ddown),
        #           file="Ddown.csv")
        # 
        # write.csv(as.data.frame(rownames(Ddown)),
        #           file="Ddown_genelist.csv")
        
        
        #give the number of "genes" with LFC > 0 (up)
        sum(filterDIET1$log2FoldChange > 0, na.rm=TRUE)
        Dup <- subset(filterDIET1, log2FoldChange > 0)
        dim(Dup)
        head(Dup)
        
        # if you need to save the results
        # write.csv(as.data.frame(Dup),
        #           file="Dup.csv")
        # 
        # write.csv(as.data.frame(rownames(Dup)),
        #           file="Dup_genelist.csv")


#################Figures DESeq2 object#########################

# transform expression levels using the regularized log transformation
rld <-rlog(ddsData, blind = FALSE)
  #why do we use blind = FALSE? If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis.
# transform expression levels using the variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(ddsData, blind = FALSE)

#do these meanSD plots on data before DESeq
notAllZero <- (rowSums(counts(ddsData))>0)
meanSdPlot(log2(counts(ddsData,normalized=FALSE)[notAllZero,]+1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

## PCA - for rld-transformed data
# Simplest PCA, not easy to interpret
plotPCA(rld, intgroup=c("age","diet"))
# Gettin' fancier with the PCA
pca.dataRLD <- plotPCA(rld, intgroup=c("age", "diet"), returnData=TRUE)
percentVar <- round(100 * attr(pca.dataRLD, "percentVar"))
ggplot(pca.dataRLD, aes(PC1, PC2, color=age, shape=diet)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() +ggtitle("RLD-transformed data")
#format the graphic 
p<-ggplot(pca.dataRLD, aes(PC1, PC2, color=age, shape=diet)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
p + theme_bw()
p + theme_classic()
p + theme_light() +theme(legend.text=element_text(size=12)) + scale_color_manual(values=c("orchid", "hotpink1", "deeppink3"))


## PCA - for vsd-transformed data
# Simplest PCA for VSD
plotPCA(vsd, intgroup=c("age","diet"))
# Gettin' fancier with the PCA
pca.dataVSD <- plotPCA(vsd, intgroup=c("age", "diet"), returnData=TRUE)
percentVar <- round(100 * attr(pca.dataVSD, "percentVar"))
ggplot(pca.dataVSD, aes(PC1, PC2, color=age, shape=diet)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+ggtitle("VSD-transformed data")
#format the graphic 
q<-ggplot(pca.dataVSD, aes(PC1, PC2, color=age, shape=diet)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+ggtitle("VSD-transformed data")
q + theme_bw()
q + theme_light()
q + theme_light() +theme(legend.text=element_text(size=12)) + scale_color_manual(values=c("orchid", "hotpink1", "deeppink3"))


## Look at these scatterplots for transformation effects ##
# Scatterplots of transformed counts for samples 1:2
par(mfrow=c(1,3))
dds3 <- estimateSizeFactors(ddsColony)
plot(log2(counts(dds3, normalized=TRUE)[,1:2]+1), pch=16, cex=0.3, main="log2 transformed")
plot(assay(rld)[,1:2],pch=16, cex=0.3, main="RLD transformed")
plot(assay(vsd)[,1:2],pch=16, cex=0.3, main="VSD transformed")


##Heat maps
library("RColorBrewer")
library("pheatmap")
select <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:50]
top500 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:500]
top100 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:100]
top3000 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:3000]
top5468 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:5468]

nt <- normTransform(ddsData)
log2.norm.counts <- assay(nt)[select,]

df <- as.data.frame(colData(ddsData)[,c("age","diet")])
  df$age = factor(df$age, levels = c("03d", "09d", "12d"))
  ageCol <- c("thistle3", "hotpink1", "deeppink3")
  names(ageCol) <- levels(df$age)
  df$diet = factor(df$diet, levels = c("0%", "25%", "50%", "75%"))
  dietCol <- c("lavender", "lightblue1", "steelblue1", "dodgerblue4")
  names(dietCol) <- levels(df$diet)
  AnnColour <- list(
    age = ageCol,
    diet = dietCol)
  
distmat <- dist(t(top100))
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
cutoff.distance <- 3  
cols <- makeColorRampPalette(c("green", "purple",    # distances 0 to 3 colored from white to red
                               "green", "purple"), # distances 3 to max(distmat) colored from green to black
                             cutoff.distance / max(distmat),
                             30)
  
#changing show_rownames from FALSE to TRUE made gene names visible. More clutter, but easy to ID blue row.
pheatmap(log2.norm.counts, cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
#the following heatmaps are less informative and interesting when using dds instead
#of rld data. Reconsider input, and if transformed first
pheatmap(assay(rld)[select,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
pheatmap(assay(rld)[top500,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
pheatmap(assay(rld)[top100,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
#without indiv bee IDs listed
pheatmap(assay(rld)[top100,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, annotation_col = df, annotation_colors = AnnColour)
pheatmap(assay(rld)[top500,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, annotation_col = df, annotation_colors = AnnColour)
pheatmap(assay(rld)[top3000,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, annotation_col = df, annotation_colors = AnnColour)
pheatmap(assay(rld)[top5468,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, annotation_col = df, annotation_colors = AnnColour)


#################Post-hoc Analysis#########################

    ###Notes: We need to run DESeq2 for Age and Diet to do the figures

  #AGE#

    ###Notes: we use here same dataset above, but including samples from zero day, and for coldata, we used only age as factor.

readcountsA <- read.csv("readcounts_age.csv", header = T, row.names = 1)
head(readcountsA)
dim(readcountsA)

coldataA <- read.csv("coldata_age.csv", header = T, row.names = 1, sep = ",")
head(coldataA)
dim(coldataA)

all(rownames(coldataA)%in%colnames(readcountsA))
all(rownames(coldataA)==colnames(readcountsA))

ddsA <- DESeqDataSetFromMatrix(countData = readcountsA,
                               colData = coldataA,
                               design= ~ Age)

ddsA

keep <- rowSums(counts(ddsA)) >= 10
ddsA <- ddsA[keep,]
ddsA

ddsA$Age <- relevel(ddsA$Age, ref = "0d")

ddsSeqA <- DESeq(ddsA)
resultsNames(ddsSeqA)

resA <- results(ddsSeqA)

#PLots for specific genes

Bimp25880<- counts(ddsSeqA['BIMP25880',], normalized = TRUE)
m1 <- list(counts = as.numeric(Bimp25880), group = as.factor(coldataA$Age))
m1 <- as.tibble(m1)
m1$group <- factor(m1$group, levels=c("0d", "03d", "09d", "12d"))
q1 <- ggplot(m1, aes(group,counts)) + geom_boxplot(aes(fill = group))
q1 <- q1 + labs(x = "age", y = "Normalized Counts ", title = "1-acyl-sn-glycerol-3-phosphate acyltransferase")
q1 <- q1 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q1 <- q1 + scale_fill_brewer(palette="PuRd")
q1


Bimp17629 <- counts(ddsSeqA['BIMP17629',], normalized = TRUE)
m2 <- list(counts = as.numeric(Bimp17629), group = as.factor(coldataA$Age))
m2 <- as.tibble(m2)
m2$group <- factor(m2$group, levels=c("0d", "03d", "09d", "12d"))
q2 <- ggplot(m2, aes(group,counts)) + geom_boxplot(aes(fill = group))
q2 <- q2 + labs(x = "age", y = "Normalized Counts ", title = "Rab-related protein 4, Rab18")
q2 <- q2 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q2 <- q2 + scale_fill_brewer(palette="PuRd")
q2

Bimp10728 <- counts(ddsSeqA['BIMP10728',], normalized = TRUE)
m3 <- list(counts = as.numeric(Bimp10728), group = as.factor(coldataA$Age))
m3 <- as.tibble(m3)
m3$group <- factor(m3$group, levels=c("0d", "03d", "09d", "12d"))
q3 <- ggplot(m3, aes(group,counts)) + geom_boxplot(aes(fill = group))
q3 <- q3 + labs(x = "age", y = "Normalized Counts ", title = "Organic anion transporting polypeptide 33Eb")
q3 <- q3 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q3 <- q3 + scale_fill_brewer(palette="PuRd")
q3

Bimp22971 <- counts(ddsSeqA['BIMP22971',], normalized = TRUE)
m4 <- list(counts = as.numeric(Bimp22971), group = as.factor(coldataA$Age))
m4 <- as.tibble(m4)
m4$group <- factor(m4$group, levels=c("0d", "03d", "09d", "12d"))
q4 <- ggplot(m4, aes(group,counts)) + geom_boxplot(aes(fill = group))
q4 <- q4 + labs(x = "age", y = "Normalized Counts ", title = "Ecdysis Triggering Hormone receptor (ETHR)")
q4 <- q4 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q4 <- q4 + scale_fill_brewer(palette="PuRd")
q4

        # if you want to put together the figures
        # figure_age <- ggarrange(q1, q2, q3,
        #                         ncol = 2, nrow = 2 )
        # figure_age


  #DIET#
    
    ###Notes: we use here same dataset above (excluding zero day), but and for coldata, we used only diet as factor.

readcountsD <- read.csv("readcounts_diet.csv", header = T, row.names = 1)
head(readcountsD)
dim(readcountsD)

coldataD <- read.csv("coldata_diet.csv", header = T, row.names = 1, sep = ",")
head(coldataD)
dim(coldataD)

all(rownames(coldataD)%in%colnames(readcountsD))
all(rownames(coldataD)==colnames(readcountsD))


ddsD <- DESeqDataSetFromMatrix(countData = readcountsD,
                               colData = coldataD,
                               design= ~ Diet)

ddsD

keep <- rowSums(counts(ddsD)) >= 10
ddsD <- ddsD[keep,]
ddsD

ddsD$Diet <- relevel(ddsD$Diet, ref = "0%")

ddsSeqD <- DESeq(ddsD)
resultsNames(ddsSeqD)

resD <- results(ddsSeqD)

#PLots for specific genes

Bimp16901<- counts(ddsSeqD['BIMP16901',], normalized = TRUE)
n1 <- list(counts = as.numeric(Bimp16901), group = as.factor(coldataD$Diet))
n1 <- as.tibble(n1)
n1$group <- factor(n1$group, levels=c("0%", "25%", "50%", "75%"))
p1 <- ggplot(n1, aes(group,counts)) + geom_boxplot(aes(fill = group))
p1 <- p1 + labs(x = "diet", y = "Normalized Counts ", title = "woc")
p1 <- p1 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p1 <- p1 + scale_fill_brewer(palette="Blues")
p1

Bimp22357 <- counts(ddsSeqD['BIMP22357',], normalized = TRUE)
n2 <- list(counts = as.numeric(Bimp22357), group = as.factor(coldataD$Diet))
n2 <- as.tibble(n2)
n2$group <- factor(n2$group, levels=c("0%", "25%", "50%", "75%"))
p2 <- ggplot(n2, aes(group,counts)) + geom_boxplot(aes(fill = group))
p2 <- p2 + labs(x = "diet", y = "Normalized Counts ", title = "Succinyl-CoA:3-ketoacid-coenzyme A transferase")
p2 <- p2 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p2 <- p2 + scale_fill_brewer(palette="Blues")
p2


Bimp20482 <- counts(ddsSeqD['BIMP20482',], normalized = TRUE)
n3 <- list(counts = as.numeric(Bimp20482), group = as.factor(coldataD$Diet))
n3 <- as.tibble(n3)
n3$group <- factor(n3$group, levels=c("0%", "25%", "50%", "75%"))
p3 <- ggplot(n3, aes(group,counts)) + geom_boxplot(aes(fill = group))
p3 <- p3 + labs(x = "diet", y = "Normalized Counts ", title = "Plenty of SH3s")
p3 <- p3 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p3 <- p3 + scale_fill_brewer(palette="Blues")
p3

Bimp17199 <- counts(ddsSeqD['BIMP17199',], normalized = TRUE)
n4 <- list(counts = as.numeric(Bimp17199), group = as.factor(coldataD$Diet))
n4 <- as.tibble(n4)
n4$group <- factor(n4$group, levels=c("0%", "25%", "50%", "75%"))
p4 <- ggplot(n4, aes(group,counts)) + geom_boxplot(aes(fill = group))
p4 <- p4 + labs(x = "diet", y = "Normalized Counts ", title = "SA")
p4 <- p4 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p4 <- p4 + scale_fill_brewer(palette="Blues")
p4


Bimp14353 <- counts(ddsSeqD['BIMP14353',], normalized = TRUE)
n5 <- list(counts = as.numeric(Bimp14353), group = as.factor(coldataD$Diet))
n5 <- as.tibble(n5)
n5$group <- factor(n5$group, levels=c("0%", "25%", "50%", "75%"))
p5 <- ggplot(n5, aes(group,counts)) + geom_boxplot(aes(fill = group))
p5 <- p5 + labs(x = "diet", y = "Normalized Counts ", title = "mei-P26")
p5 <- p5 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p5 <- p5 + scale_fill_brewer(palette="Blues")
p5


Bimp25880 <- counts(ddsSeqD['BIMP25880',], normalized = TRUE)
n6 <- list(counts = as.numeric(Bimp25880), group = as.factor(coldataD$Diet))
n6 <- as.tibble(n6)
n6$group <- factor(n6$group, levels=c("0%", "25%", "50%", "75%"))
p6 <- ggplot(n6, aes(group,counts)) + geom_boxplot(aes(fill = group))
p6 <- p6 + labs(x = "diet", y = "Normalized Counts ", title = "DmelCG3812")
p6 <- p6 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p6 <- p6 + scale_fill_brewer(palette="Blues")
p6

Bimp18385 <- counts(ddsSeqD['BIMP18385',], normalized = TRUE)
n7 <- list(counts = as.numeric(Bimp18385), group = as.factor(coldataD$Diet))
n7 <- as.tibble(n7)
n7$group <- factor(n7$group, levels=c("0%", "25%", "50%", "75%"))
p7 <- ggplot(n7, aes(group,counts)) + geom_boxplot(aes(fill = group))
p7 <- p7 + labs(x = "diet", y = "Normalized Counts ", title = "DmelCG3812")
p7 <- p7 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p7 <- p7 + scale_fill_brewer(palette="Blues")
p7

Bimp24457 <- counts(ddsSeqD['BIMP24457',], normalized = TRUE)
n8 <- list(counts = as.numeric(Bimp24457), group = as.factor(coldataD$Diet))
n8 <- as.tibble(n8)
n8$group <- factor(n8$group, levels=c("0%", "25%", "50%", "75%"))
p8 <- ggplot(n8, aes(group,counts)) + geom_boxplot(aes(fill = group))
p8 <- p8 + labs(x = "diet", y = "Normalized Counts ", title = "Mco1")
p8 <- p8 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p8 <- p8 + scale_fill_brewer(palette="Blues")
p8


Bimp17594 <- counts(ddsSeqD['BIMP17594',], normalized = TRUE)
n9 <- list(counts = as.numeric(Bimp17594), group = as.factor(coldataD$Diet))
n9 <- as.tibble(n9)
n9$group <- factor(n9$group, levels=c("0%", "25%", "50%", "75%"))
p9 <- ggplot(n9, aes(group,counts)) + geom_boxplot(aes(fill = group))
p9 <- p9 + labs(x = "diet", y = "Normalized Counts ", title = "Dmel_CG17124")
p9 <- p9 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p9 <- p9 + scale_fill_brewer(palette="Blues")
p9


Bimp10728 <- counts(ddsSeqD['BIMP10728',], normalized = TRUE)
n10 <- list(counts = as.numeric(Bimp10728), group = as.factor(coldataD$Diet))
n10 <- as.tibble(n10)
n10$group <- factor(n10$group, levels=c("0%", "25%", "50%", "75%"))
p10 <- ggplot(n10, aes(group,counts)) + geom_boxplot(aes(fill = group))
p10 <- p10 + labs(x = "diet", y = "Normalized Counts ", title = "Organic anion transporting polypeptide 33Eb")
p10 <- p10 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p10 <- p10 + scale_fill_brewer(palette="Blues")
p10

Bimp23619 <- counts(ddsSeqD['BIMP23619',], normalized = TRUE)
n11 <- list(counts = as.numeric(Bimp23619), group = as.factor(coldataD$Diet))
n11 <- as.tibble(n11)
n11$group <- factor(n11$group, levels=c("0%", "25%", "50%", "75%"))
p11 <- ggplot(n11, aes(group,counts)) + geom_boxplot(aes(fill = group))
p11 <- p11 + labs(x = "diet", y = "Normalized Counts ", title = "CG14591")
p11 <- p11 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p11 <- p11 + scale_fill_brewer(palette="Blues")
p11

figure_diet <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11,
                         ncol = 3, nrow = 4 )
figure_diet


Bimp12768 <- counts(ddsSeqD['BIMP12768',], normalized = TRUE)
n12 <- list(counts = as.numeric(Bimp12768), group = as.factor(coldataD$Diet))
n12 <- as.tibble(n12)
n12$group <- factor(n12$group, levels=c("0%", "25%", "50%", "75%"))
p12 <- ggplot(n12, aes(group,counts)) + geom_boxplot(aes(fill = group))
p12 <- p12 + labs(x = "diet", y = "Normalized Counts ", title = "Allatostatin A receptor 1")
p12 <- p12 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p12 <- p12 + scale_fill_brewer(palette="Blues")
p12


Bimp20539 <- counts(ddsSeqD['BIMP20539',], normalized = TRUE)
n13 <- list(counts = as.numeric(Bimp20539), group = as.factor(coldataD$Diet))
n13 <- as.tibble(n13)
n13$group <- factor(n12$group, levels=c("0%", "25%", "50%", "75%"))
p13 <- ggplot(n13, aes(group,counts)) + geom_boxplot(aes(fill = group))
p13 <- p13 + labs(x = "diet", y = "Normalized Counts ", title = "Capability receptor (CapaR)")
p13 <- p13 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p13 <- p13 + scale_fill_brewer(palette="Blues")
p13

