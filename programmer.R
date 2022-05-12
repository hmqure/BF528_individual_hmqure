library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)
library(dplyr)

setwd("/projectnb/bf528/users/vangogh2022/project_1/samples")
files <- list.files()

# Read sample files and normalize together. Fit PLM
samples <- ReadAffy(filenames = files)
rma_set <- affy::rma(samples)
PLM_fit <- fitPLM(samples, background=TRUE, normalize=TRUE)


#Histograms for RLE and NUSE
hist(RLE(PLM_fit,type="stats")[1,], 
     main="Distribution of Median RLE Scores",
     xlab="Score", 
     ylab="Count", 
     col="blue")

hist(NUSE(PLM_fit,type="stats")[1,], 
     main="Distribution of Median NUSE Scores",
     xlab="Score", 
     ylab="Count", 
     col="brown")


#Read Metasamples file and get batch columns
expression = expression(rma_set)
metasamples <- read.csv("/project/bf528/project_1/doc/metasamples.csv")
batch <- metasamples$normalizationcombatbatch
model <- model.matrix(~normalizationcombatmod, samples=metasamples)

#Correction with ComBat
corrected <- ComBat(dat = expression, batch = batch, model = model)
transposed <- t(corrected)
scaled <- scale(transposed, center = TRUE, scale = TRUE)
result <- t(scaled)
compobj <- prcomp(result, center = FALSE, scale = FALSE)
pcs <- samples.frame(compobj$rotation)
v <- samples.frame(summary(compobj)$importance)

#PCA variance percentages
v_pc1 <- toString(v$PC1[2]*100)
v_pc2 <- toString(v$PC2[2]*100)

#Add cancer subtype column
subtype <- cbind(pcs, subtype=metasamples$SixSubtypesClassification)

#PCA plot
pca_plot <- ggplot(subtype,aes(x=PC2, y=PC1,color=subtype)) + geom_point()
pca_plot + ggtitle("PC1 vs PC2") + 
  ylab(paste("PC2 = ",v_pc2,"%")) + xlab(paste("PC1 = ",v_pc1,"%"))