knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(tidyverse)
library(tidyr)

#Gathering csv data
mtx_path = "/projectnb/bf528/users/vangogh2022/project_1/batch_corrected_data.csv"
mtx = read_delim(mtx_path, delim = ",")
str(mtx)

colnames(mtx) = sub("_.*", "", colnames(mtx))
colnames(mtx)[1] = "Gene"

# filter based on log2(15)
min_exp= log(15,base = 2)
exp = .2 
mtxfilter1 =mtx %>% filter( (rowSums(.[-1] > min_exp)/(ncol(.)-1)) >exp)

#Get median variance
variance = function(x){
  medvar = sum((x - mean(x))^2)/(length(x) - 1)
  return(medvar)
}

# variance of genes and median
listv = apply(mtxfilter1[,-1], MARGIN = 1, FUN = variance)
medlistv = median(listv)

# minimum t-test
stat = qchisq(p = 0.01, df = ncol(mtx)-2,lower.tail = F)
df = ncol(mtx) - 2
chi = df*listv/(medlistv)

# second filter
mkfilter2 = chi>stat
pchi = 1 - pchisq(chi,ncol(mtx)-2)
mtxfilter2 = mtxfilter1[mkfilter2,]

#third filter
listvfilter3 = listv[mkfilter2]
mu = apply(mtxfilter2[,-1], MARGIN = 1, FUN = mean)
mkfilter3 = (sqrt(listvfilter3)/mu) > 0.186
mtxfilter3 = mtxfilter2[mkfilter3,]

# cluster samples passing filters
cluster = mtxfilter3

# Transpose matrix
TcolName = mtxfilter3$Gene
TrowName = names(mtxfilter3)[-1]
framecluster = as.data.frame(t(cluster[,-1]))
names(framecluster) = TcolName
row.names(framecluster) = TrowName

# Euclidean distance clustering
euclid <- dist(framecluster, method = 'euclidean')
avg_clust <- hclust(euclid, method = 'average')
plot(avg_clust,cex=0.8)
avg_tree <- cutree(avg_clust, k = 2)
plot(avg_clust, cex=0.8)
suppressPackageStartupMessages(library(dendextend))
dendro <- as.dendrogram(avg_clust)
dendro_color <- color_branches(dendro, h = 90)
plot(dendro_color)
rect.hclust(avg_clust , k = 2, border = 2:6)

# assign C3 samples to cluster
meta_data = read_delim("/project/bf528/project_1/doc/proj_metadata.csv", delim = ",")
geo_accession = rownames(framecluster)
rowName = as.data.frame(geo_accession)
classification = dplyr::left_join(x = rowName, y = meta_data[c("cit-coloncancermolecularsubtype","geo_accession")], by= "geo_accession")
C3_type = classification$`cit-coloncancermolecularsubtype`=="C3"
C3_type= C3_type[1:nrow(framecluster)]
data_heatmap <- as.matrix(t(framecluster))

# Create with Heatmap
heatmap(x = data_heatmap,ColSideColors = if_else(C3_type,"Red","blue"),cexCol = 0.8, scale = "row")
legend(x="right", legend=c("1:Lowest", "2","3","4", "5:Highest"),fill=heat.colors(5),cex = 0.8,title = "Expression")
legend(x="topright", legend=c("C3", "Other"),fill=c("Red","Blue"),cex = 0.8,title = "Colon Cancer Subtype")

#Welch test
welch = data_frame(Gene = c(),t_statistic = c(),p_value = c(), p_adj = c())
for( i in 1:ncol(framecluster)){
  t_test = t.test(framecluster[C3_type,i], framecluster[-C3_type,i], alternative = "two.sided")
  p_adj = p.adjust(t_test$p.value, method = "bonferroni", n = ncol(framecluster))
  temp = data_frame(Gene = names(framecluster)[i], t_statistic = t_test$statistic,p_value = t_test$p.value, p_adj = p_adj)
  welch = rbind(welch,temp)
}


