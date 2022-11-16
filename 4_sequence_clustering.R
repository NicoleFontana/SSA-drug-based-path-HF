#####
# Code - Section 3.4: Sequence analysis #
#####
# This code shows the sequence clustering procedure. After calculating the dissimilarity matrix 
# using the Optimal Matching distances, where the replacement costs derive from the transition rates,
# the hierarchical clustering with k = 2, ..., 8 is applied. The best hierarchical clustering 
# partition is chosen using the quality metrics and is used to initialise the PAM algorithm. 
# Once the final partition is chosen, using the visualisation tools shown in "3_sequence_analysis.R", 
# the sequences in each cluster are analysed.

#Library
library(TraMineR)
library(cluster)
library(WeightedCluster)
library(RColorBrewer)

##Combined-sequences----
load("dataset/seq_combination.Rdata")
seq_combination <- seq_combination[which(seq_combination[,"COD_REG"] %in% patients_undersampling),]
class.labels <- 0:7
class.scode <- c("No drugs", "RAS", "AA", "BB", "RAS & AA", "RAS & BB", "BB & AA", "RAS & BB & AA")

#Create the sequence object
seq.comb <- seqdef(data = seq_combination, var = 2:dim(seq_combination)[2], alphabet=class.scode, 
                   labels=class.labels,
                   xtstep = 1, id = seq_combination[,"COD_REG"], right="DEL",
                   cpal = c(brewer.pal(8, "Paired")))

#Compute the substitution-cost matrix based on transition rates 
submat <- seqsubm(seq.comb, method = "TRATE")
#Compute the Optimal Matching distances  
dist.om1 <- seqdist(seq.comb, method = "OM", indel = 1, sm = submat, full.matrix =FALSE)

##Hierarchical agglomerative clustering
clusterward <- agnes(dist.om1, diss = TRUE, method = "ward")

#Cluster from 2 to 8
wardRange <- as.clustrange(clusterward, diss=dist.om1, ncluster=8)

#Plot quality metrics
statistic <- cbind(as.data.frame(wardRange$stats[,c("ASW",  "PBC", "HC")]),"cluster" = 2:8)
statistic <- melt(statistic, id.vars = "cluster")

ggplot(data=statistic, aes(x=cluster, y=value, group=variable, color=variable)) +
  geom_line()+
  geom_point()+
  scale_x_discrete(limits=2:8)+
  labs( y = "Value", x = "Cluster") + labs(color='Metric')+
  ggtitle("Cluster quality statistics")+theme_bw()+
  theme(plot.title = element_text(),plot.title.position = "plot")

#Analysing the clusters that have been created
cluster_2 <- factor(wardRange$clustering$cluster2, labels = paste("Cluster", 1:2))
cluster_5 <- factor(wardRange$clustering$cluster5, labels = paste("Cluster", 1:5))
cluster_8 <- factor(wardRange$clustering$cluster8, labels = paste("Cluster", 1:8))

#State distribution plot
seqdplot(seq.comb, cluster_2, border=NA,with.legend=FALSE)
seqdplot(seq.comb, cluster_5, border=NA, with.legend=FALSE,use.layout=TRUE,rows=1,cols=5)
seqdplot(seq.comb, cluster_8, border=NA, with.legend=FALSE,use.layout=TRUE,rows=2,cols=4)

#Sequence frequency plot
seqfplot(seq.comb, group=cluster_2, border=NA, with.legend=FALSE)
seqfplot(seq.comb, group=cluster_5, border=NA, with.legend=FALSE,use.layout=TRUE,rows=1,cols=5)
seqfplot(seq.comb, group=cluster_8, border=NA, with.legend=FALSE,use.layout=TRUE,rows=2,cols=4)

#Mean time
seqmtplot(seq.comb,group=cluster_2,with.legend=FALSE,axes=FALSE)
seqmtplot(seq.comb,group=cluster_5,with.legend=FALSE,use.layout=TRUE,rows=1,cols=5,axes=FALSE)
seqmtplot(seq.comb,group=cluster_8,with.legend=FALSE,use.layout=TRUE,rows=2,cols=4,axes=FALSE)

##Partitioning around medoids intialize with 5 and 8 cluster of hierarchical clustering
PAM_5 <- wcKMedoids(dist.om1, k=5, initialclust=cluster_5)
PAM_8 <- wcKMedoids(dist.om1, k=8, initialclust=cluster_8)

PAM_8$stats
PAM_5$stats

##We take the partition with 8 clusters of the hierarchical clustering
save(cluster_8, file ="dataset/HC8_comb.Rdata")

##Diuretics-sequences ----
load("dataset/seq_diu.Rdata")
class.labels <- c(0,1)
class.scode <- c("No drug", "Drug")

#Create sequence object
seq.diu <- seqdef(data = seq_diu, var = 2:dim(seq_diu)[2], alphabet=class.scode,
                  labels = class.scode, xtstep = 1, id = seq_diu[,"COD_REG"], cpal = c("firebrick2","green3"))

#Compute the substitution-cost matrix based on transition rates 
submat <- seqsubm(seq.diu, method = "TRATE")
#Compute the Optimal Matching distances  
dist.om1 <- seqdist(seq.diu, method = "OM", indel = 1, sm = submat, full.matrix =FALSE)

##Hierarchical agglomerative clustering
clusterward <- agnes(dist.om1, diss = TRUE, method = "ward")

#Cluster from 2 to 8
wardRange <- as.clustrange(clusterward, diss=dist.om1, ncluster=8)

#Plot quality metrics
statistic <- cbind(as.data.frame(wardRange$stats[,c("ASW",  "PBC", "HC")]),"cluster" = 2:8)
statistic <- melt(statistic, id.vars = "cluster")

ggplot(data=statistic, aes(x=cluster, y=value, group=variable, color=variable)) +
  geom_line()+
  geom_point()+
  scale_x_discrete(limits=2:8)+
  labs( y = "Value", x = "Cluster") + labs(color='Metric')+
  ggtitle("Cluster quality statistics")+theme_bw()+
  theme(plot.title = element_text(),plot.title.position = "plot")

cluster_2 <- factor(wardRange$clustering$cluster2, labels = paste("Cluster", 1:2))

seqdplot(seq.diu, cluster_2, border=NA, with.legend=FALSE)
seqfplot(seq.diu, group=cluster_2, border=NA, with.legend=FALSE)
seqmtplot(seq.diu, group=cluster_2,with.legend=FALSE)

##Partitioning around medoids intialize with 2 cluster of hierarchical clustering
PAM2_diu <- wcKMedoids(dist.om1, k=2, initialclust=cluster_2)
PAM2_diu$stats
#Partition from PAM method is better
PAM2_diu$clustering[which(PAM2_diu$clustering==21014)] = "Cluster 1"
PAM2_diu$clustering[which(PAM2_diu$clustering==5065)] = "Cluster 2"

#Plot
seqdplot(seq.diu, group=PAM2_diu$clustering, border=NA,cols=2,with.legend=FALSE)
seqfplot(seq.diu, group=PAM2_diu$clustering, border=NA, with.legend=FALSE,use.layout=TRUE,rows=1,cols=2)
seqmtplot(seq.diu, group=PAM2_diu$clustering,with.legend=FALSE,cols=2)

##We take the partition with 2 clusters of the PAM
save(PAM2_aag, file ="dataset/PAM2_diu.Rdata")

##Antithrombotics-sequences----
load("dataset/seq_aag.Rdata")
class.labels <- c(0,1)
class.scode <- c("No drug", "Drug")
#Create sequence object
seq.aag <- seqdef(data = seq_aag, var = 2:dim(seq_aag)[2], alphabet=class.scode,
                  labels = class.scode, xtstep = 1, 
                  id = seq_aag[,"COD_REG"], cpal = c("firebrick2","green3"))

#Compute the substitution-cost matrix based on transition rates 
submat <- seqsubm(Create, method = "TRATE")
#Compute the Optimal Matching distances  
dist.om1 <- seqdist(Create, method = "OM", indel = 1, sm = submat, full.matrix =FALSE)

##Hierarchical agglomerative clustering
clusterward <- agnes(dist.om1, diss = TRUE, method = "ward")

#Cluster from 2 to 8
wardRange <- as.clustrange(clusterward, diss=dist.om1, ncluster=8)

#Plot quality metrics
statistic <- cbind(as.data.frame(wardRange$stats[,c("ASW",  "PBC", "HC")]),"cluster" = 2:8)
statistic <- melt(statistic, id.vars = "cluster")

ggplot(data=statistic, aes(x=cluster, y=value, group=variable, color=variable)) +
  geom_line()+
  geom_point()+
  scale_x_discrete(limits=2:8)+
  labs(y = "Value", x = "Cluster") + labs(color='Metric')+
  ggtitle("Cluster quality statistics")+theme_bw()+
  theme(plot.title = element_text(),plot.title.position = "plot")

cluster_2 <- factor(wardRange$clustering$cluster2, labels = paste("Cluster", 1:2))

seqdplot(Create, cluster_2, border=NA, with.legend=FALSE)
seqfplot(Create, group=cluster_2, border=NA, with.legend=FALSE)
seqmtplot(Create, group=cluster_2,with.legend=FALSE)

##Partitioning around medoids intialize with 2 cluster of hierarchical clustering
PAM2_aag <- wcKMedoids(dist.om1, k=2, initialclust=cluster_2)
PAM2_aag$stats
#Partition from PAM method is better
PAM2_aag$clustering[which(PAM2_aag$clustering==7341)] = "Cluster 1"
PAM2_aag$clustering[which(PAM2_aag$clustering==9695)] = "Cluster 2"

#Plot
seqdplot(seq.aag, group=PAM2_aag$clustering, border=NA,cols=2,with.legend=FALSE)
seqfplot(seq.aag, group=PAM2_aag$clustering, border=NA, with.legend=FALSE,use.layout=TRUE,rows=1,cols=2)
seqmtplot(seq.aag, group=PAM2_aag$clustering,with.legend=FALSE,cols=2)

##We take the partition with 2 clusters of the PAM
save(PAM2_aag, file ="dataset/PAM2_aag.Rdata")

