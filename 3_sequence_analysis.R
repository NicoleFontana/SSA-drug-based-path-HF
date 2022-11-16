#####
# Code - Section 3.3: Sequence analysis #
#####
#This code shows the analysis of the constructed sequences. For each group of sequences 
# (combined, DIU and AAG), the state distribution plot, the sequence frequency plot, the mean 
# time spent in each state and the transition rates between states are shown. Entropy and 
# turbulence and their correlation are also calculated.

#Libraries
library(TraMineR)
library(cluster)
library(WeightedCluster)
library(corrplot)
library(diagram)
library(RColorBrewer)

#Import of the patients in the dataset undersampled (since the sequence analysis is performed only on these data)
#The computation of this dataset is performed in "2_sequence_construction.R"
load("dataset/patients_final_cohort.Rdata")

##1. Combined-sequences ----
load("dataset/seq_combination.Rdata")
seq_combination <- seq_combination[which(seq_combination[,"COD_REG"] %in% patients_final_cohort),]
class.labels <- 0:7
class.scode <- c("No drugs", "RAS", "AA", "BB", "RAS & AA", "RAS & BB", "BB & AA", "RAS & BB & AA")

#Create the sequence object
seq.comb <- seqdef(data = seq_combination, var = 2:dim(seq_combination)[2], alphabet=class.scode, 
                   labels=class.labels,
                   xtstep = 1, id = seq_combination[,"COD_REG"], right="DEL",
                   cpal = c(brewer.pal(8, "Paired")))

#State distribution plot
seqdplot(seq.comb, 
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,
         border=NA,
         with.legend=FALSE, xtlab=c(1:52),
         xlab="Time [week]",
         main = "State distribution plot",
         cpal = c(brewer.pal(8, "Paired")))

#Sequence frequency plot
seqfplot(seq.comb, 
         use.layout=FALSE,
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,xtlab=c(1:52),
         with.legend = FALSE,
         #legend.prop=0.1, 
         border = NA, 
         xlab="Time [week]",
         main = "Sequence frequency plot", 
         cpal = c(brewer.pal(8, "Paired")))

#Mean time spent in each state
seqmtplot(seq.comb, 
          font.main=2,font.lab=1, 
          cex.main=1.8,cex.lab=1.3,cex.axis=1.4,
          with.legend="FALSE",
          main = "Mean time spent in each state",
          cpal = c(brewer.pal(8, "Paired")),ylim=c(0,26))
seqmeant(seq.comb, serr = TRUE)

#Entropy and Turbulence  
Entropy <- seqient(seq.comb, norm=TRUE) #within sequence entropy
Turbulence <- seqST(seq.comb, norm=TRUE)
cor(Entropy, Turbulence)

#Transition rate
comb_rate <- seqtrate(seq.comb)
round(comb_rate, 2)
colnames(comb_rate) <- c("No drugs", "RAS", "AA", "BB", "RAS & AA", "RAS & BB", "BB & AA", "RAS & BB & AA")
rownames(comb_rate) <-  c("No drugs", "RAS", "AA", "BB", "RAS & AA", "RAS & BB", "BB & AA", "RAS & BB & AA")

par(fig=c(0,1,0,0.9), new=TRUE)
corrplot(comb_rate, method="number",is.corr = FALSE,tl.cex=0.8,
         tl.col = 'black',col = c(brewer.pal(9, "BuGn"))[c(4,6:9)],col.lim = c(0, 1))
par(fig=c(0,1,0.9,0.95), new=TRUE)
title(main = "Transition rates between states",font.main=2,cex.main=1.5,adj = 0,line=3)

##2. Diuretic-sequences ----
load("dataset/seq_diu.Rdata")
class.labels <- c(0,1)
class.scode <- c("No drug", "Drug")

#Create sequence object
seq.diu <- seqdef(data = seq_diu, var = 2:dim(seq_diu)[2], alphabet=class.scode,
                  labels = class.scode, xtstep = 1, id = seq_diu[,"COD_REG"], cpal = c("firebrick2","green3"))

#State distribution plot
seqdplot(seq.diu, 
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,
         border=NA,
         with.legend=FALSE, xtlab=c(1:52),
         xlab="Time [week]",
         main = "State distribution plot",
         cpal = c("firebrick2","green3"))

#Sequence frequency plot
seqfplot(seq.diu, 
         use.layout=FALSE,
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,xtlab=c(1:52),
         with.legend = FALSE,
         #legend.prop=0.1, 
         border = NA, 
         xlab="Time [week]",
         main = "Sequence frequency plot", 
         cpal = c("firebrick2","green3"))

#Mean time spent in each state
seqmtplot(seq.diu, 
          font.main=2,font.lab=1, 
          cex.main=1.2,cex.lab=1.2,cex.axis=1.2,
          with.legend="FALSE",
          main = "Mean time spent in each state",
          cpal = c("firebrick2","green3"))
seqmeant(seq.diu, serr = TRUE)

#Transition rate
diu_rate <- seqtrate(seq.diu)
round(diu_rate, 2)
colnames(diu_rate) <- c("No drug","Drug")
rownames(diu_rate) <- c("No drug","Drug")

plotmat(t(round(diu_rate,2)), pos = c(1, 1), curve = 0.4, name = colnames(diu_rate), 
        self.cex = 0.9, self.shiftx = c(0.12, -0.12, 0.12, 0.12),relsize=1,
        lwd = 1,box.lwd = 2, cex.txt = 1, box.type = "circle", box.prop = 1)
title(main = "Transition rates between the two states",font.main=2,cex.main=1.2)

#Entropy and Turbulence
Entropy <- seqient(seq.diu, norm=TRUE) 
Turbulence <- seqST(seq.diu, norm=TRUE)
cor(Entropy, Turbulence)

##3. Antithrombotic-sequences ----
load("dataset/seq_aag.Rdata")
class.labels <- c(0,1)
class.scode <- c("No drug", "Drug")
#Create sequence object
seq.aag <- seqdef(data = seq_aag, var = 2:dim(seq_aag)[2], alphabet=class.scode,
                  labels = class.scode, xtstep = 1, 
                  id = seq_aag[,"COD_REG"], cpal = c("firebrick2","green3"))

#State distribution plot
seqdplot(seq.aag, 
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,
         border=NA,
         with.legend=FALSE, xtlab=c(1:52),
         xlab="Time [week]",
         main = "State distribution plot",
         cpal = c("firebrick2","green3"))

#Sequence frequency plot
seqfplot(seq.aag, 
         use.layout=FALSE,
         font.main=2,font.lab=1, 
         cex.main=1.2,cex.lab=0.9,xtlab=c(1:52),
         with.legend = FALSE,
         #legend.prop=0.1, 
         border = NA, 
         xlab="Time [week]",
         main = "Sequence frequency plot", 
         cpal = c("firebrick2","green3"))

#Mean time spent in each state
seqmtplot(seq.aag, 
          font.main=2,font.lab=1, 
          cex.main=1.2,cex.lab=1.2,cex.axis=1.2,
          with.legend="FALSE",
          main = "Mean time spent in each state",
          cpal = c("firebrick2","green3"))
seqmeant(seq.aag, serr = TRUE)

#Transition rate
plotmat(t(round(aag_rate,2)), pos = c(1, 1), curve = 0.4, name = colnames(aag_rate), 
        self.cex = 0.9, self.shiftx = c(0.12, -0.12, 0.12, 0.12),relsize=1,
        lwd = 1,box.lwd = 2, cex.txt = 1, box.type = "circle", box.prop = 1)
title(main = "Transition rates between the two states",font.main=2,cex.main=1.2)

#Entropy and Turbulence
Entropy <- seqient(seq.aag, norm=TRUE) 
Turbulence <- seqST(seq.aag, norm=TRUE)
cor(Entropy, Turbulence)

