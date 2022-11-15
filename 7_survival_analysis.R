#####
# Code - Section 3.6: Survival analysis #
#####
#This code presents the survival analysis using Cox's regression model. The comparison between the 
#two significant clusters of the combined-sequences reported in 3.6.1 of the manuscript is shown.

setwd("~/OneDrive - Politecnico di Milano/HF Regione Lombardia/Lavoro_nicole")

#Libraries
library(dplyr)
library(survminer)
library(tidyverse)
library(survival)
library(RColorBrewer)
library(tidyr)
library(TraMineR)

#Import the function to compute adherence
source("create_coverage_days.R")

#Load the final dataset
load("dataset/final_dataset")

##1. Stratified log-rank test ----
survdiff(Surv(timeOUT, death) ~ age_cat, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ gender, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ LOS_cat, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ MCS, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ procedures_cat, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ CABG, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ PTCA, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ ICD, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ SHOCK, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ cluster_comb, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ DIU, data = data) # Significant
survdiff(Surv(timeOUT, death) ~ AAG, data = data) # Significant

##2. Proportional-Hazard Cox model ----
mod.cox<- coxph(Surv(timeOUT, death==1)  ~ age + gender + LOS + MCS + tot_procedures + 
                           cluster + DIU + AAG, data = data)
summary(mod.cox)

# Verify Hazard assumption
test.ph <- cox.zph(mod.cox)
test.ph

# Hazard ratio 
ggforest(mod.cox, 
         data=data,
         main = "",
         fontsize=1.4) + ggtitle("Hazard Ratio (95% CI)")+theme_bw()+
         theme(plot.title = element_text(size=26,hjust=0.5,vjust = -7,face = "bold"),
         plot.title.position = "plot")

##3. Kaplan-Meier curves between groups ----
###AGE
k=4
binary_df <- with(data,data.frame(age = c(45,65,80,95),
                                  #age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M", k), 
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender","MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)

ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           legend.labs=paste("Age=",c(45,65,80,95),sep = ""),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           xlim=c(365, 2500),
           #title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ") 

###Gender 
k=2
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = c("M","F"),
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender", "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster+ DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           xlim=c(365, 2500),
           legend.labs=paste("Gender=",c("Male","Female"),sep = ""),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")

###LOS
k=3
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M",k),
                                  LOS=c(0,15,200),
                                  #LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender", "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           xlim=c(365, 2500),
           legend.labs=paste("LOS=",c(0,15,200), " days",sep = ""),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")

###MCS
k=3
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M",k),
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS=c("Good","Intermediate","Poor"),
                                  #MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender", "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           xlim=c(365, 2500),
           legend.labs=paste("MCS=",c("Good","Intermediate","Poor"),sep = ""),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")

###Procedures
k=3
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M",k),
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= c(0,1,2),
                                  MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender", "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           xlim=c(365, 2500),
           legend.labs=c("No procedures","One procedure","Two procedures"),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")

###Cluster 
k=7
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M",k),
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS = rep(levels(MCS)[1],k),
                                  #cluster_comb=rep(levels(cluster_comb)[1],k),
                                  cluster=c(levels(cluster)[1],
                                            levels(cluster)[2],
                                            levels(cluster)[3],
                                            levels(cluster)[4],
                                            levels(cluster)[5],
                                            #  levels(cluster_comb)[6],
                                            levels(cluster)[7],
                                            levels(cluster)[8]),
                                  AAG=rep("Low_Adopters",k),
                                  DIU=rep("Low_Adopters",k)
))

names <- c("gender", "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures +
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           xlim=c(365, 2500),
           #pval = TRUE,
           conf.int = TRUE, 
           legend.labs=paste("Cluster=",c(1:5,7,8),sep = ""),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2"))[c(1:6,8)], 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")

###DIU & AAG 
k=3
binary_df <- with(data,data.frame(age = rep(mean(age, na.rm = TRUE),k),
                                  gender = rep("M",k),
                                  LOS = rep(mean(LOS, na.rm = TRUE),k),
                                  tot_procedures= rep(0,k),
                                  MCS = rep(levels(MCS)[1],k),
                                  cluster=rep(levels(cluster)[1],k),
                                  DIU=c("Low_Adopters","Adopters","Low_Adopters"),
                                  AAG=c("Low_Adopters","Low_Adopters","Adopters")
                                  #AAG=rep("Low_Adopters",k),
                                  #DIU=rep("Low_Adopters",k)
))

names <- c("gender",  "MCS","AAG","DIU")
binary_df[,names] <- lapply(binary_df[,names] , factor)                
mod.cox_cluster <- coxph(Surv(timeOUT, death==1)  ~ age + gender +
                           + LOS + MCS + tot_procedures+
                           cluster + DIU + AAG, 
                         data = data)

fit <- survfit(mod.cox_cluster, newdata = binary_df)
ggsurvplot(fit, 
           data=binary_df, 
           size = 0.5,
           censor.size=0.5,
           #pval = TRUE,
           conf.int = TRUE, 
           xlim=c(365, 2500),
           legend.labs=c("Low adopters of DIU and AAG","DIU Adopters","AAG Adopters"),
           break.time.by=365,
           palette = c(brewer.pal(8, "Set2")), 
           #palette =  c(brewer.pal(8, "Paired")),
           xlab = "Time [days]",
           # title="Survival probability plot",
           ggtheme = theme_bw()+theme(plot.title = element_text(size=16),
                                      plot.title.position = "plot",
                                      legend.title=element_blank()),
           font.x = c(12),
           font.y = c(12),
           font.legend = c(12),
           legend.title=" ")


##4. Compute the adherence for patient in cluster 1 and 2 ----
cohort <- data %>% filter(cluster %in% c(1,2))

#Import dataset with all the drug purchase used for the computation of sequences
load("dataset/data_for_sequence")
data <- data_for_sequence %>% filter(event_type == "drug purchase" & is.na(qt_pharma)==FALSE & COD_REG %in% cohort$COD_REG)
data <- select(data, c("COD_REG","data_rif_ev","date_of_prescription", 
                       "end_of_prescription","qt_pharma","class_pharma","pharm", "tot_pharm"))

#Adherence to RAS
data_ras <- select(data %>% filter(class_pharma %in% "RAS"),-tot_pharm)
data_ras <- data_ras %>% group_by(COD_REG) %>% arrange(date_of_prescription,.by_group = TRUE) %>% as.data.frame()
coverage_RAS <- create_aderenza(data_ras)
aderence <- NULL
aderence <- cbind("COD_REG" = coverage_AA[,"COD_REG"],
                  "covarage" = rowSums(coverage_RAS[,2:366])) %>% as.data.frame()
aderence <- aderence %>% mutate(aderance_RAS = ifelse(covarage/365 >= 0.8, 1, 0))
cohort <- left_join(cohort, aderence[,c(1,3)], by="COD_REG")

#Adherence to AA
data_aa <- select(data %>% filter(class_pharma %in% "AA"),-tot_pharm)
data_aa <- data_aa %>% group_by(COD_REG) %>% arrange(date_of_prescription,.by_group = TRUE) %>% as.data.frame()
coverage_AA <- create_aderenza(data_aa)
aderence <- NULL
aderence <- cbind("COD_REG" = coverage_AA[,"COD_REG"],
                  "covarage" = rowSums(coverage[,2:366])) %>% as.data.frame()
aderence <- aderence %>% mutate(aderance_AA = ifelse(covarage/365 >= 0.8, 1, 0))
cohort <- left_join(cohort, aderence[,c(1,3)], by="COD_REG")

#Adherence to BB
data_bb <- select(data %>% filter(class_pharma %in% "BB"),-tot_pharm)
data_bb <- data_bb %>% group_by(COD_REG) %>% arrange(date_of_prescription,.by_group = TRUE) %>% as.data.frame()
coverage_BB <- create_aderenza(data_bb)
aderence <- NULL
aderence <- cbind("COD_REG" = coverage_BB[,"COD_REG"],
                  "covarage" = rowSums(coverage[,2:366])) %>% as.data.frame()
aderence <- aderence %>% mutate(aderance_BB = ifelse(covarage/365 >= 0.8, 1, 0))
cohort <- left_join(cohort, aderence[,c(1,3)], by="COD_REG")

#Adherence to DIU
data_diu <- select(data %>% filter(class_pharma %in% "DIU"),-tot_pharm)
data_diu <- data_diu %>% group_by(COD_REG) %>% arrange(date_of_prescription,.by_group = TRUE) %>% as.data.frame()
coverage_DIU <- create_aderenza(data_diu)
coverage_AA <- create_aderenza(data_aa)
aderence <- NULL
aderence <- cbind("COD_REG" = coverage_DIU[,"COD_REG"],
                  "covarage" = rowSums(coverage[,2:366])) %>% as.data.frame()
aderence <- aderence %>% mutate(aderance_DIU = ifelse(covarage/365 >= 0.8, 1, 0))
cohort <- left_join(cohort, aderence[,c(1,3)], by="COD_REG")

#Adherence to AAG
data_aag <- select(data %>% filter(class_pharma %in% "AAG"),-tot_pharm)
data_aag <- data_aag %>% group_by(COD_REG) %>% arrange(date_of_prescription,.by_group = TRUE) %>% as.data.frame()
coverage_AAG <- create_aderenza(data_aag)
coverage_AA <- create_aderenza(data_aa)
aderence <- NULL
aderence <- cbind("COD_REG" = coverage_AAG[,"COD_REG"],
                  "covarage" = rowSums(coverage[,2:366])) %>% as.data.frame()
aderence <- aderence %>% mutate(aderance_AAG = ifelse(covarage/365 >= 0.8, 1, 0))
cohort <- left_join(cohort, aderence[,c(1,3)], by="COD_REG")

cohort[is.na(cohort)] <- 0

##5. Comparison between patients in cluster 1 and 2 ----
cohort_cl1 <- cohort %>% filter(cluster == 1)
cohort_cl2 <- cohort %>% filter(cluster == 2)

##Test between groups - categorical variables
kruskal.test(gender ~ cluster, data = cohort)
kruskal.test(aderance_DIU ~ cluster, data = cohort)
kruskal.test(aderance_AAG ~ cluster, data = cohort)
kruskal.test(aderance_RAS ~ cluster, data = cohort)
kruskal.test(aderance_AA ~ cluster, data = cohort)
kruskal.test(MCS_cat ~ cluster, data = cohort)
kruskal.test(death ~ cluster, data = cohort)

##Test between groups - continous variables
wilcox.test(age ~ cluster, data = cohort)
wilcox.test(LOS ~ cluster, data = cohort)
wilcox.test(tot_hosp ~ cluster, data = cohort)
wilcox.test(tot_pharm ~ cluster, data = cohort)
wilcox.test(timeOUT ~ cluster, data = cohort)

##Comparison of patient characteristics in the two clusters
cohort$procedures_cat <- as.factor(cohort$procedures_cat )
summary(cohort_cl1[,c("age","LOS")])
summary(cohort_cl2[,c("age","LOS")])
c(summary(as.factor(cohort_cl1$gender)),summary(as.factor(cohort_cl2$gender)))
c(summary(as.factor(cohort_cl1$death)),summary(as.factor(cohort_cl2$death)))
c(summary(as.factor(cohort_cl1$tot_procedures)),summary(as.factor(cohort_cl2$tot_procedures)))
c(summary(as.factor(cohort_cl1$MCS)),summary(as.factor(cohort_cl2$MCS)))


##Plot the proportion of patients adhering to the 5 drug classes, in the two clusters
coll <- brewer.pal(8, "Set3")[c(4,7)]
p1 <-ggplot(data=cohort, aes(x = as.factor(cluster), fill = as.factor(aderance_RAS))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=coll) +
  labs( y = "Proportion of patients", x = "Cluster") + labs(fill='Adherence')+
  ggtitle("RAS") + theme_bw()+
  theme(plot.title = element_text())
p2<-ggplot(data=cohort, aes(x = as.factor(cluster), fill = as.factor(aderance_BB))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=coll) +
  labs( y = "Proportion of patients", x = "Cluster") + labs(fill='Adherence')+
  ggtitle("BB") + theme_bw()+
  theme(plot.title = element_text())
p3<-ggplot(data=cohort, aes(x = as.factor(cluster), fill = as.factor(aderance_DIU))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=coll) +
  labs( y = "Proportion of patients", x = "Cluster") + labs(fill='Adherence')+
  ggtitle("DIU") + theme_bw()+
  theme(plot.title = element_text())
p4<-ggplot(data=cohort, aes(x = as.factor(cluster), fill = as.factor(aderance_AAG))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=coll,labels = c("No","Yes")) +
  labs( y = "Proportion of patients", x = "Cluster") + labs(fill='Adherence')+
  ggtitle("AAG") + theme_bw()+
  theme(plot.title = element_text()) 

plot <-ggarrange(p1,p2,p3,p4, common.legend = TRUE,legend="bottom",
                 ncol = 4, nrow = 1,legend.grob=get_legend(p4, position = "bottom"))
annotate_figure(plot, top = text_grob("Proportion of patients' adherence to", size = 14))


ggplot(data=cohort, aes(x = as.factor(cluster), fill = as.factor(aderance_BB))) + 
  geom_bar(position = "fill") +
  scale_fill_manual(values=c("firebrick2","green3")) +
  labs( y = "Proportion of patients", x = "Cluster") + labs(fill='Adherence')+
  ggtitle("Patients adhering to BB") + 
  theme(text = element_text(size = 15)) 

##Comparison of sequences entropy and turbulence between groups
load("dataset/seq_combination_mortality")
class.labels <- 0:7
class.scode <- c("No drugs", "RAS", "AA", "BB", "RAS & AA", "RAS & BB", "BB & AA", "RAS & BB & AA")

#Select sequence of patients belongin to cluster 1 and 2
seq_cluster1 <- seq_combination_mortality[(seq_combination_mortality[,"COD_REG"] %in% cohort_cl1$COD_REG),]
seq_cluster1 <- seqdef(data = seq_cluster1, var = 2:dim(seq_cluster1)[2], alphabet=class.scode, 
                       labels=class.labels,
                       xtstep = 1, id = seq_cluster1[,"COD_REG"], right="DEL",
                       cpal = c(brewer.pal(8, "Paired")))
seq_cluster2 <- seq_combination_mortality[(seq_combination_mortality[,"COD_REG"] %in% cohort_cl2$COD_REG),]
seq_cluster2 <- seqdef(data = seq_cluster2, var = 2:dim(seq_cluster2)[2], alphabet=class.scode, 
                       labels=class.labels,
                       xtstep = 1, id = seq_cluster2[,"COD_REG"], right="DEL",
                       cpal = c(brewer.pal(8, "Paired")))
Entropy_1 <- seqient(seq_cluster1, norm=TRUE) 
Entropy_2 <- seqient(seq_cluster2, norm=TRUE) 
Turbulence_1 <- seqST(seq_cluster1, norm = TRUE)
Turbulence_2 <- seqST(seq_cluster2, norm = TRUE)

