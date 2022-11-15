#####
# Code - Section 2.3, 2.4 and 3.4.4: Final dataset #
#####
# This code produces the final dataset used for the predictive models. The construction of the 
# variables representing the patient's comorbidities is shown through the ICD-9-CM / CSS mapping. 
# These variables are used to calculate the Multisource Comorbidity Score for each patient. The 
# information about the clustering of sequences is added to the patient. For each patient, the 
# cluster variable will indicate belonging to the combined-sequence partition; DIU and AAG will 
# indicate whether the patient is a low-adopter/ adopter of the two drugs.

setwd("~/OneDrive - Politecnico di Milano/HF Regione Lombardia/Lavoro_nicole")

#Library
library(dplyr)

#Set seed for the dataset undersampling
set.seed(26091998)

#Load the pre-processed dataset from file "1_pre_processing.R".
load("data_pulito.Rdata") 

#Load the dataset constructed for the sequence computation "2_sequence_construction"
#(which contains updated information about hospitalisation and drug purchase)
load("dataset/data_for_sequences.Rdata") 
patients <- unique(data_for_sequences$COD_REG)

#Select, for each patients, the row that represents his reference hospitalization
data <- data %>% filter(COD_REG %in% patients & hosp == 1)

data <- select(data,c("COD_REG","death","timeOUT", "gender","age", "CCS_principal_diagnosis",
                      "metastatic","chf","dementia","renal","wtloss" ,"hemiplegia","alcohol","tumor",
                      "arrhythmia","pulmonarydz","coagulopathy","compdiabetes","anemia","electrolytes","liver",
                      "pvd","psychosis","pulmcirc","hivaids","hypertension",
                      "ICD", "SHOCK", "CABG", "PTCA"))

#Merge the information from the two datasets
data_for_sequences <- data_for_sequences %>% group_by(COD_REG) %>% arrange(hosp, .by_group = TRUE)
tot_hosp_d <- select(data_for_sequences %>% filter(row_number()==1),c("COD_REG","tot_hosp"))  

data_for_sequences <- data_for_sequences %>% group_by(COD_REG) %>% arrange(pharm, .by_group = TRUE)
tot_pharm_d <- select(data_for_sequences %>% filter(row_number()==1),c("COD_REG","tot_pharm"))  

data <- left_join(data, as.data.frame(tot_hosp_d) ,by="COD_REG")
data <- left_join(data, as.data.frame(tot_pharm_d) ,by="COD_REG")
data <- left_join(data, select(data_for_sequences %>% group_by(COD_REG) %>% filter(event_type =="hospitalization") %>% filter(row_number()==1),c(COD_REG),n_others_pharma), by="COD_REG")

##1. Create the variables which represent the comorbidities for the calculation of the MCS ----
##Using the table in appendix A, I check the CCS codes assigned to each reference 
# hospitalization and check 1 if there is a correspondence with the diagnoses for the MCS
codici <-unique(data$CCS_principal_diagnosis)#Verifico se ci sono dei codici
dim(CCS[unique(as.numeric(CCS[,1])) %in% c(1,79,38,40,44,58,156,158,5,80,81,83,95,3,202,205,211,114,115,117,248,109,54,83,139,9,96,213,58,48,100,101),])

#Separate the CCS variable, as it is made up of both the code and the diagnosis
CCS <- NULL
for (i in 1:dim(data)[1]){CCS = rbind(CCS,as.data.frame(str_split(data[i,"CCS_principal_diagnosis"],"_", n=2,simplify=TRUE)))}
data <- select(data.frame(data), -CCS_principal_diagnosis)
data <- cbind(data, "CCS" = unlist(CCS$V1), "Diagnosis" = unlist(CCS$V2))

#Create a matrix with 0 for each comorbidities
new_comorbidity <- matrix(0,dim(data)[1],16, byrow=F)
colnames(new_comorbidity) <- c("COD_REG","Tuberculosis","Parkinson", "Lymphoma","Kidney",
                               "Other_neurological","Rheumatoid", "Vascular","Cerebrovascular","Gout","Epilepsy",
                               "Ulcer", "Valvular", "Obesity","Hypothyroidism","Infarction")
new_comorbidity[,"COD_REG"] <- data[,"COD_REG"]
#Mapping the CCS in the diagnosis
for (i in 1:dim(data)[1])
{
  code = as.numeric(data[i, "CCS"])
  if(code == 1 ){new_comorbidity[i,"Tuberculosis"] = 1}
  if(code == 79){new_comorbidity[i,"Parkinson"] = 1}
  if(code %in% c(38, 40, 44, 58)){new_comorbidity[i,"Lymphoma"] = 1}
  if(code %in% c(156, 158)){new_comorbidity[i,"Kidney"] = 1}
  if(code %in% c(5, 79, 80, 81, 83, 95)){new_comorbidity[i,"Other_neurological"] = 1}
  if(code %in% c(3,202,205, 211)){new_comorbidity[i,"Rheumatoid"] = 1}
  if(code %in% c(114,115,117,248)){new_comorbidity[i,"Vascular"] = 1}
  if(code == 109){new_comorbidity[i,"Cerebrovascular"] = 1}
  if(code == 54){new_comorbidity[i,"Gout"] = 1}
  if(code == 83){new_comorbidity[i,"Epilepsy"] = 1}
  if(code == 139){new_comorbidity[i,"Ulcer"] = 1}
  if(code %in% c(9, 96,213)){new_comorbidity[i,"Valvular"] = 1}
  if(code == 58){new_comorbidity[i,"Obesity"] = 1}
  if(code == 48){new_comorbidity[i,"Hypothyroidism"] = 1}
  if(code %in% c(100,101)){new_comorbidity[i,"Infarction"] = 1}
  
}
comorbidities_tab <- left_join(data[,c(1,6:25)], as.data.frame(new_comorbidity), by="COD_REG")
#save(comorbidities_tab, file = "dataset/comorbidities_tab")

##2. Compute the multisource comorbidity score ---- 
comorbidities_tab <- comorbidities_tab %>%
           mutate(MCS = 18*metastatic +
           11*alcohol +
           10*tumor +
           10*Tuberculosis +
           8*psychosis +
           8*liver +
           #    6*anagrafe$FAR_ANXIETY +
           6*wtloss +
           6*dementia +
           #   5*anagrafe$FAR_MALIGNANCIES +
           5*Parkinson +
           5*Lymphoma +
           5*hemiplegia +
           5*coagulopathy +
           4*electrolytes+
           4*Kidney +
           #   4*anagrafe$DIALISI_RENALE +
           4*chf+
           3*Other_neurological +
           3*Rheumatoid +
           3*anemia +
           3*Cerebrovascular +
           2*compdiabetes +
           2*Vascular +
           2*Gout +
           2*Epilepsy +
           2*pulmonarydz +
           2*Ulcer +
           1*Infarction +
           1*pvd +
           1*Valvular +
           1*arrhythmia +
           1*Obesity +
           1*Hypothyroidism
  )

##3. Add necessary covariates to the dataset ----
#a. Add the MCS to each patients
data <- left_join(data, comorbidities_tab[,c("COD_REG","MCS_cont")], by="COD_REG")
data <- select(data, c(COD_REG,death,timeOUT,gender,eta,ICD,SHOCK,CABG,PTCA,tot_hosp,
                       tot_pharm,n_others_pharma,MCS_cont,tot_pharm,tot_hosp))

##Select patients in the undersampled dataset
load("dataset/patients_final_cohort.Rdata")
data <- data %>% filter(COD_REG %in% patients_final_cohort)

##b. Categorised the MCS 
data <- data %>% mutate(MCS=cut(data$MCS_cont, breaks=c(0,4, 9,42),include.lowest=TRUE))
data$MCS<- factor(data$MCS , levels=c("[0,4]","(4,9]","(9,42]"),
                  labels= c("Good","Intermediate","Poor"))

##c. Recompute the total procedures
data <- data %>% mutate(tot_procedures = as.numeric(CABG) +as.numeric( PTCA) + as.numeric(ICD) +as.numeric( SHOCK)-4)

##d. Add the LOS, the total number of day in hospital in the first year
load("dataset/seq_hosp_days.Rdata") #Import the sequence of hospitalisation of each days
tot_days_hosp <- seq_hosp_mortality_days[,"COD_REG"] %>% bind_cols(rowSums(seq_hosp_mortality_days[,2:366]))
colnames(tot_days_hosp) <- c("COD_REG","LOS")
data <- left_join(data, tot_days_hosp, by="COD_REG")
data$LOS <- as.numeric(data$LOS)

##e. Add the clustering from the combined-sequence
load("dataset/HC8_comb.Rdata")
load("dataset/patients_cluster_order.Rdata")
cluster_8<- factor(wardRange$clustering$cluster8, labels = paste("Cluster", 1:8),
                   levels=c(1,2,3,4,7,8,5,6))
cluster_data <- as.data.frame(cbind("COD_REG" = patients_cluster_order,"cluster" = cluster_8))
data <- left_join(as.data.frame(data), cluster_data, by = "COD_REG")
data$cluster <- as.factor(data$cluster)

##f. Add the DIU and AAG clustering
##DIU
load("dataset/PAM2_diu.Rdata")
load("dataset/seq_diu.Rdata")
PAM2_diu$clustering[which(PAM2_diu$clustering==21014)] = "Low_Adopters"
PAM2_diu$clustering[which(PAM2_diu$clustering==5065)] = "Adopters"
cluster_data <- as.data.frame(cbind("COD_REG" = as.numeric(seq_diu_mortality[,"COD_REG"]),
                                    "DIU" = PAM2_diu$clustering))
cluster_data$COD_REG <- as.numeric(cluster_data$COD_REG)
data <- left_join(as.data.frame(data), cluster_data, by = "COD_REG")
data$DIU[is.na(data$DIU)] <- "Low_Adopters"

##AAG
load("dataset/PAM2_aag.Rdata")
load("dataset/seq_aag.Rdata")
PAM2_aag$clustering[which(PAM2_aag$clustering==7341)] = "Low_Adopters"
PAM2_aag$clustering[which(PAM2_aag$clustering==9695)] = "Adopters"
cluster_data <- as.data.frame(cbind("COD_REG" = as.numeric(seq_aag_mortality[,"COD_REG"]),
                                    "AAG" = PAM2_aag$clustering))
cluster_data$COD_REG <- as.numeric(cluster_data$COD_REG)
data <- left_join(as.data.frame(data), cluster_data, by = "COD_REG")
data$AAG[is.na(data$AAG)] <- "Low_Adopters"

data$DIU <- factor(data$DIU, levels=c("Low_Adopters","Adopters"))
data$AAG <- factor(data$AAG, levels=c("Low_Adopters","Adopters"))

##4. Save the final dataset ----
save(data, file = "dataset/final_dataset")




