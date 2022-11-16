#####
# Code - Section 3.2: Sequence construction #
#####
# Code for section 3.2 of the manuscript shows how the sequences necessary for the analyses 
# (Combined-sequence, DIU-sequence and AAG-sequence) are constructed. The first phase consists 
# of preparing the longitudinal data to select pharmaceutical purchases in the first year of 
# observation. Then a stratified undersampling is done with respect to the survival outcome to 
# reduce the dataset. The "create_seq.R" function creates binary sequences for the five interest 
# drugs; by aggregating the sequences of RAS, BB and AA the combined sequence is obtained. It is 
# also created using "create_seq_hosp_days.R", the variable that counts the days of hospitalisation 
# in the first year of observation.

#Library
library(dplyr)

#Set seed for the dataset undersampling
set.seed(26091998)

#Import the dataset obtain from the file "1_pre_processing.R" 
load("data_pulito.Rdata")
data_rif <- data

#Import the function to construct the sequence
source("auxiliary_function/create_seq.R")
source("auxiliary_function/create_seq_hosp_days.R")

###1.Dataset preparation ----
##Exclude the reference hospitalisation of each patients
data <- data %>% filter(event_type == "drug purchase" | (event_type == "hospitalization" & hosp>1))

## Delete all those events that start after 365 days:
#purchase carried out after 365 from the reference date
#hospitalisation starting after 365 after the reference date
data <- select(data %>% filter(date_of_prescription - data_rif_ev <= 365 | 
                                 date_of_admission - data_rif_ev <= 365),-c(tot_hosp,tot_pharm,tot_events))

##Truncate the coverage days for drugs that exceed the first year of observation
data <- data %>% mutate(end_of_prescription = if_else(date_of_prescription + qt_pharma - data_rif_ev >= 365, 
                                                      data_rif_ev+365, date_of_prescription + qt_pharma))

##Truncate the hospitalisation that exceed the first year of observation
data <- data %>% mutate(date_of_discharge = if_else(date_of_discharge - data_rif_ev >= 365,  
                                                    data_rif_ev + 365, date_of_discharge))

data <- data %>% filter((event_type =="hospitalization" & date_of_admission < date_of_discharge)|
                          (event_type =="drug purchase" & date_of_prescription < end_of_prescription))

##Date check
data <- data %>% mutate(date_of_admission = if_else(date_of_admission<data_rif_ev, data_rif_ev,date_of_admission))

##Update the total number of hospitalisation and drug purchase
count_event <- data.frame(data %>% filter(event_type=="hospitalization") 
                          %>% group_by(COD_REG) %>% count(COD_REG,name="tot_hosp"))
data <- left_join(data, count_event, by="COD_REG")

##Take patients who have purchased one of the following drugs at least once: RAS, BB and AA
data$class_pharma <- ifelse(data$class_pharma %in% c("ACE","ARB"),"RAS",data$class_pharma)
ras <- select(data %>% filter(event_type=="drug purchase" & class_pharma == "RAS" & is.na(qt_pharma)==FALSE),"COD_REG")
bb <- select(data %>% filter(event_type=="drug purchase" & class_pharma == "BB" & is.na(qt_pharma)==FALSE),"COD_REG")
aa <- select(data %>% filter(event_type=="drug purchase" & class_pharma == "AA" & is.na(qt_pharma)==FALSE),"COD_REG")
patients_3 <- unique(c(ras$COD_REG, bb$COD_REG, aa$COD_REG))
data <- data %>% filter(COD_REG %in% patients_3)

##Include the first reference hospitalisation to those patients who haven't others hospitalisation
pat_drugs <- unique(select(data %>% filter(event_type=="drug purchase"), COD_REG)) 
pat_hosp <- unique(select(data %>% filter(event_type=="hospitalization"), COD_REG)) 
pat_loss <-unique(pat_drugs[!(pat_drugs$COD_REG %in% pat_hosp$COD_REG),"COD_REG"]) 
data_one_hosp <- data_rif[which(data_rif$COD_REG %in% pat_loss & data_rif$hosp == 1),]
data_one_hosp$tot_hosp <- 0 
data_one_hosp <- select(data_one_hosp, - tot_events)
data <- data %>% bind_rows(data_one_hosp)

##Update the total number of other pharma (excluding RAS, BB, AA, DIU and AAG) taken by the patient
count_pharm <- data.frame(data %>% group_by(COD_REG,.drop=FALSE) 
                          %>% filter(!is.na(class_pharma) & !(class_pharma %in% c("RAS","AA","BB","DIU","AAG"))) 
                          %>% distinct(class_pharma) %>% count(name="n_others_pharma"))

data <- left_join(select(data,-n_class_pharma), count_pharm, by="COD_REG")

##Update the total number of drugs purchase for each patient
count_event <- data.frame(data %>% filter(event_type=="drug purchase") 
                          %>% group_by(COD_REG) %>% count(COD_REG,name="tot_pharm"))
data_for_sequences <- left_join(select(data,-tot_pharm), count_event, by="COD_REG")

#Save this dataset with update information about hospitalisation and drug purchase
save(data_for_sequences, file ="dataset/data_for_sequences")

###2. Dataset undersampling, saving the patients in the final cohort ----
data <- stratified(data_for_sequences, c("death"), 0.35)
save(data$COD_REG,file="dataset/patients_final_cohort")

###3.Sequences RAS, BB, AA, DIU and AAG----
data_pharm <- data %>% select("COD_REG","event_type","data_rif_ev","date_of_prescription",
                              "end_of_prescription","qt_pharma","class_pharma","pharm")
##RAS
data_ras <- data_pharm %>% filter(event_type=="drug purchase" & class_pharma %in% c("ACE","ARB") & is.na(qt_pharma)==FALSE)
seq_ras<- create_seq(data_ras)
##BB
data_bb <- data_pharm %>% filter(event_type=="drug purchase" & class_pharma %in% c("BB") & is.na(qt_pharma)==FALSE)
seq_bb <- create_seq(data_bb)
###AA
data_aa <- data_pharm %>% filter(event_type=="drug purchase" & class_pharma %in% c("AA") & is.na(qt_pharma)==FALSE)
seq_aa <- create_seq(data_aa)
###DIU
data_diu <- data_1 %>% filter(event_type=="drug purchase" & class_pharma %in% c("DIU") & is.na(qt_pharma)==FALSE)
seq_diu <- create_seq(data_diu)
###AAG
data_aag <- data_1 %>% filter(event_type=="drug purchase" & class_pharma %in% c("AAG") & is.na(qt_pharma)==FALSE)
seq_aag <- create_seq(data_aag)

save(seq_diu, file="dataset/seq_diu.Rdata")
save(seq_aag, file="dataset/seq_aag.Rdata")

###4.Combined-sequences RAS & BB & AA ----
classe_ras <- unique(seq_ras[,"COD_REG"])
classe_bb <- unique(seq_bb[,"COD_REG"])
classe_aa <- unique(seq_aa[,"COD_REG"])

patients <- unique(c(classe_ras,classe_aa,classe_bb))
seq_combination <- matrix(NA,length(patients), 53, byrow=F)
nomi = NULL
for (i in 1:52){nomi=c(nomi, paste("week_",i,sep=""))}
colnames(seq_combination) <- c("COD_REG",nomi)
seq_combination[,1] <- patients

for (i in patients)
{
  for (j in 2:53)
  {
    if (i %in% classe_ras)
    {
      if(seq_ace[which(seq_ace[,"COD_REG"]==i), j] == 1)
        seq_combination[which(seq_combination[,"COD_REG"]==i), j] = 
          ifelse(is.na(seq_combination[which(seq_combination[,"COD_REG"]==i), j]),"RAS",paste(seq_combination[which(seq_combination[,"COD_REG"]==i),j],"& RAS"))
    }
    
    if (i %in% classe_bb)
    {
      if(seq_bb[which(seq_bb[,"COD_REG"]==i), j] == 1)
        seq_combination[which(seq_combination[,"COD_REG"]==i), j] = 
          ifelse(is.na(seq_combination[which(seq_combination[,"COD_REG"]==i), j]),"BB",paste(seq_combination[which(seq_combination[,"COD_REG"]==i),j],"& BB"))
    }
    
    if (i %in% classe_aa)
    {
      if(seq_aa[which(seq_aa[,"COD_REG"]==i), j] == 1)
        seq_combination[which(seq_combination[,"COD_REG"]==i), j] = 
          ifelse(is.na(seq_combination[which(seq_combination[,"COD_REG"]==i), j]),"AA",paste(seq_combination[which(seq_combination[,"COD_REG"]==i),j],"& AA"))
    }
  }
  
}
seq_combination[is.na(seq_combination)] <- "No drugs"
save(seq_combination, file = "dataset/seq_combination.Rdata")

###5.LOS in the first year of hospitalisation ----
data_hosp <- data %>% filter(event_type == "hospitalization")
data_hosp <- select(data_hosp, c("COD_REG","event_type","data_rif_ev","date_of_admission",
                                 "date_of_discharge","LOS","hosp","tot_hosp"))
seq_hosp_days <- create_seq_hosp_days(data_hosp)
save(seq_hosp_days, file="dataset/seq_hosp_days.Rdata")

