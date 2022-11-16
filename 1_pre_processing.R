#####
# Code - Section 1.3: Data pre-processing #
#####
#This code presents the pre-processing of the administrative dataset of the Lomabardy region.
#The covariates necessary for the analyses are created in this dataset. In particular, 
#each purchase of drugs will be assigned the relative DDD and pharmaceutical class.

#Library
library(dplyr)
library(stringr)
library(readxl)

#Import dataset
load("HF.Rdata")
data <- HF.2006_2012

#Create the column that identifies the event: hospitalization (41) or pharmaceutical purchase (30)
data <- select(data %>% mutate
               (event_type = if_else(data["tipo_prest"] == 41, "hospitalization", "drug purchase"),.after=labelOUT),
               -tipo_prest)

##Count the number of hospitalisations and drug purchases of each patient
events_type <- data %>% group_by(COD_REG,event_type)%>%count(COD_REG,event_type)
tot_hosp <- events_type[which(events_type$event_type == "hospitalization"),c(1,3)]
colnames(tot_hosp) <- c("COD_REG","tot_hosp")
tot_pharm <- events_type[which(events_type$event_type == "drug purchase"),c(1,3)]
colnames(tot_pharm) <- c("COD_REG","tot_pharm")
data <- left_join(data, tot_hosp, by="COD_REG")
data <- left_join(data, tot_pharm, by="COD_REG")

##Select patients with at least one pharmaceutical purchase
data <- data %>% filter(tot_hosp > 0 & tot_pharm > 0)

##Create the column that counts the total number of rows for each patient
data <- data %>% mutate(tot_events = tot_hosp + tot_pharm)

##Reorder each patient's rows from oldest to newest with respect to the performance date.
#If two records have the same performance date, put the hospitalisation first
data <- data %>% group_by(COD_REG) %>% arrange(data_prest,class_prest,.by_group = TRUE)

##Create follow-up variable
data <- select(data %>% mutate
               (timeOUT = as.numeric(difftime(data_studio_out,data_rif_ev,units = "days"))),
               -c(data_studio_out))

##Keep surviving patients for at least one year
data <- data %>% filter(timeOUT>365)

##Create the covariates that identify the hospitalisation diagnosis and the ATC codes of pharmaceutical purchases
data <- select(data %>% mutate 
               (CCS_principal_diagnosis = if_else(event_type == "hospitalization", class_prest, NULL),         
                 ATC_code = if_else(event_type == "drug purchase", class_prest, NULL),.after = event_type),
               -class_prest)

##Create the covariates that describe the lenght of stay in hospital for each hospitalisation,
#and the covarage days for each pharmaceutical purchase
data <- select(data %>% mutate
               (LOS = if_else(event_type == "hospitalization", qt_prest_Sum, NULL),
                 qt_pharma = if_else(event_type == "drug purchase", qt_prest_Sum, NULL),.after = ATC_code),
               -qt_prest_Sum)

##Create the covariates that indicate the date of admission and discharge from the hospitalisation
#and the date of purchase of a drug
data <- select(data %>% mutate
               (date_of_discharge = if_else(event_type == "hospitalization", data_prest, NULL),
                 date_of_admission = if_else(event_type == "hospitalization", date_of_discharge - LOS, NULL),
                 date_of_purchase = if_else(event_type == "drug purchase", data_prest, NULL),.after = qt_pharma),
               -data_prest)

##Cleaning information
data[which(data[,"event_type"] == "hospitalization"),"ASL_RESIDENZA"] <- NA
data[which(data[,"event_type"] == "drug purchase"),"strutt_id"] <- NA

##Create the variable that represents the patient's status (alive / dead) at the end of the study
data <- select(data %>% mutate (death = if_else(labelOUT == "DECEDUTO", 1, 0)),
               -labelOUT)

##Create the index of the hospitalisation and pharmaceutical purchase records
data <- data %>% group_by(COD_REG,event_type) %>% mutate(index=seq_along(tot_pharm),.after = COD_REG)
data <- data %>% mutate(hosp = ifelse(event_type=="hospitalization",index,NA),.after = COD_REG)
data <- select(data %>% mutate(pharm = ifelse(event_type=="drug purchase",index,NA),.after = hosp),-index)

##Assign the corresponding drug class to each ATC code
class_ATC <- read_excel("~/OneDrive - Politecnico di Milano/HF Regione Lombardia/Lavoro_nicole/dataset/DDD&class_pharm.xlsx",sheet="Foglio2",col_names = TRUE )
class_pharma = rep(NA,dim(data)[1])
data <- data.frame(data) #devo metterlo senno str_which ritorna dei warnings
for (i in 1:dim(class_ATC)[1])
{
  code = class_ATC[i,which(!is.na(class_ATC[i,]))]
  for (j in 2:length(code)){
    class_pharma[str_which(data[,"ATC_code"],code[[j]])]=code[[1]]
  }
}  
data <- data %>% mutate(class_pharma)

##Calculate the number of drug classes taken by each patient
count_pharm <- data.frame(data %>% filter(is.na(class_pharma) == FALSE) %>%  group_by(COD_REG) %>% distinct(class_pharma) %>% count(name="n_class_pharma"))
data <- left_join(data, count_pharm, by="COD_REG")


##Calculate the number of comorbidities at each hospitalisation event
comorbidity <- c("metastatic", "chf", "dementia", "renal", "wtloss", "hemiplegia", "alcohol", "tumor", "arrhythmia",
                 "pulmonarydz", "coagulopathy", "compdiabetes", "anemia", "electrolytes", "liver", "pvd", "psychosis",           
                 "pulmcirc", "hivaids", "hypertension")
data <- data %>% mutate(tot_comorbidity = if_else(event_type == "hospitalization", rowSums(data[,comorbidity]), NULL))
for (i in 2:dim(data)[1])
{
  if(data[i,"event_type"]=="drug purchase") {data[i,"n_comorbidity"]=data[i-1,"n_comorbidity"]}
}

##Calculate the total number of procedures the patient has undergone
procedures = c("ICD", "SHOCK","CABG", "PTCA")
data <- data %>% mutate(tot_procedures = if_else(event_type == "hospitalization", rowSums(data[,procedures]), NULL))
for (i in 2:dim(data)[1])
{
  if(data[i,"event_type"]=="drug purchase") {data[i,"tot_procedures"]=data[i-1,"tot_procedures"]}
}

##Assign the defined daily dose (or COMBO) at each drug purchase
DDD_value = read_excel("~/OneDrive - Politecnico di Milano/HF Regione Lombardia/Lavoro_nicole/dataset/DDD&class_pharm.xlsx",sheet="Foglio1")
DDD=rep(NA,dim(data)[1])
COMBO=rep(NA,dim(data)[1])

for (i in 1:dim(data)[1])
{ 
  if(data[i,"event_type"] == "drug purchase")
  {
    value = DDD_value[which(data[i,"ATC_code"] == DDD_value[,1]),2][[1]]
    if(! is.na(value)) #Controllo che non sia NA
    {
      if(value == "NO_val") COMBO[i]=1
      else 
      {
        DDD[i]=value
        COMBO[i]=0
      }
    }
  }
}

data <- data %>% mutate(DDD,COMBO)

##Change some covariates name
data <- data %>% rename(gender = SESSO,age = eta_Min)

##Reorder the covariates
reorder <- c(
  # global variables
  "COD_REG","data_rif_ev","gender","timeOUT","death", "event_type","age","tot_events","tot_pharm","tot_hosp","n_class_pharma",
  # events
  "hosp","date_of_admission","date_of_discharge","CCS_principal_diagnosis","LOS","strutt_id", # hosp block
  #comorbidities
  "tot_comorbidity","metastatic","chf","dementia","renal","wtloss" ,"hemiplegia","alcohol","tumor",
  "arrhythmia","pulmonarydz","coagulopathy","compdiabetes","anemia","electrolytes","liver",
  "pvd","psychosis","pulmcirc","hivaids","hypertension",
  #procedures
  "tot_procedures","ICD", "SHOCK", "CABG", "PTCA",
  #pharma block
  "pharm","date_of_purchase","ATC_code","class_pharma","DDD","COMBO","qt_pharma","ASL_RESIDENZA")

##Final dataset
data <- data %>% select(all_of(reorder))
save(data, file = "dataset/data_pulito.Rdata")

