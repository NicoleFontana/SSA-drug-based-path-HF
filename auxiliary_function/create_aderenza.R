create_aderenza <- function (data_class){
  #Aggiorno il numero totale di eventi dopo questi vincoli
  count_event <- data.frame(data_class %>% group_by(COD_REG) %>% count(COD_REG,name="tot_pharm_class"))
  
  data_c <- left_join(data_class, count_event, by="COD_REG") %>% as.data.frame()
  
  #Creo la matrice in cui tengo conto dei giorni di copertura (effettivi) e gli overlaps
  PRESCRIPTION <- matrix(0,dim(data_c)[1],3, byrow=F)
  colnames(PRESCRIPTION) <- c("COD_REG","days","overlaps")
  PRESCRIPTION[,1] <-data_c[,"COD_REG"]
  
  #Calcolo per ogni paziente il numero di giorni coperti dalla prescrizione
  codici <- data_c[,"COD_REG"]
  for (i in unique(codici))
  {
    for (j in 1:(count_event[count_event$COD_REG==i,"tot_pharm_class"]))
    {
      index = which(data_c[,"COD_REG"]==i)[j]
      #Calcolo i giorni di copertura, controllando le sovrapposizioni
      #Se la differenza è negativa vuol dire che l'evento i è totalmente incluso nell'evento i-1
      duration = as.numeric(data_c[index,"end_of_prescription"] - data_c[index,"date_of_prescription"])+1
      PRESCRIPTION[index,"days"] <- ifelse(duration>0, duration, 0)
      if(PRESCRIPTION[index,"days"] == 0)
      {
        data_c[index,"end_of_prescription"] <- data_c[index-1,"end_of_prescription"]
      }
      #Controllo se la data di fine dell'evento attuale sia maggiore della data di inizio del prossimo evento
      #se non è l'ultimo evento del paziente
      if(j != (count_event[count_event$COD_REG==i,2]))
      {
        #Verifico che la data della nuova prescrizione sia dopo la fine della vecchia prescrizione
        if(data_c[index,"end_of_prescription"] > data_c[index+1,"date_of_prescription"])
        {
          #Se c'è sovrapposizione, salvo i giorni di l'overlaps per nuova prescizione
          PRESCRIPTION[index+1,"overlaps"] = as.numeric(data_c[index,"end_of_prescription"] - data_c[index+1,"date_of_prescription"])+1
          #Shifto la data di inizio della nuova prescrizione al giorno dopo della fine della prescrizione precedente
          data_c[index+1,"date_of_prescription"] <- data_c[index,"end_of_prescription"]
        }
      }
      
      
    }
    
  }
  
  data_c <- cbind(data_c,as.data.frame(PRESCRIPTION)[,2:3])
  data_c$days <- ifelse(data_c$days>0,data_c$days -1, 0)
  data_c <- data_c %>% select(-qt_pharma)
  
  #Elimino le prescrizioni con 0 giorni di somministrazione e aggiorno il numero totale di eventi
  data_c <- select(data_c %>% filter(days>0), -tot_pharm_class)
  
  #Aggiorno il numero totale di eventi dopo questi vincoli e i codici
  count_event <- data.frame(data_c %>% group_by(COD_REG) %>% count(COD_REG,name="tot_pharm_class"))
  
  data_c <- left_join(data_c, count_event, by="COD_REG")
  data_c <- data_c %>% filter(tot_pharm_class>0)
  
  codici <- data_c[,"COD_REG"]
  
  #1. Parto costruendo un vettore di 365, dove ogni elemento rappresenta un giorno
  somm_day <- matrix(0,length(unique(codici)),365, byrow=F)
  rownames(somm_day)=unique(codici)
  nomi = NULL
  for (i in 1:365){nomi=c(nomi, paste("day_",i,sep=""))}
  colnames(somm_day) <- nomi
  
  for (i in unique(codici))
  {
    time = rep(0,365)
    # data primo evento farmacologico - data riferimento + 1 --> il primo giorno in cui c'è prescrizione
    ref = as.numeric(data_c[which(data_c[,"COD_REG"]==i)[1],"date_of_prescription"]-data_c[which(data_c[,"COD_REG"]==i),"data_rif_ev"][1]) + 1 
    #itero negli eventi del paziente i
    for (j in 1:count_event[count_event$COD_REG==i,"tot_pharm_class"])
    {
      duration = data_c[which(data_c[,"COD_REG"]==i)[j],"days"]
      time[ref:(ref+duration-1)] = 1
      ref = ref + duration
      if(j != count_event[count_event$COD_REG==i,"tot_pharm_class"])
      {
        no_prescription = as.numeric(data_c[which(data_c$COD_REG == i)[j+1],"date_of_prescription"]-data_c[which(data_c$COD_REG == i)[j],"end_of_prescription"])
        ref = ref + no_prescription 
      }  
    }
    somm_day[as.character(i),]=time[1:365]
  }
  
  
  #Aggiungo la colonna con i codici
  somm_day <- cbind("COD_REG" = unique(codici), somm_day)
  rownames(somm_day) <- seq(1:dim(somm_day)[1])
  return(somm_day)
}