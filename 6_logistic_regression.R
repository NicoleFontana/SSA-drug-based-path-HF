#####
# Code - Section 3.5: Logistic regressione #
#####
#This code presents the logistic regression model used to predict whether the death of 
#the patient has occurred within the follow-up time.

#Libraries
library(glmnet)
library(car)
library(dplyr)
library(caret)
library(ROCR)

#Load the final dataset from file "5_final_dataset.R"
load("dataset/final_dataset.Rdata")

#1. Compute the logistic model with all the covariates
logistic_base <- glm(death ~ age + gender + LOS + MCS + tot_procedures +
                       cluster + DIU + AAG, data = data, family = binomial('logit'))
summary(logistic_base) 

#Test collinearity
vif(logistic_base) 

#2.Plots Odds ratio
alpha = 0.05
qalpha = qnorm( 1 - alpha/2 )
##
val = rownames(coefficients(summary(logistic_base)))[-1]
coef = 1
INF <- exp(coef*logistic_base$coefficients[val] - coef*qalpha * summary(logistic_base)$coefficients[val,2])
CENT <- exp(coef*logistic_base$coefficients[val])
SUP <- exp(coef*logistic_base$coefficients[val] + coef*qalpha * summary(logistic_base)$coefficients[val,2])
df <- data.frame(yAxis = c(1:15),
                 IC_Odds_Ratio = val,
                 INF = INF,
                 x = CENT,
                 SUP = SUP)

df$IC_Odds_Ratio <- factor(df$IC_Odds_Ratio, levels=rev(df$IC_Odds_Ratio))
color <-if_else(df$SUP<1 & df$INF < 1,"green3",if_else(df$SUP>1 & df$INF<1, "black","red"))
df <- df %>% mutate(SUP_prova=c(NA,SUP[2],NA,SUP[4:15]))
df <- df %>% mutate(INF_prova=c(NA,INF[2],NA,INF[4:15]))

ggplot(df, aes(x = log(x), y =IC_Odds_Ratio), fill=x) +
  geom_point(shape="square",size = 1.5,colour=color) +
  geom_errorbar(aes(xmax = log(SUP_prova), xmin = log(INF_prova)),colour ="black",
                width=0.3,size=0.5) +
  scale_x_continuous(breaks = log(c(seq(0.4,1.2,0.2),seq(1.2,2.5,0.4))), 
                     labels = c(seq(0.4,1.2,0.2),seq(1.2,2.5,0.4)),
                     limits = log(c(0.4,2.5))) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  ggtitle("Odds Ratio (95% CI)")+theme_bw()+
  theme(text = element_text(size = 15, colour = "black")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(),plot.title.position = "plot") +
  geom_vline(aes(xintercept = log(1)),linetype="dashed")+ 
  labs( x = NULL, y=NULL)


#3. Compute the metrics using a 10-fold cross validation
train.control <- trainControl(method = "cv",number=10, savePredictions = TRUE)
data$death <- as.factor(data$death)
m_step_CV = train(
  form = death ~ age + gender + LOS + MCS + tot_procedures + cluster + DIU + AAG,
  data = data,
  trControl = train.control,
  method = "glm",
  family = "binomial"
)

#Compute accuracy for each fold
pred <- m_step_CV$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1, 0)
eachfold <- pred %>% group_by(Resample) %>% summarise_at(vars(equal), list(Accuracy = mean))              
c(mean(eachfold$Accuracy),sd(eachfold$Accuracy))

#Compute all metrics
confusionMatrix(table((m_step_CV$pred)$pred,(m_step_CV$pred)$obs), mode = "everything", positive="1")

# Compute ROC curve
prob <- predict(m_step_CV, newdata=data, type="prob")
pred <- prediction(prob[,2], data$death)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc   
