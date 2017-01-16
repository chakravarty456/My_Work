##PROJECT: PREAUTHORIZATION
# CLEAR THE ENVIRONMENT AND GET SET FOR THE ANALYSIS
rm(list=(ls(all=TRUE)))
setwd("C:/Users/Pinky/Desktop/Nivi/INSOFE/Project")

# READ THE DATA SET INTO R
data<-read.csv("PriorAuthData.csv",header=TRUE)
summary(data)
str(data)

#CHECK FOR MISSING VALUES
sum(is.na(data))

#DELETE PERSON ID AND DOCTOR ID COLUMNS AS THEY DON'T INFLUENCE THE TARGET VARIABLE
data=data[ ,c(-1,-10)]
str(data)

# ANY REPEATED OR NON-UNIQUE ENTRIES IN THE DATA MIGHT SKEW THE ANALYSIS.
# HENCE, RETRIEVING ONLY UNIQUE DATA FOR FURTHER ANALYSIS.
data <- unique(data)

#VIEW DRUG ID VS. TARGET
data1<-subset(data,select=c(Drug,Target))
table(data1$Drug, data1$Target)

##FINDINIG THE NUMBER OF CLUSTERS FOR EACH ATTRIBUTE USING K-MEANS###

#Drug ID Column
data_drug <- data1[setdiff(names(data1),"Target")]

#DrugID_WSSPlot1
#Clusters=3
wssplot1<- function(data1, nc=10, seed=1234){
  wss1<- (nrow(data_drug)-1)*sum(apply(data_drug,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss1[i] <- sum(kmeans(model.matrix(~.+0, data=data_drug), centers=i)$withinss)}
  plot(1:nc, wss1, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot1(data_drug)

#DrugSubClass_WSSPlot2
##Clusters=5
data_drugsubclass<-subset(data,select=c(DrugSubClass))

wssplot2<- function(data1, nc=15, seed=1234){
  wss2<- (nrow(data_drugsubclass)-1)*sum(apply(data_drugsubclass,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss2[i] <- sum(kmeans(model.matrix(~.+0, data=data_drugsubclass), centers=i)$withinss)}
  plot(1:nc, wss2, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot2(data_drugsubclass)

#DrugClass_WSSPlot3
#Clusters=5
data_drugclass<-subset(data,select=c(DrugClass))

wssplot3<- function(data1, nc=15, seed=1234){
  wss3<- (nrow(data_drugclass)-1)*sum(apply(data_drugclass,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss3[i] <- sum(kmeans(model.matrix(~.+0, data=data_drugclass), centers=i)$withinss)}
  plot(1:nc, wss3, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot3(data_drugsubclass)

#DrugChemicalName_WSSPlot4
#Clusters=5
data_drugchemname<-subset(data,select=c(Drug_Chemical_Name))

wssplot4<- function(data1, nc=10, seed=1234){
  wss4<- (nrow(data_drugchemname)-1)*sum(apply(data_drugchemname,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss4[i] <- sum(kmeans(model.matrix(~.+0, data=data_drugchemname), centers=i)$withinss)}
  plot(1:nc, wss4, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot4(data_drugchemname)

#GPI_WSSPlot5
#Clusters=4
data_gpi<-subset(data,select=c(GPI))

wssplot5<- function(data1, nc=10, seed=1234){
  wss5<- (nrow(data_gpi)-1)*sum(apply(data_gpi,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss5[i] <- sum(kmeans(model.matrix(~.+0, data=data_gpi), centers=i)$withinss)}
  plot(1:nc, wss5, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot5(data_gpi)

#Drug Full GPI Name_WSSPlot6
#Clusters=4
data_fullgpi<-subset(data,select=c(Drug_Full_GPI_Name))

wssplot6<- function(data1, nc=10, seed=1234){
  wss6<- (nrow(data_fullgpi)-1)*sum(apply(data_fullgpi,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss6[i] <- sum(kmeans(model.matrix(~.+0, data=data_fullgpi), centers=i)$withinss)}
  plot(1:nc, wss6, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot6(data_fullgpi)

#NDC_WSSPlot7
#Clusters = 4
data_ndc<-subset(data,select=c(NDC))

wssplot7<- function(data1, nc=10, seed=1234){
  wss7<- (nrow(data_ndc)-1)*sum(apply(data_ndc,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss7[i] <- sum(kmeans(model.matrix(~.+0, data=data_ndc), centers=i)$withinss)}
  plot(1:nc, wss7, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot7(data_ndc)

#DrugGroups_WSSPlot8
#Clusters = 8
data_druggroup<-subset(data,select=c(DrugGroup))

wssplot8<- function(data1, nc=8, seed=1234){
  wss8<- (nrow(data_druggroup)-1)*sum(apply(data_druggroup,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss8[i] <- sum(kmeans(model.matrix(~.+0, data=data_druggroup), centers=i)$withinss)}
  plot(1:nc, wss8, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot8(data_druggroup)

#RxGroupID_WSSPlot9
#Clusters = 6
data_drugrxgroup<-subset(data,select=c(RxGroupId))

wssplot9<- function(data1, nc=10, seed=1234){
  wss9<- (nrow(data_drugrxgroup)-1)*sum(apply(data_drugrxgroup,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss9[i] <- sum(kmeans(model.matrix(~.+0, data=data_drugrxgroup), centers=i)$withinss)}
  plot(1:nc, wss9, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot9(data_drugrxgroup)

#Bin_WSSPlot10
#Clusters = 4
data_bin<-subset(data,select=c(Bin))

wssplot10<- function(data1, nc=8, seed=1234){
  wss10<- (nrow(data_bin)-1)*sum(apply(data_bin,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss10[i] <- sum(kmeans(model.matrix(~.+0, data=data_bin), centers=i)$withinss)}
  plot(1:nc, wss10, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot10(data_bin)

#PCN_WSSPlot11
#Clusters = 8
data_pcn<-subset(data,select=c(PCN))

wssplot11<- function(data1, nc=8, seed=1234){
  wss11<- (nrow(data_pcn)-1)*sum(apply(data_pcn,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss11[i] <- sum(kmeans(model.matrix(~.+0, data=data_pcn), centers=i)$withinss)}
  plot(1:nc, wss11, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot11(data_pcn)

#Zip_WSSPlot12
#Clusters = 5
data_zip<-subset(data,select=c(Zip))

wssplot12<- function(data1, nc=8, seed=1234){
  wss12<- (nrow(data_zip)-1)*sum(apply(data_zip,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss12[i] <- sum(kmeans(model.matrix(~.+0, data=data_zip), centers=i)$withinss)}
  plot(1:nc, wss12, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot12(data_zip)

#State_WSSPlot13
#Clusters = 3
data_state<-subset(data,select=c(State))

wssplot13<- function(data1, nc=8, seed=1234){
  wss13<- (nrow(data_state)-1)*sum(apply(data_state,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss13[i] <- sum(kmeans(model.matrix(~.+0, data=data_state), centers=i)$withinss)}
  plot(1:nc, wss13, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot13(data_state)

#TransDate_WSSPlot14
#Clusters = 8
data_transdate<-subset(data,select=c(TransDate))

wssplot14<- function(data1, nc=10, seed=1234){
  wss14<- (nrow(data_transdate)-1)*sum(apply(data_transdate,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss14[i] <- sum(kmeans(model.matrix(~.+0, data_transdate), centers=i)$withinss)}
  plot(1:nc, wss14, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot14(data_transdate)

###########################################################################################

# CLEAR THE ENVIRONMENT AND GET SET FOR THE ANALYSIS
rm(list=(ls(all=TRUE)))
setwd("C:/Users/Pinky/Desktop/Nivi/INSOFE/Project")

# READ THE DATA SET INTO R
RawData<-read.csv("PriorAuthData.csv",header=TRUE)
summary(RawData)
str(RawData)

#PRE-PROCESSING STEPS
# CHECK IF THERE ARE ANY NULL VALUES OR MISSING VALUES
sum(is.na(RawData))
# AS THERE ARE NO MISSING VALUES, WE DO NOT HAVE TO ELIMINATE ANY ROWS OR IMPUTE DATA

# REMOVE THE USER ID AND DOCTOR ID COLUMNS AS THEY ARE NOT RELEVANT FOR OUR ANALYSIS
RawData=RawData[ ,c(-1,-10)]

# ANY REPEATED OR NON-UNIQUE ENTRIES IN THE DATA MIGHT SKEW THE ANALYSIS.
# HENCE, RETRIEVING ONLY UNIQUE DATA FOR FURTHER ANALYSIS. 
# UniqueData <- unique(unlist(data))
UniqueData <- unique(RawData)

#TAKING TARGET VARIABLE AS A DATAFRAME
TargetVar=UniqueData[ ,c(10)]
TargetVar.DF <- as.data.frame(UniqueData$Target)
colnames(TargetVar.DF)[1]="Target"
TargetVar.DF


# REMOVE THE RAWDATA DATA FRAME FROM R MEMORY
rm(RawData)


# IDENTIFY THE NUMBER OF CLUSTERS FOR DRUG
# READ THE FEATURE TO BE CLUSTERED INTO A VARIABLE
var.for.kmeans<-subset(UniqueData,select=c(Drug,Target))
# FOR ANALYSIS, PROJECT A TABLE OF DRUG AND TARGET RELATIONSHIP
table(var.for.kmeans$Drug, var.for.kmeans$Target)
# REMOVE THE TARGET VARIABLE, IT IS NOT NECESSARY FOR K MEANS
final.var.for.kmeans <- var.for.kmeans[setdiff(names(var.for.kmeans),"Target")]
# IDENTIFY THE OPTIMAL CLUSTERS REQUIRED FOR DRUG 
wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(final.var.for.kmeans)-1)*sum(apply(final.var.for.kmeans,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=i)$withinss)
  }
  plot(1:nc, wss, type="b", xlab="Number of Clusters", ylab="WithinSS")
}
# CALL THE ABOVE FUNCTION
wssplot(final.var.for.kmeans)
# FROM ABOVE PLOT, THE OPTIMAL K FOR DRUG IS IDENTIFIED AS 3
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=3)
kmeansModel

# NOW READ THE PARTITIONED DATA FROM KMEANS OUTPUT INTO A DATAFRAME
newDrug.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrug.DF)[1]="Drug Name"

# THE ABOVE FUNCTION IS INVOKED FOR ALL THE VARIABLES THAT NEED TO BE CLASSIFIED/CLUSTERED

# NOW REPEAT APPLYING K MEANS AND IDENTIFYING THE DATA FRAME FOR ALL THE VARIABLES
# DrugSubClass - OPTIMAL K = 5
final.var.for.kmeans<-subset(UniqueData,select=c(DrugSubClass))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=5)
newDrugSubClass.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrugSubClass.DF)[1]="Drug Sub Class"
ClusteredData <- c(newDrug.DF, newDrugSubClass.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# DrugClass - OPTIMAL K = 5
final.var.for.kmeans<-subset(UniqueData,select=c(DrugClass))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=5)
newDrugClass.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrugClass.DF)[1]="Drug Class"
ClusteredData <- c(ClusteredData.DF, newDrugClass.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# Drug_Chemical_Name - OPTIMAL K = 5
final.var.for.kmeans<-subset(UniqueData,select=c(Drug_Chemical_Name))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=5)
newDrug_Chemical_Name.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrug_Chemical_Name.DF)[1]="Chemical Name"
ClusteredData <- c(ClusteredData.DF, newDrug_Chemical_Name.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# GPI - OPTIMAL K = 4
final.var.for.kmeans<-subset(UniqueData,select=c(GPI))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=4)
newGPI.DF <- as.data.frame(kmeansModel$cluster)
colnames(newGPI.DF)[1]="GPI"
ClusteredData <- c(ClusteredData.DF, newGPI.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# Drug_Full_GPI_Name - OPTIMAL K = 4
final.var.for.kmeans<-subset(UniqueData,select=c(Drug_Full_GPI_Name))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=4)
newDrug_Full_GPI_Name.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrug_Full_GPI_Name.DF)[1]="Full GPI Name"
ClusteredData <- c(ClusteredData.DF, newDrug_Full_GPI_Name.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# NDC - OPTIMAL K = 4
final.var.for.kmeans<-subset(UniqueData,select=c(NDC))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=4)
newNDC.DF <- as.data.frame(kmeansModel$cluster)
colnames(newNDC.DF)[1]="NDC"
ClusteredData <- c(ClusteredData.DF, newNDC.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# DrugGroup - OPTIMAL K = 8
final.var.for.kmeans<-subset(UniqueData,select=c(DrugGroup))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=8)
newDrugGroup.DF <- as.data.frame(kmeansModel$cluster)
colnames(newDrugGroup.DF)[1]="Group"
ClusteredData <- c(ClusteredData.DF, newDrugGroup.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# RxGroupId - OPTIMAL K = 6
final.var.for.kmeans<-subset(UniqueData,select=c(RxGroupId))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=6)
newRxGroupId.DF <- as.data.frame(kmeansModel$cluster)
colnames(newRxGroupId.DF)[1]="Rx Group ID"
ClusteredData <- c(ClusteredData.DF, newRxGroupId.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# Bin - OPTIMAL K = 4
final.var.for.kmeans<-subset(UniqueData,select=c(Bin))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=4)
newBin.DF <- as.data.frame(kmeansModel$cluster)
colnames(newBin.DF)[1]="Bin"
ClusteredData <- c(ClusteredData.DF, newBin.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# PCN - OPTIMAL K = 8
final.var.for.kmeans<-subset(UniqueData,select=c(PCN))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=8)
newPCN.DF <- as.data.frame(kmeansModel$cluster)
colnames(newPCN.DF)[1]="PCN"
ClusteredData <- c(ClusteredData.DF, newPCN.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# Zip - OPTIMAL K = 5
final.var.for.kmeans<-subset(UniqueData,select=c(Zip))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=5)
newZip.DF <- as.data.frame(kmeansModel$cluster)
colnames(newZip.DF)[1]="Zip"
ClusteredData <- c(ClusteredData.DF, newZip.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# State - OPTIMAL K = 3
final.var.for.kmeans<-subset(UniqueData,select=c(State))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=3)
newState.DF <- as.data.frame(kmeansModel$cluster)
colnames(newState.DF)[1]="State"
ClusteredData <- c(ClusteredData.DF, newState.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

# TransDate - OPTIMAL K = 8
final.var.for.kmeans<-subset(UniqueData,select=c(TransDate))
kmeansModel <- kmeans(model.matrix(~.+0, data=final.var.for.kmeans), centers=8)
newTransDate.DF <- as.data.frame(kmeansModel$cluster)
colnames(newTransDate.DF)[1]="Trans Date"
ClusteredData <- c(ClusteredData.DF, newTransDate.DF)
ClusteredData.DF <- as.data.frame(ClusteredData)

##FINAL CLUSTERED DATA IN THE FORM OF DATAFRAME
ClusteredData.DF

# CLEAR ALL THE UNNCESSARY OBJECTS FROM THE ENVIRONMENT
rm(newDrug.DF, newDrugSubClass.DF, newDrugClass.DF, newDrug_Chemical_Name.DF, 
   newGPI.DF, newDrug_Full_GPI_Name.DF, newNDC.DF, newDrugGroup.DF, newRxGroupId.DF, 
   newBin.DF, newPCN.DF, newZip.DF,newState.DF,newTransDate.DF,final.var.for.kmeans, 
   kmeansModel, var.for.kmeans, wssplot)

###Working with Data Frame

#ADD TARGET VARIABLE DATAFRAME TO CLUSTERED DATAFRAME
ClusteredData <- c(ClusteredData.DF, TargetVar.DF)
ClusteredData.DF<-as.data.frame(ClusteredData)

# CONVERT THE TARGET VARIABLE INTO A FACTOR VARIABLE
ClusteredData.DF$Target <- as.factor(ClusteredData.DF$Target)
ClusteredData.DF$Drug.Name <- as.factor(ClusteredData.DF$Drug.Name)
ClusteredData.DF$Drug.Sub.Class <- as.factor(ClusteredData.DF$Drug.Sub.Class)
ClusteredData.DF$Drug.Class <- as.factor(ClusteredData.DF$Drug.Class)
ClusteredData.DF$Chemical.Name <- as.factor(ClusteredData.DF$Chemical.Name)
ClusteredData.DF$GPI <- as.factor(ClusteredData.DF$GPI)
ClusteredData.DF$Full.GPI.Name <- as.factor(ClusteredData.DF$Full.GPI.Name)
ClusteredData.DF$NDC <- as.factor(ClusteredData.DF$NDC)
ClusteredData.DF$Group <- as.factor(ClusteredData.DF$Group)
ClusteredData.DF$Rx.Group.ID <- as.factor(ClusteredData.DF$Rx.Group.ID)
ClusteredData.DF$Bin <- as.factor(ClusteredData.DF$Bin)
ClusteredData.DF$PCN <- as.factor(ClusteredData.DF$PCN)

# GENERATE TABLE TO SHOW THE EFFECT OF DEPENDANT VARIABLES ON INDEPENDENT VARIABLE
xtabs(~ Target + Drug.Name, data = ClusteredData.DF)
xtabs(~ Target + Drug.Sub.Class, data = ClusteredData.DF)
xtabs(~ Target + Drug.Class, data = ClusteredData.DF)
xtabs(~ Target + Chemical.Name, data = ClusteredData.DF)
xtabs(~ Target + GPI, data = ClusteredData.DF)
xtabs(~ Target + Full.GPI.Name, data = ClusteredData.DF)
xtabs(~ Target + NDC, data = ClusteredData.DF)
xtabs(~ Target + Group, data = ClusteredData.DF)
xtabs(~ Target + Rx.Group.ID, data = ClusteredData.DF)
xtabs(~ Target + Bin, data = ClusteredData.DF)
xtabs(~ Target + PCN, data = ClusteredData.DF)
xtabs(~ Target + Zip, data = ClusteredData.DF)
xtabs(~ Target + State, data = ClusteredData.DF)
xtabs(~ Target + Trans.Date, data = ClusteredData.DF)


#SETTING TARGET ATTRIBUTE AS
ClusteredData.DF$Target<-as.factor(as.character(ClusteredData.DF$Target))

#SPLITTING DATA INTO TRAIN AND TEST (90% TRAIN AND 10% TEST)
rows=seq(1,nrow(ClusteredData.DF),1)
set.seed(123)
trainRows=sample(rows,(90*nrow(ClusteredData.DF))/100)
trainData = ClusteredData.DF[trainRows,] 
testData = ClusteredData.DF[-trainRows,]

# PROCEED WITH MODEL GENERATION
#### MODEL 1: LOGISTIC REGRESSION
# RUN LOGISTIC REGRESSION USING ALL FEATURES
LogReg<- glm(Target ~ ., data=trainData,family=binomial)
summary(LogReg)
prob<-predict(LogReg, type="response")
pred_class<- factor(ifelse(prob> 0.5, 1, 0))
table(trainData$Target,pred_class)

# RUN LOGISTIC REGRESSION USING A SUB SET OF FEATURES
LogReg <- glm(formula = Target ~ Drug.Name + Drug.Sub.Class + Drug.Class + Chemical.Name + GPI + NDC + Group + Rx.Group.ID + Bin + PCN, family = binomial, data = trainData)
summary(LogReg)
prob<-predict(LogReg, type="response")
pred_class<- factor(ifelse(prob> 0.5, 1, 0))
table(trainData$Target,pred_class)

# RUN LOGISTIC REGRESSION AGAINST TEST DATA
test.results<- predict(LogReg,testData,type='response')
test.class<- ifelse(test.results> 0.5,1,0)
test.values<-as.data.frame(test.class)
table(testData$Target,test.class)

##ERROR METRICS ON TRAIN DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable_train = table(trainData$Target,pred_class)
accuracy = sum(diag(precTable_train))/sum(precTable_train)*100
print(accuracy)

# PRECISION
precision_train=(precTable_train[2,2])/(precTable_train[1,2]+precTable_train[2,2])*100
print(precision_train)

# RECALL
recall_train = ((precTable_train[2,2])/(precTable_train[2,1]+precTable_train[2,2]))*100
print(recall_train)

##ERROR METRICS ON TEST DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable = table(testData$Target,test.class)
accuracy = sum(diag(precTable))/sum(precTable)*100
print(accuracy)

# PRECISION
precision=(precTable[2,2])/(precTable[1,2]+precTable[2,2])*100
print(precision)

# RECALL
recall = ((precTable[2,2])/(precTable[2,1]+precTable[2,2]))*100
print(recall)

#ROC Curve
library(ROCR)
p <- predict(LogReg,trainData, type="response") 
pr <- prediction(p, trainData$Target) 
prf <- performance(pr, measure = "tpr", x.measure = "fpr") 
plot(prf) 
abline(a=0, b= 1)  
auc <- performance(pr, measure = "auc") 
auc <- auc@y.values[[1]] 
auc  # very low  

#### MODEL 2: DECISION TREES using C50 
#install.packages('C50')
library(C50)
C50 <- C5.0(Target ~ Drug.Name + Drug.Sub.Class + Drug.Class + Chemical.Name + GPI + NDC + Group + Rx.Group.ID + Bin + PCN, trials = 10, data=trainData, rules=T)
summary(C50)
C5imp(C50, pct=true)

# RUN C50 DECISION TREES AGAINST TRAIN DATA
train.results<- predict(C50,trainData,type='class')
table(trainData$Target,train.results)

# RUN C50 DECISION TREES AGAINST TEST DATA
test.results<- predict(C50,testData,type='class')
table(testData$Target,test.results)

#ERROR METRICS ON TRAIN DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable_trn = table(trainData$Target,train.results)
accuracy_train = sum(diag(precTable_trn))/sum(precTable_trn)*100
print(accuracy_train)

# PRECISION
pcEval_train=(precTable_trn[2,2])/(precTable_trn[1,2]+precTable_trn[2,2])*100
print(pcEval_train)

# RECALL
recall_train = ((precTable_trn[2,2])/(precTable_trn[2,1]+precTable_trn[2,2]))*100
print(recall_train)

#ERROR METRICS ON TEST DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable = table(testData$Target,test.results)
accuracy = sum(diag(precTable))/sum(precTable)*100
print(accuracy)

# PRECISION
pcEval=(precTable[2,2])/(precTable[1,2]+precTable[2,2])*100
print(pcEval)

# RECALL
recall = ((precTable[2,2])/(precTable[2,1]+precTable[2,2]))*100
print(recall )

#### MODEL 3: RANDOM FORESTS
#install.packages("randomForest")
library(randomForest)   
set.seed(415)
rf <- randomForest(Target ~ ., data=trainData, keep.forest=TRUE, ntree=30)
rf <- randomForest(Target ~ Drug.Name + Drug.Sub.Class + Drug.Class + Chemical.Name + GPI + NDC + Group + Rx.Group.ID + Bin + PCN, data=trainData, keep.forest=TRUE, ntree=30)
print(rf)
round(importance(rf), 2)
varImpPlot(rf)

# RUN RANDOM FOREST AGAINST TRAIN DATA
train.results<- predict(rf,trainData,type='class')
table(trainData$Target,train.results)

# RUN RANDOM FOREST AGAINST TEST DATA
test.results<- predict(rf,testData,type='class')
table(testData$Target,test.results)

#ERROR METRICS ON TRAIN DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable_train = table(trainData$Target,train.results)
accuracy_train = sum(diag(precTable_train))/sum(precTable_train)*100
print(accuracy_train)

# PRECISION
pcEval_train=(precTable_train[2,2])/(precTable_train[1,2]+precTable_train[2,2])*100
print(pcEval_train)

# RECALL
recall_train = ((precTable_train[2,2])/(precTable_train[2,1]+precTable_train[2,2]))*100
print(recall_train)

#ERROR METRICS ON TEST DATA
# CALCULATE ACCURACY, PRECISION AND RECALL
# ACCURACY
precTable = table(testData$Target,test.results)
accuracy = sum(diag(precTable))/sum(precTable)*100
print(accuracy)

# PRECISION
pcEval=(precTable[2,2])/(precTable[1,2]+precTable[2,2])*100
print(pcEval)

# RECALL
recall = ((precTable[2,2])/(precTable[2,1]+precTable[2,2]))*100
print(recall )
