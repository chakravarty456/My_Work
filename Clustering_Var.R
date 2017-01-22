###########################################################################################

# CLEAR THE ENVIRONMENT AND GET SET FOR THE ANALYSIS
rm(list=(ls(all=TRUE)))
setwd("filepath")

# READ THE DATA SET INTO R
RawData<-read.csv("file.csv",header=TRUE)
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