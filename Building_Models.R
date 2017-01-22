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