# Install ggbiplot from github
library(devtools)
install_github("vqv/ggbiplot", force=TRUE)
1

# Declare packages
library(depmixS4)
library(lubridate)
library(dplyr)
library(stats)
library(ggbiplot)
library(ggplot2)
library(zoo)
library(mice)

# Remove all object to clear environment
rm(list=ls())

# List all objects in work space
ls()

# Get data from file
setwd(getwd())
setwd("C:/Users/pnita/Downloads/")
df_full <- read.table(file="TermProjectData.txt",
                      header=TRUE, sep = ",")

####################################################
# CLEAN DATA: FILL MISSING VALUES & SEPARATE DATASET
####################################################

# Check data for missing values
sapply(df_full, function(x) sum(is.na(x)))

# Replace NA values using 'Predictive mean matching' method
temp_data <- mice(df_full,m=1, method = 'pmm', seed=500)
summary(temp_data)

# Save data
complete_data <- complete(temp_data)

# Check missing values again to make sure everything is filled
sapply(complete_data, function(x) sum(is.na(x)))

# Scale data
complete_data[c(3:9)] <- scale(complete_data[c(3:9)])

# Print head
class(complete_data)
head(complete_data)

# Extract training data (First 2 years):
train_date_start <- as.POSIXlt("16/12/2006", format="%d/%m/%Y")
train_date_end <- as.POSIXlt("16/12/2008", format="%d/%m/%Y")

train_data_full <- complete_data[as.POSIXlt(complete_data$Date, format="%d/%m/%Y") >=train_date_start &
                                   as.POSIXlt(complete_data$Date, format="%d/%m/%Y") <=train_date_end,]

# Extract testing data (last 1 year):
test_date_start <- as.POSIXlt("17/12/2008", format="%d/%m/%Y")
test_date_end <- as.POSIXlt("16/12/2009", format="%d/%m/%Y")

test_data_full <- complete_data[as.POSIXlt(complete_data$Date, format="%d/%m/%Y") >=test_date_start &
                                  as.POSIXlt(complete_data$Date, format="%d/%m/%Y") <=test_date_end,]

class(train_data_full)
head(train_data_full)

class(test_data_full)
head(test_data_full)

###############################################
# TASK 1: PERFORM PRINCIPAL COMPONENT ANALYSIS
###############################################

# Compute PCA (with scaled data)
pca_data <- train_data_full[2:nrow(train_data_full),3:9]
pca <- prcomp(pca_data, scale = TRUE, center = TRUE)

# Scree plot (to see total variance explained by each PC)
ggscreeplot(pca) + labs(title = "Scree plot") +
  theme(legend.position = "bottom") + theme_minimal()

# Print summary & Eigenvalues & loading values
summary(pca)
pca$sdev ^ 2
pca$rotation

# Biplot (display results) - 1 plot with data points, 1 plot without
ggbiplot(pca) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(title = "PCA of contributing response variables") +
  theme(legend.position = "bottom") + theme_minimal()

ggbiplot(pca, alpha = 0) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(title = "PCA of contributing response variables") +
  theme(legend.position = "bottom") + theme_minimal()

###############################################
# FILTER DATA: CHOOSE VARIABLES & TIME SERIES
###############################################

# Extract 3 chosen variables only
train_data <- train_data_full
test_data <- test_data_full

train_data <- train_data[,c(1,2,3,6)]
test_data <- test_data[,c(1,2,3,6)]

# Add weekday column: Number from 1 -> 7 = Mon -> Sun
train_data$weekDay <- wday(ymd(as.POSIXlt(train_data$Date, format="%d/%m/%Y")), week_start=1)
test_data$weekDay <- wday(ymd(as.POSIXlt(test_data$Date, format="%d/%m/%Y")), week_start=1)

# Set time window (180 minutes)
time_start <- as.POSIXlt("9:00:00", format="%H:%M:%S")
time_end <- as.POSIXlt("11:59:00", format="%H:%M:%S")

# Filter data by time window (from 9am - 11:59am | Wednesday)
train_data <- train_data[train_data$weekDay == 3 &
                           as.POSIXlt(train_data$Time, format="%H:%M:%S") >=time_start &
                           as.POSIXlt(train_data$Time, format="%H:%M:%S") <=time_end,]

test_data <- test_data[test_data$weekDay == 3 &
                         as.POSIXlt(test_data$Time, format="%H:%M:%S") >=time_start &
                         as.POSIXlt(test_data$Time, format="%H:%M:%S") <=time_end,]

# Drop last weekDay col
train_data <- train_data[,c(1:4)] 
test_data <- test_data[,c(1:4)]

class(train_data)
head(train_data)

class(test_data)
head(test_data)

###############################################
# TASK 2: HMM TRAINING & TESTING 
###############################################

############ Training Models ############

# Create new empty data frame to save Log-likelihood & BIC values
feature_train_data <- data.frame(numberOfStates=integer(),
                                 logLikelihood=double(),
                                 BIC = double(),
                                 stringsAsFactors=FALSE)

# training models starting at 4, ending at 24, by 4s
for (n_states in seq(4, 20, 4)) {
  set.seed(5)
  model <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                  data=train_data,nstates=n_states,
                  family=list(gaussian(),gaussian()),
                  ntimes=rep(c(180),times=104))
  # fitting the training model
  fitModel <- fit(model)

  # getting the loglikelihood and BIC
  log_likelihood = logLik(fitModel)/104
  BIC_value = BIC(fitModel)/104
  
  # Save values row by row
  row_values = data.frame(n_states, log_likelihood, BIC_value) 
  feature_train_data=rbind(feature_train_data,row_values)
}

# Note: states 21 to 24 would not converge

# Plotting the log likelihoods and BICs of of the trained models
ggplot(feature_train_data, aes(x = n_states)) + 
  geom_point(aes(y = log_likelihood, color='Log-likelihood')) +
  geom_point(aes(y = BIC_value, color='BIC')) +
  geom_line(aes(y = log_likelihood, color = 'Log-likelihood'), lwd=1) + 
  geom_line(aes(y = BIC_value, color = 'BIC'), lwd=1) +
  labs(title = 'Log-likelihood & BIC by number of states', x = 'nstates', y = 'Value') +
  theme_minimal()

##########################################################
### Training range of candidates 6-8
##########################################################

candidates_data <- data.frame(numberOfStates=integer(),
                                 logLikelihood=double(),
                                 BIC = double(),
                                 stringsAsFactors=FALSE)


### candidate 6 ###

set.seed(5)
model1 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                 data=train_data,nstates=6,
                 family=list(gaussian(),gaussian()),
                 ntimes=rep(c(180),times=104))

# fitting the training model
fitModel1 <- fit(model1)

# storing the model's data in a dataframe
n_states = 6
log_likelihood = logLik(fitModel1)/104
BIC_value = BIC(fitModel1)/104
row_values = data.frame(n_states, log_likelihood, BIC_value) 
candidates_data=rbind(candidates_data,row_values)


### candidate 7 ###

set.seed(5)
model2 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                 data=train_data,nstates=7,
                 family=list(gaussian(),gaussian()),
                 ntimes=rep(c(180),times=104))
# fitting the training model
fitModel2 <- fit(model2)

# storing the model's data in a dataframe
n_states = 7
log_likelihood = logLik(fitModel2)/104
BIC_value = BIC(fitModel2)/104
row_values = data.frame(n_states, log_likelihood, BIC_value) 
candidates_data=rbind(candidates_data,row_values)


### candidate 8 ### 

set.seed(5)
model3 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                data=train_data,nstates=8,
                family=list(gaussian(),gaussian()),
                ntimes=rep(c(180),times=104))
# fitting the training model
fitModel3 <- fit(model3)

# storing the model's data in a dataframe
n_states = 8
log_likelihood = logLik(fitModel3)/104
BIC_value = BIC(fitModel3)/104
row_values = data.frame(n_states, log_likelihood, BIC_value) 
candidates_data=rbind(candidates_data,row_values)

### Note: Chosen Candidates are 7, 8, 9

######################
### Testing Models ###
######################

cand_testing_data <- data.frame(numberOfStates=integer(),
                                  forback_logLik=double(),
                                  fit_logLik=double(),
                                 stringsAsFactors=FALSE)


### Candidate 6 testing ###

set.seed(5)
# Creating test model 1
test_model1 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=test_data,nstates=6,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=50))

# copying the parameters of the fitted train model 1
test_model1 <-  setpars(test_model1, getpars(fitModel1))

# forward backward
test1_forback <- forwardbackward(test_model1)

# putting the data into a data frame
n_states = 6
test_logLik = test1_forback[[6]]/50
train_logLik = candidates_data[1,2]
row_values = data.frame(n_states, test_logLik, train_logLik) 
cand_testing_data=rbind(cand_testing_data,row_values)


### Candidate 7 testing ###

set.seed(5)
# Creating test model 2
test_model2 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=test_data,nstates=7,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=50))

# copying the parameters of the fitted train model 1
test_model2 <-  setpars(test_model2, getpars(fitModel2))

# forward backward
test2_forback <- forwardbackward(test_model2)

# putting the data into a data frame
n_states = 7
test_logLik = test2_forback[[6]]/50
train_logLik = candidates_data[2,2]
row_values = data.frame(n_states, test_logLik, train_logLik) 
cand_testing_data=rbind(cand_testing_data,row_values)
cand_testing_data[2,] <- row_values


### Candidate 8 testing ###

set.seed(5)
# Creating test model 3
test_model3 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=test_data,nstates=8,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=50))

# copying the parameters of the fitted train model 2
test_model3 <-  setpars(test_model3, getpars(fitModel3))

# forward backward
test3_forback <- forwardbackward(test_model3)

# putting the data into a data frame
n_states = 8
test_logLik = test3_forback[[6]]/50
train_logLik = candidates_data[3,2]
row_values = data.frame(n_states, test_logLik, train_logLik) 
cand_testing_data=rbind(cand_testing_data,row_values)

cand_testing_data$difference <- 0

for (i in 1:3) {
  cand_testing_data[i,4] <- cand_testing_data[i,2] - cand_testing_data[i,3]
}

# Note: the best candidate is candidate 7

#################################
### TASK 3: ANOMALY DETECTION ###
#################################

### Preparing Anomalous Data sets ###

df_anom1 <- read.table(file="DataWithAnomalies1.txt",
                       header=TRUE, sep = ",")
df_anom2 <- read.table(file="DataWithAnomalies2.txt",
                       header=TRUE, sep = ",")
df_anom3 <- read.table(file="DataWithAnomalies3.txt",
                       header=TRUE, sep = ",")

# Check data for missing values
sapply(df_anom1, function(x) sum(is.na(x)))
sapply(df_anom2, function(x) sum(is.na(x)))
sapply(df_anom3, function(x) sum(is.na(x)))

# Replace NA values using 'Predictive mean matching' method
temp_anom1 <- mice(df_anom1,m=1, method = 'pmm', seed=500)
temp_anom2 <- mice(df_anom2,m=1, method = 'pmm', seed=500)
temp_anom3 <- mice(df_anom3,m=1, method = 'pmm', seed=500)

# Save data
complete_anom1 <- complete(temp_anom1)
complete_anom2 <- complete(temp_anom2)
complete_anom3 <- complete(temp_anom3)

# Check missing values again to make sure everything is filled
sapply(complete_anom1, function(x) sum(is.na(x)))
sapply(complete_anom2, function(x) sum(is.na(x)))
sapply(complete_anom3, function(x) sum(is.na(x)))

# Scale data
complete_anom1[c(3:9)] <- scale(complete_anom1[c(3:9)])
complete_anom2[c(3:9)] <- scale(complete_anom2[c(3:9)])
complete_anom3[c(3:9)] <- scale(complete_anom3[c(3:9)])

anom1_data <- complete_anom1[,c(1,2,3,6)]
anom2_data <- complete_anom2[,c(1,2,3,6)]
anom3_data <- complete_anom3[,c(1,2,3,6)]

# Add weekday column: Number from 1 -> 7 = Mon -> Sun
anom1_data$weekDay <- wday(ymd(as.POSIXlt(anom1_data$Date,
                                          format="%d/%m/%Y")), week_start=1)
anom2_data$weekDay <- wday(ymd(as.POSIXlt(anom2_data$Date,
                                          format="%d/%m/%Y")), week_start=1)
anom3_data$weekDay <- wday(ymd(as.POSIXlt(anom3_data$Date,
                                          format="%d/%m/%Y")), week_start=1)

# Set time window (180 minutes)
time_start <- as.POSIXlt("9:00:00", format="%H:%M:%S")
time_end <- as.POSIXlt("11:59:00", format="%H:%M:%S")

# Filter data by time window (from 9am - 11:59am | Wednesday)
anom1_data <- anom1_data[anom1_data$weekDay == 3 &
                           as.POSIXlt(anom1_data$Time, format="%H:%M:%S") >=time_start &
                           as.POSIXlt(anom1_data$Time, format="%H:%M:%S") <=time_end,]
anom2_data <- anom2_data[anom2_data$weekDay == 3 &
                           as.POSIXlt(anom2_data$Time, format="%H:%M:%S") >=time_start &
                           as.POSIXlt(anom2_data$Time, format="%H:%M:%S") <=time_end,]
anom3_data <- anom3_data[anom3_data$weekDay == 3 &
                           as.POSIXlt(anom3_data$Time, format="%H:%M:%S") >=time_start &
                           as.POSIXlt(anom3_data$Time, format="%H:%M:%S") <=time_end,]

# Drop last weekDay col
anom1_data <- anom1_data[,c(1:4)]
anom2_data <- anom2_data[,c(1:4)]
anom3_data <- anom3_data[,c(1:4)]


######################
### Anomaly Models ###
######################

anomaly_data <- data.frame(anomalySet=integer(),
                                forback_logLik=double(),
                                fit_logLik=double(),
                                stringsAsFactors=FALSE)


### Anomaly data set 1 ###

set.seed(5)
# testing anomaly set 1
anom_model1 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=anom1_data,nstates= 7,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=52))

# copying the parameters of the fitted train model 1
anom_model1 <-  setpars(anom_model1, getpars(fitModel2))

# forward backward
anom1_forback <- forwardbackward(anom_model1)

# putting the data into a data frame
anomaly_set = 1
anomaly_logLik = anom1_forback[[6]]/52
train_logLik = cand_testing_data[2,3]
row_values = data.frame(anomaly_set, anomaly_logLik, train_logLik) 
anomaly_data=rbind(anomaly_data,row_values)


### Anomaly data set 2 ###

set.seed(5)
# testing anomaly set 2
anom_model2 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=anom2_data,nstates= 7,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=52))

# copying the parameters of the fitted train model 1
anom_model2 <-  setpars(anom_model2, getpars(fitModel2))

# forward backward
anom2_forback <- forwardbackward(anom_model2)

# putting the data into a data frame
anomaly_set = 2
anomaly_logLik = anom2_forback[[6]]/52
train_logLik = cand_testing_data[2,3]
row_values = data.frame(anomaly_set, anomaly_logLik, train_logLik) 
anomaly_data=rbind(anomaly_data,row_values)


### Anomaly data set 3 ###

set.seed(5)
# testing anomaly set 3
anom_model3 <- depmix(response = list(Global_active_power~1,Global_intensity~1),
                      data=anom3_data,nstates= 7,
                      family=list(gaussian(),gaussian()),
                      ntimes=rep(c(180),times=52))

# copying the parameters of the fitted train model 1
anom_model3 <-  setpars(anom_model3, getpars(fitModel2))

# forward backward
anom3_forback <- forwardbackward(anom_model3)

# putting the data into a data frame
anomaly_set = 3
anomaly_logLik = anom3_forback[[6]]/52
train_logLik = cand_testing_data[2,3]
row_values = data.frame(anomaly_set, anomaly_logLik, train_logLik) 
anomaly_data=rbind(anomaly_data,row_values)

# Creates a column for the difference of the anomaly loglik and train loglik
anomaly_data$diff_subtract <- 0

for (i in 1:3) {
  anomaly_data[i,4] <- anomaly_data[i,2] - anomaly_data[i,3]
}


### conclusion ###
# Anomaly data set 3 has the highest difference of log likelihood, 
# and therefore the most anomalous. Anomaly data set 1 has the least
# log likelihood and can be considered the least anomalous data set 
# out of the 3.

