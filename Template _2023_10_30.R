#installing packages -------------------------------- 
#ignore if installed.
install.packages("bnlearn")  #to use bn algorithm
install.packages("parallel") # to parallelize algorithm runs (base package)
install.packages("foreign")  #to read spss
install.packages("readxl")   #to read excel at worksheet level
install.packages("mmr")      #to do  matrix multiplication
install.packages("data.table") #to load as data.table (part of data.frame package)
install.packages("hot.deck") #to check discrete or not
install.packages("ggplot2") #for splitting in case of continous variables


#Load all required packages-----------------------------
library(bnlearn)
library(foreign)
library(parallel)
library(readxl) 
library(mmr)  
library(data.table)
library(hot.deck)
library(ggplot2)

# preparing environment and loading data file-----------------

rm(list=ls()) #clear environment (clears variables and data)

setwd("paste folder path")           #sets working directory   
# ***important --  replace \ with /  while pasting paths

sink("Log file.txt",split=T, append = T)  #log file


raw_data<- read.spss("spss file path", to.data.frame=T, use.value.labels=F, stringsAsFactors=F) #loads spss data in R
# ***important --  replace \ with /  while pasting paths

nrow(raw_data) #check the sample size

# creating factors--------------

factor_loadings <- data.frame(read_excel("Filename.xlsx", sheet = "Factor loadings"))
factor_loadings[is.na(factor_loadings)] = 0 #assigns values 0 to NA items
var_order = factor_loadings[,1]  # ensure that attribute names are in the 1st col in factor loading sheet
raw_data = raw_data[,var_order] # Chg variable order to that of attributes in factor loading sheet and drops additional variables

names(raw_data) #check imagery variables 

for (i in names(factor_loadings)[-1])  # Normalize the loading for each factor
{
  sum = sum(factor_loadings[,i])
  factor_loadings[,i] = factor_loadings[,i]/sum
}


dim(raw_data)                 #check the dimensions of raw_data
dim(factor_loadings[,-1]) #check the dimensions of factor loading

#*** Important number of columns in raw data must match number of rows in factor_loading.

factor_data = mm(as.matrix(raw_data),as.matrix(factor_loadings[,-1])) #create factors same as Smart PLS
dim(factor_data)
write.csv(factor_data,"factor_data.csv")

#creating paths ---------------

constraints <- data.frame(read_excel("File_name.xlsx", sheet = "Permitted_Links")) #reading constraints from permitted links for bad arcs
names(constraints)[1] = 'Variables' #naming variable column
colnames(constraints) = gsub(".", "_", colnames(constraints), fixed=TRUE) # Chg column names to R standard removing spaces( ) and periods(.) to underscores(_)
constraints[1] = lapply(constraints[1], gsub, pattern = " ", replacement = "_", fixed = TRUE) # Chg variable names

badarcs = data.table() # creating empty data table for badacrs

for (i in 2:ncol(constraints))
{
  for (j in 1:nrow(constraints))
  {
    if (!is.na(constraints[j,i]) && constraints[j,i] == "N")
    {badarcs = rbind(badarcs, data.frame(from=constraints[j,'Variables'], to=names(constraints)[i]) ) }
  }
} #loading constraints to badarcs 

badarcs = data.frame(badarcs) # converting from table format to data frame format

write.csv(badarcs,"badarcs.csv")

# data prep for modelling ------------ 

#There are two cases to convert to discrete variables.
# Case 1, if there no continuous variables at imagery level .
# Case 2, if there is continuous variables at imagery level,MDS studies fall in this.

colnames(raw_data)[!(apply(raw_data,2,is.discrete))]

# above code will give name of continuous variables at imagery level, choose case accordingly. 
factor_names<-colnames(factor_data)

#Case 1 ---------------

# case 1 starts
for (i in (1:length(factor_data))){
  
  factor_data[,i] = ifelse((factor_data[,i] %in% c(0,1)),factor_data[,i],2)
  
}   
#case 1 ends

# Case 2 ------------- 

#case 2 starts
MDS<-c("fMEANINGFUL","fDIFFERENT","fSALIENT","fPOWER") 
#*** Important modify according to project;


Ex_MDS<-setdiff(factor_names,MDS) # segregates non MDS variables


for (i in MDS){
  
  factor_data[,i]<-as.numeric(ggplot2::cut_number(factor_data[,i],3))
  
}  # convert continuous to discrete MDS and Power variable

for (i in Ex_MDS){
  
  factor_data[,i]<- ifelse((factor_data[,i] %in% c(0,1)),factor_data[,i],2)
}  # convert continuous to discrete MDS and Power variable

#case 2 ends

# data back up and factorizing data ---------------

data_bkp = data.frame(factor_data) #data_backup

write.csv(colMeans(factor_data,na.rm = TRUE), "var_means.csv")
nrow(factor_data[is.na(factor_data),]) #checking for missing values
sum(!complete.cases(factor_data))

factor_data <- na.omit(factor_data) #removing missing values if any


for(i in factor_names)
{
  factor_data[,i] = as.factor(factor_data[,i])
} # adding factor levels to all variables

#modelling ---------------

cl = makeCluster(detectCores()) #  Find total # of cores in your PC or VM using detectCores()

start_time <- Sys.time()
parallel::clusterEvalQ(cl, set.seed(1))
set.seed(1)
boot50 = boot.strength(data = factor_data, R = 50, algorithm = "hc", cluster = cl,
                       algorithm.args = list( score="bde", iss = 60 , restart=10, perturb = 5 
                                              , blacklist = badarcs ))

end_time <- Sys.time()
print( end_time - start_time)    
stopCluster(cl)

write.csv(boot50,"boot50.csv")

# choosing solution and getting outputs -------------

cutoff<- seq(0,1,.05) #create cutoff from 0% to 100 % with 5 % increment

paths =c() 

for ( i  in cutoff ) {
  
  
  x <- nrow (boot50[(boot50$strength >= i) & (boot50$direction > 0.5),   ])
  
  paths <- c(paths,x)
} 

(rbind(cutoff,paths)) #shows cutoff vs no of paths - choose cutoff for strength depending upon paths

Solution1 = boot50[(boot50$strength >= 0.3) & (boot50$direction > 0.5),   ] 
#select the strength depending upon the path

avg.boot1 = averaged.network(Solution1)

write.csv(directed.arcs(avg.boot1), "avgbootDirected_Arcs.csv") #gives the output in excel

write.csv(xtabs(~.,Solution1[,1:2]),"paths to be tested.csv")

