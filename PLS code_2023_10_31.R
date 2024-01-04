# caveats for using following code:
# > results will match with Smart PLS 4, not with Smart PLS 3
# > seminr package doesn't support weighted PLS. 

# install packages for PLS ------------

install.packages('writexl') #for writing excel file
install.packages("seminr") #for PLS modelling
install.packages("dplyr") # for creating paths

# loading packages for PLS ---------

library(writexl)
library(seminr)
library(dplyr)

# defining factor structure -------------

factor_structure <-c()

for (i in 2:dim(factor_loadings)[2]){
  
  factor_structure<-c(
    constructs(composite(colnames(factor_loadings)[i],
                         factor_loadings[,1][factor_loadings[,i]!=0],
                         weights = mode_A)),
    factor_structure)
}   # creates factors with corresponding imageries


# defining paths for PLS ----------------

itr_0<-paths(Solution1[1,1],Solution1[1,2])

for (i in 2:dim(Solution1)[1]) {
  
  itr_0<- c(paths(Solution1[i,1],Solution1[i,2]),itr_0)
} # imports paths from Bayesian algorithm Solution1.

testing_paths <- relationships(itr_0)

# initial model with paths from Bayesian algo. ------------------

model_itr_1 <- estimate_pls(data = raw_data,
                            measurement_model = factor_structure, 
                            structural_model =testing_paths,
                            inner_weights = path_weighting)


path_coef1<-model_itr_1$path_coef
save_itr1<-as.data.frame(cbind(" "=row.names(path_coef1),path_coef1))

x<-as.data.frame(which(path_coef1<0,arr.ind = T))

# removing negative paths -------------

removed_paths <- as.data.frame(cbind(rownames(path_coef1)[x$row],
                       colnames(path_coef1)[x$col],
                       path_coef1[path_coef1<0]))

while (dim(x)[1]>0){
negative_paths<-relationships(as.matrix(rbind(rownames(path_coef1)[x$row],
                                              colnames(path_coef1)[x$col])))
testing_paths<-(setdiff(as.data.frame(testing_paths),
        as.data.frame(negative_paths)))

model_itr_2<- estimate_pls(data = raw_data,
                            measurement_model = factor_structure,
                            structural_model = testing_paths,
                            inner_weights = path_weighting)

path_coef1<-model_itr_2$path_coef

x<-as.data.frame(which(path_coef1<0,arr.ind = T))

removed_paths <- rbind(removed_paths,(cbind(rownames(path_coef1)[x$row],
                                        colnames(path_coef1)[x$col],
                                        path_coef1[path_coef1<0])) )
}

save_itr2<-as.data.frame(cbind(" "=row.names(path_coef1),path_coef1))

removed_negative_paths<-removed_paths
colnames(removed_negative_paths)<-c("source","target","value")


rm(removed_paths)

removed_paths<-data.frame()

# removing paths based on cut-off--------------

path_cutoff<-0.09# give a path co-eff cut-off, default is .09. 

y<-as.data.frame(cbind(which((path_coef1<path_cutoff) & (path_coef1!=0),arr.ind = T),
                       path_coef1[which((path_coef1<path_cutoff) & (path_coef1!=0),arr.ind = T)]))

colnames(y)[3]<-"value"



while (dim(y)[1]>0){
  
  lp<-y %>% group_by(col) %>% slice(which.min(value))
  
  least_path<-(as.matrix(cbind(rownames(path_coef1)[lp$row],
                                            colnames(path_coef1)[lp$col]))) 
  colnames(least_path) <- colnames(testing_paths)
  
  removed_paths<-rbind(removed_paths,cbind((least_path),lp$value))
  
  testing_paths<-(setdiff(as.data.frame(testing_paths),
                      as.data.frame(least_path)))
  
  model_itr_3 <- estimate_pls(data = raw_data,
                              measurement_model = factor_structure, 
                              structural_model =testing_paths,
                              inner_weights = path_weighting)
  
  path_coef1<-model_itr_3$path_coef
  
  y<-as.data.frame(cbind(which((path_coef1<path_cutoff) & (path_coef1!=0),arr.ind = T),
                         path_coef1[which((path_coef1<path_cutoff) & (path_coef1!=0),arr.ind = T)]))
  colnames(y)[3]<-"value"
  
} # iteratively removes least path, from each variable until all paths are above cut-off.  

colnames(removed_paths)[3]<-"value"
removed_paths<-rbind(removed_negative_paths,removed_paths)
removed_paths[,3]<-as.numeric(removed_paths[,3])
tested_values <- as.data.frame(unclass(xtabs(value~source+target,data=removed_paths)))
tested_values<-cbind(" "=row.names(tested_values),tested_values)

save_itr3<- as.data.frame(cbind(" "=row.names(path_coef1),path_coef1))


#saving the outputs ----------------
write_xlsx(list(itr1=save_itr1,
                itr2=save_itr2,
                itr3=save_itr3,
                testing_sheet=tested_values),
           "output.xlsx") 
# creates excel sheet with 
# > iteration 1 initial model
# > iteration 2 with negative paths removed
# > iteration 3 final model co-eff
# > paths below cut-off co-eff (testing sheet) and negative paths

