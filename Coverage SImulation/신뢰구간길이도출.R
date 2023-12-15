install.packages("openxlsx")
library("openxlsx")

result_obj <- readRDS("coverage_result_DGP3(multi-t(3)).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()
#17
result_obj[[17]]$description
simulation1<-data.frame(result_obj[[17]]$obj)
colnames(simulation1)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation1$quant_ci <- 2* 1.96*simulation1$sigma_quant/sqrt(500)
simulation1$expec_ci <- 2* 1.96*simulation1$sigma_expec/sqrt(500)
write.xlsx(simulation1, file="simulation1.csv")

#23
result_obj[[23]]$description
simulation2<-data.frame(result_obj[[23]]$obj)
colnames(simulation2)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation2$quant_ci <- 2* 1.96*simulation2$sigma_quant/sqrt(500)
simulation2$expec_ci <- 2* 1.96*simulation2$sigma_expec/sqrt(500)
write.xlsx(simulation2, file="simulation2.csv")

#29
result_obj[[29]]$description
simulation3<-data.frame(result_obj[[29]]$obj)
colnames(simulation3)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation3$quant_ci <- 2* 1.96*simulation3$sigma_quant/sqrt(500)
simulation3$expec_ci <- 2* 1.96*simulation3$sigma_expec/sqrt(500)
write.xlsx(simulation3, file="simulation3.csv")

#t(5)
result_obj <- readRDS("coverage_result_DGP3(multi-t(5)).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()
#17
result_obj[[17]]$description
simulation4<-data.frame(result_obj[[17]]$obj)
colnames(simulation4)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation4$quant_ci <- 2* 1.96*simulation4$sigma_quant/sqrt(500)
simulation4$expec_ci <- 2* 1.96*simulation4$sigma_expec/sqrt(500)
write.xlsx(simulation4, file="simulation4.csv")

#23
result_obj[[23]]$description
simulation5<-data.frame(result_obj[[23]]$obj)
colnames(simulation5)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation5$quant_ci <- 2* 1.96*simulation5$sigma_quant/sqrt(500)
simulation5$expec_ci <- 2* 1.96*simulation5$sigma_expec/sqrt(500)
write.xlsx(simulation5, file="simulation5.csv")

#29
result_obj[[29]]$description
simulation6<-data.frame(result_obj[[29]]$obj)
colnames(simulation6)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation6$quant_ci <- 2* 1.96*simulation6$sigma_quant/sqrt(500)
simulation6$expec_ci <- 2* 1.96*simulation6$sigma_expec/sqrt(500)
write.xlsx(simulation6, file="simulation6.csv")


#t(7)
result_obj <- readRDS("coverage_result_DGP3(multi-t(7)).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()
#17
result_obj[[17]]$description
simulation7<-data.frame(result_obj[[17]]$obj)
colnames(simulation7)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation7$quant_ci <- 2* 1.96*simulation7$sigma_quant/sqrt(500)
simulation7$expec_ci <- 2* 1.96*simulation7$sigma_expec/sqrt(500)
write.xlsx(simulation7, file="simulation7.csv")

#23
result_obj[[23]]$description
simulation8<-data.frame(result_obj[[23]]$obj)
colnames(simulation8)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation8$quant_ci <- 2* 1.96*simulation8$sigma_quant/sqrt(500)
simulation8$expec_ci <- 2* 1.96*simulation8$sigma_expec/sqrt(500)
write.xlsx(simulation8, file="simulation8.csv")

#29
result_obj[[29]]$description
simulation9<-data.frame(result_obj[[29]]$obj)
colnames(simulation9)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation9$quant_ci <- 2* 1.96*simulation9$sigma_quant/sqrt(500)
simulation9$expec_ci <- 2* 1.96*simulation9$sigma_expec/sqrt(500)
write.xlsx(simulation9, file="simulation9.csv")

result_obj <- readRDS("coverage_result_DGP3(multi-t(7)).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()
#17
result_obj[[17]]$description
simulation7<-data.frame(result_obj[[17]]$obj)
colnames(simulation7)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation7$quant_ci <- 2* 1.96*simulation7$sigma_quant/sqrt(500)
simulation7$expec_ci <- 2* 1.96*simulation7$sigma_expec/sqrt(500)
simulation7
write.xlsx(simulation7, file="simulation7.csv")
#23
result_obj[[23]]$description
simulation8<-data.frame(result_obj[[23]]$obj)
colnames(simulation8)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation8$quant_ci <- 2* 1.96*simulation8$sigma_quant/sqrt(500)
simulation8$expec_ci <- 2* 1.96*simulation8$sigma_expec/sqrt(500)
write.xlsx(simulation8, file="simulation8.csv")

#29
result_obj[[29]]$description
simulation9<-data.frame(result_obj[[29]]$obj)
colnames(simulation9)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation9$quant_ci <- 2* 1.96*simulation9$sigma_quant/sqrt(500)
simulation9$expec_ci <- 2* 1.96*simulation9$sigma_expec/sqrt(500)
simulation9
write.xlsx(simulation9, file="simulation9.csv")

#normal
result_obj <- readRDS("coverage_result_DGP1(normal).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()

#17
result_obj[[17]]$description
simulation10<-data.frame(result_obj[[17]]$obj)
colnames(simulation10)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation10$quant_ci <- 2* 1.96*simulation10$sigma_quant/sqrt(500)
simulation10$expec_ci <- 2* 1.96*simulation10$sigma_expec/sqrt(500)
write.xlsx(simulation10, file="simulation10.csv")

#23
result_obj[[23]]$description
simulation11<-data.frame(result_obj[[23]]$obj)
colnames(simulation11)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation11$quant_ci <- 2* 1.96*simulation11$sigma_quant/sqrt(500)
simulation11$expec_ci <- 2* 1.96*simulation11$sigma_expec/sqrt(500)
write.xlsx(simulation11, file="simulation11.csv")

#29
result_obj[[29]]$description
simulation12<-data.frame(result_obj[[29]]$obj)
colnames(simulation12)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation12$quant_ci <- 2* 1.96*simulation12$sigma_quant/sqrt(500)
simulation12$expec_ci <- 2* 1.96*simulation12$sigma_expec/sqrt(500)
write.xlsx(simulation12, file="simulation12.csv")


#rocket Type
result_obj <- readRDS("coverage_result_DGP2(multi-Rocket type).RDS")
n <- length(result_obj)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()

#11 n=500 ,tau = 0.1 , rho =0.3
result_obj[[2]]$description
simulation14<-data.frame(result_obj[[2]]$obj)
colnames(simulation14)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation14$quant_ci <- 2* 1.96*simulation14$sigma_quant/sqrt(500)
simulation14$expec_ci <- 2* 1.96*simulation14$sigma_expec/sqrt(500)
write.xlsx(simulation14, file="simulation14.csv")



#14n=500 ,tau = 0.5 , rho =0.3
result_obj[[14]]$description
simulation13<-data.frame(result_obj[[14]]$obj)
colnames(simulation13)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation13$quant_ci <- 2* 1.96*simulation13$sigma_quant/sqrt(500)
simulation13$expec_ci <- 2* 1.96*simulation13$sigma_expec/sqrt(500)
write.xlsx(simulation13, file="simulation13.csv")


#11 n=500 ,tau = 0.9 , rho =0.3
result_obj[[26]]$description
simulation15<-data.frame(result_obj[[26]]$obj)
colnames(simulation15)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

simulation15$quant_ci <- 2* 1.96*simulation15$sigma_quant/sqrt(500)
simulation15$expec_ci <- 2* 1.96*simulation15$sigma_expec/sqrt(500)
write.xlsx(simulation15, file="simulation15.csv")
