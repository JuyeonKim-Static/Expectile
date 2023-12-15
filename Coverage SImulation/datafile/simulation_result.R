ibrary("openxlsx")
result_obj <- readRDS("coverage_result_DGP1(normal).RDS")
n <- length(result_obj)
n
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
result_vec <- c()
for (i in 1:n){
  simulation1<-data.frame(result_obj[[i]]$obj)
  colnames(simulation1)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")
  result_vec <- rbind(result_vec, 
                      c(apply(simulation1, 2, mean), 
                        apply(simulation1, 2, sd)))
  # write.xlsx(simulation1, sheet = "sheet1",file="simulation1.xlsx")
}
colnames(result_vec) <- c("quant_corr_mean","sigma_quant_mean","coverage_quant_mean" ,"expec_corr_mean" ,"sigma_expec_mean","coverage_expec_mean", "quant_corr_sd","sigma_quant_sd","coverage_quant_sd" ,"expec_corr_sd" ,"sigma_expec_sd","coverage_expec_sd")
write.csv(result_vec, "ddd.csv")

result_obj2 <- readRDS("coverage_result_DGP3(multi-t(3)).RDS")
n <- length(result_obj2)
for(i in 1:n){
  cat(result_obj[[i]]$description, "\n")
}
simulation2_t3<-data.frame(result_obj2[[1]]$obj)
colnames(simulation2_t3)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

write.xlsx(simulation2_t3, file="simulation2_t3.xlsx")


result_obj3 <- readRDS("coverage_result_DGP3(multi-t(5)).RDS")
n <- length(result_obj3)
for(i in 1:n){
  cat(result_obj3[[i]]$description, "\n")
}
simulation2_t5<-data.frame(result_obj3[[1]]$obj)
colnames(simulation2_t5)<-c( "quant_corr","sigma_quant","coverage_quant" ,"expec_corr" ,"sigma_expec","coverage_expec")

write.xlsx(simulation2_t5, file="simulation2_t5.xlsx")

