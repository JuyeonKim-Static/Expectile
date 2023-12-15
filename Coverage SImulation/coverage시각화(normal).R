install.packages("readxl")
library(readxl)
library(ggplot2)

#normal값에 대한 coverge시각화

#n=100_rho=0.3
df<- read_excel("coverage시각화.xlsx" , sheet = "n=100_rho=0.3") 
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 100 , rho = 0.3")

#n=500_rho=0.3
df<- read_excel("coverage시각화.xlsx" , sheet = "n=500_rho=0.3")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 500 , rho = 0.3")


#n=1000_rho=0.3
df<- read_excel("coverage시각화.xlsx" , sheet = "n=1000_rho=0.3")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 1000 , rho = 0.3")

#n=100_rho=0.5
df<- read_excel("coverage시각화.xlsx" , sheet = "n=100_rho=0.5")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 100 , rho = 0.5")


#n=500_rho=0.5
df<- read_excel("coverage시각화.xlsx" , sheet = "n=500_rho=0.5")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 500 , rho = 0.5")


#n=1000_rho=0.5
df<- read_excel("coverage시각화.xlsx" , sheet = "n=1000_rho=0.5")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 1000 , rho = 0.5")


#n=100_rho=0.7
df<- read_excel("coverage시각화.xlsx" , sheet = "n=100_rho=0.7")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 100 , rho = 0.7")


#n=500_rho=0.7
df<- read_excel("coverage시각화.xlsx" , sheet = "n=500_rho=0.7")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 500 , rho = 0.7")

#n=1000_rho=0.7
df<- read_excel("coverage시각화.xlsx" , sheet = "n=1000_rho=0.7")
Model=c(rep("quantile",9) ,rep("expectile",9))

ggplot(df, aes(x=tau, y=coverage_mean, colour= Model)) + 
  geom_line(size = 1) +   # line graph를 추가합니다. 
  geom_point(size=3, fill="white") +
  scale_color_manual(values = c("blue","red"))+
  labs(title ="n = 1000 , rho = 0.7"
