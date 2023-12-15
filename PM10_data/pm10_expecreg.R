library(tidyverse)
library(quantreg)

expec.reg <- function(x, y, tau, M = 1e+03){
  if(length(x) != length(y)){
    stop("Lengths of x and y are different")
  }
  n <- length(x)
  weight.new <- rep(tau, n); weight <- vector(length = n)
  iter <- 1
  while(sum((weight-weight.new)^2)!=0){
    weight <- weight.new
    #print(paste(iter, "th iteration"))
    lmfit <- lm(y ~ x, weights = weight)
    weight.new <- tau * (y > lmfit$fitted.values) + (1-tau) * (y <= lmfit$fitted.values)
    iter <- iter + 1 
    if(iter == M){
      warning("procedure does not converge"); 
      break; 
    }
  }
  # if(!is.na(coeff[1])){
  coeff <- lmfit$coefficients
  names(coeff) <- NULL
  # }
  return(coeff)
}

expec.corr <- function(x, y, tau, M = 1e+03){
  x.new <- x[!is.na(x) & !is.na(y)]
  y.new <- y[!is.na(x) & !is.na(y)]
  x <- x.new
  y <- y.new
  
  if(length(x) != length(y)){
    stop("Lengths of x and y are different")
  }
  
  theta_yx <- expec.reg(x, y, tau = tau, M = M)
  beta_yx  <- theta_yx[2]
  alpha_yx <- theta_yx[1]
  theta_xy <- expec.reg(y, x, tau = tau, M = M)
  beta_xy  <- theta_xy[2]
  alpha_xy <- theta_xy[1]
  return( sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0)) )
}

quant.corr <- function(x, y, tau){
  x.new <- x[!is.na(x) & !is.na(y)]
  y.new <- y[!is.na(x) & !is.na(y)]
  x <- x.new
  y <- y.new
  
  if(length(x) != length(y)){
    stop("Lengths of x and y are different")
  }
  
  beta_yx <- rq(y ~ x, tau = tau)$coeff[2]
  beta_xy <- rq(x ~ y, tau = tau)$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  return(quant_corr_tau)
}

pm10 <- read.csv("Daily/2009.csv", fileEncoding = 'euc-kr')
pm10
for(year in 2010:2022){
  pm10_temp <- read.csv(paste("Daily/", year, ".csv", sep = ""), fileEncoding = 'euc-kr')
  pm10 <- rbind(pm10, pm10_temp)
}
head(pm10)

colnames(pm10) <- c("LocNo", "LocName", "Time", "PM10")

pm10 <- as_tibble(pm10)
pm10
df_pm10 <- pm10 %>% separate(Time, c("Year", "Month", "Date"), "\\.")
df_pm10
LocNo_vec <- df_pm10 %>% select(LocNo) %>% unique()
LocNo_vec
LocNo_100  <- df_pm10 %>% filter(LocNo ==  100) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_108 <- df_pm10 %>% filter(LocNo == 108) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_119 <- df_pm10 %>% filter(LocNo == 119) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_135 <- df_pm10 %>% filter(LocNo == 135) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_136 <- df_pm10 %>% filter(LocNo == 136) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_143 <- df_pm10 %>% filter(LocNo == 143) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_146 <- df_pm10 %>% filter(LocNo == 146) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_152 <- df_pm10 %>% filter(LocNo == 152) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_156 <- df_pm10 %>% filter(LocNo == 156) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_160 <- df_pm10 %>% filter(LocNo == 160) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_175 <- df_pm10 %>% filter(LocNo == 175) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_185 <- df_pm10 %>% filter(LocNo == 185) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_192 <- df_pm10 %>% filter(LocNo == 192) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_201 <- df_pm10 %>% filter(LocNo == 201) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))
LocNo_232 <- df_pm10 %>% filter(LocNo == 232) %>% select(Year, Month, Date, PM10) %>% mutate(logPM10 = log(PM10))



#일자별을 기준으로 데이터프레임을 바꿈.
LocNo_all <- LocNo_100 %>%
  merge(LocNo_108, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_100", "_108")) %>%
  merge(LocNo_119, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_135, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_119", "_135")) %>%
  merge(LocNo_136, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_143, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_136", "_143")) %>%
  merge(LocNo_146, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_152, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_146", "_152")) %>%
  merge(LocNo_156, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_160, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_156", "_160")) %>%
  merge(LocNo_175, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_185, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_175", "_185")) %>%
  merge(LocNo_192, by = c("Year", "Month", "Date"), all = TRUE) %>%
  merge(LocNo_201, by = c("Year", "Month", "Date"), all = TRUE, suffix = c("_192", "_201"))


# 결과 확인
LocNo_all
LocNo_all <- LocNo_all %>%
  merge(LocNo_232, by = c("Year", "Month", "Date"), all = TRUE)
LocNo_all <- rename(LocNo_all, PM10_232 = PM10, logPM10_232 = logPM10)
head(LocNo_all)

expec_corr_100 <- expec_corr_108 <- expec_corr_119 <- expec_corr_135 <- expec_corr_136 <- expec_corr_143 <- vector(length = 9)
expec_corr_146 <- expec_corr_152 <- expec_corr_156 <- expec_corr_160 <- expec_corr_175 <- expec_corr_185 <- vector(length = 9)
expec_corr_192 <- expec_corr_201 <- expec_corr_232 <- vector(length = 9)
quant_corr_100 <- quant_corr_108 <- quant_corr_119 <- quant_corr_135 <- quant_corr_136 <- quant_corr_143 <- vector(length = 9)
quant_corr_146 <- quant_corr_152 <- quant_corr_156 <- quant_corr_160 <- quant_corr_175 <- quant_corr_185 <- vector(length = 9)
quant_corr_192 <- quant_corr_201 <- quant_corr_232 <- vector(length = 9)


mean_pm10 <- df_pm10 %>% group_by(LocNo, Month, Date) %>% summarize(mean = mean(PM10), logmean = mean(log(PM10)))

anom_100 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo ==  "100"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_100 - mean, anom.log = logPM10_100 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_108 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "108"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_108 - mean, anom.log = logPM10_108 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_119 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "119"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_119- mean, anom.log = logPM10_119 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_135 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "135"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_135 - mean, anom.log = logPM10_135 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_136 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "136"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_136 - mean, anom.log = logPM10_136 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_143 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "143"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_143 - mean, anom.log = logPM10_143 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_146 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "146"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_146 - mean, anom.log = logPM10_146 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_152 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "152"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_152 - mean, anom.log = logPM10_152 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_156 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "156"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_156 - mean, anom.log = logPM10_156 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_160 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "160"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_160 - mean, anom.log = logPM10_160 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_175 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "175"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_175 - mean, anom.log = logPM10_175 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_185 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "185"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_185 - mean, anom.log = logPM10_185 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_192 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "192"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_192 - mean, anom.log = logPM10_192 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_201 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "201"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_201 - mean, anom.log = logPM10_201 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_232 <- LocNo_all %>% merge(mean_pm10 %>% filter(LocNo == "232"), by = c("Month", "Date")) %>% 
  mutate(anom = PM10_232 - mean, anom.log = logPM10_232 - logmean) %>%
  arrange(Year, Month, Date) %>% select(Year, Month, Date, anom, anom.log)
anom_232
for(i in 1:9){
  expec_corr_100[i] <- expec.corr(anom_100$anom, anom_108$anom, tau = 0.1*i) #강원
  expec_corr_108[i] <- expec.corr(anom_108$anom, anom_108$anom, tau = 0.1*i) # 서울
  expec_corr_119[i] <- expec.corr(anom_119$anom, anom_108$anom, tau = 0.1*i) #경기
  expec_corr_135[i] <- expec.corr(anom_135$anom, anom_108$anom, tau = 0.1*i) # 충청북도
  expec_corr_136[i] <- expec.corr(anom_136$anom, anom_108$anom, tau = 0.1*i) # 경상북도
  expec_corr_143[i] <- expec.corr(anom_143$anom, anom_108$anom, tau = 0.1*i) # 대구
  expec_corr_146[i] <- expec.corr(anom_146$anom, anom_108$anom, tau = 0.1*i) # 전라북도
  expec_corr_152[i] <- expec.corr(anom_152$anom, anom_108$anom, tau = 0.1*i) # 울산
  expec_corr_156[i] <- expec.corr(anom_156$anom, anom_108$anom, tau = 0.1*i) # 광주
  expec_corr_160[i] <- expec.corr(anom_160$anom, anom_108$anom, tau = 0.1*i) # 부산
  expec_corr_175[i] <- expec.corr(anom_175$anom, anom_108$anom, tau = 0.1*i) # 전라남도
  expec_corr_185[i] <- expec.corr(anom_185$anom, anom_108$anom, tau = 0.1*i) # 제주도
  expec_corr_192[i] <- expec.corr(anom_192$anom, anom_108$anom, tau = 0.1*i) # 경상남도
  expec_corr_201[i] <- expec.corr(anom_201$anom, anom_108$anom, tau = 0.1*i) # 인천
  expec_corr_232[i] <- expec.corr(anom_232$anom, anom_108$anom, tau = 0.1*i) # 충청남도
  quant_corr_100[i] <- quant.corr(anom_100$anom, anom_108$anom, tau = 0.1*i) #
  quant_corr_108[i] <- quant.corr(anom_108$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_119[i] <- quant.corr(anom_119$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_135[i] <- quant.corr(anom_135$anom, anom_108$anom, tau = 0.1*i) #
  quant_corr_136[i] <- quant.corr(anom_136$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_143[i] <- quant.corr(anom_143$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_146[i] <- quant.corr(anom_146$anom, anom_108$anom, tau = 0.1*i) #
  quant_corr_152[i] <- quant.corr(anom_152$anom, anom_108$anom, tau = 0.1*i) #
  quant_corr_156[i] <- quant.corr(anom_156$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_160[i] <- quant.corr(anom_160$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_175[i] <- quant.corr(anom_175$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_185[i] <- quant.corr(anom_185$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_192[i] <- quant.corr(anom_192$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_201[i] <- quant.corr(anom_201$anom, anom_108$anom, tau = 0.1*i) # 
  quant_corr_232[i] <- quant.corr(anom_232$anom, anom_108$anom, tau = 0.1*i) # 
}
expec_corr_df <- rbind(expec_corr_100,
                       expec_corr_108,
                       expec_corr_119,
                       expec_corr_135, 
                       expec_corr_136,
                       expec_corr_143,
                       expec_corr_146, 
                       expec_corr_152, 
                       expec_corr_156,
                       expec_corr_160, 
                       expec_corr_175, 
                       expec_corr_185, 
                       expec_corr_192,
                       expec_corr_201,
                       expec_corr_232)

quant_corr_df <- rbind(quant_corr_100,
                       quant_corr_108,
                       quant_corr_119,
                       quant_corr_135, 
                       quant_corr_136,
                       quant_corr_143,
                       quant_corr_146, 
                       quant_corr_152, 
                       quant_corr_156,
                       quant_corr_160, 
                       quant_corr_175, 
                       quant_corr_185, 
                       quant_corr_192,
                       quant_corr_201,
                       quant_corr_232)


column_names <- seq(0.1, 0.9, by = 0.1)
colnames(expec_corr_df) <- column_names
row_names  <- c( "Gangwon", "Seoul","Gyeonggi", "Chungcheongbuk", "Gyeongsangbuk", "Daegu",
                 "Jeollabuk", "Ulsan", "Gwangju", "Busan", "Jeollanam", "Jeju",
                 "Gyeongsangnam", "Incheon","Chungcheongnam")
rownames(expec_corr_df) <-row_names
expec_corr_df
write.csv(expec_corr_df, file = "expec_corr_seoul.csv")

rownames(quant_corr_df) <-row_names
colnames(quant_corr_df) <- column_names
quant_corr_df
write.csv(quant_corr_df, file = "quant_corr_seoul.csv")



plot(c(1:9)/10,  expec_corr_100, type = "o", col = 2, ylim = c(0.35, 1))
lines(c(1:9)/10, expec_corr_108, type = "o", col = 1)
lines(c(1:9)/10, expec_corr_119, type = "o", col = 3)
lines(c(1:9)/10, expec_corr_135, type = "o", col = 4)
lines(c(1:9)/10, expec_corr_136, type = "o", col = 5)
lines(c(1:9)/10, expec_corr_143, type = "o", col = 6)
lines(c(1:9)/10, expec_corr_146, type = "o", col = 7)
lines(c(1:9)/10, expec_corr_152, type = "o", col = 8)
lines(c(1:9)/10, expec_corr_156, type = "o", col = 9)
lines(c(1:9)/10, expec_corr_160, type = "o", col = 10)
lines(c(1:9)/10, expec_corr_175, type = "o", col = 11)
lines(c(1:9)/10, expec_corr_185, type = "o", col = 12)
lines(c(1:9)/10, expec_corr_192, type = "o", col = 13)
lines(c(1:9)/10, expec_corr_201, type = "o", col = 14)

plot(c(1:9)/10,  quant_corr_100, type = "o", lty = 2, col = 2, ylim = c(0.35, 1))
lines(c(1:9)/10, quant_corr_108, type = "o", lty = 2,col = 1)
lines(c(1:9)/10, quant_corr_119, type = "o", lty = 2,col = 3)
lines(c(1:9)/10, quant_corr_135, type = "o", lty = 2, col = 4)
lines(c(1:9)/10, quant_corr_136, type = "o", lty = 2,col = 5)
lines(c(1:9)/10, quant_corr_143, type = "o", lty = 2, col = 6)
lines(c(1:9)/10, quant_corr_146, type = "o", lty = 2,col = 7)
lines(c(1:9)/10, quant_corr_152, type = "o", lty = 2,col = 8)
lines(c(1:9)/10, quant_corr_156, type = "o", lty = 2,col = 9)
lines(c(1:9)/10, quant_corr_160, type = "o", lty = 2,col = 10)
lines(c(1:9)/10, quant_corr_175, type = "o", lty = 2,col = 11)
lines(c(1:9)/10, quant_corr_185, type = "o", lty = 2,col = 12)
lines(c(1:9)/10, quant_corr_192, type = "o", lty = 2,col = 13)
lines(c(1:9)/10, quant_corr_201, type = "o", lty = 2,col = 14)



labels <- c("Seoul", "Gangwon", "Gyeonggi", "Chungcheongbuk", "Gyeongsangbuk", "Daegu",
            "Jeollabuk", "Ulsan", "Gwangju", "Busan", "Jeollanam", "Jeju",
            "Gyeongsangnam", "Ganghwa")

# 범례 추가
legend("bottomleft" ,legend = labels, lty = 1, col = 1:14, ncol = 1, cex = 0.8)

library(xts)

par(mfrow = c(1, 1))

anom_108_date <- (anom_108 %>% transmute(date = as.POSIXct(paste0(Year, "-", Month, "-", Date))))$date
xts(anom_108$anom, order.by = anom_108_date) %>% plot()

anom_232_date <- (anom_232 %>% transmute(date = as.POSIXct(paste0(Year, "-", Month, "-", Date))))$date
xts(anom_232$anom, order.by = anom_232_date) %>% lines(col = "blue")

