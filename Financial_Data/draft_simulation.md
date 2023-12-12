rolling_window_simulation
================
Juyeon Kim
2023-12-12

``` r
#Select df
df1<-read.csv("S&P_DAX(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
df2<-read.csv("S&P_FTSE(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
df3<-read.csv("S&P_CAC(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
df4<-read.csv("S&P_KOSPI(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
```

``` r
#expectile regreesion function
expec.reg <- function(x, y, tau){
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
    if(iter == 500){
      warning("procedure does not converge"); 
      break; 
    }
    # if(!is.na(coeff[1])){
    coeff <- lmfit$coefficients
    names(coeff) <- NULL
    # }
  }
  return(coeff)
}
```

\##Estimate CI

``` r
library(quantreg)
```

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

``` r
library(zoo)
```

    ## 
    ## Attaching package: 'zoo'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.Date, as.Date.numeric

``` r
library(latex2exp)
library(hdrcde)
```

    ## This is hdrcde 3.4

``` r
tau <-seq( from = 0.05, to = 0.95, by=0.01 )  
df<-df1 #for to change df
x<-df[,2]
y<- df[,3]

df_quantile <-c()
df_expectile<-c()
sigma_tau_quant<- c()

for (a in 1: length(tau)){
  alpha_yx<-  rq(y ~ x, tau = tau[a])$coeff[1]
  alpha_xy<-  rq(x ~ y, tau = tau[a])$coeff[1]
  beta_yx <-  rq(y ~ x, tau = tau[a])$coeff[2]
  beta_xy <-  rq(x ~ y, tau = tau[a])$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  df_quantile[a] <- quant_corr_tau}

n <- length(x)
#밀도함수 추정
#y|x분포
f_1 =vector(length = n)
for(i in 1:n){
  new_cde1 <- cde(x, y, x.margin = x[i]) #, y.margin = y_margin)
  y_margin <- new_cde1$y
  res <- y_margin - (alpha_yx + beta_yx * x[i])
  where_y <- which.min(sapply(res, function(x){ifelse(x<0, Inf, x)}))
  if(where_y == 1){
    f_1[i] <- new_cde1$z[where_y]
  }else{
    f_1[i]<- (abs(res[where_y-1]) * new_cde1$z[where_y] + abs(res[where_y]) * new_cde1$z[where_y-1]) / 
      (abs(res[where_y-1]) + abs(res[where_y]))
  }
}

#x|y분포
f_2 =vector(length = n)
for (i in 1:n){
  new_cde2 <- cde(y, x, x.margin = y[i]) # , y.margin = x_margin)
  x_margin <- new_cde2$y
  res <- x_margin - (alpha_xy + beta_xy * y[i])
  where_x <- which.min(sapply(res, function(x){ifelse(x<0, Inf, x)}))
  if(where_x == 1){
    f_2[i] <- new_cde2$z[where_x]
  }else{
    f_2[i]<- (abs(res[where_x-1]) * new_cde2$z[where_x] + abs(res[where_x]) * new_cde2$z[where_x-1]) / 
      (abs(res[where_x-1]) + abs(res[where_x]))
  }
}
#R1
a1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
for (i in c(1:n)){
  a1[1,1,i] <- 1
  a1[1,2,i] <- x[i]
  a1[2,1,i] <- x[i]
  a1[2,2,i] <- (x[i])^2
  a1[ , ,i] <- f_1[i] * a1[, ,i]
}
R1 <- cbind(apply(a1, c(1,2), mean), matrix(0, nrow = 2, ncol = 2))
#R2
b1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
for (i in c(1:n)){
  b1[1,1,i] <- 1
  b1[1,2,i] <- y[i]
  b1[2,1,i] <- y[i]
  b1[2,2,i] <- (y[i])^2
  b1[ , ,i] <- f_2[i] * b1[ , ,i]
}

R2 <- cbind(matrix(0, nrow = 2, ncol = 2), apply(b1, c(1,2), mean))

#M행렬
M <- rbind(R1, R2)
#H행렬
H <- matrix(0, ncol = 4, nrow = 4)

#Quantile
for (a in 1: length(tau)){
  alpha_yx<-  rq(y ~ x, tau = tau[a])$coeff[1]
  alpha_xy<-  rq(x ~ y, tau = tau[a])$coeff[1]
  beta_yx <-  rq(y ~ x, tau = tau[a])$coeff[2]
  beta_xy <-  rq(x ~ y, tau = tau[a])$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  df_quantile[a] <- quant_corr_tau
  for (i in 1:n){
    x1 <- c(1,x[i])* (tau[a] - ifelse( (y[i]- (alpha_yx + beta_yx * x[i]))<0 ,1 ,0))
    y1 <- c(1,y[i])* (tau[a] - ifelse( (x[i]- (alpha_xy + beta_xy * y[i]))< 0 ,1 ,0))
    d<-matrix(append(x1,y1) , nrow =1)
    H <- H + (t(d) %*% d)
  }
  H <- H / n
  #G행렬
  G <- 1/2*matrix( c(0, ifelse(beta_yx * beta_xy > 0, sign(beta_yx) * sqrt(beta_xy /beta_yx), 0), 
                     0, ifelse(beta_yx * beta_xy > 0, sign(beta_xy) * sqrt(beta_yx /beta_xy), 0)))
  v <- solve(M, G)
  sigma_tau_quant[a]  <- sqrt(t(v) %*% H %*% v)
}

#Expectile
sigma_tau_expec <-c()
for ( b in 1 :length(tau)){
  theta_yx <- expec.reg(x, y, tau = tau[b])
  beta_yx  <- theta_yx[2]
  alpha_yx <- theta_yx[1]
  theta_xy <- expec.reg(y, x, tau = tau[b])
  beta_xy  <- theta_xy[2]
  alpha_xy <- theta_xy[1]
  
  df_expectile[b] <-  sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  rho_tau_expec <-  sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  e_yx <- y - alpha_yx - beta_yx* x
  e_xy <- x - alpha_xy - beta_xy* y
  w_yx <- tau[b] * (e_yx > 0) + (1-tau[b]) * (e_yx <= 0)
  w_xy <- tau[b] * (e_xy > 0) + (1-tau[b]) * (e_xy <= 0)
  
  g <- rho_tau_expec/2*matrix(c(1/beta_yx, 1/beta_xy), ncol = 1)
  A <- matrix(c(0,0,1,0,0,0,0,1), ncol = 4)
  det_yx <- mean(w_yx)*mean(w_yx*x^2) - (mean(w_yx*x))^2
  det_xy <- mean(w_xy)*mean(w_xy*y^2) - (mean(w_xy*y))^2
  W.inv <- matrix(c(c(mean(w_yx*x^2), -mean(w_yx*x), 0, 0, -mean(w_yx*x), mean(w_yx), 0, 0)/det_yx,
                    c(0, 0, mean(w_xy*y^2), -mean(w_xy*y), 
                      0, 0, -mean(w_xy*y), mean(w_xy))/det_xy),ncol = 4)/2
  
  V <- cov(2*cbind(e_yx*w_yx, x*e_yx*w_yx, e_xy*w_xy, y*e_xy*w_xy))
  sigma_tau_expec[b] <- sqrt(t(g) %*% A %*% W.inv %*% V %*% W.inv %*% t(A) %*% g)
}

quant_upper <- df_quantile + 1.96*sigma_tau_quant/sqrt(n)
quant_lower <- df_quantile - 1.96*sigma_tau_quant/sqrt(n)
expec_upper <- df_expectile + 1.96*sigma_tau_expec/sqrt(n)
expec_lower <-df_expectile - 1.96*sigma_tau_expec/sqrt(n)

pearson <-cor(x,y)
#plot
plot(tau ,df_quantile , type = 'l', col = 'red' , main = "S&P500 and DAX(Quantile VS Expectile )",ylim = c(0.4,0.8),ylab = "p")
polygon(c(tau, rev(tau)), c(quant_upper, rev(quant_lower)), col = "#FFBBBB66", border = NA)
polygon(c(tau, rev(tau)), c(expec_upper, rev(expec_lower)), col = "#BBAAFF66", border = NA)
lines(tau ,df_expectile , type = 'l', col = 'blue')
abline(h = pearson,col="black",lty=2)
legend("topright", legend=c("quantile","expectile") , fill=c("red","blue"))
```

![draft_simulation_files/figure-gfm/unnamed-chunk-3-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/a44dce2344c50957dd4ae6c3b42a0c6f9b9c994d/Financial_Data/image/draft_DAX_corrplot.png))<!-- -->

\##Comparison Between Expectile and Quantile

``` r
#Quantile
library(quantreg)
library(zoo)
library(latex2exp)

df<-df1 #change to other df
df$date<- as.Date(df1$날짜 ,format = "%Y.%m.%d")
x <- df[,2]
y <- df[,3]
f_pearson_rolling_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width = 500  ,
  function(x) cor(x[,1],x[,2]), 
  by.column= FALSE,
  align = "center")

rq(y ~ x, tau = 0.1)$coeff
```

    ## (Intercept)           x 
    ## -0.01188326  0.66406402

``` r
rq(x ~ y, tau = 0.1)$coeff
```

    ##  (Intercept)            y 
    ## -0.009757307  0.545717826

``` r
quant_reg_0.1 <- function(x,y){
    beta_yx <- rq(y ~ x, tau = 0.1)$coeff[2]
    beta_xy <- rq(x ~ y, tau = 0.1)$coeff[2]
    quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
    return(quant_corr_tau)
    }

f_quant0.1_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width = 500 ,
  function(x) quant_reg_0.1(x[,2],x[,1]), 
  by.column=FALSE,
  align = "center")

quant_reg_0.9 <- function(x,y){
  beta_yx <- rq(y ~ x, tau = 0.9)$coeff[2]
  beta_xy <- rq(x ~ y, tau = 0.9)$coeff[2]
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  return(quant_corr_tau)
}


f_quant0.9_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width = 500 ,
  function(x) quant_reg_0.9(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

global_crisis_start_date <- as.Date("2007-08-01")
global_crisis_end_date <- as.Date("2009-03-31")
corona_crisis_start_date <- as.Date("2020-03-09")

#2022 plot
year_seq <- seq(as.Date("2002-02-11"), as.Date("2022-12-31"), by = "1 year")
plot(f_pearson_rolling_win, col="black", lwd=1, type="l", main="S&P 500 and DAX(2001~2022)_Quantile",
     xaxt='n', xlab ="",ylab = TeX("$rho_{tau}"),
     ylim= c(0,0.8))
abline(v = global_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = global_crisis_end_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = corona_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)
lines(f_quant0.1_win, col="red", lwd=1 )
lines(f_quant0.9_win, col="blue", lwd=1)
axis(1, at = year_seq, labels = format(year_seq, "%Y"), cex.axis = 0.6)
text(x = as.Date("2008-01-01") + 150 , y = 0.35, "Global \n financial \n Crisis", col = "black", cex = 0.8)
text(x = as.Date("2019-11-17") + 300 , y = 0.7 , "Corona \n Crisis", col = "black", cex = 0.8)

legend("bottomright",
       c(TeX("\\hat{\\rho}"), TeX("\\hat{$rho_{0.1}}"), TeX("\\hat{$rho_{0.9}}")),
       col = c("black", "red", "blue"),
       lty = c("solid", "solid", "solid"),
       ncol = 3, bty = "n")
```

![draft_simulation_files/figure-gfm/unnamed-chunk-4-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/a44dce2344c50957dd4ae6c3b42a0c6f9b9c994d/Financial_Data/image/draft_DAX_q.png)<!-- -->



``` r
#Expectile
df$date<- as.Date(df$날짜,format = "%Y.%m.%d")
x <- zoo(df[,2])
y <- zoo(df[,3])

f_pearson_rolling_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width =500 ,
  function(x) cor(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

expec_reg_0.1 <- function(x,y){
  theta_yx <- expec.reg(x, y, tau = 0.05)
  beta_yx  <- theta_yx[2]
  alpha_yx <- theta_yx[1]
  theta_xy <- expec.reg(y, x, tau = 0.05)
  beta_xy  <- theta_xy[2]
  alpha_xy <- theta_xy[1]
  expec_corr_tau  <-  sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
}

f_expec_0.1_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width =500 ,
  function(x) expec_reg_0.1(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

expec_reg_0.9 <- function(x,y){
  theta_yx <- expec.reg(x, y, tau = 0.95)
  beta_yx  <- theta_yx[2]
  alpha_yx <- theta_yx[1]
  theta_xy <- expec.reg(y, x, tau = 0.95)
  beta_xy  <- theta_xy[2]
  alpha_xy <- theta_xy[1]
  expec_corr_tau  <-  sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
}

f_expec_0.9_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width =500 ,
  function(x) expec_reg_0.9(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

global_crisis_start_date <- as.Date("2007-08-01")
global_crisis_end_date <- as.Date("2009-03-31")
corona_crisis_start_date <- as.Date("2020-03-09")

#2022 plot
corona_crisis_start_date <- as.Date("2019-11-17")
year_seq <- seq(as.Date("2002-02-11"), as.Date("2022-12-31"), by = "1 year")
plot(f_pearson_rolling_win, col="black", lwd=1, type="l", main="S&P 500 and DAX(2001~2022)_Expectile",
     xaxt='n', xlab ="",ylab = TeX("$rho_{tau}"),
     ylim= c(0,0.8))
abline(v = global_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = global_crisis_end_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = corona_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)
lines(f_expec_0.1_win, col="red", lwd=1 )
lines(f_expec_0.9_win, col="blue", lwd=1)
axis(1, at = year_seq, labels = format(year_seq, "%Y"), cex.axis = 0.6)
text(x = as.Date("2008-01-01") + 150 , y = 0.35, "Global \n financial \n Crisis", col = "black", cex = 0.8)
text(x = as.Date("2019-11-17") + 300 , y = 0.7 , "Corona \n Crisis", col = "black", cex = 0.8)

legend("bottomright",
       c(TeX("\\hat{\\rho}"), TeX("\\hat{$rho_{0.1}}"), TeX("\\hat{$rho_{0.9}}")),
       col = c("black", "red", "blue"),
       lty = c("solid", "solid", "solid"),
       ncol = 3, bty = "n")
```

![https://github.com/JuyeonKim-Static/Expectile/blob/a44dce2344c50957dd4ae6c3b42a0c6f9b9c994d/Financial_Data/image/draft_DAX_e.png](https://github.com/JuyeonKim-Static/Expectile/blob/93eec41f5420faaa4d3266e1f1f1475a1e112a21/Financial_Data/image/draft_DAX_e.png))<!-- -->


##Confidence Interval for Quantile and Expectile

```{r}
df<- df1
tau<- 0.9
df$date<- as.Date(df$날짜,format = "%Y.%m.%d")
x <- df[,2]
y <- df[,3]

quant_reg_0.9 <- function(x,y){
  alpha_yx <-rq(y ~ x, tau = 0.9)$coeff[1]
  beta_yx <- rq(y ~ x, tau = 0.9)$coeff[2]
  alpha_xy <-rq(x ~ y, tau = 0.9)$coeff[1]
  beta_xy <- rq(x ~ y, tau = 0.9)$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  #밀도함수추정 for long time
  start_time <- Sys.time()
  n <- 500
  f_1 =vector(length = n)
  for(i in 1:n){
    new_cde1 <- cde(x, y, x.margin = x[i]) #, y.margin = y_margin)
    y_margin <- new_cde1$y
    res <- y_margin - (alpha_yx + beta_yx * x[i])
    where_y <- which.min(sapply(res, function(x){ifelse(x<0, Inf, x)}))
    if(where_y == 1){
      f_1[i] <- new_cde1$z[where_y]
    }else{
      f_1[i]<- (abs(res[where_y-1]) * new_cde1$z[where_y] + abs(res[where_y]) * new_cde1$z[where_y-1]) / 
        (abs(res[where_y-1]) + abs(res[where_y]))
    }
  }
  #x|y분포
  f_2 =vector(length = n)
  for (i in 1:n){
    new_cde2 <- cde(y, x, x.margin = y[i]) # , y.margin = x_margin)
    x_margin <- new_cde2$y
    res <- x_margin - (alpha_xy + beta_xy * y[i])
    where_x <- which.min(sapply(res, function(x){ifelse(x<0, Inf, x)}))
    if(where_x == 1){
      f_2[i] <- new_cde2$z[where_x]
    }else{
      f_2[i]<- (abs(res[where_x-1]) * new_cde2$z[where_x] + abs(res[where_x]) * new_cde2$z[where_x-1]) / 
        (abs(res[where_x-1]) + abs(res[where_x]))
    }
  }
  #R1
  a1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
  for (i in c(1:n)){
    a1[1,1,i] <- 1
    a1[1,2,i] <- x[i]
    a1[2,1,i] <- x[i]
    a1[2,2,i] <- (x[i])^2
    a1[ , ,i] <- f_1[i] * a1[, ,i]
  }
  R1 <- cbind(apply(a1, c(1,2), mean), matrix(0, nrow = 2, ncol = 2))
  #R2
  b1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
  for (i in c(1:n)){
    b1[1,1,i] <- 1
    b1[1,2,i] <- y[i]
    b1[2,1,i] <- y[i]
    b1[2,2,i] <- (y[i])^2
    b1[ , ,i] <- f_2[i] * b1[ , ,i]
  }
  
  R2 <- cbind(matrix(0, nrow = 2, ncol = 2), apply(b1, c(1,2), mean))
  
  #M행렬
  M <- rbind(R1, R2)
  
  #H행렬
  H <- matrix(0, ncol = 4, nrow = 4)
  for (i in 1:n){
    x1 <- c(1,x[i])* (tau - ifelse( (y[i]- (alpha_yx + beta_yx * x[i]))<0 ,1 ,0))
    y1 <- c(1,y[i])* (tau - ifelse( (x[i]- (alpha_xy + beta_xy * y[i]))< 0 ,1 ,0))
    d<-matrix(append(x1,y1) , nrow =1)
    H <- H + (t(d) %*% d)
  }
  H <- H / n
  #G행렬
  G <- 1/2*matrix( c(0, ifelse(beta_yx * beta_xy > 0, sign(beta_yx) * sqrt(beta_xy /beta_yx), 0), 
                     0, ifelse(beta_yx * beta_xy > 0, sign(beta_xy) * sqrt(beta_yx /beta_xy), 0)))
  v <- solve(M, G)
  sigma_tau_quant<- sqrt(t(v) %*% H %*% v)
  upper <- quant_corr_tau + 1.96*sigma_tau_quant/sqrt(500)
  lower <- quant_corr_tau - 1.96*sigma_tau_quant/sqrt(500)
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  print(paste("Elapsed time:", elapsed_time))
  return(list(quant_corr_tau = quant_corr_tau ,sigma_tau_quant = sigma_tau_quant ,
              upper = upper , lower = lower))
}

f_quant0.9_win<- rollapply(
  zoo(cbind(x,y),df$date), 
  width = 500 ,
  function(x) quant_reg_0.9(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

f_quant_0.9_df <- fortify.zoo(f_quant0.9_win, index = TRUE)

#expectile
expec_reg_0.9 <- function(x,y){
  theta_yx <- expec.reg(x, y, tau = 0.9)
  beta_yx  <- theta_yx[2]
  alpha_yx <- theta_yx[1]
  theta_xy <- expec.reg(y, x, tau = 0.9)
  beta_xy  <- theta_xy[2]
  alpha_xy <- theta_xy[1]
  rho_tau_expec <-  sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  e_yx <- y - alpha_yx - beta_yx* x
  e_xy <- x - alpha_xy - beta_xy* y
  w_yx <- tau * (e_yx > 0) + (1-tau) * (e_yx <= 0)
  w_xy <- tau * (e_xy > 0) + (1-tau) * (e_xy <= 0)
  g <- rho_tau_expec/2*matrix(c(1/beta_yx, 1/beta_xy), ncol = 1)
  A <- matrix(c(0,0,1,0,0,0,0,1), ncol = 4)
  det_yx <- mean(w_yx)*mean(w_yx*x^2) - (mean(w_yx*x))^2
  det_xy <- mean(w_xy)*mean(w_xy*y^2) - (mean(w_xy*y))^2
  W.inv <- matrix(c(c(mean(w_yx*x^2), -mean(w_yx*x), 0, 0, -mean(w_yx*x), mean(w_yx), 0, 0)/det_yx,
                    c(0, 0, mean(w_xy*y^2), -mean(w_xy*y), 
                      0, 0, -mean(w_xy*y), mean(w_xy))/det_xy),ncol = 4)/2
  V <- cov(2*cbind(e_yx*w_yx, x*e_yx*w_yx, e_xy*w_xy, y*e_xy*w_xy))
  sigma_tau_expec <- sqrt(t(g) %*% A %*% W.inv %*% V %*% W.inv %*% t(A) %*% g)
  upper <- rho_tau_expec +1.96*sigma_tau_expec/sqrt(500)
  lower <- rho_tau_expec -1.96*sigma_tau_expec/sqrt(500)
  return(list( rho_tau_expec = rho_tau_expec, sigma_tau_expec = sigma_tau_expec, 
               upper = upper ,lower = lower))
}


f_expec_0.9_win <- rollapply(
  zoo(cbind(x,y),df$date), 
  width =500 ,
  function(x) expec_reg_0.9(x[,1],x[,2]), 
  by.column=FALSE,
  align = "center")

f_expec_0.9_df <- fortify.zoo(f_expec_0.9_win, index = TRUE)

#plot
year_seq <- seq(as.Date("2002-02-11"), as.Date("2022-12-31"), by = "1 year")
global_crisis_start_date <- as.Date("2007-08-01")
global_crisis_end_date <- as.Date("2009-03-31")
corona_crisis_start_date <- as.Date("2019-11-17")
plot(f_quant_0.9_df$Index, f_quant_0.9_df$quant_corr_tau , col="red", lwd=1, type="l", main="S&P 500 and DAX(2001~2022)",
     xlab = "Date",xaxt='n', ylim= c(0,1), ylab = TeX("$rho_{tau}"))
abline(v = global_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = global_crisis_end_date, col = "gray", lty = "dashed", lwd = 1)
abline(v = corona_crisis_start_date, col = "gray", lty = "dashed", lwd = 1)

lines(f_quant_0.9_df$Index, f_quant_0.9_df$upper, col="red", lwd=1 )
lines(f_quant_0.9_df$Index, f_quant_0.9_df$lower, col="red", lwd=1)
polygon(c(f_quant_0.9_df$Index, rev(f_quant_0.9_df$Index)),
        c(f_quant_0.9_df$upper, rev(f_quant_0.9_df$lower)), col = "#FFBBBB66", border = NA)

axis(1, at = year_seq, labels = format(year_seq, "%Y"), cex.axis = 0.5)
text(x = as.Date("2008-01-01") + 150 , y = 0.3, "Global \n financial \n Crisis", col = "black", cex = 0.8)
text(x = as.Date("2019-11-17") + 300 , y = 0.9 , "Corona \n Crisis", col = "black", cex = 0.8)


lines(f_expec_0.9_df$Index,f_expec_0.9_df$rho_tau_expec , col="blue", lwd=1)
lines(f_expec_0.9_df$Index, f_expec_0.9_df$upper, col="blue", lwd=1 )
lines(f_expec_0.9_df$Index, f_expec_0.9_df$lower, col="blue", lwd=1)
polygon(c(f_expec_0.9_df$Index, rev(f_expec_0.9_df$Index)),
        c(f_expec_0.9_df$upper, rev(f_expec_0.9_df$lower)), col = "#BBAAFF66", border = NA)
legend("bottomleft", lty = c("solid", "solid"),
       legend=c("quantile","expectile") , col =c("red","blue"),bty = "n") 

```



