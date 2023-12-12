rolling_window_simulation
================
Juyeon Kim
2023-12-12

\##Expectile

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

# df<-read.csv("S&P_DAX(2001~2017).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_FTSE(2001~2017).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_CAC(2001~2017).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_KOSPI(2001~2017).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# 
df<-read.csv("S&P_DAX(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_FTSE(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_CAC(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")
# df<-read.csv("S&P_KOSPI(2001~2022).csv" ,fileEncoding = "CP949", encoding = "UTF-8")

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

![](draft_simulation_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
