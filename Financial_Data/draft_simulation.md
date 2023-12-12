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
df_quantile
```

    ##  [1] 0.6035571 0.6025315 0.6022479 0.6031188 0.6016500 0.6019897 0.6021963
    ##  [8] 0.5974986 0.5980950 0.5937924 0.5905833 0.5911337 0.5854480 0.5802684
    ## [15] 0.5817596 0.5826730 0.5762683 0.5731373 0.5756738 0.5742107 0.5741160
    ## [22] 0.5737452 0.5729492 0.5724724 0.5682433 0.5648355 0.5686390 0.5655157
    ## [29] 0.5618645 0.5617534 0.5598111 0.5596437 0.5593840 0.5568013 0.5510576
    ## [36] 0.5465188 0.5426736 0.5410316 0.5408783 0.5380300 0.5335811 0.5313786
    ## [43] 0.5301722 0.5254211 0.5259198 0.5269224 0.5280863 0.5259276 0.5267648
    ## [50] 0.5241808 0.5246571 0.5225618 0.5222029 0.5226813 0.5218643 0.5202940
    ## [57] 0.5188323 0.5178481 0.5202939 0.5227579 0.5209720 0.5204598 0.5166590
    ## [64] 0.5148198 0.5187161 0.5221628 0.5235193 0.5209853 0.5207279 0.5169244
    ## [71] 0.5152612 0.5170681 0.5196516 0.5227242 0.5234163 0.5246917 0.5215508
    ## [78] 0.5167015 0.5217422 0.5159843 0.5185953 0.5192726 0.5280540 0.5324044
    ## [85] 0.5300865 0.5335456 0.5391399 0.5473004 0.5473976 0.5480357 0.5430832

``` r
sigma_tau_quant
```

    ##  [1] 1.630299 1.679967 1.745031 1.825924 1.876933 1.945829 2.005025 2.000799
    ##  [9] 2.050195 2.090848 2.121131 2.200563 2.224086 2.243381 2.267615 2.297474
    ## [17] 2.319303 2.358989 2.375646 2.391010 2.451050 2.472587 2.499688 2.520696
    ## [25] 2.546732 2.567270 2.575030 2.588983 2.626832 2.636299 2.645148 2.661830
    ## [33] 2.668977 2.678235 2.683869 2.689657 2.696509 2.698057 2.699750 2.700925
    ## [41] 2.702123 2.704194 2.705579 2.704093 2.705364 2.705950 2.707113 2.706257
    ## [49] 2.707564 2.709080 2.706128 2.706092 2.705398 2.702452 2.690154 2.691179
    ## [57] 2.691536 2.687296 2.681505 2.678562 2.672609 2.661717 2.651362 2.647193
    ## [65] 2.647898 2.644387 2.645440 2.636904 2.613102 2.603310 2.592063 2.586901
    ## [73] 2.565196 2.552577 2.543190 2.541743 2.527211 2.516184 2.467488 2.446868
    ## [81] 2.436081 2.374623 2.326017 2.297816 2.283735 2.272770 2.222496 2.177950
    ## [89] 2.136451 2.081382 2.043256

``` r
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
sigma_tau_expec
```

    ##  [1] 2.175583 2.048133 1.964418 1.877755 1.843076 1.798614 1.745687 1.709124
    ##  [9] 1.658058 1.629800 1.604097 1.572714 1.531642 1.506308 1.489284 1.471786
    ## [17] 1.449264 1.431018 1.415381 1.402937 1.390770 1.371788 1.359284 1.348962
    ## [25] 1.340554 1.323693 1.317474 1.311867 1.307579 1.302058 1.295260 1.289882
    ## [33] 1.287341 1.285139 1.282134 1.280461 1.280075 1.279581 1.279842 1.280890
    ## [41] 1.281799 1.282977 1.284510 1.286361 1.288506 1.291187 1.294155 1.297525
    ## [49] 1.301232 1.305063 1.309234 1.314255 1.319284 1.324651 1.331644 1.338260
    ## [57] 1.345556 1.355442 1.361776 1.370089 1.380838 1.390644 1.400223 1.414819
    ## [65] 1.427876 1.439324 1.451330 1.465758 1.479823 1.496029 1.511576 1.529617
    ## [73] 1.547012 1.569229 1.584551 1.615725 1.637290 1.656044 1.688243 1.718840
    ## [81] 1.755152 1.790718 1.826223 1.877580 1.909655 1.947351 2.002062 2.071384
    ## [89] 2.151364 2.255377 2.374678

``` r
df_expectile
```

    ##  [1] 0.5847309 0.5847475 0.5843296 0.5837431 0.5830207 0.5822885 0.5815613
    ##  [8] 0.5807661 0.5799746 0.5792406 0.5785678 0.5778607 0.5771914 0.5765382
    ## [15] 0.5758829 0.5752419 0.5746233 0.5740300 0.5734295 0.5728551 0.5722875
    ## [22] 0.5717309 0.5711986 0.5706871 0.5701946 0.5697106 0.5692198 0.5687382
    ## [29] 0.5682588 0.5677816 0.5673164 0.5668509 0.5663904 0.5659333 0.5654768
    ## [36] 0.5650316 0.5645934 0.5641633 0.5637414 0.5633276 0.5629233 0.5625316
    ## [43] 0.5621496 0.5617774 0.5614165 0.5610666 0.5607278 0.5604004 0.5600840
    ## [50] 0.5597792 0.5594857 0.5592077 0.5589455 0.5586992 0.5584642 0.5582424
    ## [57] 0.5580330 0.5578357 0.5576467 0.5574748 0.5573141 0.5571667 0.5570357
    ## [64] 0.5569215 0.5568201 0.5567363 0.5566669 0.5566139 0.5565767 0.5565768
    ## [71] 0.5565973 0.5566443 0.5567165 0.5568060 0.5569180 0.5570566 0.5572288
    ## [78] 0.5574297 0.5576636 0.5579449 0.5583210 0.5587950 0.5593817 0.5601416
    ## [85] 0.5610251 0.5620162 0.5631746 0.5646017 0.5661790 0.5680440 0.5699520

``` r
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



