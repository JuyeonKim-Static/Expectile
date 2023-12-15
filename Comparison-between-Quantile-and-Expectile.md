Comparison between Quantile and Expectile
================
Juyeon Kim
2023-12-12

\#Normal distribution

``` r
library(mvtnorm)
library(quantreg)
```

    ## Loading required package: SparseM

    ## 
    ## Attaching package: 'SparseM'

    ## The following object is masked from 'package:base':
    ## 
    ##     backsolve

``` r
rho <- 0.5
mvn.dat <- as.data.frame(rmvnorm(100000, sigma = matrix(c(1, rho, rho, 1), ncol = 2)))
names(mvn.dat) <- c("x", "y")
plot(mvn.dat$x, mvn.dat$y) 
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# quantile regression 
#install.packages("quantreg")

tau <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
quantile_p= c()
for (i in 1: length(tau)){
  beta_yx <- rq(y ~ x, tau = tau[i], data = mvn.dat)$coeff[2]
  beta_xy <- rq(x ~ y, tau = tau[i], data = mvn.dat)$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  quantile_p[i] = quant_corr_tau
}

p = cor(mvn.dat$x, mvn.dat$y)
p
```

    ## [1] 0.4972591

``` r
# expectile regression
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
      print(paste("procedure does not converge")); 
      break; 
    }
    coeff <- lmfit$coefficients[2]
    names(coeff) <- NULL
  }
  return (coeff)
}

expectile_p <- c()
for (i in 1: length(tau)){
  beta_yx <- expec.reg(mvn.dat$x, mvn.dat$y, tau = tau[i])
  beta_xy <- expec.reg(mvn.dat$y, mvn.dat$x, tau = tau[i])
  names(beta_yx) <- names(beta_xy) <- NULL
  expec_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  expectile_p[i] = expec_corr_tau
}


#plot
library(latex2exp)
plot(tau , quantile_p , type = 'o', col = 'red' , main = "(a) Bivariate normal", 
     xlab = TeX('$\\Tau$'), ylab = TeX('$\\hat{p_\\Tau}$'),lty=6 , pch=16)
abline(h = p , lty = 3 ,col="black" )
par(new=TRUE)
plot(tau , expectile_p, type = 'o', col = 'blue', axes = FALSE ,xlab = TeX('$\\Tau$') ,ylab = TeX('$\\hat{p_\\Tau}$'),pch = 17,lty=6 )
legend("topright", legend=c( TeX('quantile_$\\hat{p_\\Tau}$') , TeX('expectile_$\\hat{p_\\Tau}$')),col=c("red" , "blue"), lty = c(6, 6),pch = c(16,17) )
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->
\#Rocket Type

``` r
library(mvtnorm)
library(quantreg) 
c <- -1.645
rho <- 0.5
mvn.dat <- as.data.frame(rmvnorm(100000, sigma = matrix(c(1, rho, rho, 1), ncol = 2)))
names(mvn.dat) <- c("x", "y")
I <- ifelse( (mvn.dat$x<= c) & (mvn.dat$y<= c),1 ,0)
mvn.dat1 <- I* rnorm(n =100000, mean = 0, sd = 1)
mvn.dat$x <- mvn.dat$x + mvn.dat1                        
mvn.dat$y <- mvn.dat$y + mvn.dat1
plot(mvn.dat$x, mvn.dat$y)
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
tau <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
quantile_p<- c()
for (i in 1: length(tau)){
  beta_yx <- rq(y ~ x, tau = tau[i], data = mvn.dat)$coeff[2]
  beta_xy <- rq(x ~ y, tau = tau[i], data = mvn.dat)$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  quantile_p[i] = quant_corr_tau
}

p = cor(mvn.dat$x, mvn.dat$y)
# expectile regression
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
      print(paste("procedure does not converge")); 
      break; 
    }
    coeff <- lmfit$coefficients[2]
    names(coeff) <- NULL
  }
  return (coeff)
}

expectile_p= c()
for (i in 1: length(tau)){
  beta_yx <- expec.reg(mvn.dat$x, mvn.dat$y, tau = tau[i])
  beta_xy <- expec.reg(mvn.dat$y, mvn.dat$x, tau = tau[i])
  names(beta_yx) <- names(beta_xy) <- NULL
  expec_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  expectile_p[i] = expec_corr_tau
}

#plot
library(latex2exp)
plot(tau , quantile_p , type = 'o', col = 'red' , main = "(b) Rocket type", 
     xlab = TeX('$\\Tau$'), ylab = TeX('$\\hat{p_\\Tau}$'),lty=6 , pch=16)
abline(h = p , lty = 3 ,col="black" )
par(new=TRUE)
plot(tau , expectile_p, type = 'o', col = 'blue', axes = FALSE ,xlab = TeX('$\\Tau$') ,ylab = TeX('$\\hat{p_\\Tau}$'),pch = 17,lty=6 )
legend("topright", legend=c( TeX('quantile_$\\hat{p_\\Tau}$') , TeX('expectile_$\\hat{p_\\Tau}$')),col=c("red" , "blue"), lty = c(6, 6),pch = c(16,17) )
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

\#T-distribution

``` r
library(mvtnorm)
library(quantreg)
# df =3
rho <- 0.5
mvn.dat <- as.data.frame(rmvt(10000, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = 3 ) )
names(mvn.dat) <- c("x", "y")
plot(mvn.dat$x, mvn.dat$y) 
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
tau <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
quantile_p<-c()
for (i in 1: length(tau)){
  beta_yx <- rq(y ~ x, tau = tau[i], data = mvn.dat)$coeff[2]
  beta_xy <- rq(x ~ y, tau = tau[i], data = mvn.dat)$coeff[2]
  names(beta_yx) <- names(beta_xy) <- NULL
  quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  quantile_p[i] = quant_corr_tau
}


#자유도에 따른 분포
nu <- c(3,5,7,10,100)
df_quantile<- matrix(c(1 : (length(nu)*length(tau))) ,nrow = length(nu))

for ( a in 1 : length(nu) ){
  mvn.dat <- as.data.frame(rmvt(1000,sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = nu[a]))
  names(mvn.dat) <- c("x", "y")
  for ( b in 1: length(tau)){
    beta_yx <- rq(y ~ x, tau = tau[b], data = mvn.dat)$coeff[2]
    beta_xy <- rq(x ~ y, tau = tau[b], data = mvn.dat)$coeff[2]
    names(beta_yx) <- names(beta_xy) <- NULL
    quant_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
    df_quantile[a,b] <- quant_corr_tau
  }}


#자유도에 따른 값 변화 양상
plot(tau ,df_quantile[1,], type = 'o', col = 'red' , main = "Comparing quantile correlations based on degrees of freedom",
     ylim = c(0,1) , xlab = TeX('$\\Tau$'), ylab = TeX('$\\hat{p_\\Tau}$'))
lines(tau ,df_quantile[2,], type = 'o', col = 'blue')
lines(tau ,df_quantile[3,], type = 'o', col = 'green')
lines(tau ,df_quantile[4,], type = 'o', col = 'orange')
lines(tau ,df_quantile[5,], type = 'o', col = 'black')
legend("topright", legend=c("3","5","7","10","100") , fill=c("red","blue" , "green" ,"orange" ,"black"))
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
df_expectile<- matrix(c(1 : (length(nu)*length(tau))) ,nrow = length(nu))
for ( a in 1 : length(nu) ){
  mvn.dat <- as.data.frame(rmvt(1000,sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = nu[a]))
  names(mvn.dat) <- c("x", "y")
  for ( b in 1: length(tau)){
  beta_yx <- expec.reg(mvn.dat$x, mvn.dat$y, tau = tau[i])
  beta_xy <- expec.reg(mvn.dat$y, mvn.dat$x, tau = tau[i])
  names(beta_yx) <- names(beta_xy) <- NULL
  expec_corr_tau <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
  df_expectile[a,b] <- expec_corr_tau
  }
}
  
plot(tau ,df_expectile[1,], type = 'o', col = 'red' , main = "Comparing expectile correlations based on degrees of freedom",
     ylim = c(0,1) , xlab = TeX('$\\Tau$'), ylab = TeX('$\\hat{p_\\Tau}$'))
lines(tau ,df_expectile[2,], type = 'o', col = 'blue')
lines(tau ,df_expectile[3,], type = 'o', col = 'green')
lines(tau ,df_expectile[4,], type = 'o', col = 'orange')
lines(tau ,df_expectile[5,], type = 'o', col = 'black')
legend("topright", legend=c("3","5","7","10","100") , fill=c("red","blue" , "green" ,"orange" ,"black"))
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
#plot
plot(tau , quantile_p , type = 'o', col = 'red' , main = "(c)T-distribution(df = 3)",
     xlab = TeX('$\\Tau$'), ylab = TeX('$\\hat{p_\\Tau}$'))
abline(h = p , col="black")
par(new=TRUE)
plot(tau , df_expectile[1,] , type = 'o', col = 'blue', axes = FALSE ,ylab = "")
axis(side = 4, col = "blue", col.axis = "blue")
legend("bottomright", legend=c("quantile","expectile") , fill=c("red","blue"))
```

![](real-Comparison-between-Quantile-and-Expectile_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->
