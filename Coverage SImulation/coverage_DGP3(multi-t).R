library(hdrcde)
library(mvtnorm)
library(quantreg)
library(foreach)
library(doParallel)

getwd()
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

rho_vec <- c(0.3, 0.5, 0.7) 
tau_vec <- c(1:9)/10 
set.seed(230427)# artibrarily chosen; should be changed at every simulation
# df = 3 -> seed 230325
# df = 5 -> seed 230326
# df = 7 ->seed  230427
simul_num <- 1
result_obj <- list()

for(i in 1:3){
  for(j in 1:9){
    for(df_num in 4){
      rho <- rho_vec[i] # should be changed \in {0.3, 0.5, 0.7}
      tau <- tau_vec[j] #should be changed \in {0.1, 0.2, ..., 0.9}
      df <- df_vec[df_num]
      
      if(j %% 2 == 0){
        mvn.dat <- as.data.frame(mvtnorm::rmvt(100, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        mvn.dat <- as.data.frame(mvtnorm::rmvt(500, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        mvn.dat <- as.data.frame(mvtnorm::rmvt(1000, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        next()
      }
      
      ##### n = 100 ##### 
      
      cores <- detectCores() - 1 # to avoid overload
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      time1 <- Sys.time() 
      
      n <- 100
      coverage_tau <- foreach(c = 1:1000, .combine = rbind, .packages = c("hdrcde", "mvtnorm", "quantreg")) %dopar% {
        mvn.dat <- as.data.frame(mvtnorm::rmvt(n, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        names(mvn.dat) <- c("x", "y")
        mvn_x <- mvn.dat$x
        mvn_y <- mvn.dat$y
        
        ## quantile correlation
        theta_yx <- rq(y ~ x, tau = tau, data = mvn.dat)$coeff
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- rq(x ~ y, tau = tau, data = mvn.dat)$coeff
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_quant <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        # y|x분포
        f_1 =vector(length = n)
        for(i in 1:n){
          # y_margin <-seq(from = range(mvn_y)[1], to = range(mvn_y)[2],length = 100)
          new_cde1 <- cde(mvn_x, mvn_y, x.margin = mvn_x[i]) #, y.margin = y_margin)
          y_margin <- new_cde1$y
          res <- y_margin - (alpha_yx + beta_yx * mvn_x[i])
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
          # x_margin <-seq(from = range(mvn_x)[1], to = range(mvn_x)[2],length = 100)
          new_cde2 <- cde(mvn_y, mvn_x, x.margin = mvn_y[i]) # , y.margin = x_margin)
          x_margin <- new_cde2$y
          res <- x_margin - (alpha_xy + beta_xy * mvn.dat$y[i])
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
          a1[1,2,i] <- mvn.dat$x[i]
          a1[2,1,i] <- mvn.dat$x[i]
          a1[2,2,i] <- (mvn.dat$x[i])^2
          a1[ , ,i] <- f_1[i] * a1[, ,i]
        }
        # t1<-0;t2<-0;t3<-0;t4<-0;
        # for (i in c(1:n)){
        #   t1<-  t1 + (f_1[i] * a1[1,1,i])
        #   t2<-  t2 + (f_1[i] * a1[1,2,i])
        #   t3<-  t3 + (f_1[i] * a1[2,1,i])
        #   t4<-  t4 + (f_1[i] * a1[2,1,i])
        #   }
        R1 <- cbind(apply(a1, c(1,2), mean), matrix(0, nrow = 2, ncol = 2))
        #R2
        b1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
        for (i in c(1:n)){
          b1[1,1,i] <- 1
          b1[1,2,i] <- mvn.dat$y[i]
          b1[2,1,i] <- mvn.dat$y[i]
          b1[2,2,i] <- (mvn.dat$y[i])^2
          b1[ , ,i] <- f_2[i] * b1[ , ,i]
        }
        # t5<-0;t6<-0;t7<-0;t8<-0;
        # for (i in c(1:n)){
        #   t5<- t5 + (f_2[i] * b1[1,1,i])
        #   t6<- t6 + (f_2[i] * b1[1,2,i])
        #   t7<- t7 + (f_2[i] * b1[2,1,i])
        #   t8<- t8 + (f_2[i] * b1[2,2,i])
        # }
        R2 <- cbind(matrix(0, nrow = 2, ncol = 2), apply(b1, c(1,2), mean))
        
        #M행혛
        # M <- matrix(c(t1,t3,0,0,t2,t4,0,0,0,0,t5,t7,0,0,t6,t8),ncol =4)
        # M <- M / n
        M <- rbind(R1, R2)
        
        #H행렬
        H <- matrix(0, ncol = 4, nrow = 4)
        for (i in 1:n){
          x1 <- c(1,mvn_x[i])* (tau - ifelse( (mvn_y[i]- (alpha_yx + beta_yx * mvn_x[i]))<0 ,1 ,0))
          y1 <- c(1,mvn_y[i])* (tau - ifelse( (mvn_x[i]- (alpha_xy + beta_xy * mvn_y[i]))< 0 ,1 ,0))
          d<-matrix(append(x1,y1) , nrow =1)
          H <- H + (t(d) %*% d)
        }
        H <- H / n
        #G행렬
        G <- 1/2*matrix( c(0, ifelse(beta_yx * beta_xy > 0, sign(beta_yx) * sqrt(beta_xy /beta_yx), 0), 
                           0, ifelse(beta_yx * beta_xy > 0, sign(beta_xy) * sqrt(beta_yx /beta_xy), 0)))
        # t(G) %*% solve(M) %*% H %*% solve(M) %*% G
        v <- solve(M, G)
        sigma_tau_quant <- sqrt(t(v) %*% H %*% v)
        coverage_tau_quant <- (rho_tau_quant - 1.96*sigma_tau_quant/sqrt(n) <= rho) &
          (rho_tau_quant + 1.96*sigma_tau_quant/sqrt(n) >= rho)
        
        ## expectile correlation
        theta_yx <- expec.reg(mvn_x, mvn_y, tau = tau)
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- expec.reg(mvn_y, mvn_x, tau = tau)
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_expec <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        
        e_yx <- mvn_y - alpha_yx - beta_yx*mvn_x
        e_xy <- mvn_x - alpha_xy - beta_xy*mvn_y
        w_yx <- tau * (e_yx > 0) + (1-tau) * (e_yx <= 0)
        w_xy <- tau * (e_xy > 0) + (1-tau) * (e_xy <= 0)
        
        g <- rho_tau_expec/2*matrix(c(1/beta_yx, 1/beta_xy), ncol = 1)
        A <- matrix(c(0,0,1,0,0,0,0,1), ncol = 4)
        det_yx <- mean(w_yx)*mean(w_yx*mvn_x^2) - (mean(w_yx*mvn_x))^2
        det_xy <- mean(w_xy)*mean(w_xy*mvn_y^2) - (mean(w_xy*mvn_y))^2
        W.inv <- matrix(c(c(mean(w_yx*mvn_x^2), -mean(w_yx*mvn_x), 0, 0, -mean(w_yx*mvn_x), mean(w_yx), 0, 0)/det_yx,
                          c(0, 0, mean(w_xy*mvn_y^2), -mean(w_xy*mvn_y), 
                            0, 0, -mean(w_xy*mvn_y), mean(w_xy))/det_xy),ncol = 4)/2
        
        V <- cov(2*cbind(e_yx*w_yx, mvn_x*e_yx*w_yx, e_xy*w_xy, mvn_y*e_xy*w_xy))
        sigma_tau_expec <- sqrt(t(g) %*% A %*% W.inv %*% V %*% W.inv %*% t(A) %*% g)
        coverage_tau_expec <- (rho_tau_expec - 1.96*sigma_tau_expec/sqrt(n) <= rho) & 
          (rho_tau_expec + 1.96*sigma_tau_expec/sqrt(n) >= rho)
        
        # merge both results
        c(rho_tau_quant, sigma_tau_quant, coverage_tau_quant, 
          rho_tau_expec, sigma_tau_expec, coverage_tau_expec)
      }
      
      stopCluster(cl)
      
      time2 <- Sys.time()
      
      
      names(coverage_tau) <- c("quant_corr", "sigma_quant", "coverage_quant", 
                               "expec_corr", "sigma_expec", "coverage_expec")
      
      cat ("Simulation: DGP3, multivariate t with df ", df, "\n", 
           "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
           "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
           "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
           difftime(time2, time1, unit = "mins"), "mins", "\n")
      
      result_obj[[simul_num]] <- list(obj = coverage_tau, 
                                      description = paste("Simulation: DGP3, multivariate t with df ", df, "\n", 
                                                          "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
                                                          "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
                                                          "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
                                                          difftime(time2, time1, unit = "mins"), "mins", 
                                                          sep = " "))
      simul_num <- simul_num + 1
      saveRDS(result_obj, paste0("/Users/gimjuyeon/Documents/학부연구생/연구/coverage_result_DGP3(multi-t(", df, ")).RDS"))
      
      ##### n = 500 ##### 
      
      cores <- detectCores() - 1 # to avoid overload
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      time1 <- Sys.time() 
      
      n <- 500
      coverage_tau <- foreach(c = 1:1000, .combine = rbind, .packages = c("hdrcde", "mvtnorm", "quantreg")) %dopar% {
        mvn.dat <- as.data.frame(mvtnorm::rmvt(n, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        names(mvn.dat) <- c("x", "y")
        mvn_x <- mvn.dat$x
        mvn_y <- mvn.dat$y
        
        ## quantile correlation
        theta_yx <- rq(y ~ x, tau = tau, data = mvn.dat)$coeff
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- rq(x ~ y, tau = tau, data = mvn.dat)$coeff
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_quant <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        # y|x분포
        f_1 =vector(length = n)
        for(i in 1:n){
          # y_margin <-seq(from = range(mvn_y)[1], to = range(mvn_y)[2],length = 100)
          new_cde1 <- cde(mvn_x, mvn_y, x.margin = mvn_x[i]) #, y.margin = y_margin)
          y_margin <- new_cde1$y
          res <- y_margin - (alpha_yx + beta_yx * mvn_x[i])
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
          # x_margin <-seq(from = range(mvn_x)[1], to = range(mvn_x)[2],length = 100)
          new_cde2 <- cde(mvn_y, mvn_x, x.margin = mvn_y[i]) # , y.margin = x_margin)
          x_margin <- new_cde2$y
          res <- x_margin - (alpha_xy + beta_xy * mvn.dat$y[i])
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
          a1[1,2,i] <- mvn.dat$x[i]
          a1[2,1,i] <- mvn.dat$x[i]
          a1[2,2,i] <- (mvn.dat$x[i])^2
          a1[ , ,i] <- f_1[i] * a1[, ,i]
        }
        # t1<-0;t2<-0;t3<-0;t4<-0;
        # for (i in c(1:n)){
        #   t1<-  t1 + (f_1[i] * a1[1,1,i])
        #   t2<-  t2 + (f_1[i] * a1[1,2,i])
        #   t3<-  t3 + (f_1[i] * a1[2,1,i])
        #   t4<-  t4 + (f_1[i] * a1[2,1,i])
        #   }
        R1 <- cbind(apply(a1, c(1,2), mean), matrix(0, nrow = 2, ncol = 2))
        #R2
        b1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
        for (i in c(1:n)){
          b1[1,1,i] <- 1
          b1[1,2,i] <- mvn.dat$y[i]
          b1[2,1,i] <- mvn.dat$y[i]
          b1[2,2,i] <- (mvn.dat$y[i])^2
          b1[ , ,i] <- f_2[i] * b1[ , ,i]
        }
        # t5<-0;t6<-0;t7<-0;t8<-0;
        # for (i in c(1:n)){
        #   t5<- t5 + (f_2[i] * b1[1,1,i])
        #   t6<- t6 + (f_2[i] * b1[1,2,i])
        #   t7<- t7 + (f_2[i] * b1[2,1,i])
        #   t8<- t8 + (f_2[i] * b1[2,2,i])
        # }
        R2 <- cbind(matrix(0, nrow = 2, ncol = 2), apply(b1, c(1,2), mean))
        
        #M행혛
        # M <- matrix(c(t1,t3,0,0,t2,t4,0,0,0,0,t5,t7,0,0,t6,t8),ncol =4)
        # M <- M / n
        M <- rbind(R1, R2)
        
        #H행렬
        H <- matrix(0, ncol = 4, nrow = 4)
        for (i in 1:n){
          x1 <- c(1,mvn_x[i])* (tau - ifelse( (mvn_y[i]- (alpha_yx + beta_yx * mvn_x[i]))<0 ,1 ,0))
          y1 <- c(1,mvn_y[i])* (tau - ifelse( (mvn_x[i]- (alpha_xy + beta_xy * mvn_y[i]))< 0 ,1 ,0))
          d<-matrix(append(x1,y1) , nrow =1)
          H <- H + (t(d) %*% d)
        }
        H <- H / n
        #G행렬
        G <- 1/2*matrix( c(0, ifelse(beta_yx * beta_xy > 0, sign(beta_yx) * sqrt(beta_xy /beta_yx), 0), 
                           0, ifelse(beta_yx * beta_xy > 0, sign(beta_xy) * sqrt(beta_yx /beta_xy), 0)))
        # t(G) %*% solve(M) %*% H %*% solve(M) %*% G
        v <- solve(M, G)
        sigma_tau_quant <- sqrt(t(v) %*% H %*% v)
        coverage_tau_quant <- (rho_tau_quant - 1.96*sigma_tau_quant/sqrt(n) <= rho) &
          (rho_tau_quant + 1.96*sigma_tau_quant/sqrt(n) >= rho)
        
        ## expectile correlation
        theta_yx <- expec.reg(mvn_x, mvn_y, tau = tau)
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- expec.reg(mvn_y, mvn_x, tau = tau)
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_expec <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        
        e_yx <- mvn_y - alpha_yx - beta_yx*mvn_x
        e_xy <- mvn_x - alpha_xy - beta_xy*mvn_y
        w_yx <- tau * (e_yx > 0) + (1-tau) * (e_yx <= 0)
        w_xy <- tau * (e_xy > 0) + (1-tau) * (e_xy <= 0)
        
        g <- rho_tau_expec/2*matrix(c(1/beta_yx, 1/beta_xy), ncol = 1)
        A <- matrix(c(0,0,1,0,0,0,0,1), ncol = 4)
        det_yx <- mean(w_yx)*mean(w_yx*mvn_x^2) - (mean(w_yx*mvn_x))^2
        det_xy <- mean(w_xy)*mean(w_xy*mvn_y^2) - (mean(w_xy*mvn_y))^2
        W.inv <- matrix(c(c(mean(w_yx*mvn_x^2), -mean(w_yx*mvn_x), 0, 0, -mean(w_yx*mvn_x), mean(w_yx), 0, 0)/det_yx,
                          c(0, 0, mean(w_xy*mvn_y^2), -mean(w_xy*mvn_y), 
                            0, 0, -mean(w_xy*mvn_y), mean(w_xy))/det_xy),ncol = 4)/2
        
        V <- cov(2*cbind(e_yx*w_yx, mvn_x*e_yx*w_yx, e_xy*w_xy, mvn_y*e_xy*w_xy))
        sigma_tau_expec <- sqrt(t(g) %*% A %*% W.inv %*% V %*% W.inv %*% t(A) %*% g)
        coverage_tau_expec <- (rho_tau_expec - 1.96*sigma_tau_expec/sqrt(n) <= rho) & 
          (rho_tau_expec + 1.96*sigma_tau_expec/sqrt(n) >= rho)
        
        # merge both results
        c(rho_tau_quant, sigma_tau_quant, coverage_tau_quant, 
          rho_tau_expec, sigma_tau_expec, coverage_tau_expec)
      }
      
      stopCluster(cl)
      
      time2 <- Sys.time()
      
      
      names(coverage_tau) <- c("quant_corr", "sigma_quant", "coverage_quant", 
                               "expec_corr", "sigma_expec", "coverage_expec")
      
      cat ("Simulation: DGP3, multivariate t with df ", df, "\n", 
           "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
           "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
           "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
           difftime(time2, time1, unit = "mins"), "mins", "\n")
      
      result_obj[[simul_num]] <- list(obj = coverage_tau, 
                                      description = paste("Simulation: DGP3, multivariate t with df ", df, "\n", 
                                                          "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
                                                          "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
                                                          "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
                                                          difftime(time2, time1, unit = "mins"), "mins", 
                                                          sep = " "))
      simul_num <- simul_num + 1
      saveRDS(result_obj, paste0("/Users/gimjuyeon/Documents/학부연구생/연구/coverage_result_DGP3(multi-t(", df, ")).RDS"))
      
      ##### n = 1000 ##### 
      
      cores <- detectCores() - 1 # to avoid overload
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      
      time1 <- Sys.time() 
      
      n <- 1000
      coverage_tau <- foreach(c = 1:1000, .combine = rbind, .packages = c("hdrcde", "mvtnorm", "quantreg")) %dopar% {
        mvn.dat <- as.data.frame(mvtnorm::rmvt(n, sigma = matrix(c(1, rho, rho, 1), ncol = 2), df = df))
        names(mvn.dat) <- c("x", "y")
        mvn_x <- mvn.dat$x
        mvn_y <- mvn.dat$y
        
        ## quantile correlation
        theta_yx <- rq(y ~ x, tau = tau, data = mvn.dat)$coeff
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- rq(x ~ y, tau = tau, data = mvn.dat)$coeff
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_quant <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        # y|x분포
        f_1 =vector(length = n)
        for(i in 1:n){
          # y_margin <-seq(from = range(mvn_y)[1], to = range(mvn_y)[2],length = 100)
          new_cde1 <- cde(mvn_x, mvn_y, x.margin = mvn_x[i]) #, y.margin = y_margin)
          y_margin <- new_cde1$y
          res <- y_margin - (alpha_yx + beta_yx * mvn_x[i])
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
          # x_margin <-seq(from = range(mvn_x)[1], to = range(mvn_x)[2],length = 100)
          new_cde2 <- cde(mvn_y, mvn_x, x.margin = mvn_y[i]) # , y.margin = x_margin)
          x_margin <- new_cde2$y
          res <- x_margin - (alpha_xy + beta_xy * mvn.dat$y[i])
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
          a1[1,2,i] <- mvn.dat$x[i]
          a1[2,1,i] <- mvn.dat$x[i]
          a1[2,2,i] <- (mvn.dat$x[i])^2
          a1[ , ,i] <- f_1[i] * a1[, ,i]
        }
        # t1<-0;t2<-0;t3<-0;t4<-0;
        # for (i in c(1:n)){
        #   t1<-  t1 + (f_1[i] * a1[1,1,i])
        #   t2<-  t2 + (f_1[i] * a1[1,2,i])
        #   t3<-  t3 + (f_1[i] * a1[2,1,i])
        #   t4<-  t4 + (f_1[i] * a1[2,1,i])
        #   }
        R1 <- cbind(apply(a1, c(1,2), mean), matrix(0, nrow = 2, ncol = 2))
        #R2
        b1 <-array(0, dim =c(2,2,n)) #비어있는 행렬 생성
        for (i in c(1:n)){
          b1[1,1,i] <- 1
          b1[1,2,i] <- mvn.dat$y[i]
          b1[2,1,i] <- mvn.dat$y[i]
          b1[2,2,i] <- (mvn.dat$y[i])^2
          b1[ , ,i] <- f_2[i] * b1[ , ,i]
        }
        # t5<-0;t6<-0;t7<-0;t8<-0;
        # for (i in c(1:n)){
        #   t5<- t5 + (f_2[i] * b1[1,1,i])
        #   t6<- t6 + (f_2[i] * b1[1,2,i])
        #   t7<- t7 + (f_2[i] * b1[2,1,i])
        #   t8<- t8 + (f_2[i] * b1[2,2,i])
        # }
        R2 <- cbind(matrix(0, nrow = 2, ncol = 2), apply(b1, c(1,2), mean))
        
        #M행혛
        # M <- matrix(c(t1,t3,0,0,t2,t4,0,0,0,0,t5,t7,0,0,t6,t8),ncol =4)
        # M <- M / n
        M <- rbind(R1, R2)
        
        #H행렬
        H <- matrix(0, ncol = 4, nrow = 4)
        for (i in 1:n){
          x1 <- c(1,mvn_x[i])* (tau - ifelse( (mvn_y[i]- (alpha_yx + beta_yx * mvn_x[i]))<0 ,1 ,0))
          y1 <- c(1,mvn_y[i])* (tau - ifelse( (mvn_x[i]- (alpha_xy + beta_xy * mvn_y[i]))< 0 ,1 ,0))
          d<-matrix(append(x1,y1) , nrow =1)
          H <- H + (t(d) %*% d)
        }
        H <- H / n
        #G행렬
        G <- 1/2*matrix( c(0, ifelse(beta_yx * beta_xy > 0, sign(beta_yx) * sqrt(beta_xy /beta_yx), 0), 
                           0, ifelse(beta_yx * beta_xy > 0, sign(beta_xy) * sqrt(beta_yx /beta_xy), 0)))
        # t(G) %*% solve(M) %*% H %*% solve(M) %*% G
        v <- solve(M, G)
        sigma_tau_quant <- sqrt(t(v) %*% H %*% v)
        coverage_tau_quant <- (rho_tau_quant - 1.96*sigma_tau_quant/sqrt(n) <= rho) &
          (rho_tau_quant + 1.96*sigma_tau_quant/sqrt(n) >= rho)
        
        ## expectile correlation
        theta_yx <- expec.reg(mvn_x, mvn_y, tau = tau)
        beta_yx  <- theta_yx[2]
        alpha_yx <- theta_yx[1]
        theta_xy <- expec.reg(mvn_y, mvn_x, tau = tau)
        beta_xy  <- theta_xy[2]
        alpha_xy <- theta_xy[1]
        rho_tau_expec <- sign(beta_yx) * sqrt(beta_yx * beta_xy * (beta_yx * beta_xy > 0))
        
        e_yx <- mvn_y - alpha_yx - beta_yx*mvn_x
        e_xy <- mvn_x - alpha_xy - beta_xy*mvn_y
        w_yx <- tau * (e_yx > 0) + (1-tau) * (e_yx <= 0)
        w_xy <- tau * (e_xy > 0) + (1-tau) * (e_xy <= 0)
        
        g <- rho_tau_expec/2*matrix(c(1/beta_yx, 1/beta_xy), ncol = 1)
        A <- matrix(c(0,0,1,0,0,0,0,1), ncol = 4)
        det_yx <- mean(w_yx)*mean(w_yx*mvn_x^2) - (mean(w_yx*mvn_x))^2
        det_xy <- mean(w_xy)*mean(w_xy*mvn_y^2) - (mean(w_xy*mvn_y))^2
        W.inv <- matrix(c(c(mean(w_yx*mvn_x^2), -mean(w_yx*mvn_x), 0, 0, -mean(w_yx*mvn_x), mean(w_yx), 0, 0)/det_yx,
                          c(0, 0, mean(w_xy*mvn_y^2), -mean(w_xy*mvn_y), 
                            0, 0, -mean(w_xy*mvn_y), mean(w_xy))/det_xy),ncol = 4)/2
        
        V <- cov(2*cbind(e_yx*w_yx, mvn_x*e_yx*w_yx, e_xy*w_xy, mvn_y*e_xy*w_xy))
        sigma_tau_expec <- sqrt(t(g) %*% A %*% W.inv %*% V %*% W.inv %*% t(A) %*% g)
        coverage_tau_expec <- (rho_tau_expec - 1.96*sigma_tau_expec/sqrt(n) <= rho) & 
          (rho_tau_expec + 1.96*sigma_tau_expec/sqrt(n) >= rho)
        
        # merge both results
        c(rho_tau_quant, sigma_tau_quant, coverage_tau_quant, 
          rho_tau_expec, sigma_tau_expec, coverage_tau_expec)
      }
      
      stopCluster(cl)
      
      time2 <- Sys.time()
      
      
      names(coverage_tau) <- c("quant_corr", "sigma_quant", "coverage_quant", 
                               "expec_corr", "sigma_expec", "coverage_expec")
      
      cat ("Simulation: DGP3, multivariate t with df ", df, "\n", 
           "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
           "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
           "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
           difftime(time2, time1, unit = "mins"), "mins", "\n")
      
      result_obj[[simul_num]] <- list(obj = coverage_tau, 
                                      description = paste("Simulation: DGP3, multivariate t with df ", df, "\n", 
                                                          "n: ", n, "tau:",tau ,"rho: " , rho ,"\n",
                                                          "coverage (quantile corr):", sum(coverage_tau[,3])/1000 ,"\n",
                                                          "coverage (expectile corr):", sum(coverage_tau[,6])/1000 ,"\n",
                                                          difftime(time2, time1, unit = "mins"), "mins", 
                                                          sep = " "))
      simul_num <- simul_num + 1
      saveRDS(result_obj, paste0("/Users/gimjuyeon/Documents/학부연구생/연구/coverage_result_DGP3(multi-t(", df, ")).RDS"))
      
      Sys.sleep(100)
    }
  }
}

