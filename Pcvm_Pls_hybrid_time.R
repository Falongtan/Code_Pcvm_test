# All computations were performed on a workstation equipped with an Intel Core Ultra 7 265K CPU (3.90 GHz) and 48 GB of RAM, running R version 4.5.1. 
# To facilitate efficient computation, simulation tasks were parallelized across 8 cores.

# H11
if(T){
  PCvM1_linear_only <- function(p, n, a, pho) {
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening)
    library(foreach)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }   
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2 
    cl <- makeCluster(cores) 
    registerDoParallel(cl, cores=cores)
    
    results <- foreach(k = 1:s, .combine='rbind', 
                       .packages = c('MASS','glmnet','LassoSIR','VariableScreening'),
                       .export = c("pvalue_integ_Brown")) %dopar% 
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    
        
        #data splitting  
        index_x <- sample(1:n, floor(n/2), replace = FALSE)  
        index_x <- sort(index_x)
        
        x1 <- x[-index_x,] 
        y1 <- as.matrix(y[-index_x,]) 
        x2 <- x[index_x,]
        y2 <- as.matrix(y[index_x,])
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        
        # residuals and beta_0 projections based on lasso and data x1 y1
        lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
        lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- which(is.finite(lasso_beta1) & lasso_beta1 != 0)
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          U1 <- y1 - mean(y1)                                                      
          beta1_hat <- lasso_beta1
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[, index_beta1_non0, drop=FALSE]
          sec_model1 <- glm(y1~x1_sec, family = gaussian)
          sec_beta1 <- unname(sec_model1$coefficients)[-1]
          sec_beta1[!is.finite(sec_beta1)] <- 0
          lasso_beta1[index_beta1_non0] <- sec_beta1
          beta1_hat <- lasso_beta1
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1-pred1                                                           
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        # residuals and beta0-projection based on lasso and data x2 y2
        lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
        lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- which(is.finite(lasso_beta2) & lasso_beta2 != 0)
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          U2 <- y2 - mean(y2)                                                      
          beta2_hat <- lasso_beta2
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0, drop=FALSE]                              
          sec_model2 <- glm(y2~x2_sec, family = gaussian)
          sec_beta2 <- unname(sec_model2$coefficients)[-1]
          sec_beta2[!is.finite(sec_beta2)] <- 0
          lasso_beta2[index_beta2_non0] <- sec_beta2
          beta2_hat <- lasso_beta2
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2                                                           
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        # projections based on screening and lassosir for x2 U2
        screen_num2 <- floor(n2/log(n2))
        if(p <= screen_num2){
          sir_U2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
        }else{
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
          x2_resi_screen <- x2[,index_U2, drop=FALSE]
          sir_U2 <- LassoSIR(x2_resi_screen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
          sdr_Upro2[index_U2,] <- sir_Ubeta2
          sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))
        }
        
        # Schmidt orthogonalization
        pro2 <- cbind(beta2_pro, sdr_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) {
              v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii]
            }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) {
            next
          }
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        
        # ====================================================================
        # PCvM1
        # ====================================================================
        pro2 <- pro2_orth                                                          
        pro2_num <- ncol(pro2)                                                     
        pval_matrix1_Pcvm <- matrix(nrow=1, ncol=pro2_num)                         
        
        for(q in seq_len(pro2_num)){
          x1_pro <- x1 %*% pro2[,q]                                                  
          x1_Indictor <- ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro), 1, 0)
          hat_A1 <- rbind(t(x1_pro), rep(1,n1))  
          Gamma1_inv <- array(0, dim=c(2, 2, n1))
          for (i in 1:n1) { Gamma1_inv[,,i] = MASS::ginv(hat_A1 %*% (t(hat_A1) * x1_Indictor[i,]))}                                                                        
          Integral_1 <- ((rep(1,2) %*% t(U1)) * hat_A1) %*% t(x1_Indictor)              
          
          martingle1_sec <- diag(0,n1)                                               
          for (l in 1:n1){
            martingle1_sec[l,] <- (t(hat_A1[,l]) %*% Gamma1_inv[,,l] %*% Integral_1[,l]) %*% x1_Indictor[l,] 
          }
          martingle_sta1 <- (1/sqrt(n1))*t(U1) %*% x1_Indictor - (1/sqrt(n1)) * colSums(martingle1_sec)  
          
          ordx1_pro <- sort(x1_pro)                                                
          t_1 <- ordx1_pro[floor(0.99 * n1)]                                       
          sigma1_square <- mean(U1^2)   
          F_1 <- (mean(x1_pro <= t_1))^2
          PcvM_martingle1 <- (1/(sigma1_square*F_1)) * mean((x1_pro <= t_1) * (t(martingle_sta1)^2))  
          pval_matrix1_Pcvm[,q] <- pvalue_integ_Brown(PcvM_martingle1)
        }
        pval_cauchy1_Pcvm <- 1 - pcauchy(mean(tan((0.5 - pval_matrix1_Pcvm[1:pro2_num])*pi))) 
        
        pval_cauchy1_Pcvm
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    pow <- mean(results <= 0.05, na.rm=TRUE)
    
    return(list(martingale_cauchy1_power = pow))
  }
  
  Pcvm_linear_run <- function(p, n, a, pho) { 
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening)
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }   
    beta0 <- c(rep(1,5),rep(0,p-5))
    beta1 <- c(rep(1,10),rep(0,p-10)) 
    
    cores <- detectCores(logical=F)-2 
    cl <- makeCluster(cores) 
    registerDoParallel(cl, cores=cores)
    
    pcvm_linear_power <- foreach(k = 1:s, .combine='rbind', 
                                 .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych'),
                                 .export = c("pvalue_integ_Brown")) %dopar% 
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    
        
        index_x <- sample(1:n, floor(n/2),replace = FALSE)  
        index_x <- sort(index_x)
        
        x1 <- x[-index_x,] 
        y1 <- as.matrix(y[-index_x,]) 
        
        x2 <- x[index_x,]
        y2 <- as.matrix(y[index_x,])
        
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        
        # lasso x1 y1
        lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
        lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- which(is.finite(lasso_beta1) & lasso_beta1 != 0)
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          U1 <- y1 - mean(y1)                                                      
          beta1_hat <- lasso_beta1
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[, index_beta1_non0, drop=FALSE]
          sec_model1 <- glm(y1~x1_sec, family = gaussian)
          sec_beta1 <- unname(sec_model1$coefficients)[-1]
          sec_beta1[!is.finite(sec_beta1)] <- 0
          lasso_beta1[index_beta1_non0] <- sec_beta1
          beta1_hat <- lasso_beta1
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1-pred1                                                           
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        screen_num1 <- floor(n1/log(n1))
        if(p <= screen_num1){
          sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sdr_Upro1 <- sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
        }else{
          rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
          index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
          x1_resi_screen <- x1[,index_U1, drop=FALSE]
          sir_U1 <- LassoSIR(x1_resi_screen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sdr_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
          sdr_Upro1[index_U1,] <- sir_Ubeta1
          sdr_Upro1 <- sdr_Upro1/sqrt(colSums(sdr_Upro1^2))
        } 
        
        # lasso x2 y2
        lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
        lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- which(is.finite(lasso_beta2) & lasso_beta2 != 0)
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          U2 <- y2 - mean(y2)                                                      
          beta2_hat <- lasso_beta2
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0, drop=FALSE]                              
          sec_model2 <- glm(y2~x2_sec, family = gaussian)
          sec_beta2 <- unname(sec_model2$coefficients)[-1]
          sec_beta2[!is.finite(sec_beta2)] <- 0
          lasso_beta2[index_beta2_non0] <- sec_beta2
          beta2_hat <- lasso_beta2
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2                                                           
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        screen_num2 <- floor(n2/log(n2))
        if(p <= screen_num2){
          sir_U2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
        }else{
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
          x2_resi_screen <- x2[,index_U2, drop=FALSE]
          sir_U2 <- LassoSIR(x2_resi_screen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
          sdr_Upro2[index_U2,] <- sir_Ubeta2
          sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))
        }
        
        # Schmidt orthogonalization
        pro1 <- cbind(beta1_pro, sdr_Upro1)
        pro1_orth <- matrix(numeric(0), nrow = nrow(pro1), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro1))) {
          v1 <- pro1[, jj]
          if (ncol(pro1_orth) > 0) {
            for (ii in seq_len(ncol(pro1_orth))) {
              v1 <- v1 - sum(v1 * pro1_orth[, ii]) * pro1_orth[, ii]
            }
          }
          v1_norm <- sqrt(sum(v1^2))
          if (!is.finite(v1_norm) || v1_norm < tol) next
          pro1_orth <- cbind( pro1_orth, v1 / v1_norm)
        }
        
        pro2 <- cbind(beta2_pro, sdr_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) {
              v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii]
            }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) next
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        
        # ====================================================================
        # PCvM Test 1
        # ====================================================================
        pro2 <- pro2_orth                                                          
        pro2_num <- ncol(pro2)                                                     
        pval_matrix1_Pcvm <- matrix(nrow=1, ncol=pro2_num)                         
        for(q in seq_len(pro2_num)){
          x1_pro <- x1%*%pro2[,q]                                                  
          x1_Indictor<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          hat_A1 <- rbind(t(x1_pro), rep(1,n1))  
          Gamma1_inv <- array(0, dim=c(2, 2, n1))
          for (i in 1:n1) { 
            Gamma1_inv[,,i] = ginv(hat_A1 %*% (t(hat_A1) * x1_Indictor[i,])) 
          }                                                                        
          Integral_1 <- ((rep(1,2) %*% t(U1)) * hat_A1) %*% t(x1_Indictor)              
          
          martingle1_sec<-diag(0,n1)                                               
          for (l in 1:n1){ 
            martingle1_sec[l,]<- (t(hat_A1[,l]) %*% Gamma1_inv[,,l] %*% Integral_1[,l]) %*% x1_Indictor[l,] 
          }
          martingle_sta1<-(1/sqrt(n1))*t(U1) %*% x1_Indictor - (1/sqrt(n1)) * colSums(martingle1_sec)  
          ordx1_pro <- sort(x1_pro)                                                
          t_1 <- ordx1_pro[floor(0.99 * n1)]                                       
          sigma1_square <- mean(U1^2)   
          F_1 <- (mean(x1_pro <= t_1))^2
          PcvM_martingle1 <- (1/(sigma1_square*F_1)) * mean((x1_pro <= t_1) * (t(martingle_sta1)^2))  
          pval_matrix1_Pcvm[,q] <- pvalue_integ_Brown(PcvM_martingle1)
        }
        
        # ====================================================================
        # PCvM Test 2
        # ====================================================================
        pro1 <- pro1_orth                                                          
        pro1_num <- ncol(pro1)                                                     
        pval_matrix2_Pcvm <- matrix(nrow=1, ncol=pro1_num)                         
        for(qq in seq_len(pro1_num)){
          x2_pro <- x2%*%pro1[,qq]                                                 
          x2_Indictor<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
          hat_A2 <- rbind(t(x2_pro), rep(1,n2))  
          Gamma2_inv <- array(0, dim=c(2, 2, n2))
          for (ii in 1:n2) { 
            Gamma2_inv[,,ii] = ginv(hat_A2 %*% (t(hat_A2) * x2_Indictor[ii,])) 
          }                                                                        
          Integral_2 <- ((rep(1,2) %*% t(U2)) * hat_A2) %*% t(x2_Indictor)              
          martingle2_sec<-diag(0,n2)                                               
          for (ll in 1:n2){ 
            martingle2_sec[ll,]<- (t(hat_A2[,ll]) %*% Gamma2_inv[,,ll] %*% Integral_2[,ll]) %*% x2_Indictor[ll,] 
          }
          martingle_sta2<-(1/sqrt(n2))*t(U2) %*% x2_Indictor - (1/sqrt(n2)) * colSums(martingle2_sec) 
          ordx2_pro <- sort(x2_pro)                                                
          t_2 <- ordx2_pro[floor(0.99 * n2)]                                       
          sigma2_square <- mean(U2^2)                                              
          F_2 <- (mean(x2_pro <= t_2))^2
          PcvM_martingle2 <- (1/(sigma2_square*F_2)) * mean((x2_pro <= t_2) * (t(martingle_sta2)^2))  
          pval_matrix2_Pcvm[,qq] <- pvalue_integ_Brown(PcvM_martingle2)
        }
        
        pval_matrix_Pcvm <- cbind(pval_matrix1_Pcvm, pval_matrix2_Pcvm)
        pval_num <- pro1_num + pro2_num
        pval_cauchy1_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_Pcvm[1:pro2_num])*pi)))  
        pval_cauchy2_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_Pcvm[1:pro1_num])*pi)))  
        pval_cauchy_Pcvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_Pcvm[1:pval_num])*pi)))   
        
        c(pval_cauchy1_Pcvm, pval_cauchy2_Pcvm, pval_cauchy_Pcvm)
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    return(list(martingale_cauchy1_power = mean(pcvm_linear_power[,1]<=0.05,na.rm=TRUE), 
                martingale_cauchy2_power = mean(pcvm_linear_power[,2]<=0.05,na.rm=TRUE), 
                martingale_cauchy_power = mean(pcvm_linear_power[,3]<=0.05,na.rm=TRUE)))
  }
  
  Hybrid_linear_run <- function(p, n, a, pho) { 
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening)
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }   
    beta0 <- c(rep(1,5),rep(0,p-5))
    beta1 <- c(rep(1,10),rep(0,p-10)) 
    
    cores <- detectCores(logical=F)-2 
    cl <- makeCluster(cores) 
    registerDoParallel(cl, cores=cores)
    
    hybrid_linear_power <- foreach(k = 1:s, .combine='c',
                                   .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych'),
                                   .export = c("pvalue_integ_Brown")) %dopar% 
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    
        
        index_x <- sample(1:n, floor(n/2),replace = FALSE)  
        index_x <- sort(index_x)
        
        x1 <- x[-index_x,] 
        y1 <- as.matrix(y[-index_x,]) 
        x2 <- x[index_x,]
        y2 <- as.matrix(y[index_x,])
        n1 <- nrow(x1); n2 <- nrow(x2)
        
        # lasso x1 y1
        lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
        lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- which(is.finite(lasso_beta1) & lasso_beta1 != 0)
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          U1 <- y1 - mean(y1)                                                      
          beta1_hat <- lasso_beta1
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[, index_beta1_non0, drop=FALSE]
          sec_model1 <- glm(y1~x1_sec, family = gaussian)
          sec_beta1 <- unname(sec_model1$coefficients)[-1]
          sec_beta1[!is.finite(sec_beta1)] <- 0
          lasso_beta1[index_beta1_non0] <- sec_beta1
          beta1_hat <- lasso_beta1
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1-pred1                                                           
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        screen_num1 <- floor(n1/log(n1))
        if(p <= screen_num1){
          sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sdr_Upro1 <- sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
        }else{
          rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
          index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
          x1_resi_screen <- x1[,index_U1, drop=FALSE]
          sir_U1 <- LassoSIR(x1_resi_screen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sdr_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
          sdr_Upro1[index_U1,] <- sir_Ubeta1
          sdr_Upro1 <- sdr_Upro1/sqrt(colSums(sdr_Upro1^2))
        } 
        
        # if(p <= screen_num1){
        #   sir_y1 <- LassoSIR(x1, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        #   sir_ybeta1<- sir_y1$beta
        #   sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
        # }else{
        #   rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                      
        #   index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        
        #   x1_yscreen <- x1[,index_y1]
        #   sir_y1 <- LassoSIR(x1_yscreen, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        #   sir_ybeta1<- sir_y1$beta
        #   sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
        #   sdr_ypro1[index_y1,] <- sir_ybeta1
        #   sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
        # }
        
        # lasso x2 y2
        lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
        lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- which(is.finite(lasso_beta2) & lasso_beta2 != 0)
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          U2 <- y2 - mean(y2)                                                      
          beta2_hat <- lasso_beta2
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0, drop=FALSE]                              
          sec_model2 <- glm(y2~x2_sec, family = gaussian)
          sec_beta2 <- unname(sec_model2$coefficients)[-1]
          sec_beta2[!is.finite(sec_beta2)] <- 0
          lasso_beta2[index_beta2_non0] <- sec_beta2
          beta2_hat <- lasso_beta2
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2                                                           
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        screen_num2 <- floor(n2/log(n2))
        if(p <= screen_num2){
          sir_U2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
        }else{
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
          x2_resi_screen <- x2[,index_U2, drop=FALSE]
          sir_U2 <- LassoSIR(x2_resi_screen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
          sir_Ubeta2<- sir_U2$beta
          sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
          sdr_Upro2[index_U2,] <- sir_Ubeta2
          sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))
        }
        
        # if(p <= screen_num2){
        #   sir_y2 <- LassoSIR(x2, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        #   sir_ybeta2<- sir_y2$beta
        #   sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                      
        # }else{ 
        #   rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                      
        #   index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        
        #   x2_yscreen <- x2[,index_y2]
        #   sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        #   sir_ybeta2 <- sir_ylasso2$beta
        #   sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
        #   sdr_ypro2[index_y2,] <- sir_ybeta2
        #   sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                        
        # } 
        
        pro1 <- cbind(beta1_pro, sdr_Upro1)
        pro1_orth <- matrix(numeric(0), nrow = nrow(pro1), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro1))) {
          v1 <- pro1[, jj]
          if (ncol(pro1_orth) > 0) {
            for (ii in seq_len(ncol(pro1_orth))) {
              v1 <- v1 - sum(v1 * pro1_orth[, ii]) * pro1_orth[, ii]
            }
          }
          v1_norm <- sqrt(sum(v1^2))
          if (!is.finite(v1_norm) || v1_norm < tol) next
          pro1_orth <- cbind( pro1_orth, v1 / v1_norm)
        }
        
        pro2 <- cbind(beta2_pro, sdr_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) {
              v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii]
            }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) next
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        
        # ====================================================================
        # PCvM1
        # ====================================================================
        pro2_pcvm <- pro2_orth                                                          
        pro2_num <- ncol(pro2_pcvm)                                                     
        pval_matrix1_Pcvm <- matrix(nrow=1, ncol=pro2_num)                         
        for(q in seq_len(pro2_num)){
          x1_pro <- x1%*%pro2_pcvm[,q]                                                  
          x1_Indictor<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          hat_A1 <- rbind(t(x1_pro), rep(1,n1))  
          Gamma1_inv <- array(0, dim=c(2, 2, n1))
          for (i in 1:n1) { 
            Gamma1_inv[,,i] = ginv(hat_A1 %*% (t(hat_A1) * x1_Indictor[i,])) 
          }                                                                        
          Integral_1 <- ((rep(1,2)%*%t(U1))*hat_A1)%*%t(x1_Indictor)             
          
          martingle1_sec<-diag(0,n1)                                               
          for (l in 1:n1){ 
            martingle1_sec[l,]<- (t(hat_A1[,l]) %*% Gamma1_inv[,,l] %*% Integral_1[,l]) %*% x1_Indictor[l,] 
          }
          martingle_sta1<-(1/sqrt(n1))*t(U1) %*% x1_Indictor - (1/sqrt(n1)) * colSums(martingle1_sec)  
          ordx1_pro <- sort(x1_pro)                                                
          t_1 <- ordx1_pro[floor(0.99 * n1)]                                       
          sigma1_square <- mean(U1^2)   
          F_1 <- (mean(x1_pro <= t_1))^2
          PcvM_martingle1 <- (1/(sigma1_square*F_1)) * mean((x1_pro <= t_1) * (t(martingle_sta1)^2))  
          pval_matrix1_Pcvm[,q] <- pvalue_integ_Brown(PcvM_martingle1)
        }
        
        # ====================================================================
        # PCvM2
        # ====================================================================
        pro1_pcvm <- pro1_orth                                                          
        pro1_num <- ncol(pro1_pcvm)                                                     
        pval_matrix2_Pcvm <- matrix(nrow=1, ncol=pro1_num)                         
        for(qq in seq_len(pro1_num)){
          x2_pro <- x2%*%pro1_pcvm[,qq]                                                 
          x2_Indictor<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
          hat_A2 <- rbind(t(x2_pro), rep(1,n2))  
          Gamma2_inv <- array(0, dim=c(2, 2, n2))
          for (ii in 1:n2) { 
            Gamma2_inv[,,ii] = ginv(hat_A2 %*% (t(hat_A2) * x2_Indictor[ii,])) 
          }                                                                        
          Integral_2 <- ((rep(1,2)%*%t(U2))*hat_A2)%*%t(x2_Indictor)              
          
          martingle2_sec<-diag(0,n2)                                               
          for (ll in 1:n2){ 
            martingle2_sec[ll,]<- (t(hat_A2[,ll]) %*% Gamma2_inv[,,ll] %*% Integral_2[,ll]) %*% x2_Indictor[ll,] 
          }
          martingle_sta2<-(1/sqrt(n2))*t(U2) %*% x2_Indictor - (1/sqrt(n2)) * colSums(martingle2_sec) 
          ordx2_pro <- sort(x2_pro)                                                
          t_2 <- ordx2_pro[floor(0.99 * n2)]                                       
          sigma2_square <- mean(U2^2)                                              
          F_2 <- (mean(x2_pro <= t_2))^2
          PcvM_martingle2 <- (1/(sigma2_square*F_2)) * mean((x2_pro <= t_2) * (t(martingle_sta2)^2))  
          pval_matrix2_Pcvm[,qq] <- pvalue_integ_Brown(PcvM_martingle2)
        }
        pval_matrix_Pcvm <- cbind(pval_matrix1_Pcvm, pval_matrix2_Pcvm)
        
        # PLS computations
        pro2_pls_num <- ncol(pro2_orth)                                              
        h1 <- n1^(-2/9)                                                          
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_pls_num) 
        errormat1 <- U1%*%t(U1)                                                  
        for(q in seq_len(pro2_pls_num)){
          x_pro1 <- x1%*%pro2_orth[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
        }
        
        pro1_pls_num <- ncol(pro1_orth)                                              
        h2 <- n2^(-2/9)                                                          
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_pls_num) 
        errormat2 <- U2%*%t(U2)                                                  
        for(l in seq_len(pro1_pls_num)){
          x_pro2 <- x2%*%pro1_orth[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2      
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                    
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
        }
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
        
        # Hybrid PCvM + PLS 
        pval_matrix_hybrid <- cbind(pval_matrix_Pcvm, pval_matrix_PLS)
        
        1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid)*pi), na.rm=TRUE))
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    return(list(hybrid_cauchy_power = mean(hybrid_linear_power <= 0.05, na.rm=TRUE)))
  }
  
  PLS_lm_parallel_run <- function(p, n, a, pho){
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening)
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    c_h <- 1
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      sigma <- toeplitz(pho^(0:(p-1)))
    }   
    beta0 <- c(rep(1,5),rep(0,p-5))
    beta1 <- c(rep(1,10),rep(0,p-10))
    
    cores <- detectCores(logical=F)-2
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    
    PLS_power <- foreach(k = 1:s, .combine='rbind', .errorhandling = "remove",
                         .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych')) %dopar%  
      {
          set.seed(k)
          # result <- tryCatch({ 
          x <- mvrnorm(n, mu, sigma) 
          y <- x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)  
          
          # 2. Data splitting  
          index_x <- sample(1:n, floor(n/2), replace = FALSE)  
          index_x <- sort(index_x)
          x1 <- x[-index_x, ]; y1 <- as.matrix(y[-index_x, ]) 
          x2 <- x[index_x, ];  y2 <- as.matrix(y[index_x, ])
          n1 <- nrow(x1); n2 <- nrow(x2)
          
          # 3. Residuals based on Lasso and data x1 y1
          lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
          lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
          index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
          if(length(index_beta1_non0) == 0){
            pred1 <- mean(y1)
          }else{
            x1_sec <- x1[,index_beta1_non0, drop=FALSE]
            sec_model1 <- glm(y1~x1_sec, family = gaussian)
            pred1 <- matrix(unname(predict(sec_model1, newx = x1_sec, type="response")), ncol = 1)
          }
          U1 <- y1 - pred1
          
          # Projections based on screening and lassosir for x1 U1
          screen_num1 <- floor(n1/log(n1))
          sdr_Upro1 <- tryCatch({
            if(p <= screen_num1){
              sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_U1$beta/sqrt(colSums(sir_U1$beta^2))
            }else{
              rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")
              index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]
              sir_U1 <- LassoSIR(x1[,index_U1, drop=FALSE], U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_U1$beta))
              sir_Upro1[index_U1,] <- sir_U1$beta
              sir_Upro1/sqrt(colSums(sir_Upro1^2))
            }
          }, error=function(e) matrix(NA, p, 1))
          
          # Projections based on screening and lassosir for x1 y1
          sdr_ypro1 <- tryCatch({
            if(p <= screen_num1){
              sir_y1 <- LassoSIR(x1, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_y1$beta/sqrt(colSums(sir_y1$beta^2))
            }else{
              rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")
              index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]
              sir_y1 <- LassoSIR(x1[,index_y1, drop=FALSE], y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_y1$beta))
              sir_ypro1[index_y1,] <- sir_y1$beta
              sir_ypro1/sqrt(colSums(sir_ypro1^2))
            }
          }, error=function(e) matrix(NA, p, 1))
          
          # 4. Residuals based on lasso and data x2 y2
          lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
          lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
          index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2) != 0]
          if(length(index_beta2_non0) == 0){
            pred2 <- mean(y2)
          }else{
            x2_sec <- x2[,index_beta2_non0, drop=FALSE]
            sec_model2 <- glm(y2~x2_sec, family = gaussian)
            pred2 <- matrix(unname(predict(sec_model2, newx = x2_sec, type="response")), ncol = 1)
          }
          U2 <- y2 - pred2
          
          # Projections based on screening and lassosir for x2 U2
          screen_num2 <- floor(n2/log(n2))
          sdr_Upro2 <- tryCatch({
            if(p <= screen_num2){
              sir_U2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_U2$beta/sqrt(colSums(sir_U2$beta^2))
            }else{
              rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")
              index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]
              sir_U2 <- LassoSIR(x2[,index_U2, drop=FALSE], U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_U2$beta))
              sir_Upro2[index_U2,] <- sir_U2$beta
              sir_Upro2/sqrt(colSums(sir_Upro2^2))
            }
          }, error=function(e) matrix(NA, p, 1))
          
          # Projections based on screening and lassosir for x2 y2
          sdr_ypro2 <- tryCatch({
            if(p <= screen_num2){
              sir_y2 <- LassoSIR(x2, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_y2$beta/sqrt(colSums(sir_y2$beta^2))
            }else{
              rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")
              index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]
              sir_y2 <- LassoSIR(x2[,index_y2, drop=FALSE], y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
              sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_y2$beta))
              sir_ypro2[index_y2,] <- sir_y2$beta
              sir_ypro2/sqrt(colSums(sir_ypro2^2))
            }
          }, error=function(e) matrix(NA, p, 1))
          
          # 5. Construct PLS test statistic based on x1 y1 and projections from x2 y2
          pro2 <- cbind(sdr_Upro2, sdr_ypro2)
          pro2 <- tryCatch({
            qr_obj2 <- qr(pro2)
            if( qr_obj2$rank == ncol(pro2) ){ pro2 }else{ pro2[, qr_obj2$pivot[1:qr_obj2$rank], drop = FALSE] }
          }, error =function(e) NA)
          pro2 <- as.matrix(pro2)
          pro2_num <- ncol(pro2)
          h1 <- c_h*(n1^(-2/9))
          
          pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) 
          errormat1 <- U1%*%t(U1)
          for(q in seq_len(pro2_num)){
            if(any(is.na(pro2[,q])) ){
              pval_matrix1[,q] <- NA
            }else{
              x_pro1 <- x1%*%pro2[,q]
              x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1
              indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
              kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1
              Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
              pval_matrix1[,q] <- 1-pnorm(Tn1)
            }
          }
          
          # 6. Construct PLS test statistics based on x2 y2 and projections from x1 y1
          pro1 <- cbind(sdr_Upro1, sdr_ypro1)
          pro1 <- tryCatch({
            qr_obj1 <- qr(pro1)
            if( qr_obj1$rank == ncol(pro1) ){ pro1 }else{ pro1[, qr_obj1$pivot[1:qr_obj1$rank], drop = FALSE] }
          }, error =function(e) NA)
          pro1 <- as.matrix(pro1)
          pro1_num <- ncol(pro1)
          h2 <- c_h*(n2^(-2/9))
          
          pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) 
          errormat2 <- U2%*%t(U2)
          for(l in seq_len(pro1_num)){
            if(any(is.na(pro1[,l])) ){
              pval_matrix2[,l] <- NA
            }else{
              x_pro2 <- x2%*%pro1[,l]
              x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2
              indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
              kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2
              Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
              pval_matrix2[,l] <- 1-pnorm(Tn2)
            }
          }
 
          pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
          pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
          pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
          pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
          
          pval_min1 <- min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
          pval_min2 <- min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
          pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
          pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi)))
          
          c(pval_fisher1, pval_fisher2, pval_fisher, pval_min1, pval_min2, pval_min)
          
        # }, error = function(e){ c(NA, NA, NA, NA, NA, NA) })
        # return(result)
      }
    stopCluster(cl)
    
    # Power calculation 
    PLS_fisher1_power <- mean(PLS_power[,1]<=0.05,na.rm=TRUE)
    PLS_fisher2_power <- mean(PLS_power[,2]<=0.05,na.rm=TRUE)
    PLS_fisher_power <- mean(PLS_power[,3]<=0.05,na.rm=TRUE)
    PLS_min1_power <- mean(PLS_power[,4],na.rm=TRUE)
    PLS_min2_power <- mean(PLS_power[,5],na.rm=TRUE)
    PLS_min_power <- mean(PLS_power[,6]<=0.05,na.rm=TRUE)
    
    return(list(PLS_fisher1_power=PLS_fisher1_power,
                PLS_fisher2_power=PLS_fisher2_power,
                PLS_fisher_power=PLS_fisher_power,
                PLS_min1_power=PLS_min1_power,
                PLS_min2_power=PLS_min2_power,
                PLS_min_power=PLS_min_power))
  }
  
  grp_rp_lm_parallel_run <- function(p, n, a, pho, target_method = c("RP", "GRP")){
    target_method <- match.arg(target_method)
    library(MASS)
    library(matrixcalc)
    library(RPtests)
    library(GRPtests)
    library(foreach) 
    library(parallel) 
    library(iterators)
    library(doParallel)
    
    s <- 1000   
    mu <- rep(0, p)
    sigma <- if(pho == 0) diag(rep(1,p)) else toeplitz(pho^(0:(p-1)))
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2
    cl <- makeCluster(cores) 
    registerDoParallel(cl, cores=cores)
    rp_grp_power <- foreach(k = 1:s, .combine='rbind', .errorhandling = "remove",
                            .packages = c('MASS','matrixcalc','RPtests','GRPtests')) %dopar%  
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        y <- x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)
        if (target_method == "RP") {
          return(RPtest(x, y, test="nonlin", B=49L, nperms=2, resid_type = "Lasso"))
        } else {
          return(GRPtest(x, y, fam = "gaussian", nsplits = 1, RP_function = NULL, output_all=FALSE, penalize=TRUE))
        }
      }
    stopCluster(cl)
    
    return(list(Power = mean(rp_grp_power <= 0.05, na.rm=TRUE)))
  }
 
  # The choice for the bandwidth for martingale-based test  # 1-Lii  n^{-1/5} 
  bandwidth_choice <- function(x, y){
    n <- nrow(x)
    p <- ncol(x)
    c_h <- seq(0.25, 3, 0.25)
    len_h <- length(c_h)
    h <- c_h*n^(-1/5)
    CV <- rep(0,len_h)                                                          # martingale with beta0 and alpha as projections
    
    for (i in 1:len_h) {
      Ker_inter <- (x %*% rep(1,n) - rep(1,n) %*% t(x))/h[i]                    # kernel function matrix
      #kernel <- exp(-0.5*Ker_inter^2)                                          # Gaussian kernel
      
      Ker_indictor <- ifelse(abs(Ker_inter) <= 1, 1, 0) 
      kernel <- (3/4)*(1-Ker_inter^2)*Ker_indictor
      
      Lii <- diag(kernel) / (rowSums(kernel) + n^(-12))                         # denom <- pmax(1 - Lii, 1e-8)
      # denom <- pmax(1 - mean(Lii), 1e-8)                                        # denom <- pmax(1 - Lii, 1e-8)
      denom <- pmax(1 - Lii, 1e-8)                                             # gcv
      
      R1 <- mean(((y[,1] - (kernel%*%y[,1])/(rowSums(kernel) + n^(-12))) / denom)^2)
      R2 <- mean(((y[,2] - (kernel%*%y[,2])/(rowSums(kernel) + n^(-12))) / denom)^2)
      R3 <- mean(((y[,3] - (kernel%*%y[,3])/(rowSums(kernel) + n^(-12))) / denom)^2)
      CV[i] <- R1+R2+R3 
    }
    index <- which.min(CV)
    if (length(index) == 0) {
      return(NA_real_)
    }
    # cat(sprintf("Selected c_h = %.1f \n", c_h[index]))
    h0 <- c_h[index] * n^(-1/5)              # The bandwidth for the martingale with beta1 and alpha as projections
    return(c(h0))
  }
 
  # The p-value of int_0^1 B(t)^2 dt where B(t) is the standard brown motion
  pvalue_integ_Brown <- function(x){
    p0 <- c(1, 1, 0.9994, 0.9945, 0.9824, 0.9642, 0.9417, 0.9169, 0.891, 0.8648, 0.839, 0.8138, 0.7894, 0.7659,
            0.7434, 0.7218, 0.7012, 0.6814, 0.6626, 0.6445, 0.6273, 0.6108, 0.5949, 0.5798, 0.5652, 0.5513,
            0.5378, 0.5249, 0.5125, 0.5006, 0.489, 0.4779, 0.4672, 0.4568, 0.4468, 0.4371, 0.4278, 0.4187,
            0.4099, 0.4014, 0.3931, 0.3851, 0.3773, 0.3697, 0.3623, 0.3552, 0.3482, 0.3414, 0.3348, 0.3284,
            0.3222, 0.3161, 0.3101, 0.3043, 0.2987, 0.2932, 0.2878, 0.2825, 0.2774, 0.2724, 0.2675, 0.2627,
            0.258, 0.2534, 0.2489, 0.2446, 0.2403, 0.2361, 0.232, 0.228, 0.224, 0.2202, 0.2164, 0.2127, 0.2091,
            0.2056, 0.2021, 0.1987, 0.1953, 0.1921, 0.1889, 0.1857, 0.1826, 0.1796, 0.1767, 0.1738, 0.1709,
            0.1681, 0.1654, 0.1627, 0.16, 0.1574, 0.1549, 0.1524, 0.1499, 0.1475, 0.1451, 0.1428, 0.1405, 0.1383,
            0.1361, 0.1339, 0.1318, 0.1297, 0.1277, 0.1257, 0.1237, 0.1218, 0.1198, 0.118, 0.1161, 0.1143, 0.1125,
            0.1108, 0.1091, 0.1074, 0.1057, 0.1041, 0.1025, 0.1009, 0.0994, 0.0978, 0.0963, 0.0949, 0.0934, 0.092,
            0.0906, 0.0892, 0.0878, 0.0865, 0.0852, 0.0839, 0.0826, 0.0814, 0.0802, 0.0789, 0.0778, 0.0766, 0.0754,
            0.0743, 0.0732, 0.0721, 0.071, 0.07, 0.0689, 0.0679, 0.0669, 0.0659, 0.0649, 0.0639, 0.063, 0.0621, 0.0611,
            0.0602, 0.0593, 0.0585, 0.0576, 0.0568, 0.0559, 0.0551, 0.0543, 0.0535, 0.0527, 0.0519, 0.0512, 0.0504, 0.0497,
            0.049, 0.0482, 0.0475, 0.0469, 0.0462, 0.0455, 0.0448, 0.0442, 0.0435, 0.0429, 0.0423, 0.0417, 0.0411, 0.0405,
            0.0399, 0.0393, 0.0388, 0.0382, 0.0376, 0.0371, 0.0366, 0.036, 0.0355, 0.035, 0.0345, 0.034, 0.0335, 0.033, 0.0326,
            0.0321, 0.0317, 0.0312, 0.0308, 0.0303, 0.0299, 0.0295, 0.029, 0.0286, 0.0282, 0.0278, 0.0274, 0.027, 0.0266, 0.0263,
            0.0259, 0.0255, 0.0252, 0.0248, 0.0245, 0.0241, 0.0238, 0.0234, 0.0231, 0.0228, 0.0225, 0.0221, 0.0218, 0.0215,
            0.0212, 0.0209, 0.0206, 0.0203, 0.0201, 0.0198, 0.0195, 0.0192, 0.019, 0.0187, 0.0184, 0.0182, 0.0179, 0.0177,
            0.0174, 0.0172, 0.0169, 0.0167, 0.0165, 0.0162, 0.016, 0.0158, 0.0156, 0.0153, 0.0151, 0.0149, 0.0147, 0.0145,
            0.0143, 0.0141, 0.0139, 0.0137, 0.0135, 0.0133, 0.0132, 0.013, 0.0128, 0.0126, 0.0124, 0.0123, 0.0121, 0.0119,
            0.0118, 0.0116, 0.0114, 0.0113, 0.0111, 0.011, 0.0108, 0.0107, 0.0105, 0.0104, 0.0102, 0.0101, 0.01, 0.0098,
            0.0097, 0.0096, 0.0094, 0.0093, 0.0092, 0.009, 0.0089, 0.0088, 0.0087, 0.0086, 0.0084, 0.0083, 0.0082, 0.0081,
            0.008, 0.0079, 0.0078, 0.0077, 0.0076, 0.0075, 0.0074, 0.0073, 0.0072, 0.0071, 0.007, 0.0069, 0.0068, 0.0067,
            0.0066, 0.0065, 0.0064, 0.0063, 0.0062, 0.0062, 0.0061, 0.006, 0.0059, 0.0058, 0.0057, 0.0057, 0.0056, 0.0055,
            0.0054, 0.0054, 0.0053, 0.0052, 0.0052, 0.0051, 0.005, 0.0049, 0.0049, 0.0048, 0.0047, 0.0047, 0.0046, 0.0046,
            0.0045, 0.0044, 0.0044, 0.0043, 0.0043, 0.0042, 0.0041, 0.0041, 0.004, 0.004, 0.0039, 0.0039, 0.0038, 0.0038,
            0.0037, 0.0037, 0.0036, 0.0036, 0.0035, 0.0035, 0.0034, 0.0034, 0.0033, 0.0033, 0.0032, 0.0032, 0.0032, 0.0031,
            0.0031, 0.003, 0.003, 0.003, 0.0029, 0.0029, 0.0028, 0.0028, 0.0028, 0.0027, 0.0027, 0.0026, 0.0026, 0.0026,
            0.0025, 0.0025, 0.0025, 0.0024, 0.0024, 0.0024, 0.0023, 0.0023, 0.0023, 0.0023, 0.0022, 0.0022, 0.0022, 0.0021,
            0.0021, 0.0021, 0.0021, 0.002, 0.002, 0.002, 0.0019)
    x0 <- seq(0,3.99,by=0.01)
    if(is.na(x)){
      p <- NA
    }else if(x > 3.99 & x<8){
      p <- -0.0019*x/4.01 + 8*0.0019/4.01
    }else if (x >= 8){
      p <- 0
    }else if(x==0){
      p <- 1
    }else{
      i <- sum(x0 < x)
      p <- (p0[i+1]-p0[i])*x/(x0[i+1]-x0[i])+(x0[i+1]*p0[i]-p0[i+1]*x0[i])/(x0[i+1]-x0[i])
    }
    return(p)
  }
  
}

# H21
if(T){
  PCvM1_logit_only <- function(p, n, a, pho){
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening)
    library(foreach)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2 
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    
    results <- foreach(k = 1:s, .combine='rbind', 
                       .packages = c('MASS', 'glmnet','LassoSIR','VariableScreening'),
                       .export = c("pvalue_integ_Brown", "bandwidth_choice")) %dopar%
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                        
        pr <- 1/(1 + exp(-z))
        y <- matrix(rbinom(n, 1, pr),ncol = 1)
        
        index_x <- sample(1:n, floor(n/2), replace = FALSE)
        index_x <- sort(index_x)
        x1 <- x[-index_x,]
        x2 <- x[index_x,]
        y1 <- as.matrix(y[-index_x,])
        y2 <- as.matrix(y[index_x,])
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        
 
        lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
        beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
          U1 <- y1 - pred1                                                         
          beta1_hat <- beta1_lasso
          intercept1 <- coef(lasso_model1, s = "lambda.min")[1]                    
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[,index_beta1_non0]                                          
          sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
          beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
          beta1_hat <- beta1_lasso
          intercept1 <- unname(sec_model1$coefficients)[1]
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1 - pred1                                                         
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
        beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso)!=0]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
          U2 <- y2-pred2                                                           
          beta2_hat <- beta2_lasso
          intercept2 <- coef(lasso_model2, s = "lambda.min")[1]
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0]
          sec_model2 <- glm(y2~x2_sec,family = binomial(link = "logit"))
          beta2_lasso[index_beta2_non0] <- unname(sec_model2$coefficients)[-1]
          beta2_hat <- beta2_lasso
          intercept2 <- unname(sec_model2$coefficients)[1]
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        screen_num2 <- floor(n2/log(n2)/2)
        sir_Upro2 <- tryCatch({
          if(p <= screen_num2){
            sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2 <- sir_Ulasso2$beta
            sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
          }else{
            rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
            index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
            x2_Uscreen <- x2[,index_U2]
            sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2 <- sir_Ulasso2$beta
            sir_Upro2_temp <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
            sir_Upro2_temp[index_U2,] <- sir_Ubeta2
            sir_Upro2_temp/sqrt(colSums(sir_Upro2_temp^2))
          }}, error =function(e){ matrix(NA,p,1) })
        
        # Schmidt orthogonalization
        pro2 <- cbind(beta2_pro, sir_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) {
              v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii]
            }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) {
            next
          }
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        sir_Upro2_orth <- pro2_orth[, -1, drop = FALSE]
        
        # ====================================================================
        # PCvM1
        # ====================================================================
        pro2_num <- ncol(sir_Upro2_orth)                                                                                     
        pval_matrix1_PCvm <- matrix(nrow=1, ncol=pro2_num+1)                                                               
        deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
        
        for(q in seq_len(pro2_num)){
          x1_pro <- x1 %*% sir_Upro2_orth[,q]
          band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
          h1 <- bandwidth_choice(x1_pro, band_y1)                                                                            
          
          Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                                  
          Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
          kernel_1 <- (3/4)*(1 - Ker_inter1^2)*Ker_indictor1
          
          Indictor_1 <- ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          
          A1_1 <- ((kernel_1 %*% deri_link1) * x1_pro) / (kernel_1 %*% U1^2 + n1^(-12))                                                
          A1_2 <- kernel_1 %*% ((x1 %*% beta1_hat) * deri_link1) / (kernel_1 %*% U1^2 + n1^(-12))                                     
          A1_3 <- (kernel_1 %*% deri_link1) / (kernel_1 %*% U1^2 + n1^(-12))                                                         
          hat_A1 <- rbind(t(A1_1), t(A1_2), t(A1_3))
          
          Gamma_inv1 <- array(0, dim=c(3, 3, n1))
          for (i in 1:n1) {
            Gamma_inv1[,,i] = MASS::ginv(((rep(1,3) %*% t(U1)^2)*hat_A1) %*% (t(hat_A1) * (as.matrix(Indictor_1[i,]) %*% rep(1,3))))
          }
          integral_1 <- ((rep(1,3) %*% t(U1))*hat_A1) %*% t(Indictor_1)
          
          martingle_sec1 <- diag(0,n1)  
          for (l in 1:n1){
            martingle_sec1[l,] <- (U1[l]^2 * t(hat_A1[,l])) %*% Gamma_inv1[,,l] %*% integral_1[,l] %*% Indictor_1[l,]
          }
          martingle_sta1 <- (1/sqrt(n1))*t(U1) %*% Indictor_1 - (1/sqrt(n1)) * colSums(martingle_sec1)   
          
          ordx1_pro <- sort(x1_pro)  
          t_1 <- ordx1_pro[floor(0.99 * n1)]
          psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
          PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        
          pval_matrix1_PCvm[,q] <- pvalue_integ_Brown(PCvM_martingle1)
        }
        
        # beta2_pro 
        x1_betapro <- x1 %*% beta2_pro
        beta_Indictor1 <- ifelse(x1_betapro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_betapro), 1, 0)
        hat_beta_A1 <- rbind(t(x1_betapro), matrix(1,1,n1))
        
        Gamma_beta_inv1 <- array(0, dim=c(2, 2, n1))
        for (i in 1:n1) {
          Gamma_beta_inv1[,,i] = MASS::ginv(((rep(1,2) %*% t(U1)^2)*hat_beta_A1) %*% (t(hat_beta_A1) * (as.matrix(beta_Indictor1[i,]) %*% rep(1,2))))
        }
        integral_beta_1 <- ((rep(1,2) %*% t(U1))*hat_beta_A1) %*% t(beta_Indictor1)                      
        
        martingle_beta_sec1 <- diag(0,n1)                                                                                      
        for (l in 1:n1){
          martingle_beta_sec1[l,] <- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,]
        }
        martingle_beta_sta1 <- (1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      
        
        ordx1_betapro <- sort(x1_betapro)  
        beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
        psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
        PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  
        pval_matrix1_PCvm[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1)
        
        pval_cauchy1_PCvm <- 1 - pcauchy(mean(tan((0.5 - pval_matrix1_PCvm[1:(pro2_num+1)])*pi)))
        
        pval_cauchy1_PCvm
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    pow <- mean(results <= 0.05, na.rm=TRUE)
    
    return(list(martingale_cauchy1_power = pow))
  }

  Pcvm_logit_run <- function(p, n, a, pho){
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(harmonicmeanp)
    library(VariableScreening)
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2 
    random_port <- sample(10000:60000, 1)
    cl <- makeCluster(cores, port = random_port)
    registerDoParallel(cl, cores=cores)
    
    pcvm_power <- foreach(k = 1:s, .combine='rbind', 
                          .packages = c('MASS','harmonicmeanp', 'glmnet','LassoSIR','VariableScreening','psych'),
                          .export = c("pvalue_integ_Brown",'bandwidth_choice')) %dopar%
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                              
        pr <- 1/(1 + exp(-z))
        y <- matrix(rbinom(n, 1, pr),ncol = 1)
        
        index_x<-sample(1:n,floor(n/2),replace = FALSE)
        index_x<-sort(index_x)
        
        x1 <- x[-index_x,]
        x2 <- x[index_x,]
        y1 <- as.matrix(y[-index_x,])
        y2 <- as.matrix(y[index_x,])
        
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        
        lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
        beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
          U1 <- y1 - pred1                                                         
          beta1_hat <- beta1_lasso
          intercept1 <- coef(lasso_model1, s = "lambda.min")[1]                    
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[,index_beta1_non0]                                          
          sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
          beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
          beta1_hat <- beta1_lasso
          intercept1 <- unname(sec_model1$coefficients)[1]
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1 - pred1                                                         
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        screen_num1 <- floor(n1/log(n1))
        sir_Upro1 <- tryCatch({
          if(p <= screen_num1){
            sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta1<- sir_U1$beta
            sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
          }else{
            rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
            index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
            x1_Uscreen <- x1[,index_U1]
            sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta1<- sir_Ulasso1$beta
            sir_Upro1_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
            sir_Upro1_t[index_U1,] <- sir_Ubeta1
            sir_Upro1_t/sqrt(colSums(sir_Upro1_t^2))
          }}, error =function(e){ matrix(NA,p,1) })
        
        lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
        beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso)!=0]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
          U2 <- y2-pred2                                                           
          beta2_hat <- beta2_lasso
          intercept2 <- coef(lasso_model2, s = "lambda.min")[1]
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0]
          sec_model2 <- glm(y2~x2_sec,family = binomial(link = "logit"))
          beta2_lasso[index_beta2_non0] <- unname(sec_model2$coefficients)[-1]
          beta2_hat <- beta2_lasso
          intercept2 <- unname(sec_model2$coefficients)[1]
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        screen_num2 <- floor(n2/log(n2))
        sir_Upro2 <- tryCatch({
          if(p <= screen_num2){
            sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2<- sir_Ulasso2$beta
            sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
          }else{
            rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
            index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
            x2_Uscreen <- x2[,index_U2]
            sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2 <- sir_Ulasso2$beta
            sir_Upro2_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
            sir_Upro2_t[index_U2,] <- sir_Ubeta2
            sir_Upro2_t/sqrt(colSums(sir_Upro2_t^2))
          }}, error =function(e){matrix(NA,p,1) })
        
        pro1 <- cbind(beta1_pro, sir_Upro1)                                  
        pro1_orth <- matrix(numeric(0), nrow = nrow(pro1), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro1))) {
          v1 <- pro1[, jj]
          if (ncol(pro1_orth) > 0) {
            for (ii in seq_len(ncol(pro1_orth))) { v1 <- v1 - sum(v1 * pro1_orth[, ii]) * pro1_orth[, ii] }
          }
          v1_norm <- sqrt(sum(v1^2))
          if (!is.finite(v1_norm) || v1_norm < tol) next
          pro1_orth <- cbind( pro1_orth, v1 / v1_norm)
        }
        sir_Upro1_orth <- pro1_orth[, -1, drop = FALSE]
        
        pro2 <- cbind(beta2_pro, sir_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) { v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii] }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) next
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        sir_Upro2_orth <- pro2_orth[, -1, drop = FALSE]
        
        pro2_num <- ncol(sir_Upro2_orth)                                                                                     
        pval_matrix1_PCvm <- matrix(nrow=1, ncol=pro2_num+1)                                                               
        deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
        for(q in seq_len(pro2_num)){
          x1_pro <- x1%*%sir_Upro2_orth[,q]
          band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
          h1 <- bandwidth_choice(x1_pro, band_y1)                                                                            
          Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                                  
          Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
          kernel_1 <- (3/4)*(1-Ker_inter1^2)*Ker_indictor1
          Indictor_1<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          A1_1 <- ((kernel_1%*%deri_link1)*x1_pro)/(kernel_1%*%U1^2+n1^(-12))                                                
          A1_2 <- kernel_1%*%((x1 %*% beta1_hat)*deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                     
          A1_3 <- (kernel_1%*%deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                                         
          hat_A1 <- rbind(t(A1_1), t(A1_2), t(A1_3))
          Gamma_inv1 <- array(0, dim=c(3, 3, n1))
          for (i in 1:n1) { Gamma_inv1[,,i] = ginv(((rep(1,3)%*%t(U1)^2)*hat_A1)%*%(t(hat_A1) * (as.matrix(Indictor_1[i,]) %*%rep(1,3)))) }
          integral_1 <- ((rep(1,3)%*%t(U1))*hat_A1)%*%t(Indictor_1)
          martingle_sec1<-diag(0,n1)  
          for (l in 1:n1){ martingle_sec1[l,]<- (U1[l]^2 * t(hat_A1[,l])) %*% Gamma_inv1[,,l] %*% integral_1[,l] %*% Indictor_1[l,] }
          martingle_sta1<-(1/sqrt(n1))*t(U1) %*% Indictor_1 - (1/sqrt(n1)) * colSums(martingle_sec1)   
          ordx1_pro <- sort(x1_pro)  
          t_1 <- ordx1_pro[floor(0.99 * n1)]
          psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
          PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        
          pval_matrix1_PCvm[,q] <- pvalue_integ_Brown(PCvM_martingle1)
        }
        x1_betapro <- x1%*%beta2_pro
        beta_Indictor1<-ifelse(x1_betapro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_betapro),1,0)
        hat_beta_A1 <- rbind(t(x1_betapro), matrix(1,1,n1))
        Gamma_beta_inv1 <- array(0, dim=c(2, 2, n1))
        for (i in 1:n1) { Gamma_beta_inv1[,,i] = ginv(((rep(1,2)%*%t(U1)^2)*hat_beta_A1)%*%(t(hat_beta_A1) * (as.matrix(beta_Indictor1[i,]) %*%rep(1,2)))) }
        integral_beta_1 <- ((rep(1,2)%*%t(U1))*hat_beta_A1)%*%t(beta_Indictor1)                      
        martingle_beta_sec1<-diag(0,n1)                                                                                      
        for (l in 1:n1){ martingle_beta_sec1[l,]<- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,] }
        martingle_beta_sta1<-(1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      
        ordx1_betapro <- sort(x1_betapro)  
        beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
        psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
        PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  
        pval_matrix1_PCvm[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1)
        
        pro1_num <- ncol(sir_Upro1_orth)                                                                                     
        pval_matrix2_PCvm <- matrix(nrow=1, ncol=pro1_num+1)                                                               
        deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
        for(qq in seq_len(pro1_num)){
          x2_pro <- x2%*%sir_Upro1_orth[,qq]
          band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
          h2 <- bandwidth_choice(x2_pro, band_y2)                                                                            
          Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                                                  
          Ker_indictor2 <- ifelse(abs(Ker_inter2) <= 1, 1, 0)
          kernel_2 <- (3/4)*(1-Ker_inter2^2)*Ker_indictor2
          Indictor_2<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
          A2_1 <- ((kernel_2%*%deri_link2)*x2_pro)/(kernel_2%*%U2^2+n2^(-12))                                                
          A2_2 <- kernel_2%*%((x2 %*% beta2_hat)*deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                     
          A2_3 <- (kernel_2%*%deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                                         
          hat_A2 <- rbind(t(A2_1), t(A2_2), t(A2_3))
          Gamma_inv2 <- array(0, dim=c(3, 3, n2))
          for (ii in 1:n2){ Gamma_inv2[,,ii] = ginv(((rep(1,3)%*%t(U2)^2)*hat_A2)%*%(t(hat_A2) * (as.matrix(Indictor_2[ii,]) %*%rep(1,3)))) }
          integral_2 <- ((rep(1,3)%*%t(U2))*hat_A2)%*%t(Indictor_2)
          martingle_sec2<-diag(0,n2)                                                                                       
          for (ll in 1:n2){ martingle_sec2[ll,]<- (U2[ll]^2 * t(hat_A2[,ll])) %*% Gamma_inv2[,,ll] %*% integral_2[,ll] %*% Indictor_2[ll,] }
          martingle_sta2<-(1/sqrt(n2))*t(U2) %*% Indictor_2 - (1/sqrt(n2)) * colSums(martingle_sec2)  
          ordx2_pro <- sort(x2_pro)                                                                                        
          t_2 <- ordx2_pro[floor(0.99 * n2)]
          psi_2 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2))
          PCvM_martingle2 <- (1/psi_2^2) * mean(U2^2 * (x2_pro <= t_2) * (t(martingle_sta2)^2))       
          pval_matrix2_PCvm[,qq] <- pvalue_integ_Brown(PCvM_martingle2)
        }
        x2_betapro <- x2%*%beta1_pro
        beta_Indictor2<-ifelse(x2_betapro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_betapro),1,0)
        hat_beta_A2 <- rbind(t(x2_betapro), matrix(1,1,n2))
        Gamma_beta_inv2 <- array(0, dim=c(2, 2, n2))
        for (ii in 1:n2) { Gamma_beta_inv2[,,ii] = ginv(((rep(1,2)%*%t(U2)^2)*hat_beta_A2)%*%(t(hat_beta_A2) * (as.matrix(beta_Indictor2[ii,]) %*%rep(1,2)))) }
        integral_beta_2 <- ((rep(1,2)%*%t(U2))*hat_beta_A2)%*%t(beta_Indictor2)
        martingle_beta_sec2<-diag(0,n2)                                                                                    
        for (ll in 1:n2){ martingle_beta_sec2[ll,]<- (U2[ll]^2 * t(hat_beta_A2[,ll])) %*% Gamma_beta_inv2[,,ll] %*% integral_beta_2[,ll] %*% beta_Indictor2[ll,] }
        martingle_beta_sta2<-(1/sqrt(n2))*t(U2) %*% beta_Indictor2 - (1/sqrt(n2)) * colSums(martingle_beta_sec2)   
        ordx2_betapro <- sort(x2_betapro)                                                                                  
        beta_t2 <- ordx2_betapro[floor(0.99 * n2)]
        psi_beta_2 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2))
        PCvM_beta_martingle2 <- (1/psi_beta_2^2) * mean(U2^2 * (x2_betapro <= beta_t2) * (t(martingle_beta_sta2)^2))   
        pval_matrix2_PCvm[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2)
        
        pval_matrix_PCvm <- cbind(pval_matrix1_PCvm, pval_matrix2_PCvm)
        pval_num <- ncol(pval_matrix_PCvm)
        pval_cauchy1_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm[1:(pro2_num+1)])*pi)))         
        pval_cauchy2_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm[1:(pro1_num+1)])*pi)))         
        pval_cauchy_PCvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm[1:pval_num])*pi)))              
        
        c(pval_cauchy1_PCvm, pval_cauchy2_PCvm, pval_cauchy_PCvm)
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    return(list(martingale_cauchy1_power = mean(pcvm_power[,1] <= 0.05,na.rm=TRUE),
                martingale_cauchy2_power = mean(pcvm_power[,2] <= 0.05,na.rm=TRUE),
                martingale_cauchy_power  = mean(pcvm_power[,3] <= 0.05,na.rm=TRUE)))
  }
 
  Hybrid_logit_run <- function(p, n, a, pho){
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(harmonicmeanp)
    library(VariableScreening)
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    mu <- rep(0, p)
    if(pho == 0){
      sigma <- diag(rep(1,p))
    }else{
      v <- pho^(0:(p-1))
      sigma <- toeplitz(v)
    }
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2 
    random_port <- sample(10000:60000, 1)
    cl <- makeCluster(cores, port = random_port)
    registerDoParallel(cl, cores=cores)
    
    hybrid_power <- foreach(k = 1:s, .combine='c',#.combine='rbind', 
                            .packages = c('MASS','harmonicmeanp', 'glmnet','LassoSIR','VariableScreening','psych'),
                            .export = c("pvalue_integ_Brown",'bandwidth_choice')) %dopar%
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                              
        pr <- 1/(1 + exp(-z))
        y <- matrix(rbinom(n, 1, pr),ncol = 1)
        
        index_x<-sample(1:n,floor(n/2),replace = FALSE)
        index_x<-sort(index_x)
        
        x1 <- x[-index_x,]
        x2 <- x[index_x,]
        y1 <- as.matrix(y[-index_x,])
        y2 <- as.matrix(y[index_x,])
        
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        
        # Model 1
        lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
        beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
        len_beta1_non0 <- length(index_beta1_non0)
        if(len_beta1_non0 == 0){
          pred1 <- mean(y1)
          U1 <- y1 - pred1                                                         
          beta1_hat <- beta1_lasso
          intercept1 <- coef(lasso_model1, s = "lambda.min")[1]                    
          beta1_pro <- rep(1,p)/sqrt(p)
        }else{
          x1_sec <- x1[,index_beta1_non0]                                          
          sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
          beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
          beta1_hat <- beta1_lasso
          intercept1 <- unname(sec_model1$coefficients)[1]
          pred1 = predict(sec_model1, type="response")
          pred1 = matrix(unname(pred1), ncol = 1)
          U1 <- y1 - pred1                                                         
          beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
          beta1_pro <- as.matrix(beta1_pro)
        }
        
        screen_num1 <- floor(n1/log(n1))
        sir_Upro1 <- tryCatch({
          if(p <= screen_num1){
            sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta1<- sir_U1$beta
            sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
          }else{
            rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
            index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
            x1_Uscreen <- x1[,index_U1]
            sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta1<- sir_Ulasso1$beta
            sir_Upro1_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
            sir_Upro1_t[index_U1,] <- sir_Ubeta1
            sir_Upro1_t/sqrt(colSums(sir_Upro1_t^2))
          }}, error =function(e){ matrix(NA,p,1) })
        
        # sir_ypro1 <- tryCatch({
        #   if(p <= screen_num1){
        #     sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
        #     sir_ybeta1<- sir_ylasso1$beta
        #     sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
        #   }else{
        #     rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                      
        #     index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        
        #     x1_yscreen <- x1[,index_y1]
        #     sir_ylasso1 <- LassoSIR(x1_yscreen, y1, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
        #     sir_ybeta1<- sir_ylasso1$beta
        #     sir_ypro1_t <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
        #     sir_ypro1_t[index_y1,] <- sir_ybeta1
        #     sir_ypro1_t/sqrt(colSums(sir_ypro1_t^2))
        #   }}, error =function(e){ matrix(NA,p,1) })
        
        # Model 2
        lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
        beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso)!=0]
        len_beta2_non0 <- length(index_beta2_non0)
        if(len_beta2_non0 == 0){
          pred2 <- mean(y2)
          U2 <- y2-pred2                                                           
          beta2_hat <- beta2_lasso
          intercept2 <- coef(lasso_model2, s = "lambda.min")[1]
          beta2_pro <- rep(1,p)/sqrt(p)
        }else{
          x2_sec <- x2[,index_beta2_non0]
          sec_model2 <- glm(y2~x2_sec,family = binomial(link = "logit"))
          beta2_lasso[index_beta2_non0] <- unname(sec_model2$coefficients)[-1]
          beta2_hat <- beta2_lasso
          intercept2 <- unname(sec_model2$coefficients)[1]
          pred2 = predict(sec_model2, type="response")
          pred2 = matrix(unname(pred2), ncol = 1)
          U2 <- y2-pred2
          beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
          beta2_pro <- as.matrix(beta2_pro)
        }
        
        screen_num2 <- floor(n2/log(n2))
        sir_Upro2 <- tryCatch({
          if(p <= screen_num2){
            sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2<- sir_Ulasso2$beta
            sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
          }else{
            rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
            index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
            x2_Uscreen <- x2[,index_U2]
            sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ubeta2 <- sir_Ulasso2$beta
            sir_Upro2_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
            sir_Upro2_t[index_U2,] <- sir_Ubeta2
            sir_Upro2_t/sqrt(colSums(sir_Upro2_t^2))
          }}, error =function(e){matrix(NA,p,1) })
        
        # sir_ypro2 <- tryCatch({
        #   if(p <= screen_num2){
        #     sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
        #     sir_ybeta2<- sir_ylasso2$beta
        #     sir_ybeta2/sqrt(colSums(sir_ybeta2^2))
        #   }else{
        #     rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                      
        #     index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        
        #     x2_yscreen <- x2[,index_y2]
        #     sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
        #     sir_ybeta2 <- sir_ylasso2$beta
        #     sir_ypro2_t <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
        #     sir_ypro2_t[index_y2,] <- sir_ybeta2
        #     sir_ypro2_t/sqrt(colSums(sir_ypro2_t^2))
        #   }}, error =function(e){matrix(NA,p,1) })
        
        pro1 <- cbind(beta1_pro, sir_Upro1)                                  
        pro1_orth <- matrix(numeric(0), nrow = nrow(pro1), ncol = 0)
        tol <- 1e-12
        for (jj in seq_len(ncol(pro1))) {
          v1 <- pro1[, jj]
          if (ncol(pro1_orth) > 0) {
            for (ii in seq_len(ncol(pro1_orth))) { v1 <- v1 - sum(v1 * pro1_orth[, ii]) * pro1_orth[, ii] }
          }
          v1_norm <- sqrt(sum(v1^2))
          if (!is.finite(v1_norm) || v1_norm < tol) next
          pro1_orth <- cbind( pro1_orth, v1 / v1_norm)
        }
        sir_Upro1_orth <- pro1_orth[, -1, drop = FALSE]
        
        pro2 <- cbind(beta2_pro, sir_Upro2)
        pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
        for (jj in seq_len(ncol(pro2))) {
          v2 <- pro2[, jj]
          if (ncol(pro2_orth) > 0) {
            for (ii in seq_len(ncol(pro2_orth))) { v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii] }
          }
          v2_norm <- sqrt(sum(v2^2))
          if (!is.finite(v2_norm) || v2_norm < tol) next
          pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        }
        sir_Upro2_orth <- pro2_orth[, -1, drop = FALSE]
        
        # PCvM Computing
        pro2_num <- ncol(sir_Upro2_orth)                                                                                     
        pval_matrix1_PCvm <- matrix(nrow=1, ncol=pro2_num+1)                                                               
        deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
        for(q in seq_len(pro2_num)){
          x1_pro <- x1%*%sir_Upro2_orth[,q]
          band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
          h1 <- bandwidth_choice(x1_pro, band_y1)                                                                            
          Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                                  
          Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
          kernel_1 <- (3/4)*(1-Ker_inter1^2)*Ker_indictor1
          Indictor_1<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          A1_1 <- ((kernel_1%*%deri_link1)*x1_pro)/(kernel_1%*%U1^2+n1^(-12))                                                
          A1_2 <- kernel_1%*%((x1 %*% beta1_hat)*deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                     
          A1_3 <- (kernel_1%*%deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                                         
          hat_A1 <- rbind(t(A1_1), t(A1_2), t(A1_3))
          Gamma_inv1 <- array(0, dim=c(3, 3, n1))
          for (i in 1:n1) { Gamma_inv1[,,i] = ginv(((rep(1,3)%*%t(U1)^2)*hat_A1)%*%(t(hat_A1) * (as.matrix(Indictor_1[i,]) %*%rep(1,3)))) }
          integral_1 <- ((rep(1,3)%*%t(U1))*hat_A1)%*%t(Indictor_1)
          martingle_sec1<-diag(0,n1)  
          for (l in 1:n1){ martingle_sec1[l,]<- (U1[l]^2 * t(hat_A1[,l])) %*% Gamma_inv1[,,l] %*% integral_1[,l] %*% Indictor_1[l,] }
          martingle_sta1<-(1/sqrt(n1))*t(U1) %*% Indictor_1 - (1/sqrt(n1)) * colSums(martingle_sec1)   
          ordx1_pro <- sort(x1_pro)  
          t_1 <- ordx1_pro[floor(0.99 * n1)]
          psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
          PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        
          pval_matrix1_PCvm[,q] <- pvalue_integ_Brown(PCvM_martingle1)
        }
        x1_betapro <- x1%*%beta2_pro
        beta_Indictor1<-ifelse(x1_betapro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_betapro),1,0)
        hat_beta_A1 <- rbind(t(x1_betapro), matrix(1,1,n1))
        Gamma_beta_inv1 <- array(0, dim=c(2, 2, n1))
        for (i in 1:n1) { Gamma_beta_inv1[,,i] = ginv(((rep(1,2)%*%t(U1)^2)*hat_beta_A1)%*%(t(hat_beta_A1) * (as.matrix(beta_Indictor1[i,]) %*%rep(1,2)))) }
        integral_beta_1 <- ((rep(1,2)%*%t(U1))*hat_beta_A1)%*%t(beta_Indictor1)                      
        martingle_beta_sec1<-diag(0,n1)                                                                                      
        for (l in 1:n1){ martingle_beta_sec1[l,]<- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,] }
        martingle_beta_sta1<-(1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      
        ordx1_betapro <- sort(x1_betapro)  
        beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
        psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
        PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  
        pval_matrix1_PCvm[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1)
        
        pro1_num <- ncol(sir_Upro1_orth)                                                                                     
        pval_matrix2_PCvm <- matrix(nrow=1, ncol=pro1_num+1)                                                               
        deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
        for(qq in seq_len(pro1_num)){
          x2_pro <- x2%*%sir_Upro1_orth[,qq]
          band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
          h2 <- bandwidth_choice(x2_pro, band_y2)                                                                            
          Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                                                  
          Ker_indictor2 <- ifelse(abs(Ker_inter2) <= 1, 1, 0)
          kernel_2 <- (3/4)*(1-Ker_inter2^2)*Ker_indictor2
          Indictor_2<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
          A2_1 <- ((kernel_2%*%deri_link2)*x2_pro)/(kernel_2%*%U2^2+n2^(-12))                                                
          A2_2 <- kernel_2%*%((x2 %*% beta2_hat)*deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                     
          A2_3 <- (kernel_2%*%deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                                         
          hat_A2 <- rbind(t(A2_1), t(A2_2), t(A2_3))
          Gamma_inv2 <- array(0, dim=c(3, 3, n2))
          for (ii in 1:n2){ Gamma_inv2[,,ii] = ginv(((rep(1,3)%*%t(U2)^2)*hat_A2)%*%(t(hat_A2) * (as.matrix(Indictor_2[ii,]) %*%rep(1,3)))) }
          integral_2 <- ((rep(1,3)%*%t(U2))*hat_A2)%*%t(Indictor_2)
          martingle_sec2<-diag(0,n2)                                                                                       
          for (ll in 1:n2){ martingle_sec2[ll,]<- (U2[ll]^2 * t(hat_A2[,ll])) %*% Gamma_inv2[,,ll] %*% integral_2[,ll] %*% Indictor_2[ll,] }
          martingle_sta2<-(1/sqrt(n2))*t(U2) %*% Indictor_2 - (1/sqrt(n2)) * colSums(martingle_sec2)  
          ordx2_pro <- sort(x2_pro)                                                                                        
          t_2 <- ordx2_pro[floor(0.99 * n2)]
          psi_2 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2))
          PCvM_martingle2 <- (1/psi_2^2) * mean(U2^2 * (x2_pro <= t_2) * (t(martingle_sta2)^2))       
          pval_matrix2_PCvm[,qq] <- pvalue_integ_Brown(PCvM_martingle2)
        }
        x2_betapro <- x2%*%beta1_pro
        beta_Indictor2<-ifelse(x2_betapro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_betapro),1,0)
        hat_beta_A2 <- rbind(t(x2_betapro), matrix(1,1,n2))
        Gamma_beta_inv2 <- array(0, dim=c(2, 2, n2))
        for (ii in 1:n2) { Gamma_beta_inv2[,,ii] = ginv(((rep(1,2)%*%t(U2)^2)*hat_beta_A2)%*%(t(hat_beta_A2) * (as.matrix(beta_Indictor2[ii,]) %*%rep(1,2)))) }
        integral_beta_2 <- ((rep(1,2)%*%t(U2))*hat_beta_A2)%*%t(beta_Indictor2)
        martingle_beta_sec2<-diag(0,n2)                                                                                    
        for (ll in 1:n2){ martingle_beta_sec2[ll,]<- (U2[ll]^2 * t(hat_beta_A2[,ll])) %*% Gamma_beta_inv2[,,ll] %*% integral_beta_2[,ll] %*% beta_Indictor2[ll,] }
        martingle_beta_sta2<-(1/sqrt(n2))*t(U2) %*% beta_Indictor2 - (1/sqrt(n2)) * colSums(martingle_beta_sec2)   
        ordx2_betapro <- sort(x2_betapro)                                                                                  
        beta_t2 <- ordx2_betapro[floor(0.99 * n2)]
        psi_beta_2 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2))
        PCvM_beta_martingle2 <- (1/psi_beta_2^2) * mean(U2^2 * (x2_betapro <= beta_t2) * (t(martingle_beta_sta2)^2))   
        pval_matrix2_PCvm[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2)
        
        pval_matrix_PCvm <- cbind(pval_matrix1_PCvm, pval_matrix2_PCvm)
        
        # PLS computations for Hybrid
        pro2_num_pls <- ncol(pro2_orth)                                                                                    
        h1 <- n1^(-2/9)                                                                                                    
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) 
        errormat1 <- U1%*%t(U1)                                                                                            
        for(q in seq_len(pro2_num_pls)){
          x_pro1 <- x1%*%pro2_orth[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
        }
        
        pro1_num_pls <- ncol(pro1_orth)                                                                                    
        h2 <- n2^(-2/9)                                                                                                    
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) 
        errormat2 <- U2%*%t(U2)  
        for(l in seq_len(pro1_num_pls)){
          x_pro2 <- x2%*%pro1_orth[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
        }
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
        
        # Hybrid PCvM + PLS
        pval_matrix_hybrid <- cbind(pval_matrix_PCvm, pval_matrix_PLS)
        pval_num_h <- ncol(pval_matrix_hybrid)
        # pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
        1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid)*pi), na.rm=TRUE))
        
      }
    
    stopImplicitCluster()
    stopCluster(cl)
    
    # return(list(pval_fisher1_power = mean(hybrid_power[,1] <= 0.05,na.rm=TRUE),
    #             pval_fisher2_power = mean(hybrid_power[,2] <= 0.05,na.rm=TRUE),
    #             pval_fisher_cauchy_power = mean(hybrid_power[,3] <= 0.05,na.rm=TRUE),
    #             hybrid_cauchy_power = mean(hybrid_power[,4] <= 0.05,na.rm=TRUE)))
    return(list(hybrid_cauchy_power = mean(hybrid_power <= 0.05, na.rm=TRUE)))
  }
  
  PLS_logit_parallel_run <- function(p, n, a, pho){
    library(MASS)
    library(glmnet)
    library(LassoSIR)
    library(VariableScreening) 
    library(psych)
    library(foreach)
    library(parallel)
    library(iterators)
    library(doParallel)
    
    s <- 1000 
    c_h <- 1
    mu <- rep(0, p)
    sigma <- if(pho == 0) diag(rep(1,p)) else toeplitz(pho^(0:(p-1)))
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    PLS_power <- foreach(k = 1:s, .combine='rbind', .errorhandling = "remove",
                         .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych')) %dopar%  
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                               
        pr <- 1/(1 + exp(-z)) 
        y <- matrix(rbinom(n, 1, pr),ncol = 1)
        
        index_x<-sort(sample(1:n,floor(n/2),replace = FALSE))  
        x1 <- x[-index_x,]; y1 <- as.matrix(y[-index_x,])
        x2 <- x[index_x,];  y2 <- as.matrix(y[index_x,]) 
        n1 <- nrow(x1); n2 <- nrow(x2)
        
        # x1 estimation
        lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
        beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
        index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
        if(length(index_beta1_non0) == 0){
          pred1 <- mean(y1)
        }else{
          x1_sec <- x1[,index_beta1_non0, drop=FALSE]                                          
          sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
          pred1 = matrix(unname(predict(sec_model1, newx = x1_sec, type="response")), ncol = 1) 
        }
        U1 <- y1 - pred1                                                             
        
        screen_num1 <- floor(n1/log(n1))
        sdr_Upro1 <- tryCatch({
          if(p <= screen_num1){
            sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_U1$beta/sqrt(colSums(sir_U1$beta^2))
          }else{
            rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
            index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
            sir_Ulasso1 <- LassoSIR(x1[,index_U1,drop=FALSE], U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ulasso1$beta))
            sir_Upro1[index_U1,] <- sir_Ulasso1$beta
            sir_Upro1/sqrt(colSums(sir_Upro1^2))
          }
        }, error =function(e){ NA })
        
        sdr_ypro1 <- tryCatch({
          if(p <= screen_num1){
            sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
            sir_ylasso1$beta/sqrt(colSums(sir_ylasso1$beta^2))
          }else{
            rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")
            index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]
            sir_ylasso1 <- LassoSIR(x1[,index_y1,drop=FALSE], y1, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
            sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ylasso1$beta))
            sir_ypro1[index_y1,] <- sir_ylasso1$beta
            sir_ypro1/sqrt(colSums(sir_ypro1^2))
          }
        }, error =function(e){ NA })

        # x2 estimation
        lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
        beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
        index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso) != 0]
        if(length(index_beta2_non0) == 0){
          pred2 <- mean(y2)
        }else{
          x2_sec <- x2[,index_beta2_non0, drop=FALSE]                                          
          sec_model2 <- glm(y2~x2_sec, family = binomial(link = "logit"))
          pred2 = matrix(unname(predict(sec_model2, newx = x2_sec, type="response")), ncol = 1) 
        }
        U2 <- y2 - pred2                                                             
        
        screen_num2 <- floor(n2/log(n2))
        sdr_Upro2 <- tryCatch({
          if(p <= screen_num2){
            sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Ulasso2$beta/sqrt(colSums(sir_Ulasso2$beta^2))
          }else{
            rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
            index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
            sir_Ulasso2 <- LassoSIR(x2[,index_U2, drop=FALSE], U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
            sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ulasso2$beta))
            sir_Upro2[index_U2,] <- sir_Ulasso2$beta
            sir_Upro2/sqrt(colSums(sir_Upro2^2))
          }
        }, error =function(e){ NA })
        
        sdr_ypro2 <- tryCatch({
          if(p <= screen_num2){
            sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
            sir_ylasso2$beta/sqrt(colSums(sir_ylasso2$beta^2))
          }else{
            rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")
            index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]
            sir_ylasso2 <- LassoSIR(x2[,index_y2, drop=FALSE], y2, H=2, choosing.d="automatic", solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
            sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ylasso2$beta))
            sir_ypro2[index_y2,] <- sir_ylasso2$beta
            sir_ypro2/sqrt(colSums(sir_ypro2^2))
          }
        }, error =function(e){ NA })

        # PLS 1
        pro2 <- cbind(sdr_Upro2, sdr_ypro2)                                          
        pro2 <- tryCatch({
          qr_obj2 <- qr(pro2)
          if( qr_obj2$rank == ncol(pro2) ) pro2 else pro2[, qr_obj2$pivot[1:qr_obj2$rank], drop = FALSE]
        }, error =function(e) NA)
        pro2 <- as.matrix(pro2)
        pro2_num <- ncol(pro2)                                                      
        h1 <- c_h*(n1^(-2/9))                                                       
        
        pval_matrix1 <- matrix(nrow=1, ncol=pro2_num) 
        errormat1 <- U1%*%t(U1)  
        for(q in 1:pro2_num){
          if(any(is.na(pro2[,q])) ){ pval_matrix1[,q] <- NA }else{
            x_pro1 <- x1%*%pro2[,q]
            x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   
            indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
            kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                
            Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
            pval_matrix1[,q] <- 1-pnorm(Tn1)
          }
        }
        
        # PLS 2
        pro1 <- cbind(sdr_Upro1, sdr_ypro1)                                          
        pro1 <- tryCatch({
          qr_obj1 <- qr(pro1)
          if( qr_obj1$rank == ncol(pro1) ) pro1 else pro1[, qr_obj1$pivot[1:qr_obj1$rank], drop = FALSE]
        }, error =function(e) NA)
        pro1 <- as.matrix(pro1)
        pro1_num <- ncol(pro1)                                                      
        h2 <- c_h*(n2^(-2/9))                                                       
        
        pval_matrix2 <- matrix(nrow=1, ncol=pro1_num) 
        errormat2 <- U2%*%t(U2)  
        for(l in 1:pro1_num){
          if(any(is.na(pro1[,l])) ){ pval_matrix2[,l] <- NA }else{
            x_pro2 <- x2%*%pro1[,l]
            x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2   
            indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
            kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                
            Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
            pval_matrix2[,l] <- 1-pnorm(Tn2)
          }
        }
        
        pval_fisher1 <- 1- pchisq(-2*sum(log(pval_matrix1[1:pro2_num])),df=2*pro2_num)
        pval_fisher2 <- 1- pchisq(-2*sum(log(pval_matrix2[1:pro1_num])),df=2*pro1_num)
        pval_fisher_m <- cbind(pval_fisher1, pval_fisher2)
        pval_fisher <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
        
        pval_min1 <- min(pval_matrix1)<=1-(0.95)^(1/pro2_num) 
        pval_min2 <- min(pval_matrix2)<=1-(0.95)^(1/pro1_num) 
        pval_min_m <- cbind(1 - (1 - min(pval_matrix1))^pro2_num, 1 - (1 - min(pval_matrix2))^pro1_num) 
        pval_min <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
        
        c(pval_fisher, pval_min)
      }
    stopCluster(cl)
    
    PLS_fisher_power <- mean(PLS_power[,1]<=0.05,na.rm=TRUE)
    PLS_min_power <- mean(PLS_power[,2]<=0.05,na.rm=TRUE)
    
    return(list(PLS_fisher_power=PLS_fisher_power, PLS_min_power=PLS_min_power))
  }
 
  grp_logit_parallel_run <- function(p, n, a, pho){
    library(MASS); library(matrixcalc); library(GRPtests)
    library(foreach); library(parallel); library(iterators); library(doParallel)
    
    s <- 1000
    mu <- rep(0, p)
    sigma <- if(pho == 0) diag(rep(1,p)) else toeplitz(pho^(0:(p-1)))
    beta0 <- c(rep(1,5),rep(0,p-5))
    
    cores <- detectCores(logical=F)-2
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores=cores)
    grp_power <- foreach(k = 1:s, .combine='rbind', .errorhandling = "remove",
                         .packages = c('MASS','matrixcalc','GRPtests')) %dopar%
      {
        set.seed(k)
        x <- mvrnorm(n, mu, sigma)
        z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                              
        pr <- 1/(1 + exp(-z))
        y <- matrix(rbinom(n, 1, pr),ncol = 1)
        GRPtest(x, y, fam = "binomial", nsplits = 1, RP_function = NULL, output_all=FALSE, penalize=TRUE)
      }
    stopCluster(cl)
    
    return(list(Power = mean(grp_power<=0.05,na.rm=TRUE)))
  }
 
}


# ====== Parameter Settings ======
n_1 <- 300
n_2 <- 600
p_vec <- c(50, 100, 300, 600, 900, 1200)   
a_vec <- c(0, 1)                            
pho_vec <- c(0)                          
benchmark_results <- data.frame()

closeAllConnections()  

for (kk in seq_along(pho_vec)) {
  for (ii in seq_along(a_vec)) {
    for (jj in seq_along(p_vec)) {
      
      curr_p <- p_vec[jj]
      curr_a <- a_vec[ii]
      curr_pho <- pho_vec[kk] 
      curr_seed <- 12345 + jj
      
      cat(sprintf("\n>>> Parameter combination: p = %d, a = %.1f, pho = %.1f <<<\n", curr_p, curr_a, curr_pho))
      
      # ===================== Linear test =====================
      
      cat("[Linear] Running PCvM1 (Only) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_pcvm1 <- tryCatch(PCvM1_linear_only(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_lin_pcvm1 <- PCvM1_linear_only(curr_p, n_1, curr_a, curr_pho)
      t_lin_pcvm1 <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_pcvm1 <- if(is.list(res_lin_pcvm1)) res_lin_pcvm1$martingale_cauchy1_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_pcvm1, pow_lin_pcvm1))
      
      cat("[Linear] Running PCvM (Full) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_pcvm <- tryCatch(Pcvm_linear_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_lin_pcvm <- Pcvm_linear_run(curr_p, n_1, curr_a, curr_pho)
      t_lin_pcvm <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_pcvm <- if(is.list(res_lin_pcvm[1])) res_lin_pcvm$martingale_cauchy_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_pcvm, pow_lin_pcvm))

      # # --- Hybrid ---
      cat("[Linear] Running Hybrid ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_hyb <- tryCatch(Hybrid_linear_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_lin_hyb <- Hybrid_linear_run(curr_p, n_1, curr_a, curr_pho)
      t_lin_hyb <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_hyb <- if(is.list(res_lin_hyb[1])) res_lin_hyb$hybrid_cauchy_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_hyb, pow_lin_hyb))

      cat("[Linear] Running CPB(PLS) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_pls <- tryCatch(PLS_lm_parallel_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_lin_pls <- PLS_lm_parallel_run(curr_p, n_1, curr_a, curr_pho)
      t_lin_pls <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_pls <- if(is.list(res_lin_pls[1])) res_lin_pls$PLS_fisher_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_pls, pow_lin_pls))

      cat("[Linear] Running RPtest ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_rp <- tryCatch(grp_rp_lm_parallel_run(curr_p, n, curr_a, curr_pho, target_method="RP"), error=function(e) NA)
      res_lin_rp <- grp_rp_lm_parallel_run(curr_p, n_1, curr_a, curr_pho, target_method="RP")
      t_lin_rp <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_rp <- if(is.list(res_lin_rp[1])) res_lin_rp$Power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_rp, pow_lin_rp))

      cat("[Linear] Running GRPtest ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_lin_grp <- tryCatch(grp_rp_lm_parallel_run(curr_p, n, curr_a, curr_pho, target_method="GRP"), error=function(e) NA)
      res_lin_grp <- grp_rp_lm_parallel_run(curr_p, n_1, curr_a, curr_pho, target_method="GRP")
      t_lin_grp <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_lin_grp <- if(is.list(res_lin_grp[1])) res_lin_grp$Power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_lin_grp, pow_lin_grp))
      
      # # ===================== Logit test =====================

      cat("[Logit] Running PCvM1 (Only) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_log_pcvm1 <- tryCatch(PCvM1_logit_only(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_log_pcvm1 <- PCvM1_logit_only(curr_p, n_2, curr_a, curr_pho)
      t_log_pcvm1 <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_log_pcvm1 <- if(is.list(res_log_pcvm1[1])) res_log_pcvm1$martingale_cauchy1_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_log_pcvm1, pow_log_pcvm1))

      cat("[Logit] Running PCvM (Full) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_log_pcvm <- tryCatch(Pcvm_logit_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_log_pcvm <- Pcvm_logit_run(curr_p, n_2, curr_a, curr_pho)
      t_log_pcvm <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_log_pcvm <- if(is.list(res_log_pcvm[1])) res_log_pcvm$martingale_cauchy_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_log_pcvm, pow_log_pcvm))

      # --- Hybrid ---
      cat("[Logit] Running Hybrid ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_log_hyb <- tryCatch(Hybrid_logit_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_log_hyb <- Hybrid_logit_run(curr_p, n_2, curr_a, curr_pho)
      t_log_hyb <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_log_hyb <- if(is.list(res_log_hyb[1])) res_log_hyb$hybrid_cauchy_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_log_hyb, pow_log_hyb))

      cat("[Logit] Running CPB(PLS) ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_log_pls <- tryCatch(PLS_logit_parallel_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_log_pls <- PLS_logit_parallel_run(curr_p, n_2, curr_a, curr_pho)
      t_log_pls <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_log_pls <- if(is.list(res_log_pls[1])) res_log_pls$PLS_fisher_power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_log_pls, pow_log_pls))

      cat("[Logit] Running GRPtest ... ")
      t0 <- Sys.time()
      set.seed(curr_seed)
      # res_log_grp <- tryCatch(grp_logit_parallel_run(curr_p, n, curr_a, curr_pho), error=function(e) NA)
      res_log_grp <- grp_logit_parallel_run(curr_p, n_2, curr_a, curr_pho)
      t_log_grp <- as.numeric(difftime(Sys.time(), t0, units="secs"))
      pow_log_grp <- if(is.list(res_log_grp[1])) res_log_grp$Power else NA
      cat(sprintf("Time: %.2f s | Power: %.3f\n", t_log_grp, pow_log_grp))
      # 
      # ===================== Save results =====================
      row_res <- data.frame(
        p = curr_p, a = curr_a, pho = curr_pho,
        # Linear (Time & Power)
        Lin_Time_PCvM1 = t_lin_pcvm1, Lin_Pow_PCvM1 = pow_lin_pcvm1,  
        Lin_Time_PCvM  = t_lin_pcvm,  Lin_Pow_PCvM  = pow_lin_pcvm,
        Lin_Time_Hyb   = t_lin_hyb,   Lin_Pow_Hyb   = pow_lin_hyb ,
        Lin_Time_PLS   = t_lin_pls,   Lin_Pow_PLS   = pow_lin_pls,
        Lin_Time_RP    = t_lin_rp,    Lin_Pow_RP    = pow_lin_rp,
        Lin_Time_GRP   = t_lin_grp,   Lin_Pow_GRP   = pow_lin_grp ,
        # Logit (Time & Power)
        Log_Time_PCvM1 = t_log_pcvm1, Log_Pow_PCvM1 = pow_log_pcvm1, 
        Log_Time_PCvM  = t_log_pcvm,  Log_Pow_PCvM  = pow_log_pcvm,
        Log_Time_Hyb   = t_log_hyb,   Log_Pow_Hyb   = pow_log_hyb,
        Log_Time_PLS   = t_log_pls,   Log_Pow_PLS   = pow_log_pls,
        Log_Time_GRP   = t_log_grp,   Log_Pow_GRP   = pow_log_grp
      )
      
      benchmark_results <- rbind(benchmark_results, row_res)
      
      write.csv(benchmark_results, "Results_time.csv", row.names = FALSE)
    }
  }
}
 


