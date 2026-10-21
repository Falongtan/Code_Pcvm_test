#The hybrid test, the martingale-based test with orthogonal projections, and local-smoothing test for testing high dimensional linear and logistical regression models
#projections with schmidt orthogonalization 
Pcvm_Pls_hybrid_logit_h <- function(p,n,a,pho){
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
  
  cores <- detectCores(logical = FALSE) - 2
  random_port <- sample(10000:60000, 1)
  cl <- makeCluster(cores, port = random_port)
  registerDoParallel(cl, cores=cores)
  
  hybrid_power <- foreach(k = 1:s, .combine='rbind',
                          .packages = c('MASS','harmonicmeanp', 'glmnet','LassoSIR','VariableScreening','psych'),
                          .export = c("pvalue_integ_Brown",'bandwidth_choice')) %dopar%
    {
      x <- mvrnorm(n, mu, sigma)
      z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                                # H21
      
      pr <- 1/(1 + exp(-z))
      y <- matrix(rbinom(n, 1, pr),ncol = 1)
      
      #data splitting
      index_x<-sample(1:n,floor(n/2),replace = FALSE)
      index_x<-sort(index_x)
      
      x1 <- x[-index_x,]
      x2 <- x[index_x,]
      
      y1 <- as.matrix(y[-index_x,])
      y2 <- as.matrix(y[index_x,])
      
      n1 <- nrow(x1)
      n2 <- nrow(x2)
      
      # residuals and beta_0 projections based on lasso and data x1 y1
      lasso_model1 <- cv.glmnet(x1, y1, family="binomial", intercept = T)
      beta1_lasso <- coef(lasso_model1, s = "lambda.min")[-1]
      index_beta1_non0 <- seq(1:p)[as.numeric(beta1_lasso) != 0]
      len_beta1_non0 <- length(index_beta1_non0)
      if(len_beta1_non0 == 0){
        pred1 <- mean(y1)
        U1 <- y1 - pred1                                                         # residuals based on x1 y1
        beta1_hat <- beta1_lasso
        intercept1 <- coef(lasso_model1, s = "lambda.min")[1]                    # intercept term
        beta1_pro <- rep(1,p)/sqrt(p)
      }else{
        x1_sec <- x1[,index_beta1_non0]                                          # second estimation
        sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
        beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
        beta1_hat <- beta1_lasso
        intercept1 <- unname(sec_model1$coefficients)[1]
        pred1 = predict(sec_model1, type="response")
        pred1 = matrix(unname(pred1), ncol = 1)
        U1 <- y1 - pred1                                                         # residuals based on x1 y1
        beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
        beta1_pro <- as.matrix(beta1_pro)
      }
      
      # projections based on screening and lassosir for x1 U1
      screen_num1 <- floor(n1/log(n1))
      sir_Upro1 <- tryCatch({
        if(p <= screen_num1){
          sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic",
                             solution.path=FALSE, categorical=FALSE, nfolds=5,
                             screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
        }else{
          rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
          index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
          x1_Uscreen <- x1[,index_U1]
          sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta1<- sir_Ulasso1$beta
          sir_Upro1_temp <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
          sir_Upro1_temp[index_U1,] <- sir_Ubeta1
          sir_Upro1_temp/sqrt(colSums(sir_Upro1_temp^2))
        }}, error =function(e){ matrix(NA,p,1) })
      
      # projections based on screening and lassosir for x1 y1
      sir_ypro1 <- tryCatch({
        if(p <= screen_num1){
          sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta1<- sir_ylasso1$beta
          sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
        }else{
          rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                      
          index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        
          x1_yscreen <- x1[,index_y1]
          sir_ylasso1 <- LassoSIR(x1_yscreen, y1, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta1<- sir_ylasso1$beta
          sir_ypro1_temp <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
          sir_ypro1_temp[index_y1,] <- sir_ybeta1
          sir_ypro1_temp/sqrt(colSums(sir_ypro1_temp^2))
        }
      }, error =function(e){ matrix(NA,p,1) })
      
      # residuals and beta0-projection based on lasso and data x2 y2
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
      
      #sir projections based on screening and x2 U2
      screen_num2 <- floor(n2/log(n2))
      sir_Upro2 <- tryCatch({
        if(p <= screen_num2){
          sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta2<- sir_Ulasso2$beta
          sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
        }else{
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
          x2_Uscreen <- x2[,index_U2]
          sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta2 <- sir_Ulasso2$beta
          sir_Upro2_temp <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
          sir_Upro2_temp[index_U2,] <- sir_Ubeta2
          sir_Upro2_temp/sqrt(colSums(sir_Upro2_temp^2))
        }}, error =function(e){matrix(NA,p,1) })
      
      #sir projections based on screening and x2 y2
      sir_ypro2 <- tryCatch({
        if(p <= screen_num2){
          sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta2<- sir_ylasso2$beta
          sir_ybeta2/sqrt(colSums(sir_ybeta2^2))
        }else{
          rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                      
          index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        
          x2_yscreen <- x2[,index_y2]
          sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta2 <- sir_ylasso2$beta
          sir_ypro2_temp <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
          sir_ypro2_temp[index_y2,] <- sir_ybeta2
          sir_ypro2_temp/sqrt(colSums(sir_ypro2_temp^2))
        }}, error =function(e){matrix(NA,p,1) })
      
      # Schmidt orthogonalization
      pro1 <- cbind(beta1_pro, sir_Upro1)                                  
      pro1_orth <- matrix(numeric(0), nrow = nrow(pro1), ncol = 0)
      tol <- 1e-12
      deleted1 <- integer(0)
      kept1 <- integer(0)
      for (jj in seq_len(ncol(pro1))) {
        v1 <- pro1[, jj]
        if (ncol(pro1_orth) > 0) {
          for (ii in seq_len(ncol(pro1_orth))) {
            v1 <- v1 - sum(v1 * pro1_orth[, ii]) * pro1_orth[, ii]
          }
        }
        v1_norm <- sqrt(sum(v1^2))
        if (!is.finite(v1_norm) || v1_norm < tol) {
          deleted1 <- c(deleted1, jj)
          next
        }
        pro1_orth <- cbind( pro1_orth, v1 / v1_norm)
        kept1 <- c(kept1, jj)
      }
      sir_Upro1_orth <- pro1_orth[, -1, drop = FALSE]
      
      pro2 <- cbind(beta2_pro, sir_Upro2)
      pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
      deleted2 <- integer(0)
      kept2 <- integer(0)
      for (jj in seq_len(ncol(pro2))) {
        v2 <- pro2[, jj]
        if (ncol(pro2_orth) > 0) {
          for (ii in seq_len(ncol(pro2_orth))) {
            v2 <- v2 - sum(v2 * pro2_orth[, ii]) * pro2_orth[, ii]
          }
        }
        v2_norm <- sqrt(sum(v2^2))
        if (!is.finite(v2_norm) || v2_norm < tol) {
          deleted2 <- c(deleted2, jj)
          next
        }
        pro2_orth <- cbind( pro2_orth, v2 / v2_norm)
        kept2 <- c(kept2, jj)
      }
      sir_Upro2_orth <- pro2_orth[, -1, drop = FALSE]
      
      # =========================================================================
      # 1. PCvM Test (Bandwidth Sensitivity Analysis: 0.5, 1.0, 2.0)
      # =========================================================================
      pro2_num <- ncol(sir_Upro2_orth)                                                                                
      pval_matrix1_PCvm_05 <- matrix(nrow=1, ncol=pro2_num+1)                                                            
      pval_matrix1_PCvm_10 <- matrix(nrow=1, ncol=pro2_num+1)                                                            
      pval_matrix1_PCvm_20 <- matrix(nrow=1, ncol=pro2_num+1)  
      
      deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
      
      for(q in seq_len(pro2_num)){
        x1_pro <- x1%*%sir_Upro2_orth[,q]
        band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
        h1_base <- bandwidth_choice(x1_pro, band_y1)                                                            
        
        for (m_idx in 1:3) {
          mult <- c(0.5, 1.0, 2.0)[m_idx]
          h1 <- mult * h1_base
          
          Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                  
          Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
          kernel_1 <- (3/4)*(1-Ker_inter1^2)*Ker_indictor1
          
          Indictor_1<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
          
          A1_1 <- ((kernel_1%*%deri_link1)*x1_pro)/(kernel_1%*%U1^2+n1^(-12))                                
          A1_2 <- kernel_1%*%((x1 %*% beta1_hat)*deri_link1)/(kernel_1%*%U1^2+n1^(-12))                      
          A1_3 <- (kernel_1%*%deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                         
          hat_A1 <- rbind(t(A1_1), t(A1_2), t(A1_3))
          
          Gamma_inv1 <- array(0, dim=c(3, 3, n1))
          for (i in 1:n1) {
            Gamma_inv1[,,i] = ginv(((rep(1,3)%*%t(U1)^2)*hat_A1)%*%(t(hat_A1) * (as.matrix(Indictor_1[i,]) %*%rep(1,3))))
          }
          integral_1 <- ((rep(1,3)%*%t(U1))*hat_A1)%*%t(Indictor_1)
          
          martingle_sec1<-diag(0,n1)  
          for (l in 1:n1){
            martingle_sec1[l,]<- (U1[l]^2 * t(hat_A1[,l])) %*% Gamma_inv1[,,l] %*% integral_1[,l] %*% Indictor_1[l,]
          }
          martingle_sta1<-(1/sqrt(n1))*t(U1) %*% Indictor_1 - (1/sqrt(n1)) * colSums(martingle_sec1)   
          
          ordx1_pro <- sort(x1_pro)  
          t_1 <- ordx1_pro[floor(0.99 * n1)]
          psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
          PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        
          pval <- pvalue_integ_Brown(PCvM_martingle1)
          
          if(m_idx==1) pval_matrix1_PCvm_05[,q] <- pval
          if(m_idx==2) pval_matrix1_PCvm_10[,q] <- pval
          if(m_idx==3) pval_matrix1_PCvm_20[,q] <- pval
        }
      }
      
      # beta projection PCvM: Calculated once (no bandwidth involved)
      x1_betapro <- x1%*%beta2_pro
      beta_Indictor1<-ifelse(x1_betapro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_betapro),1,0)
      
      hat_beta_A1 <- rbind(t(x1_betapro), matrix(1,1,n1))
      Gamma_beta_inv1 <- array(0, dim=c(2, 2, n1))
      for (i in 1:n1) {
        Gamma_beta_inv1[,,i] = ginv(((rep(1,2)%*%t(U1)^2)*hat_beta_A1)%*%(t(hat_beta_A1) * (as.matrix(beta_Indictor1[i,]) %*%rep(1,2))) )
      }
      integral_beta_1 <- ((rep(1,2)%*%t(U1))*hat_beta_A1)%*%t(beta_Indictor1)                      
      
      martingle_beta_sec1<-diag(0,n1)                                                              
      for (l in 1:n1){
        martingle_beta_sec1[l,]<- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,]
      }
      martingle_beta_sta1<-(1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      
      
      ordx1_betapro <- sort(x1_betapro)  
      beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
      psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
      PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  
      pval_beta1 <- pvalue_integ_Brown(PCvM_beta_martingle1)
      
      pval_matrix1_PCvm_05[pro2_num+1] <- pval_beta1
      pval_matrix1_PCvm_10[pro2_num+1] <- pval_beta1
      pval_matrix1_PCvm_20[pro2_num+1] <- pval_beta1
      
      
      # For Data part 2
      pro1_num <- ncol(sir_Upro1_orth)                                                                             
      pval_matrix2_PCvm_05 <- matrix(nrow=1, ncol=pro1_num+1)                                                         
      pval_matrix2_PCvm_10 <- matrix(nrow=1, ncol=pro1_num+1)                                                         
      pval_matrix2_PCvm_20 <- matrix(nrow=1, ncol=pro1_num+1)  
      
      deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
      
      for(qq in seq_len(pro1_num)){
        x2_pro <- x2%*%sir_Upro1_orth[,qq]
        band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
        h2_base <- bandwidth_choice(x2_pro, band_y2)                                                            
        
        for (m_idx in 1:3) {
          mult <- c(0.5, 1.0, 2.0)[m_idx]
          h2 <- mult * h2_base
          
          Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                                  
          Ker_indictor2 <- ifelse(abs(Ker_inter2) <= 1, 1, 0)
          kernel_2 <- (3/4)*(1-Ker_inter2^2)*Ker_indictor2
          
          Indictor_2<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
          
          A2_1 <- ((kernel_2%*%deri_link2)*x2_pro)/(kernel_2%*%U2^2+n2^(-12))                                
          A2_2 <- kernel_2%*%((x2 %*% beta2_hat)*deri_link2)/(kernel_2%*%U2^2+n2^(-12))                      
          A2_3 <- (kernel_2%*%deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                         
          hat_A2 <- rbind(t(A2_1), t(A2_2), t(A2_3))
          
          Gamma_inv2 <- array(0, dim=c(3, 3, n2))
          for (ii in 1:n2){
            Gamma_inv2[,,ii] = ginv(((rep(1,3)%*%t(U2)^2)*hat_A2)%*%(t(hat_A2) * (as.matrix(Indictor_2[ii,]) %*%rep(1,3))))
          }
          integral_2 <- ((rep(1,3)%*%t(U2))*hat_A2)%*%t(Indictor_2)
          
          martingle_sec2<-diag(0,n2)                                                                         
          for (ll in 1:n2){
            martingle_sec2[ll,]<- (U2[ll]^2 * t(hat_A2[,ll])) %*% Gamma_inv2[,,ll] %*% integral_2[,ll] %*% Indictor_2[ll,]
          }
          martingle_sta2<-(1/sqrt(n2))*t(U2) %*% Indictor_2 - (1/sqrt(n2)) * colSums(martingle_sec2)  
          
          ordx2_pro <- sort(x2_pro)                                                                          
          t_2 <- ordx2_pro[floor(0.99 * n2)]
          psi_2 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2))
          PCvM_martingle2 <- (1/psi_2^2) * mean(U2^2 * (x2_pro <= t_2) * (t(martingle_sta2)^2))       
          pval <- pvalue_integ_Brown(PCvM_martingle2)
          
          if(m_idx==1) pval_matrix2_PCvm_05[,qq] <- pval
          if(m_idx==2) pval_matrix2_PCvm_10[,qq] <- pval
          if(m_idx==3) pval_matrix2_PCvm_20[,qq] <- pval
        }
      }
      
      # beta projection PCvM
      x2_betapro <- x2%*%beta1_pro
      beta_Indictor2<-ifelse(x2_betapro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_betapro),1,0)
      
      hat_beta_A2 <- rbind(t(x2_betapro), matrix(1,1,n2))
      Gamma_beta_inv2 <- array(0, dim=c(2, 2, n2))
      for (ii in 1:n2) {
        Gamma_beta_inv2[,,ii] = ginv(((rep(1,2)%*%t(U2)^2)*hat_beta_A2)%*%(t(hat_beta_A2) * (as.matrix(beta_Indictor2[ii,]) %*%rep(1,2)))  )
      }
      integral_beta_2 <- ((rep(1,2)%*%t(U2))*hat_beta_A2)%*%t(beta_Indictor2)
      
      martingle_beta_sec2<-diag(0,n2)                                                                    
      for (ll in 1:n2){
        martingle_beta_sec2[ll,]<- (U2[ll]^2 * t(hat_beta_A2[,ll])) %*% Gamma_beta_inv2[,,ll] %*% integral_beta_2[,ll] %*% beta_Indictor2[ll,]
      }
      martingle_beta_sta2<-(1/sqrt(n2))*t(U2) %*% beta_Indictor2 - (1/sqrt(n2)) * colSums(martingle_beta_sec2)   
      
      ordx2_betapro <- sort(x2_betapro)                                                                  
      beta_t2 <- ordx2_betapro[floor(0.99 * n2)]
      psi_beta_2 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2))
      PCvM_beta_martingle2 <- (1/psi_beta_2^2) * mean(U2^2 * (x2_betapro <= beta_t2) * (t(martingle_beta_sta2)^2))   
      pval_beta2 <- pvalue_integ_Brown(PCvM_beta_martingle2)
      
      pval_matrix2_PCvm_05[pro1_num+1] <- pval_beta2
      pval_matrix2_PCvm_10[pro1_num+1] <- pval_beta2
      pval_matrix2_PCvm_20[pro1_num+1] <- pval_beta2
      
      
      # PCvM Aggregation
      pval_matrix_PCvm_05 <- cbind(pval_matrix1_PCvm_05, pval_matrix2_PCvm_05)
      pval_matrix_PCvm_10 <- cbind(pval_matrix1_PCvm_10, pval_matrix2_PCvm_10)
      pval_matrix_PCvm_20 <- cbind(pval_matrix1_PCvm_20, pval_matrix2_PCvm_20)
      
      pval_num <- ncol(pval_matrix_PCvm_05)
      
      pval_cauchy1_PCvm_05 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_05[1:(pro2_num+1)])*pi)))          
      pval_cauchy2_PCvm_05 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_05[1:(pro1_num+1)])*pi)))          
      pval_cauchy_PCvm_05  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_05[1:pval_num])*pi)))               
      
      pval_cauchy1_PCvm_10 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_10[1:(pro2_num+1)])*pi)))          
      pval_cauchy2_PCvm_10 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_10[1:(pro1_num+1)])*pi)))          
      pval_cauchy_PCvm_10  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_10[1:pval_num])*pi)))               
      
      pval_cauchy1_PCvm_20 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_20[1:(pro2_num+1)])*pi)))          
      pval_cauchy2_PCvm_20 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_20[1:(pro1_num+1)])*pi)))          
      pval_cauchy_PCvm_20  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_20[1:pval_num])*pi)))               
      
      # =========================================================================
      # 2. Hybrid & PLS Test (Fixed default bandwidth)
      # =========================================================================
      if(T){
        pro2_num_pls <- ncol(pro2_orth)                                                              
        h1 <- n1^(-2/9)                                                                              
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) 
        errormat1 <- U1%*%t(U1)                                                                      
        for(q in 1:pro2_num_pls){
          x_pro1 <- x1%*%pro2_orth[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval1 <- 1-pnorm(Tn1)
          pval_matrix1_PLS[,q] <- pval1
        }
        
        pro1_num_pls <- ncol(pro1_orth)                                                              
        h2 <- n2^(-2/9)                                                                              
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) 
        errormat2 <- U2%*%t(U2)  
        for(l in 1:pro1_num_pls){
          x_pro2 <- x2%*%pro1_orth[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval2 <- 1-pnorm(Tn2)
          pval_matrix2_PLS[,l] <- pval2
        }
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
      }
      
      # Hybrid Combinations
      pval_matrix_hybrid_05 <- cbind(pval_matrix_PCvm_05, pval_matrix_PLS)
      pval_cauchy_hybrid_05 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_05[1:ncol(pval_matrix_hybrid_05)])*pi)))
      
      pval_matrix_hybrid_10 <- cbind(pval_matrix_PCvm_10, pval_matrix_PLS)
      pval_cauchy_hybrid_10 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_10[1:ncol(pval_matrix_hybrid_10)])*pi)))
      
      pval_matrix_hybrid_20 <- cbind(pval_matrix_PCvm_20, pval_matrix_PLS)
      pval_cauchy_hybrid_20 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_20[1:ncol(pval_matrix_hybrid_20)])*pi)))
      
      c(
        pval_cauchy1_PCvm_05, pval_cauchy2_PCvm_05, pval_cauchy_PCvm_05, pval_cauchy_hybrid_05, 
        pval_cauchy1_PCvm_10, pval_cauchy2_PCvm_10, pval_cauchy_PCvm_10, pval_cauchy_hybrid_10,  
        pval_cauchy1_PCvm_20, pval_cauchy2_PCvm_20, pval_cauchy_PCvm_20, pval_cauchy_hybrid_20
      )
    }
  stopImplicitCluster()
  stopCluster(cl)
  
  # Power calculation 
  return(list(
    PCvM1_05_power = mean(hybrid_power[,1] <= 0.05, na.rm=TRUE),
    PCvM2_05_power = mean(hybrid_power[,2] <= 0.05, na.rm=TRUE),
    PCvM_05_power  = mean(hybrid_power[,3] <= 0.05, na.rm=TRUE),
    Hybrid_05_power= mean(hybrid_power[,4] <= 0.05, na.rm=TRUE),
    
    PCvM1_10_power = mean(hybrid_power[,5] <= 0.05, na.rm=TRUE),
    PCvM2_10_power = mean(hybrid_power[,6] <= 0.05, na.rm=TRUE),
    PCvM_10_power  = mean(hybrid_power[,7] <= 0.05, na.rm=TRUE),
    Hybrid_10_power= mean(hybrid_power[,8] <= 0.05, na.rm=TRUE),
    
    PCvM1_20_power = mean(hybrid_power[,9] <= 0.05, na.rm=TRUE),
    PCvM2_20_power = mean(hybrid_power[,10] <= 0.05, na.rm=TRUE),
    PCvM_20_power  = mean(hybrid_power[,11] <= 0.05, na.rm=TRUE),
    Hybrid_20_power= mean(hybrid_power[,12] <= 0.05, na.rm=TRUE)
  ))
}
 
source("utils.R") # bandwidth_choice  pvalue_integ_Brown
 
# Bandwidth sensitivity analysis
if(T){ 
  # comparison between GRP test,RP test and PLS test with bandwidth sensitivity
  a <- c(0,.5, 1)          
  n <- c(600)       
  p <-  c(50,100,300,600,900,1200)    
  c_h <- 1                     
  pho <-  c(0, .4, .8)   
  r <- 12    
  compar_lm_result <- matrix(0, nrow = r * length(pho) * length(n) * length(a), ncol = length(p) + 4) 
  for (kk in 1:length(pho)){
    for (ll in 1:length(n)){      
      for (ii in 1:length(a)){
        for (jj in 1:length(p)){
          t3 <- Sys.time()  
          result_pls <- Pcvm_Pls_hybrid_logit_h(p[jj], n[ll], a[ii], pho[kk])     
          base_row <- r * length(a) * length(n) * (kk - 1) + r * length(a) * (ll - 1) + ii 
          # 0.5*h
          compar_lm_result[base_row + 0*length(a), jj] <- result_pls$PCvM1_05_power
          compar_lm_result[base_row + 1*length(a), jj] <- result_pls$PCvM2_05_power
          compar_lm_result[base_row + 2*length(a), jj] <- result_pls$PCvM_05_power
          compar_lm_result[base_row + 3*length(a), jj] <- result_pls$Hybrid_05_power
          # 1.0*h (Base)
          compar_lm_result[base_row + 4*length(a), jj] <- result_pls$PCvM1_10_power
          compar_lm_result[base_row + 5*length(a), jj] <- result_pls$PCvM2_10_power
          compar_lm_result[base_row + 6*length(a), jj] <- result_pls$PCvM_10_power
          compar_lm_result[base_row + 7*length(a), jj] <- result_pls$Hybrid_10_power
          # 2.0*h
          compar_lm_result[base_row + 8*length(a), jj] <- result_pls$PCvM1_20_power
          compar_lm_result[base_row + 9*length(a), jj] <- result_pls$PCvM2_20_power
          compar_lm_result[base_row + 10*length(a), jj] <- result_pls$PCvM_20_power
          compar_lm_result[base_row + 11*length(a), jj] <- result_pls$Hybrid_20_power 
          
          t4 <- Sys.time() 
          cat("p = ", p[jj], ", n = ", n[ll], ", a = ", a[ii], ", pho = ", pho[kk], 
              "| PCvM1_power (0.5*h_0) = ", result_pls$PCvM1_05_power,
              ", PCvM1_power (1*h_0)  = ", result_pls$PCvM1_10_power,
              ", PCvM1_power (2*h_0)  = ", result_pls$PCvM1_20_power,
              ", PCvM_power (0.5*h_0) = ", result_pls$PCvM_05_power,
              ", PCvM_power (1*h_0)  = ", result_pls$PCvM_10_power,
              ", PCvM_power (2*h_0)  = ", result_pls$PCvM_20_power,
              ", Hybrid_power (0.5*h_0)  = ", result_pls$Hybrid_05_power,
              ", Hybrid_power (1*h_0)  = ", result_pls$Hybrid_10_power,
              ", Hybrid_power (2*h_0)  = ", result_pls$Hybrid_20_power, 
              " | time:", as.numeric(difftime(t4, t3, units="secs")), "s \n")
          
          for(step in 0:(r-1)){
            compar_lm_result[base_row + step*length(a), (length(p)+1):(length(p)+4)] <- c(n[ll], a[ii], pho[kk], c_h)
          }
          
        }
      }
    }
  }
  write.table(compar_lm_result, "H21_h_sensitivity.txt", sep = "  ", row.names = FALSE, col.names = FALSE)
 
}

