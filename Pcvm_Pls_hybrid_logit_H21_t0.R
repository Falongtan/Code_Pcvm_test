# H=5 t0 0.99, 0.95, 0.90 
Pcvm_Pls_hybrid_logit_t0 <- function(p,n,a,pho){
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
  #library(RPtests)
  #library(GRPtests)
  # library(tictoc)   #library packages
  
  s <- 1000
  mu <- rep(0, p)
  if(pho == 0){
    sigma <- diag(rep(1,p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }
  beta0 <- c(rep(1,5),rep(0,p-5))
  #beta1 <- c(rep(1,10),rep(0,p-10))
  
  #parallel
  # tic()
  cores <- detectCores(logical=F)-2
  random_port <- sample(10000:60000, 1)
  cl <- makeCluster(cores, port = random_port)
  registerDoParallel(cl, cores=cores)
  hybrid_power <- foreach(k = 1:s, .combine='rbind',
                          .packages = c('MASS','harmonicmeanp', 'glmnet','LassoSIR','VariableScreening','psych'),
                          .export = c("pvalue_integ_Brown",'bandwidth_choice')) %dopar%
    {
      set.seed(k)
      # result <- tryCatch({
      # source("pvalue_integ_Brown.R")
      # source("bandwidth_choice.R")
      x <- mvrnorm(n, mu, sigma)
      z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                                # H21
      # z <- x %*% beta0 + a * (x[,1]*x[,2] + x[,2]*x[,3] + x[,3]*x[,4] + x[,4]*x[,5])    # H22
      #      z <- x %*% beta0 + a * 2 * cos(0.6 * x %*% beta0 * pi)                                # H23 
      
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
        #U1 <- sec_model1$residuals
        beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
        beta1_pro <- as.matrix(beta1_pro)
      }
      
      # projections based on screening and lassosir for x1 U1
      screen_num1 <- floor(n1/log(n1))
      sir_Upro1 <- tryCatch({
        if(p <= screen_num1){
          # sir projection without screening
          sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic",
                             solution.path=FALSE, categorical=FALSE, nfolds=5,
                             screening=FALSE)
          sir_Ubeta1<- sir_U1$beta
          sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
        }else{
          #screening
          rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      # screening method based on distance correlation
          index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        # chosen index of X
          x1_Uscreen <- x1[,index_U1]
          sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta1<- sir_Ulasso1$beta
          sir_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
          sir_Upro1[index_U1,] <- sir_Ubeta1
          sir_Upro1/sqrt(colSums(sir_Upro1^2))
        }}, error =function(e){ matrix(NA,p,1) })
      
      # projections based on screening and lassosir for x1 y1
      sir_ypro1 <- tryCatch({
        if(p <= screen_num1){
          # sir projection without screening
          sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta1<- sir_ylasso1$beta
          sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
        }else{
          #screening
          rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                      # screening method based on distance correlation
          index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        # chosen index of X
          x1_yscreen <- x1[,index_y1]
          sir_ylasso1 <- LassoSIR(x1_yscreen, y1, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta1<- sir_ylasso1$beta
          sir_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
          sir_ypro1[index_y1,] <- sir_ybeta1
          sir_ypro1/sqrt(colSums(sir_ypro1^2))
        }
      }, error =function(e){ matrix(NA,p,1) })
      
      # residuals and beta0-projection based on lasso and data x2 y2
      # residuals and beta0-projections based on x2 y2
      lasso_model2 <- cv.glmnet(x2, y2, family="binomial", intercept = T)
      beta2_lasso <- coef(lasso_model2, s = "lambda.min")[-1]
      index_beta2_non0 <- seq(1:p)[as.numeric(beta2_lasso)!=0]
      len_beta2_non0 <- length(index_beta2_non0)
      if(len_beta2_non0 == 0){
        pred2 <- mean(y2)
        U2 <- y2-pred2                                                            # residuals based on x2 y2
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
        #U2 <- sec_model2$residuals
        beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
        beta2_pro <- as.matrix(beta2_pro)
      }
      
      #sir projections based on screening and x2 U2
      screen_num2 <- floor(n2/log(n2))
      sir_Upro2 <- tryCatch({
        if(p <= screen_num2){
          # sir projection without screen
          sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta2<- sir_Ulasso2$beta
          sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
        }else{
          #screening
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      # screening method based on distance correlation
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        # chosen index of X
          x2_Uscreen <- x2[,index_U2]
          sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic",
                                  solution.path=FALSE, categorical=FALSE, nfolds=5,
                                  screening=FALSE)
          sir_Ubeta2 <- sir_Ulasso2$beta
          sir_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
          sir_Upro2[index_U2,] <- sir_Ubeta2
          sir_Upro2/sqrt(colSums(sir_Upro2^2))
        }}, error =function(e){matrix(NA,p,1) })
      
      #sir projections based on screening and x2 y2
      sir_ypro2 <- tryCatch({
        if(p <= screen_num2){
          # sir projection without screen
          sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta2<- sir_ylasso2$beta
          sir_ybeta2/sqrt(colSums(sir_ybeta2^2))
        }else{
          #screening
          rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                      # screening method based on distance correlation
          index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        # chosen index of X
          x2_yscreen <- x2[,index_y2]
          sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=2, choosing.d="automatic",
                                  solution.path=FALSE, categorical=TRUE, nfolds=5,
                                  screening=TRUE)
          sir_ybeta2 <- sir_ylasso2$beta
          sir_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
          sir_ypro2[index_y2,] <- sir_ybeta2
          sir_ypro2/sqrt(colSums(sir_ypro2^2))
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
      
      #construct martingale-based test statistic based on x1 y1 and projections sir_Upro2_orth based on x2 y2
      pro2_num <- ncol(sir_Upro2_orth)                                                                                # the number of SIR-U projections
      pval_matrix1_PCvm_99 <- matrix(nrow=1, ncol=pro2_num+1)                                                            
      pval_matrix1_PCvm_95 <- matrix(nrow=1, ncol=pro2_num+1)                                                            
      pval_matrix1_PCvm_90 <- matrix(nrow=1, ncol=pro2_num+1)
      
      deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
      
      for(q in seq_len(pro2_num)){
        x1_pro <- x1%*%sir_Upro2_orth[,q]
        band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
        h1 <- bandwidth_choice(x1_pro, band_y1)                                                            # choose the bandwidth
        
        Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                  # kernel function matrix
        Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
        kernel_1 <- (3/4)*(1-Ker_inter1^2)*Ker_indictor1
        
        Indictor_1<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
        
        A1_1 <- ((kernel_1%*%deri_link1)*x1_pro)/(kernel_1%*%U1^2+n1^(-12))                                # the 1st compoenent of A1
        A1_2 <- kernel_1%*%((x1 %*% beta1_hat)*deri_link1)/(kernel_1%*%U1^2+n1^(-12))                      # the 2st compoenent of A1
        A1_3 <- (kernel_1%*%deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                         # the 3st compoenent of A1
        hat_A1 <- rbind(t(A1_1), t(A1_2), t(A1_3))
        
        # gamma matrix
        Gamma_inv1 <- array(0, dim=c(3, 3, n1))
        for (i in 1:n1) {
          Gamma_inv1[,,i] = ginv(((rep(1,3)%*%t(U1)^2)*hat_A1)%*%(t(hat_A1) * (as.matrix(Indictor_1[i,]) %*%rep(1,3))))
        }
        integral_1 <- ((rep(1,3)%*%t(U1))*hat_A1)%*%t(Indictor_1)
        
        martingle_sec1<-diag(0,n1)  # the second term in the martingale transformation
        for (l in 1:n1){
          martingle_sec1[l,]<- (U1[l]^2 * t(hat_A1[,l])) %*% Gamma_inv1[,,l] %*% integral_1[,l] %*% Indictor_1[l,]
        }
        martingle_sta1<-(1/sqrt(n1))*t(U1) %*% Indictor_1 - (1/sqrt(n1)) * colSums(martingle_sec1)   # the test statistic of the martingale transformation
        
        ordx1_pro <- sort(x1_pro)  # order x%*%projection
        
        # 0.99
        t_1_99 <- ordx1_pro[floor(0.99 * n1)]
        psi_1_99 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1_99))
        PCvM_martingle1_99 <- (1/psi_1_99^2) * mean(U1^2 * (x1_pro <= t_1_99) * (t(martingle_sta1)^2))  
        pval_matrix1_PCvm_99[,q] <- pvalue_integ_Brown(PCvM_martingle1_99)
        
        # 0.95
        t_1_95 <- ordx1_pro[floor(0.95 * n1)]
        psi_1_95 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1_95))
        PCvM_martingle1_95 <- (1/psi_1_95^2) * mean(U1^2 * (x1_pro <= t_1_95) * (t(martingle_sta1)^2))  
        pval_matrix1_PCvm_95[,q] <- pvalue_integ_Brown(PCvM_martingle1_95)
        
        # 0.90
        t_1_90 <- ordx1_pro[floor(0.90 * n1)]
        psi_1_90 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1_90))
        PCvM_martingle1_90 <- (1/psi_1_90^2) * mean(U1^2 * (x1_pro <= t_1_90) * (t(martingle_sta1)^2))  
        pval_matrix1_PCvm_90[,q] <- pvalue_integ_Brown(PCvM_martingle1_90)
      }
      
      #construct martingale-based test statistic based on x1 y1 and projections beta2_pro based on x2 y2
      x1_betapro <- x1%*%beta2_pro
      beta_Indictor1<-ifelse(x1_betapro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_betapro),1,0)
      
      hat_beta_A1 <- rbind(t(x1_betapro), matrix(1,1,n1))
      # gamma matrix
      Gamma_beta_inv1 <- array(0, dim=c(2, 2, n1))
      for (i in 1:n1) {
        Gamma_beta_inv1[,,i] = ginv(((rep(1,2)%*%t(U1)^2)*hat_beta_A1)%*%(t(hat_beta_A1) * (as.matrix(beta_Indictor1[i,]) %*%rep(1,2))))
      }
      integral_beta_1 <- ((rep(1,2)%*%t(U1))*hat_beta_A1)%*%t(beta_Indictor1)                      # integral part in the martingale transformation
      
      martingle_beta_sec1<-diag(0,n1)                                                              # The second term in the martingale transformation
      for (l in 1:n1){
        martingle_beta_sec1[l,]<- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,]
      }
      martingle_beta_sta1<-(1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      # the test statistic of the martingale transformation
      
      ordx1_betapro <- sort(x1_betapro)  # order x%*%projection
      
      # 0.99
      beta_t1_99 <- ordx1_betapro[floor(0.99 * n1)]
      psi_beta_1_99 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1_99))
      PCvM_beta_martingle1_99 <- (1/psi_beta_1_99^2) * mean(U1^2 * (x1_betapro <= beta_t1_99) * (t(martingle_beta_sta1)^2))   
      pval_matrix1_PCvm_99[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1_99)
      
      # 0.95
      beta_t1_95 <- ordx1_betapro[floor(0.95 * n1)]
      psi_beta_1_95 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1_95))
      PCvM_beta_martingle1_95 <- (1/psi_beta_1_95^2) * mean(U1^2 * (x1_betapro <= beta_t1_95) * (t(martingle_beta_sta1)^2))   
      pval_matrix1_PCvm_95[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1_95)
      
      # 0.90
      beta_t1_90 <- ordx1_betapro[floor(0.90 * n1)]
      psi_beta_1_90 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1_90))
      PCvM_beta_martingle1_90 <- (1/psi_beta_1_90^2) * mean(U1^2 * (x1_betapro <= beta_t1_90) * (t(martingle_beta_sta1)^2))   
      pval_matrix1_PCvm_90[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1_90)
      
      #construct martingale-based test statistics based on x2 y2 and projections sir_Upro1_orth based on x1 y1
      pro1_num <- ncol(sir_Upro1_orth)                                                                             # the number of projections
      pval_matrix2_PCvm_99 <- matrix(nrow=1, ncol=pro1_num+1)                                                            
      pval_matrix2_PCvm_95 <- matrix(nrow=1, ncol=pro1_num+1)                                                            
      pval_matrix2_PCvm_90 <- matrix(nrow=1, ncol=pro1_num+1)
      deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
      
      for(qq in seq_len(pro1_num)){
        x2_pro <- x2%*%sir_Upro1_orth[,qq]
        band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
        h2 <- bandwidth_choice(x2_pro, band_y2)                                                            # choose the bandwidth
        
        Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                                  # kernel function matrix
        Ker_indictor2 <- ifelse(abs(Ker_inter2) <= 1, 1, 0)
        kernel_2 <- (3/4)*(1-Ker_inter2^2)*Ker_indictor2
        
        Indictor_2<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
        
        A2_1 <- ((kernel_2%*%deri_link2)*x2_pro)/(kernel_2%*%U2^2+n2^(-12))                                # the 1st compoenent of A1
        A2_2 <- kernel_2%*%((x2 %*% beta2_hat)*deri_link2)/(kernel_2%*%U2^2+n2^(-12))                      # the 2st compoenent of A1
        A2_3 <- (kernel_2%*%deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                         # the 3st compoenent of A1
        hat_A2 <- rbind(t(A2_1), t(A2_2), t(A2_3))
        
        #gamma matrix
        Gamma_inv2 <- array(0, dim=c(3, 3, n2))
        for (ii in 1:n2){
          Gamma_inv2[,,ii] = ginv(((rep(1,3)%*%t(U2)^2)*hat_A2)%*%(t(hat_A2) * (as.matrix(Indictor_2[ii,]) %*%rep(1,3))))
        }
        integral_2 <- ((rep(1,3)%*%t(U2))*hat_A2)%*%t(Indictor_2)
        
        martingle_sec2<-diag(0,n2)                                                                         # the second term in the martingale transformation
        for (ll in 1:n2){
          martingle_sec2[ll,]<- (U2[ll]^2 * t(hat_A2[,ll])) %*% Gamma_inv2[,,ll] %*% integral_2[,ll] %*% Indictor_2[ll,]
        }
        martingle_sta2<-(1/sqrt(n2))*t(U2) %*% Indictor_2 - (1/sqrt(n2)) * colSums(martingle_sec2)  # the test statistic of the martingale transformation
        
        ordx2_pro <- sort(x2_pro)                                                                          # order x%*%projection
        
        # 0.99
        t_2_99 <- ordx2_pro[floor(0.99 * n2)]
        psi_2_99 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2_99))
        PCvM_martingle2_99 <- (1/psi_2_99^2) * mean(U2^2 * (x2_pro <= t_2_99) * (t(martingle_sta2)^2))       
        pval_matrix2_PCvm_99[,qq] <- pvalue_integ_Brown(PCvM_martingle2_99)
        
        # 0.95
        t_2_95 <- ordx2_pro[floor(0.95 * n2)]
        psi_2_95 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2_95))
        PCvM_martingle2_95 <- (1/psi_2_95^2) * mean(U2^2 * (x2_pro <= t_2_95) * (t(martingle_sta2)^2))       
        pval_matrix2_PCvm_95[,qq] <- pvalue_integ_Brown(PCvM_martingle2_95)
        
        # 0.90
        t_2_90 <- ordx2_pro[floor(0.90 * n2)]
        psi_2_90 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2_90))
        PCvM_martingle2_90 <- (1/psi_2_90^2) * mean(U2^2 * (x2_pro <= t_2_90) * (t(martingle_sta2)^2))       
        pval_matrix2_PCvm_90[,qq] <- pvalue_integ_Brown(PCvM_martingle2_90)
      }
      
      #construct martingale-based test statistics based on x2 y2 and projections beta1_pro based on x1 y1
      x2_betapro <- x2%*%beta1_pro
      beta_Indictor2<-ifelse(x2_betapro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_betapro),1,0)
      
      hat_beta_A2 <- rbind(t(x2_betapro), matrix(1,1,n2))
      
      # gamma matrix
      Gamma_beta_inv2 <- array(0, dim=c(2, 2, n2))
      for (ii in 1:n2) {
        Gamma_beta_inv2[,,ii] = ginv(((rep(1,2)%*%t(U2)^2)*hat_beta_A2)%*%(t(hat_beta_A2) * (as.matrix(beta_Indictor2[ii,]) %*%rep(1,2))))
      }
      integral_beta_2 <- ((rep(1,2)%*%t(U2))*hat_beta_A2)%*%t(beta_Indictor2)
      
      martingle_beta_sec2<-diag(0,n2)                                                                    # The second term in the martingale transformation
      for (ll in 1:n2){
        martingle_beta_sec2[ll,]<- (U2[ll]^2 * t(hat_beta_A2[,ll])) %*% Gamma_beta_inv2[,,ll] %*% integral_beta_2[,ll] %*% beta_Indictor2[ll,]
      }
      martingle_beta_sta2<-(1/sqrt(n2))*t(U2) %*% beta_Indictor2 - (1/sqrt(n2)) * colSums(martingle_beta_sec2)   # The test statistic of the martingale transformation
      
      ordx2_betapro <- sort(x2_betapro)                                                                  # order x%*%proj
      
      # 0.99
      beta_t2_99 <- ordx2_betapro[floor(0.99 * n2)]
      psi_beta_2_99 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2_99))
      PCvM_beta_martingle2_99 <- (1/psi_beta_2_99^2) * mean(U2^2 * (x2_betapro <= beta_t2_99) * (t(martingle_beta_sta2)^2))   
      pval_matrix2_PCvm_99[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2_99)
      
      # 0.95
      beta_t2_95 <- ordx2_betapro[floor(0.95 * n2)]
      psi_beta_2_95 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2_95))
      PCvM_beta_martingle2_95 <- (1/psi_beta_2_95^2) * mean(U2^2 * (x2_betapro <= beta_t2_95) * (t(martingle_beta_sta2)^2))   
      pval_matrix2_PCvm_95[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2_95)
      
      # 0.90
      beta_t2_90 <- ordx2_betapro[floor(0.90 * n2)]
      psi_beta_2_90 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2_90))
      PCvM_beta_martingle2_90 <- (1/psi_beta_2_90^2) * mean(U2^2 * (x2_betapro <= beta_t2_90) * (t(martingle_beta_sta2)^2))   
      pval_matrix2_PCvm_90[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2_90)
      
      pval_num <- pro1_num + 1 + pro2_num + 1 
      # 0.99
      pval_matrix_PCvm_99 <- cbind(pval_matrix1_PCvm_99, pval_matrix2_PCvm_99)
      pval_cauchy1_PCvm_99 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_99[1:(pro2_num+1)])*pi)))              
      pval_cauchy2_PCvm_99 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_99[1:(pro1_num+1)])*pi)))              
      pval_cauchy_PCvm_99  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_99[1:pval_num])*pi)))                
      
      # 0.95
      pval_matrix_PCvm_95 <- cbind(pval_matrix1_PCvm_95, pval_matrix2_PCvm_95)
      pval_cauchy1_PCvm_95 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_95[1:(pro2_num+1)])*pi)))              
      pval_cauchy2_PCvm_95 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_95[1:(pro1_num+1)])*pi)))              
      pval_cauchy_PCvm_95  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_95[1:pval_num])*pi)))  
      
      # 0.90
      pval_matrix_PCvm_90 <- cbind(pval_matrix1_PCvm_90, pval_matrix2_PCvm_90)
      pval_cauchy1_PCvm_90 <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm_90[1:(pro2_num+1)])*pi)))              
      pval_cauchy2_PCvm_90 <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm_90[1:(pro1_num+1)])*pi)))              
      pval_cauchy_PCvm_90  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm_90[1:pval_num])*pi)))
      
      # construct PLS test statistic
      if(T){
        #construct PLS test statistic based on x1 y1 and projections based on x2 y2
        pro2_pls <- cbind(sir_Upro2_orth, beta2_pro)
        pro2_num_pls <- ncol(pro2_pls)                                                              # the number of projections
        h1 <- n1^(-2/9)                                                                             # bandwidth
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) # p-value matrix
        errormat1 <- U1%*%t(U1)                                                                     # residual matrix based on x1 y1
        for(q in 1:pro2_num_pls){
          x_pro1 <- x1%*%pro2_pls[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel
          #PLS test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval1 <- 1-pnorm(Tn1)
          pval_matrix1_PLS[,q] <- pval1
        }
        
        #construct PLS test statistics based on x2 y2 and projections based on x1 y1
        pro1_pls <- cbind(sir_Upro1_orth, beta1_pro)
        pro1_num_pls <- ncol(pro1_pls)                                                              # the number of projections
        h2 <- n2^(-2/9)                                                                             # bandwidth
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) # p-value matrix
        errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
        for(l in 1:pro1_num_pls){
          x_pro2 <- x2%*%pro1_pls[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2   # kernel function matrix
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                 # Epanechnikov kernel
          # PLS test statistics
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval2 <- 1-pnorm(Tn2)
          pval_matrix2_PLS[,l] <- pval2
        }
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
        pval_num_pls <- pro1_num_pls + pro2_num_pls
      }
      
      pval_matrix_hybrid_99 <- cbind(pval_matrix_PCvm_99, pval_matrix_PLS)
      pval_cauchy_hybrid_99 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_99[1:ncol(pval_matrix_hybrid_99)])*pi)))
      
      pval_matrix_hybrid_95 <- cbind(pval_matrix_PCvm_95, pval_matrix_PLS)
      pval_cauchy_hybrid_95 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_95[1:ncol(pval_matrix_hybrid_95)])*pi)))
      
      pval_matrix_hybrid_90 <- cbind(pval_matrix_PCvm_90, pval_matrix_PLS)
      pval_cauchy_hybrid_90 <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid_90[1:ncol(pval_matrix_hybrid_90)])*pi)))
      
      # construct PLS test statistic
      if(T){
        #construct PLS test statistic based on x1 y1 and projections based on x2 y2
        pro2_pls <- cbind(sir_Upro2, sir_ypro2)                                                               # projections based on x2 y2
        pro2_num_pls <- ncol(pro2_pls)                                                                        # the number of projections
        h1 <- n1^(-2/9)                                                                                       # bandwidth
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) # p-value matrix
        errormat1 <- U1%*%t(U1)                                                                               # residual matrix based on x1 y1
        for(q in 1:pro2_num_pls){
          x_pro1 <- x1%*%pro2_pls[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1            # kernel function matrix
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                          # Epanechnikov kernel
          #PLS test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval1 <- 1-pnorm(Tn1)
          pval_matrix1_PLS[,q] <- pval1
        }
        
        #construct PLS test statistics based on x2 y2 and projections based on x1 y1
        pro1_pls <- cbind(sir_Upro1, sir_ypro1)                                                               # projections based on x1 y1
        pro1_num_pls <- ncol(pro1_pls)                                                                        # the number of projections
        h2 <- n2^(-2/9)                                                                                       # bandwidth
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) # p-value matrix
        errormat2 <- U2%*%t(U2)                                                                               #residual matrix based on x2 y2
        for(l in 1:pro1_num_pls){
          x_pro2 <- x2%*%pro1_pls[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        # kernel function matrix
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      # Epanechnikov kernel
          # PLS test statistics
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval2 <- 1-pnorm(Tn2)
          pval_matrix2_PLS[,l] <- pval2
        }
        
        pval_fisher1_PLS <- 1- pchisq(-2*sum(log(pval_matrix1_PLS[1:pro2_num_pls])),df=2*pro2_num_pls)
        pval_fisher2_PLS <- 1- pchisq(-2*sum(log(pval_matrix2_PLS[1:pro1_num_pls])),df=2*pro1_num_pls)
        pval_fisher_m <- cbind(pval_fisher1_PLS, pval_fisher2_PLS)
        pval_fisher_cauchy <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))
      }
      
      result <- c(
        pval_cauchy1_PCvm_99, pval_cauchy2_PCvm_99, pval_cauchy_PCvm_99, pval_cauchy_hybrid_99, 
        pval_cauchy1_PCvm_95, pval_cauchy2_PCvm_95, pval_cauchy_PCvm_95, pval_cauchy_hybrid_95,  
        pval_cauchy1_PCvm_90, pval_cauchy2_PCvm_90, pval_cauchy_PCvm_90, pval_cauchy_hybrid_90,  
        pval_fisher1_PLS, pval_fisher2_PLS, pval_fisher_cauchy  
      )
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)
  
  # 0.99
  PCvM1_99_power <- mean(hybrid_power[,1] <= 0.05,na.rm=TRUE)
  PCvM2_99_power <- mean(hybrid_power[,2] <= 0.05,na.rm=TRUE)
  PCvM_99_power  <- mean(hybrid_power[,3] <= 0.05,na.rm=TRUE)
  Hybrid_99_power<- mean(hybrid_power[,4] <= 0.05,na.rm=TRUE)
  
  # 0.95
  PCvM1_95_power <- mean(hybrid_power[,5] <= 0.05,na.rm=TRUE)
  PCvM2_95_power <- mean(hybrid_power[,6] <= 0.05,na.rm=TRUE)
  PCvM_95_power  <- mean(hybrid_power[,7] <= 0.05,na.rm=TRUE)
  Hybrid_95_power<- mean(hybrid_power[,8] <= 0.05,na.rm=TRUE)
  
  # 0.90
  PCvM1_90_power <- mean(hybrid_power[,9] <= 0.05,na.rm=TRUE)
  PCvM2_90_power <- mean(hybrid_power[,10] <= 0.05,na.rm=TRUE)
  PCvM_90_power  <- mean(hybrid_power[,11] <= 0.05,na.rm=TRUE)
  Hybrid_90_power<- mean(hybrid_power[,12] <= 0.05,na.rm=TRUE)
  
  # PLS
  PLS1_power <- mean(hybrid_power[,13] <= 0.05,na.rm=TRUE)
  PLS2_power <- mean(hybrid_power[,14] <= 0.05,na.rm=TRUE)
  PLS_power  <- mean(hybrid_power[,15] <= 0.05,na.rm=TRUE)
  
  return(list(
    PCvM1_99_power = PCvM1_99_power, PCvM2_99_power = PCvM2_99_power, PCvM_99_power = PCvM_99_power,
    PCvM1_95_power = PCvM1_95_power, PCvM2_95_power = PCvM2_95_power, PCvM_95_power = PCvM_95_power,
    PCvM1_90_power = PCvM1_90_power, PCvM2_90_power = PCvM2_90_power, PCvM_90_power = PCvM_90_power,
    PLS1_power = PLS1_power, PLS2_power = PLS2_power, PLS_power = PLS_power,
    Hybrid_99_power = Hybrid_99_power, Hybrid_95_power = Hybrid_95_power, Hybrid_90_power = Hybrid_90_power
  ))
}
 
# grp tests for testing logistic regression models
grp_logit_parallel<- function(p,n,a,pho){
  library(MASS)
  library(matrixcalc)
  library(RPtests)
  library(GRPtests)
  # library(tictoc)
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
  #parallel
  # tic()
  cores <- detectCores(logical=F)-2
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores=cores)
  grp_power <- foreach(k = 1:s, .combine='rbind',
                       # .packages = c('MASS','matrixcalc','RPtests','GRPtests','foreach','parallel','iterators','doParallel','tictoc'),
                       .packages = c('MASS','matrixcalc','RPtests','GRPtests')) %dopar%
    {
      set.seed(k)
      x <- mvrnorm(n, mu, sigma)
      z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                             # H21
       # z <- x %*% beta0 + a * (x[,1]*x[,2] + x[,2]*x[,3] + x[,3]*x[,4] + x[,4]*x[,5])   # H22
#      z <- x %*% beta0 + a * 2 * cos(0.6 * x %*% beta0 * pi)                                # H23
      pr <- 1/(1 + exp(-z))
      y <- matrix(rbinom(n, 1, pr),ncol = 1)
      pval_grp <- GRPtest(x, y, fam = "binomial", nsplits = 1,RP_function = NULL,output_all=FALSE,penalize=TRUE)
      result <- c(pval_grp)
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)

  # toc()
  #power
  grp_power <- mean(grp_power<=0.05,na.rm=TRUE)
  return(grp_power)
}

# The choice for the bandwidth for martingale-based test
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

 
a <- c(0, .5, 1)
n <- 600                      
p <- c(50, 100, 300, 600, 900, 1200)   
c_h <- 1                      
pho <- c(0, .4, .8)           
r <- 16
compar_lm_result <- matrix(0, nrow = r*length(pho)*length(a), ncol = length(p)+4)
for (kk in 1:length(pho)){
  for (ii in 1:length(a)){
    for (jj in 1:length(p)){
      t3 <- Sys.time()
      
      set.seed(12345 + jj)
      result_pls <- Pcvm_Pls_hybrid_logit_t0(p[jj], n, a[ii], pho[kk])
      
      set.seed(12345 + jj)
      result_grp_rp <- grp_logit_parallel(p[jj], n, a[ii], pho[kk]) 
      
      base_row <- r * length(a) * (kk - 1) + ii
      
      # PCvM
      compar_lm_result[base_row + 0*length(a), jj] <- result_pls$PCvM1_99_power
      compar_lm_result[base_row + 1*length(a), jj] <- result_pls$PCvM2_99_power
      compar_lm_result[base_row + 2*length(a), jj] <- result_pls$PCvM_99_power
      compar_lm_result[base_row + 3*length(a), jj] <- result_pls$PCvM1_95_power
      compar_lm_result[base_row + 4*length(a), jj] <- result_pls$PCvM2_95_power
      compar_lm_result[base_row + 5*length(a), jj] <- result_pls$PCvM_95_power
      compar_lm_result[base_row + 6*length(a), jj] <- result_pls$PCvM1_90_power
      compar_lm_result[base_row + 7*length(a), jj] <- result_pls$PCvM2_90_power
      compar_lm_result[base_row + 8*length(a), jj] <- result_pls$PCvM_90_power
      # PLS
      compar_lm_result[base_row + 9*length(a), jj]  <- result_pls$PLS1_power
      compar_lm_result[base_row + 10*length(a), jj] <- result_pls$PLS2_power
      compar_lm_result[base_row + 11*length(a), jj] <- result_pls$PLS_power
      # Hybrid
      compar_lm_result[base_row + 12*length(a), jj] <- result_pls$Hybrid_99_power
      compar_lm_result[base_row + 13*length(a), jj] <- result_pls$Hybrid_95_power
      compar_lm_result[base_row + 14*length(a), jj] <- result_pls$Hybrid_90_power
      # GRP
      compar_lm_result[base_row + 15*length(a), jj] <- result_grp_rp
      
      for(step in 0:15){
        compar_lm_result[base_row + step*length(a), (length(p)+1):(length(p)+4)] <- c(n, a[ii], pho[kk], c_h)
      }
      
      t4 <- Sys.time()
      
      cat(sprintf("p=%d, n=%d, a=%.1f, pho=%.1f | time=%.2fs", p[jj], n, a[ii], pho[kk], as.numeric(difftime(t4, t3, units="secs"))))
      cat(sprintf(", PCvM(99/95/90)   = %.3f / %.3f / %.3f", result_pls$PCvM_99_power, result_pls$PCvM_95_power, result_pls$PCvM_90_power))
      cat(sprintf(", Hybrid(99/95/90) = %.3f / %.3f / %.3f", result_pls$Hybrid_99_power, result_pls$Hybrid_95_power, result_pls$Hybrid_90_power))
      cat(sprintf(", PLS=%.3f, GRP=%.3f\n", result_pls$PLS_power, result_grp_rp))
    } 
  }
}
write.table(compar_lm_result, "H21_t0.txt", sep = "  ", row.names = FALSE, col.names = FALSE)