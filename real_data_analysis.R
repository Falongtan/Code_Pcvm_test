# The hybrid test, the martingale-based test, and local-smoothing test for testing high dimensional linear and logistical regression models

# linear + test
Pcvm_Pls_hybrid_linear <- function(p,n,a,pho){ 
  library(MASS)
  library(glmnet)
  library(LassoSIR)
  library(VariableScreening)
  library(psych)
  library(foreach)
  library(parallel)
  library(iterators)
  library(doParallel)
  #library packages
  s <- 1000 
  mu <- rep(0, p)
  if(pho == 0){
    sigma <- diag(rep(1,p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }   
  beta0 <- c(rep(1,5),rep(0,p-5))
  beta1 <- c(rep(1,10),rep(0,p-10)) # /sqrt(10)
  #parallel
  # tic()
  cores <- detectCores(logical=F)-2
  cl <- makeCluster(cores) 
  registerDoParallel(cl, cores=cores)
  hybrid_linear_power <- foreach(k = 1:s, .combine='rbind', 
                                 .packages = c('MASS','glmnet','LassoSIR','VariableScreening','psych'),
                                 .export = c("pvalue_integ_Brown")) %dopar% 
    {
      set.seed(k)
      x <- mvrnorm(n, mu, sigma)
      y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    # H11
      # y <- x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n)            # H12
      # y <- x %*% beta0 + a *exp(0.25 * x%*%beta1) + rnorm(n)                          # H13
      # y <- x %*% beta0 + a *exp(0.5 * x%*%beta1) + rnorm(n)                          # H13
      
      #data splitting  
      index_x<-sample(1:n, floor(n/2),replace = FALSE)  
      index_x<-sort(index_x)
      
      x1 <- x[-index_x,] 
      y1 <- as.matrix(y[-index_x,]) 
      
      x2 <- x[index_x,]
      y2 <- as.matrix(y[index_x,])
      
      n1 <- nrow(x1)
      n2 <- nrow(x2)
      
      # residuals and beta_0 projections based on lasso and data x1 y1
      lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
      lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
      # index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
      index_beta1_non0 <- which(is.finite(lasso_beta1) & lasso_beta1 != 0)
      len_beta1_non0 <- length(index_beta1_non0)
      if(len_beta1_non0 == 0){
        U1 <- y1 - mean(y1)                                                        # residual based on x1 y1
        beta1_hat <- lasso_beta1
        beta1_pro <- rep(1,p)/sqrt(p)
      }else{
        # x1_sec <- x1[,index_beta1_non0]                                               # second estimation
        x1_sec <- x1[, index_beta1_non0, drop=FALSE]
        sec_model1 <- glm(y1~x1_sec, family = gaussian)
        sec_beta1 <- unname(sec_model1$coefficients)[-1]
        sec_beta1[!is.finite(sec_beta1)] <- 0
        lasso_beta1[index_beta1_non0] <- sec_beta1
        beta1_hat <- lasso_beta1
        pred1 = predict(sec_model1, type="response")
        pred1 = matrix(unname(pred1), ncol = 1)
        U1 <- y1-pred1                                                             # residual based on x1 y1
        #U1 <- sec_model1$residuals
        beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
        beta1_pro <- as.matrix(beta1_pro)
      }
      
      # construct projections
      # projections based on screening and lassosir for x1 U1
      screen_num1 <- floor(n1/log(n1))
      if(p <= screen_num1){
        # projection based on lassosir without screening
        sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, 
                           categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta1<- sir_U1$beta
        sdr_Upro1 <- sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
      }else{
        #screening
        rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
        index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
        x1_resi_screen <- x1[,index_U1, drop=FALSE]
        # projection based on lassosir after screening
        sir_U1 <- LassoSIR(x1_resi_screen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                           nfolds=5, screening=FALSE)
        sir_Ubeta1<- sir_U1$beta
        sdr_Upro1 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
        sdr_Upro1[index_U1,] <- sir_Ubeta1
        sdr_Upro1 <- sdr_Upro1/sqrt(colSums(sdr_Upro1^2))
      } 
      
      # projections based on screening and lassosir for x1 y1
      if(p <= screen_num1){
        # laosssir projection without screening
        sir_y1 <- LassoSIR(x1, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                           nfolds=5, screening=FALSE)
        sir_ybeta1<- sir_y1$beta
        sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
      }else{
        #screening
        rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                     # screening method 
        index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
        x1_yscreen <- x1[,index_y1]
        
        # projection based on lassosir after screening
        sir_y1 <- LassoSIR(x1_yscreen, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                           nfolds=5, screening=FALSE)
        sir_ybeta1<- sir_y1$beta
        sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
        sdr_ypro1[index_y1,] <- sir_ybeta1
        sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
      }
      
      # residuals and beta0-projection based on lasso and data x2 y2
      lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
      lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
      # index_beta2_non0 <- seq(1:p)[as.numeric(lasso_beta2) != 0]
      index_beta2_non0 <- which(is.finite(lasso_beta2) & lasso_beta2 != 0)
      len_beta2_non0 <- length(index_beta2_non0)
      if(len_beta2_non0 == 0){
        U2 <- y2 - mean(y2)                                                        # residual based on x2 y2
        beta2_hat <- lasso_beta2
        beta2_pro <- rep(1,p)/sqrt(p)
      }else{
        x2_sec <- x2[,index_beta2_non0, drop=FALSE]                                               # second estimation
        sec_model2 <- glm(y2~x2_sec, family = gaussian)
        sec_beta2 <- unname(sec_model2$coefficients)[-1]
        sec_beta2[!is.finite(sec_beta2)] <- 0
        lasso_beta2[index_beta2_non0] <- sec_beta2
        beta2_hat <- lasso_beta2
        pred2 = predict(sec_model2, type="response")
        pred2 = matrix(unname(pred2), ncol = 1)
        U2 <- y2-pred2                                                             # residual based on x2 y2
        #U2 <- sec_model2$residuals
        beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
        beta2_pro <- as.matrix(beta2_pro)
      }
      
      
      # projections based on screening and lassosir for x2 U2
      screen_num2 <- floor(n2/log(n2))
      if(p <= screen_num2){
        # projection based on lassosir without screening
        sir_U2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, 
                           categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta2<- sir_U2$beta
        sdr_Upro2 <- sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
      }else{
        #screening
        rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
        index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
        x2_resi_screen <- x2[,index_U2, drop=FALSE]
        # projection based on lassosir after screening
        sir_U2 <- LassoSIR(x2_resi_screen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                           nfolds=5, screening=FALSE)
        sir_Ubeta2<- sir_U2$beta
        sdr_Upro2 <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
        sdr_Upro2[index_U2,] <- sir_Ubeta2
        sdr_Upro2 <- sdr_Upro2/sqrt(colSums(sdr_Upro2^2))
      }
      
      #sir projections based on screening and lassosir for x2 y2
      if(p <= screen_num2){
        # Lassosir projection without screen
        sir_y2 <- LassoSIR(x2, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                           nfolds=5, screening=FALSE)
        sir_ybeta2<- sir_y2$beta
        sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                        # sdr projection based on lassosir and x2 y2
      }else{ 
        #screening
        rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                        # screening method based on distance correlation
        index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                          # chosen index of X
        x2_yscreen <- x2[,index_y2]
        
        # projection based on lassosir after screening
        sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, 
                                nfolds=5, screening=FALSE)
        sir_ybeta2 <- sir_ylasso2$beta
        sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
        sdr_ypro2[index_y2,] <- sir_ybeta2
        sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                          # sdr projection based on lassosir and x2 y2
      } 
      
      # Schmidt orthogonalization
      pro1 <- cbind(beta1_pro, sdr_Upro1)
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
      
      pro2 <- cbind(beta2_pro, sdr_Upro2)
      pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
      # tol <- 1e-12
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
      
      #construct martingale-based test based on x1 y1 and projections based on x2 y2
      pro2 <- pro2_orth                                                            # projections based on x2 y2
      pro2_num <- ncol(pro2)                                                       # the number of projections
      pval_matrix1_Pcvm <- matrix(nrow=1, ncol=pro2_num)                           # p-value matrix of martingale-based test based on x1 y1
      
      for(q in 1:pro2_num){
        x1_pro <- x1%*%pro2[,q]                                                    # projected x1 with pro2 
        x1_Indictor<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
        
        hat_A1 <- rbind(t(x1_pro), rep(1,n1))  
        Gamma1_inv <- array(0, dim=c(2, 2, n1))
        for (i in 1:n1) {
          Gamma1_inv[,,i] = ginv((hat_A1%*%diag(x1_Indictor[i,]))%*%t(hat_A1))
        }                                                                          # The inverse of the matrix gamma
        Integral_1 <- hat_A1%*%diag(as.vector(U1))%*%t(x1_Indictor)                # The integral term in the martingale transformation 
        
        #The test statistic of the martingale transformation
        martingle1_sec<-diag(0,n1)                                                 # the second term in the martingale transformation
        for (l in 1:n1){
          martingle1_sec[l,]<- (t(hat_A1[,l]) %*% Gamma1_inv[,,l] %*% Integral_1[,l]) %*% x1_Indictor[l,] 
        }
        martingle_sta1<-(1/sqrt(n1))*t(U1) %*% x1_Indictor - (1/sqrt(n1)) * colSums(martingle1_sec)  # The test statistic of martingale transformation 
        
        ordx1_pro <- sort(x1_pro)                                                   # order x1_pro
        t_1 <- ordx1_pro[floor(0.99 * n1)]                                          # 99% quantile of x1_pro 
        sigma1_square <- mean(U1^2)   # estimation of variance of error
        F_1 <- (mean(x1_pro <= t_1))^2
        PcvM_martingle1 <- (1/(sigma1_square*F_1)) * mean((x1_pro <= t_1) * (t(martingle_sta1)^2))  #The CvM test statistic based on martingale transformation 
        pval_matrix1_Pcvm[,q] <- pvalue_integ_Brown(PcvM_martingle1)
      }
      
      #construct martingale-based test based on x2 y2 and projections based on x1 y1                                        
      pro1 <- pro1_orth                                                             # projections based on x1 y1
      pro1_num <- ncol(pro1)                                                        # the number of projections
      pval_matrix2_Pcvm <- matrix(nrow=1, ncol=pro1_num)                            # p-value matrix based on x2 y2
      
      for(qq in 1:pro1_num){
        x2_pro <- x2%*%pro1[,qq]                                                    # projected x1 with pro2 
        x2_Indictor<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
        
        hat_A2 <- rbind(t(x2_pro), rep(1,n2))  
        Gamma2_inv <- array(0, dim=c(2, 2, n2))
        for (ii in 1:n2) {
          Gamma2_inv[,,ii] = ginv((hat_A2%*%diag(x2_Indictor[ii,]))%*%t(hat_A2))
        }                                                                          # The inverse of the matrix gamma
        Integral_2 <- hat_A2%*%diag(as.vector(U2))%*%t(x2_Indictor)                # The integral term in the martingale transformation 
        
        #The test statistic of the martingale transformation
        martingle2_sec<-diag(0,n2)                                                 # The second term in the martingale transformation
        for (ll in 1:n2){
          martingle2_sec[ll,]<- (t(hat_A2[,ll]) %*% Gamma2_inv[,,ll] %*% Integral_2[,ll]) %*% x2_Indictor[ll,] 
        }
        martingle_sta2<-(1/sqrt(n2))*t(U2) %*% x2_Indictor - (1/sqrt(n2)) * colSums(martingle2_sec)  # The test statistic of martingale transformation
        
        ordx2_pro <- sort(x2_pro)                                                  # order x2_pro
        t_2 <- ordx2_pro[floor(0.99 * n2)]                                         # 99% quantile of x2_pro 
        sigma2_square <- mean(U2^2)                                                # estimation of variance of error
        F_2 <- (mean(x2_pro <= t_2))^2
        PcvM_martingle2 <- (1/(sigma2_square*F_2)) * mean((x2_pro <= t_2) * (t(martingle_sta2)^2))  #The CvM test statistic based on martingale transformation 
        pval_matrix2_Pcvm[,qq] <- pvalue_integ_Brown(PcvM_martingle2)
      }
      pval_matrix_Pcvm <- cbind(pval_matrix1_Pcvm, pval_matrix2_Pcvm)
      pval_num <- pro1_num + pro2_num
      
      pval_cauchy1_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_Pcvm[1:pro2_num])*pi)))  # cauchy combination of martingale-based test based on x1 y1
      pval_cauchy2_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_Pcvm[1:pro1_num])*pi)))  # cauchy combination of martingale-based test based on x1 y1
      pval_cauchy_Pcvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_Pcvm[1:pval_num])*pi)))   # cauchy combination of martingale-based test based on x y
      
      # construct hybrid test statistic
      if(T){
        # construct hybrid test statistic with martingale test statistic and local smoothing test based on x1 y1 and projections based on x2 y2
        pro2_pls_num <- ncol(pro2_orth)                                               # the number of projections based on x2 y2
        h1 <- n1^(-2/9)                                                               # bandwidth 
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_pls_num) # p-value matrix
        errormat1 <- U1%*%t(U1)                                                       #residual matrix based on x1 y1
        for(q in 1:pro2_pls_num){
          x_pro1 <- x1%*%pro2_orth[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
          #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                            # Gaussian kernel
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel 
          #test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
        }
        
        #construct test statistics based on x2 y2 and projections based on x1 y1                                       
        pro1_pls_num <- ncol(pro1_orth)                                                  # the number of projections based on x1 y1
        h2 <- n2^(-2/9)                                                                  # bandwidth 
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_pls_num) # p-value matrix
        errormat2 <- U2%*%t(U2)                                                          #residual matrix based on x2 y2
        for(l in 1:pro1_pls_num){
          x_pro2 <- x2%*%pro1_orth[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2      # kernel function matrix
          #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                               # Gaussian kernel
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                    # Epanechnikov kernel 
          #test statistics
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
        }
        
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
        pval_pls_num <- pro1_pls_num + pro2_pls_num
      }
      
      # pval_cauchy_hybrid 
      pval_matrix_hybrid <- cbind(pval_matrix_Pcvm, pval_matrix_PLS)
      pval_num <- ncol(pval_matrix_hybrid)
      pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num])*pi)))
      
      # compare with local smoothing test (Tan et al. 2025)
      if(T){
        #construct local smoothing test statistic based on x1 y1 and projections based on x2 y2
        pro2_pls <- cbind(sdr_Upro2, sdr_ypro2)                                           # projections based on x2 y2
        pro2_pls_num <- ncol(pro2_pls)                                                    # the number of projections
        h1 <- n1^(-2/9)                                                                   # bandwidth 
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_pls_num) # p-value matrix
        errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
        for(q in 1:pro2_pls_num){
          x_pro1 <- x1%*%pro2_pls[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1       # kernel function matrix
          #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                                # Gaussian kernel
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                     # Epanechnikov kernel 
          #test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
        }
        
        #construct test statistics based on x2 y2 and projections based on x1 y1
        pro1_pls <- cbind(sdr_Upro1, sdr_ypro1)                                          # projections based on x1 y1
        pro1_pls_num <- ncol(pro1_pls)                                                   # the number of projections
        h2 <- n2^(-2/9)                                                                  # bandwidth 
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_pls_num) # p-value matrix
        errormat2 <- U2%*%t(U2)                                                          #residual matrix based on x2 y2
        for(l in 1:pro1_pls_num){
          x_pro2 <- x2%*%pro1_pls[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2      # kernel function matrix
          #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                               # Gaussian kernel
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                    # Epanechnikov kernel 
          #test statistics
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
        }
        
        pval_fisher1_PLS  <- 1- pchisq(-2*sum(log(pval_matrix1_PLS[1:pro2_pls_num])),df=2*pro2_pls_num)
        pval_fisher2_PLS  <- 1- pchisq(-2*sum(log(pval_matrix2_PLS[1:pro1_pls_num])),df=2*pro1_pls_num)
        pval_fisher_m <- cbind(pval_fisher1_PLS, pval_fisher2_PLS)
        pval_fisher_cauchy <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
      }
      
      result <- c(pval_cauchy1_Pcvm, pval_cauchy2_Pcvm, pval_cauchy_Pcvm,
                  pval_fisher1_PLS, pval_fisher2_PLS, pval_fisher_cauchy,
                  pval_cauchy_hybrid) 
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)
  
  # toc()
  #power
  
  martingale_cauchy1_power <- mean(hybrid_linear_power[,1]<=0.05,na.rm=TRUE)
  martingale_cauchy2_power <- mean(hybrid_linear_power[,2]<=0.05,na.rm=TRUE)
  martingale_cauchy_power <- mean(hybrid_linear_power[,3]<=0.05,na.rm=TRUE)
  pval_fisher1_power <- mean(hybrid_linear_power[,4]<=0.05,na.rm=TRUE)
  pval_fisher2_power <- mean(hybrid_linear_power[,5]<=0.05,na.rm=TRUE)
  pval_fisher_power <- mean(hybrid_linear_power[,6]<=0.05,na.rm=TRUE)
  
  hybrid_cauchy_power <- mean(hybrid_linear_power[,7]<=0.05,na.rm=TRUE) 
  
  return(list(martingale_cauchy1_power = martingale_cauchy1_power, 
              martingale_cauchy2_power = martingale_cauchy2_power, 
              martingale_cauchy_power = martingale_cauchy_power,
              pval_fisher1_power = pval_fisher1_power,
              pval_fisher2_power = pval_fisher2_power,
              pval_fisher_power = pval_fisher_power,
              hybrid_cauchy_power = hybrid_cauchy_power)
  )
}
Pcvm_Pls_hybrid_linear_test <- function(x,y){ 
  library(MASS)
  library(glmnet)
  library(LassoSIR)
  library(VariableScreening)
  library(psych) 
  
  n <- nrow(x)
  p <- ncol(x)
  # data splitting  
  index_x<-sample(1:n, floor(n/2),replace = FALSE)  
  index_x<-sort(index_x)
  
  x1 <- as.matrix(x[-index_x, , drop=FALSE]) 
  y1 <- as.matrix(y[-index_x]) 
  
  x2 <- as.matrix(x[index_x, , drop=FALSE])
  y2 <- as.matrix(y[index_x])
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  # residuals and beta_0 projections based on lasso and data x1 y1
  lasso_model1 <- cv.glmnet(x1, y1, family="gaussian", intercept = T)
  lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
  index_beta1_non0 <- which(is.finite(lasso_beta1) & lasso_beta1 != 0)
  len_beta1_non0 <- length(index_beta1_non0)
  if(len_beta1_non0 == 0){
    U1 <- y1 - mean(y1)                                                      # residual based on x1 y1
    beta1_hat <- lasso_beta1
    beta1_pro <- matrix(rep(1,p)/sqrt(p), ncol=1)
  }else{
    x1_sec <- x1[,index_beta1_non0, drop=FALSE]                              # second estimation
    sec_model1 <- glm(y1~x1_sec, family = gaussian)
    sec_beta1 <- unname(sec_model1$coefficients)[-1]
    lasso_beta1[index_beta1_non0] <- sec_beta1
    beta1_hat <- lasso_beta1
    pred1 = predict(sec_model1, newx = x1_sec, type="response")
    pred1 = matrix(unname(pred1), ncol = 1)
    U1 <- y1-pred1                                                           # residual based on x1 y1
    beta1_pro <- matrix(beta1_hat/sqrt(sum(beta1_hat^2)), ncol=1)
    
    # Collinearity protection for real data
    if(sum(is.na(sec_beta1))>0){ 
      library(caret)
      lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
      cor_matrix <- cor(x1_sec, use = "pairwise.complete.obs")
      high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
      if(length(high_cor_vars) > 0){
        index_beta1_non0 <- index_beta1_non0[-high_cor_vars]
        x1_sec <- x1[,index_beta1_non0, drop=FALSE ]  
      }
      sec_model1 <- glm(y1~x1_sec, family = gaussian)
      sec_beta1 <- unname(sec_model1$coefficients)[-1]
      sec_beta1[!is.finite(sec_beta1)] <- 0
      lasso_beta1[index_beta1_non0] <- sec_beta1
      beta1_hat <- lasso_beta1
      pred1 = predict(sec_model1, newx = x1_sec, type="response")
      pred1 = matrix(unname(pred1), ncol = 1)
      U1 <- y1-pred1                                                           
      beta1_pro <- matrix(beta1_hat/sqrt(sum(beta1_hat^2)), ncol=1)
    }
  } 
  
  # projections based on screening and lassosir for x1 U1
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
  
  # projections based on screening and lassosir for x1 y1
  if(p <= screen_num1){
    sir_y1 <- LassoSIR(x1, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sdr_ypro1 <- sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
  }else{
    rank_y1 <- screenIID(X=x1, Y=y1, method = "DC-SIS")                      
    index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        
    x1_yscreen <- x1[,index_y1, drop=FALSE]
    sir_y1 <- LassoSIR(x1_yscreen, y1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
    sir_ybeta1<- sir_y1$beta
    sdr_ypro1 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
    sdr_ypro1[index_y1,] <- sir_ybeta1
    sdr_ypro1 <- sdr_ypro1/sqrt(colSums(sdr_ypro1^2))
  }
  
  # residuals and beta0-projection based on lasso and data x2 y2
  lasso_model2 <- cv.glmnet(x2, y2, family="gaussian", intercept = T)
  lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
  index_beta2_non0 <- which(is.finite(lasso_beta2) & lasso_beta2 != 0)
  len_beta2_non0 <- length(index_beta2_non0)
  if(len_beta2_non0 == 0){
    U2 <- y2 - mean(y2)                                                      
    beta2_hat <- lasso_beta2
    beta2_pro <- matrix(rep(1,p)/sqrt(p), ncol=1)
  }else{
    x2_sec <- x2[,index_beta2_non0, drop=FALSE]                              
    sec_model2 <- glm(y2~x2_sec, family = gaussian)
    sec_beta2 <- unname(sec_model2$coefficients)[-1]
    sec_beta2[!is.finite(sec_beta2)] <- 0
    lasso_beta2[index_beta2_non0] <- sec_beta2
    beta2_hat <- lasso_beta2
    pred2 = predict(sec_model2, newx = x2_sec, type="response")
    pred2 = matrix(unname(pred2), ncol = 1)
    U2 <- y2-pred2                                                           
    beta2_pro <- matrix(beta2_hat/sqrt(sum(beta2_hat^2)), ncol=1)
    
    # Collinearity protection for real data
    if(sum(is.na(sec_beta2))>0){  
      library(caret)
      lasso_beta2 <- coef(lasso_model2, s = "lambda.min")[-1]
      cor_matrix <- cor(x2_sec, use = "pairwise.complete.obs")
      high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9, verbose = TRUE)
      if(length(high_cor_vars) > 0){
        index_beta2_non0 <- index_beta2_non0[-high_cor_vars]
        x2_sec <- x2[,index_beta2_non0, drop=FALSE ]  
      }
      sec_model2 <- glm(y2~x2_sec, family = gaussian)
      sec_beta2 <- unname(sec_model2$coefficients)[-1]
      sec_beta2[!is.finite(sec_beta2)] <- 0
      lasso_beta2[index_beta2_non0] <- sec_beta2
      beta2_hat <- lasso_beta2
      pred2 = predict(sec_model2, newx = x2_sec, type="response")
      pred2 = matrix(unname(pred2), ncol = 1)
      U2 <- y2-pred2                                                           
      beta2_pro <- matrix(beta2_hat/sqrt(sum(beta2_hat^2)), ncol=1)
    }
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
  
  # sir projections based on screening and lassosir for x2 y2
  if(p <= screen_num2){
    sir_y2 <- LassoSIR(x2, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
    sir_ybeta2<- sir_y2$beta
    sdr_ypro2 <- sir_ybeta2/sqrt(colSums(sir_ybeta2^2))                      
  }else{ 
    rank_y2 <- screenIID(X=x2, Y=y2, method = "DC-SIS")                      
    index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        
    x2_yscreen <- x2[,index_y2, drop=FALSE]
    sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
    sir_ybeta2 <- sir_ylasso2$beta
    sdr_ypro2 <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
    sdr_ypro2[index_y2,] <- sir_ybeta2
    sdr_ypro2 <- sdr_ypro2/sqrt(colSums(sdr_ypro2^2))                        
  } 
  
  # =========================================================================
  # Schmidt Orthogonalization 
  # =========================================================================
  pro1 <- cbind(beta1_pro, sdr_Upro1)
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
  
  pro2 <- cbind(beta2_pro, sdr_Upro2)
  pro2_orth <- matrix(numeric(0), nrow = nrow(pro2), ncol = 0)
  # tol <- 1e-12
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
  
  # ==============================================================
  # PCvM Test
  # ==============================================================
  #construct martingale-based test based on x1 y1 and projections based on x2 y2
  pro2 <- pro2_orth                                                            # projections based on x2 y2
  pro2_num <- ncol(pro2)                                                       # the number of projections
  pval_matrix1_Pcvm <- matrix(nrow=1, ncol=pro2_num)                           # p-value matrix of martingale-based test based on x1 y1
  
  for(q in 1:pro2_num){
    x1_pro <- x1%*%pro2[,q]                                                    # projected x1 with pro2 
    x1_Indictor<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
    
    hat_A1 <- rbind(t(x1_pro), rep(1,n1))  
    Gamma1_inv <- array(0, dim=c(2, 2, n1))
    for (i in 1:n1) {
      Gamma1_inv[,,i] = ginv((hat_A1%*%diag(x1_Indictor[i,]))%*%t(hat_A1))
    }                                                                          # The inverse of the matrix gamma
    Integral_1 <- hat_A1%*%diag(as.vector(U1))%*%t(x1_Indictor)                # The integral term in the martingale transformation 
    
    martingle1_sec<-diag(0,n1)                                                 # the second term in the martingale transformation
    for (l in 1:n1){
      martingle1_sec[l,]<- (t(hat_A1[,l]) %*% Gamma1_inv[,,l] %*% Integral_1[,l]) %*% x1_Indictor[l,] 
    }
    martingle_sta1<-(1/sqrt(n1))*t(U1) %*% x1_Indictor - (1/sqrt(n1)) * colSums(martingle1_sec)  # The test statistic of martingale transformation 
    
    ordx1_pro <- sort(x1_pro)                                                  # order x1_pro
    t_1 <- ordx1_pro[floor(0.99 * n1)]                                         # 99% quantile of x1_pro 
    sigma1_square <- mean(U1^2)   # estimation of variance of error
    F_1 <- (mean(x1_pro <= t_1))^2
    PcvM_martingle1 <- (1/(sigma1_square*F_1)) * mean((x1_pro <= t_1) * (t(martingle_sta1)^2))  #The CvM test statistic based on martingale transformation 
    pval_matrix1_Pcvm[,q] <- pvalue_integ_Brown(PcvM_martingle1)
  }
  
  #construct martingale-based test based on x2 y2 and projections based on x1 y1                                       
  pro1 <- pro1_orth                                                            # projections based on x1 y1
  pro1_num <- ncol(pro1)                                                       # the number of projections
  pval_matrix2_Pcvm <- matrix(nrow=1, ncol=pro1_num)                           # p-value matrix based on x2 y2
  
  for(qq in 1:pro1_num){
    x2_pro <- x2%*%pro1[,qq]                                                   # projected x1 with pro2 
    x2_Indictor<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
    
    hat_A2 <- rbind(t(x2_pro), rep(1,n2))  
    Gamma2_inv <- array(0, dim=c(2, 2, n2))
    for (ii in 1:n2) {
      Gamma2_inv[,,ii] = ginv((hat_A2%*%diag(x2_Indictor[ii,]))%*%t(hat_A2))
    }                                                                          # The inverse of the matrix gamma
    Integral_2 <- hat_A2%*%diag(as.vector(U2))%*%t(x2_Indictor)                # The integral term in the martingale transformation 
    
    martingle2_sec<-diag(0,n2)                                                 # The second term in the martingale transformation
    for (ll in 1:n2){
      martingle2_sec[ll,]<- (t(hat_A2[,ll]) %*% Gamma2_inv[,,ll] %*% Integral_2[,ll]) %*% x2_Indictor[ll,] 
    }
    martingle_sta2<-(1/sqrt(n2))*t(U2) %*% x2_Indictor - (1/sqrt(n2)) * colSums(martingle2_sec)  # The test statistic of martingale transformation
    
    ordx2_pro <- sort(x2_pro)                                                  # order x2_pro
    t_2 <- ordx2_pro[floor(0.99 * n2)]                                         # 99% quantile of x2_pro 
    sigma2_square <- mean(U2^2)                                                # estimation of variance of error
    F_2 <- (mean(x2_pro <= t_2))^2
    PcvM_martingle2 <- (1/(sigma2_square*F_2)) * mean((x2_pro <= t_2) * (t(martingle_sta2)^2))  #The CvM test statistic based on martingale transformation 
    pval_matrix2_Pcvm[,qq] <- pvalue_integ_Brown(PcvM_martingle2)
  }
  pval_matrix_Pcvm <- cbind(pval_matrix1_Pcvm, pval_matrix2_Pcvm)
  pval_num <- pro1_num + pro2_num
  
  pval_cauchy1_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_Pcvm[1:pro2_num])*pi)))  # cauchy combination of martingale-based test based on x1 y1
  pval_cauchy2_Pcvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_Pcvm[1:pro1_num])*pi)))  # cauchy combination of martingale-based test based on x1 y1
  pval_cauchy_Pcvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_Pcvm[1:pval_num])*pi)))   # cauchy combination of martingale-based test based on x y
  
  # ==============================================================
  # Hybrid and PLS
  # ==============================================================
  if(T){
    # construct hybrid test statistic with martingale test statistic and local smoothing test based on x1 y1 and projections based on x2 y2
    pro2_pls_num <- ncol(pro2_orth)                                               # the number of projections based on x2 y2
    h1 <- n1^(-2/9)                                                               # bandwidth 
    
    pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_pls_num) # p-value matrix
    errormat1 <- U1%*%t(U1)                                                       #residual matrix based on x1 y1
    for(q in 1:pro2_pls_num){
      x_pro1 <- x1%*%pro2_orth[,q]
      x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
      indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
      kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel 
      #test statistics
      Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
      pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
    }
    
    #construct test statistics based on x2 y2 and projections based on x1 y1                                       
    pro1_pls_num <- ncol(pro1_orth)                                                # the number of projections based on x1 y1
    h2 <- n2^(-2/9)                                                                # bandwidth 
    
    pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_pls_num) # p-value matrix
    errormat2 <- U2%*%t(U2)                                                        #residual matrix based on x2 y2
    for(l in 1:pro1_pls_num){
      x_pro2 <- x2%*%pro1_orth[,l]
      x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2      # kernel function matrix
      indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
      kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                    # Epanechnikov kernel 
      #test statistics
      Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
      pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
    }
    
    pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
    pval_pls_num <- pro1_pls_num + pro2_pls_num
    
    pval_cauchy1_PLS <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PLS[1:pro2_pls_num])*pi)))  # cauchy combination based on x1 y1
    pval_cauchy2_PLS <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PLS[1:pro1_pls_num])*pi)))  # cauchy combination based on x1 y1
    pval_cauchy_PLS  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PLS[1:pval_pls_num])*pi)))   # cauchy combination based on x y
    
    pval_min1_PLS <- 1 - (1 - min(pval_matrix1_PLS))^pro2_pls_num
    pval_min2_PLS <- 1 - (1 - min(pval_matrix2_PLS))^pro1_pls_num
    pval_min_m <- cbind(pval_min1_PLS, pval_min2_PLS) 
    pval_min_cauchy_PLS <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi))) 
  }
  
  # pval_cauchy_hybrid 
  pval_matrix_hybrid <- cbind(pval_matrix_Pcvm, pval_matrix_PLS)
  pval_num_h <- ncol(pval_matrix_hybrid)
  pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
  
  # ==============================================================
  # Standalone PLS Test (Tan et al. 2025)
  # ==============================================================
  if(T){
    #construct local smoothing test statistic based on x1 y1 and projections based on x2 y2
    pro2_pls <- cbind(sdr_Upro2, sdr_ypro2)                                           # projections based on x2 y2
    pro2_pls_num <- ncol(pro2_pls)                                                    # the number of projections
    h1 <- n1^(-2/9)                                                                   # bandwidth 
    
    pval_matrix1_PLS_std <- matrix(nrow=1, ncol=pro2_pls_num) # p-value matrix
    errormat1 <- U1%*%t(U1)  #residual matrix based on x1 y1
    for(q in 1:pro2_pls_num){
      x_pro1 <- x1%*%pro2_pls[,q]
      x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1       # kernel function matrix
      indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0) 
      kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                     # Epanechnikov kernel 
      #test statistics
      Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
      pval_matrix1_PLS_std[,q] <- 1-pnorm(Tn1)
    }
    
    #construct test statistics based on x2 y2 and projections based on x1 y1
    pro1_pls <- cbind(sdr_Upro1, sdr_ypro1)                                          # projections based on x1 y1
    pro1_pls_num <- ncol(pro1_pls)                                                   # the number of projections
    h2 <- n2^(-2/9)                                                                  # bandwidth 
    
    pval_matrix2_PLS_std <- matrix(nrow=1, ncol=pro1_pls_num) # p-value matrix
    errormat2 <- U2%*%t(U2)                                                          #residual matrix based on x2 y2
    for(l in 1:pro1_pls_num){
      x_pro2 <- x2%*%pro1_pls[,l]
      x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2      # kernel function matrix
      indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0) 
      kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                    # Epanechnikov kernel 
      #test statistics
      Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
      pval_matrix2_PLS_std[,l] <- 1-pnorm(Tn2)
    }
    
    pval_fisher1_PLS  <- 1- pchisq(-2*sum(log(pval_matrix1_PLS_std[1:pro2_pls_num])),df=2*pro2_pls_num)
    pval_fisher2_PLS  <- 1- pchisq(-2*sum(log(pval_matrix2_PLS_std[1:pro1_pls_num])),df=2*pro1_pls_num)
    pval_fisher_m <- cbind(pval_fisher1_PLS, pval_fisher2_PLS)
    pval_fisher_cauchy_PLS <- 1- pcauchy(mean(tan((0.5-pval_fisher_m)*pi)))  
  }
  
  result <- list(pval_cauchy1_Pcvm = pval_cauchy1_Pcvm, 
                 pval_cauchy2_Pcvm = pval_cauchy2_Pcvm, 
                 pval_cauchy_Pcvm = pval_cauchy_Pcvm,
                 pval_cauchy1_PLS = pval_cauchy1_PLS, 
                 pval_cauchy2_PLS = pval_cauchy2_PLS, 
                 pval_cauchy_PLS = pval_cauchy_PLS,
                 pval_min1_PLS = pval_min1_PLS, 
                 pval_min2_PLS = pval_min2_PLS, 
                 pval_min_cauchy_PLS = pval_min_cauchy_PLS,
                 pval_fisher1_PLS = pval_fisher1_PLS, 
                 pval_fisher2_PLS = pval_fisher2_PLS, 
                 pval_fisher_cauchy_PLS = pval_fisher_cauchy_PLS,
                 pval_cauchy_hybrid = pval_cauchy_hybrid) 
  
  return(result) 
}
# logit + test
Pcvm_Pls_hybrid_logit <- function(p,n,a,pho){
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
  #library(tictoc)    #library packages
  
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
      # z <- x %*% beta0 + a * (x[,1]*x[,2] + x[,2]*x[,3] + x[,3]*x[,4] + x[,4]*x[,5])        # H22
      # z <- x %*% beta0 + a * 2 * cos(0.6 * x %*% beta0 * pi)                                # H23 
      
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
          rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                     # screening method based on distance correlation
          index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                       # chosen index of X
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
          rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                     # screening method based on distance correlation
          index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                       # chosen index of X
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
          rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                     # screening method based on distance correlation
          index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                       # chosen index of X
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
          rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                     # screening method based on distance correlation
          index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                       # chosen index of X
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
      pro2_num <- ncol(sir_Upro2_orth)                                                                  # the number of SIR-U projections
      pval_matrix1_PCvm <- matrix(nrow=1, ncol=pro2_num+1)                                              # p-value matrix
      
      #deri_link1 <- exp(-x1 %*% beta1_hat - intercept1)/(1+exp(-x1 %*% beta1_hat - intercept1))^2      # The derivative of the link funciont \mu(x)
      deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
      
      for(q in seq_len(pro2_num)){
        x1_pro <- x1%*%sir_Upro2_orth[,q]
        band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
        h1 <- bandwidth_choice(x1_pro, band_y1)                                            # choose the bandwidth
        
        Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                  # kernel function matrix
        #kernel_1 <- exp(-0.5*Ker_inter1^2)                                                # Gaussian kernel
        Ker_indictor1 <- ifelse(abs(Ker_inter1) <= 1, 1, 0)
        kernel_1 <- (3/4)*(1-Ker_inter1^2)*Ker_indictor1
        
        Indictor_1<-ifelse(x1_pro %*% rep(1,n1) <= rep(1,n1)%*%t(x1_pro),1,0)
        
        A1_1 <- ((kernel_1%*%deri_link1)*x1_pro)/(kernel_1%*%U1^2+n1^(-12))                          # the 1st compoenent of A1
        A1_2 <- kernel_1%*%((x1 %*% beta1_hat)*deri_link1)/(kernel_1%*%U1^2+n1^(-12))                # the 2st compoenent of A1
        A1_3 <- (kernel_1%*%deri_link1)/(kernel_1%*%U1^2+n1^(-12))                                   # the 3st compoenent of A1
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
        t_1 <- ordx1_pro[floor(0.99 * n1)]
        psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
        PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        # the CvM test statistic based on martingale transformation
        pval_matrix1_PCvm[,q] <- pvalue_integ_Brown(PCvM_martingle1)
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
      integral_beta_1 <- ((rep(1,2)%*%t(U1))*hat_beta_A1)%*%t(beta_Indictor1)                     # integral part in the martingale transformation
      
      martingle_beta_sec1<-diag(0,n1)                                                             # The second term in the martingale transformation
      for (l in 1:n1){
        martingle_beta_sec1[l,]<- (U1[l]^2 * t(hat_beta_A1[,l])) %*% Gamma_beta_inv1[,,l] %*% integral_beta_1[,l] %*% beta_Indictor1[l,]
      }
      martingle_beta_sta1<-(1/sqrt(n1))*t(U1) %*% beta_Indictor1 - (1/sqrt(n1)) * colSums(martingle_beta_sec1)      # the test statistic of the martingale transformation
      
      ordx1_betapro <- sort(x1_betapro)  # order x%*%projection
      beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
      psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
      PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  # the CvM test statistic based on martingale transformation
      pval_matrix1_PCvm[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1)
      
      
      #construct martingale-based test statistics based on x2 y2 and projections sir_Upro1_orth based on x1 y1
      pro1_num <- ncol(sir_Upro1_orth)                                                                     # the number of projections
      pval_matrix2_PCvm <- matrix(nrow=1, ncol=pro1_num+1)                                                 # p-value matrix
      
      #deri_link2 <- exp(-x2 %*% beta2_hat - intercept2)/(1+exp(-x2 %*% beta2_hat - intercept2))^2         # The derivative of the link funciont \mu(x)
      deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
      
      for(qq in seq_len(pro1_num)){
        x2_pro <- x2%*%sir_Upro1_orth[,qq]
        band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
        h2 <- bandwidth_choice(x2_pro, band_y2)                                                    # choose the bandwidth
        
        Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                          # kernel function matrix
        #kernel_2 <- exp(-0.5*Ker_inter2^2)                                                        # Gaussian kernel
        
        Ker_indictor2 <- ifelse(abs(Ker_inter2) <= 1, 1, 0)
        kernel_2 <- (3/4)*(1-Ker_inter2^2)*Ker_indictor2
        
        Indictor_2<-ifelse(x2_pro %*% rep(1,n2) <= rep(1,n2)%*%t(x2_pro),1,0)
        
        A2_1 <- ((kernel_2%*%deri_link2)*x2_pro)/(kernel_2%*%U2^2+n2^(-12))                           # the 1st compoenent of A1
        A2_2 <- kernel_2%*%((x2 %*% beta2_hat)*deri_link2)/(kernel_2%*%U2^2+n2^(-12))                 # the 2st compoenent of A1
        A2_3 <- (kernel_2%*%deri_link2)/(kernel_2%*%U2^2+n2^(-12))                                    # the 3st compoenent of A1
        hat_A2 <- rbind(t(A2_1), t(A2_2), t(A2_3))
        
        #gamma matrix
        Gamma_inv2 <- array(0, dim=c(3, 3, n2))
        for (ii in 1:n2){
          Gamma_inv2[,,ii] = ginv(((rep(1,3)%*%t(U2)^2)*hat_A2)%*%(t(hat_A2) * (as.matrix(Indictor_2[ii,]) %*%rep(1,3))))
        }
        integral_2 <- ((rep(1,3)%*%t(U2))*hat_A2)%*%t(Indictor_2)
        
        martingle_sec2<-diag(0,n2)                                                                  # the second term in the martingale transformation
        for (ll in 1:n2){
          martingle_sec2[ll,]<- (U2[ll]^2 * t(hat_A2[,ll])) %*% Gamma_inv2[,,ll] %*% integral_2[,ll] %*% Indictor_2[ll,]
        }
        martingle_sta2<-(1/sqrt(n2))*t(U2) %*% Indictor_2 - (1/sqrt(n2)) * colSums(martingle_sec2)  # the test statistic of the martingale transformation
        
        ordx2_pro <- sort(x2_pro)                                                                   # order x%*%projection
        t_2 <- ordx2_pro[floor(0.99 * n2)]
        psi_2 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2))
        PCvM_martingle2 <- (1/psi_2^2) * mean(U2^2 * (x2_pro <= t_2) * (t(martingle_sta2)^2))       # The CvM test statistic based on martingale transformation
        pval_matrix2_PCvm[,qq] <- pvalue_integ_Brown(PCvM_martingle2)
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
      
      martingle_beta_sec2<-diag(0,n2)                                                    # The second term in the martingale transformation
      for (ll in 1:n2){
        martingle_beta_sec2[ll,]<- (U2[ll]^2 * t(hat_beta_A2[,ll])) %*% Gamma_beta_inv2[,,ll] %*% integral_beta_2[,ll] %*% beta_Indictor2[ll,]
      }
      martingle_beta_sta2<-(1/sqrt(n2))*t(U2) %*% beta_Indictor2 - (1/sqrt(n2)) * colSums(martingle_beta_sec2)   # The test statistic of the martingale transformation
      
      ordx2_betapro <- sort(x2_betapro)                                                  # order x%*%proj
      beta_t2 <- ordx2_betapro[floor(0.99 * n2)]
      psi_beta_2 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2))
      PCvM_beta_martingle2 <- (1/psi_beta_2^2) * mean(U2^2 * (x2_betapro <= beta_t2) * (t(martingle_beta_sta2)^2))   # The CvM test statistic based on martingale transformation
      pval_matrix2_PCvm[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2)
      
      pval_matrix_PCvm <- cbind(pval_matrix1_PCvm, pval_matrix2_PCvm)
      pval_num <- ncol(pval_matrix_PCvm)
      pval_cauchy1_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm[1:(pro2_num+1)])*pi)))          # cauchy combination based on x1 y1
      pval_cauchy2_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm[1:(pro1_num+1)])*pi)))          # cauchy combination based on x2 y2
      pval_cauchy_PCvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm[1:pval_num])*pi)))               # cauchy combination based on x y
      
      # construct hybrid test statistic
      if(T){
        #construct hybrid test based Martingale test statistics and PLS test statistics based on x1 y1 and projections based on x2 y2
        pro2_num_pls <- ncol(pro2_orth)                                                        # the number of projections based on x2 y2
        h1 <- n1^(-2/9)                                                                        # bandwidth
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) # p-value matrix
        errormat1 <- U1%*%t(U1)                                                                # residual matrix based on x1 y1
        for(q in 1:pro2_num_pls){
          x_pro1 <- x1%*%pro2_orth[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
          #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                            # Gaussian kernel
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel
          #PLS test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval1 <- 1-pnorm(Tn1)
          pval_matrix1_PLS[,q] <- pval1
        }
        
        #construct hybrid test based Martingale test statistics and PLS test statistics based on x2 y2 and projections based on x1 y1
        pro1_num_pls <- ncol(pro1_orth)                                                    # the number of projections based on x1 y1
        h2 <- n2^(-2/9)                                                                    # bandwidth
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) # p-value matrix
        errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
        for(l in 1:pro1_num_pls){
          x_pro2 <- x2%*%pro1_orth[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        # kernel function matrix
          #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                                 # Gaussian kernel
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      # Epanechnikov kernel
          # PLS test statistics
          Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
          pval2 <- 1-pnorm(Tn2)
          pval_matrix2_PLS[,l] <- pval2
        }
        pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
        pval_num_pls <- pro1_num_pls + pro2_num_pls
      }
      
      # pval_cauchy_hybrid
      pval_matrix_hybrid <- cbind(pval_matrix_PCvm, pval_matrix_PLS)
      pval_num_h <- ncol(pval_matrix_hybrid)
      pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
      
      # construct PLS test statistic (Tan et al. 2025) for comparison
      if(T){
        #construct PLS test statistic based on x1 y1 and projections based on x2 y2
        pro2_pls <- cbind(sir_Upro2, sir_ypro2)                                               # projections based on x2 y2
        pro2_num_pls <- ncol(pro2_pls)                                                        # the number of projections
        h1 <- n1^(-2/9)                                                                       # bandwidth
        
        pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_num_pls) # p-value matrix
        errormat1 <- U1%*%t(U1)                                                               # residual matrix based on x1 y1
        for(q in 1:pro2_num_pls){
          x_pro1 <- x1%*%pro2_pls[,q]
          x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1           # kernel function matrix
          #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                                    # Gaussian kernel
          indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
          kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                         # Epanechnikov kernel
          #PLS test statistics
          Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
          pval1 <- 1-pnorm(Tn1)
          pval_matrix1_PLS[,q] <- pval1
        }
        
        #construct PLS test statistics based on x2 y2 and projections based on x1 y1
        pro1_pls <- cbind(sir_Upro1, sir_ypro1)                                           # projections based on x1 y1
        pro1_num_pls <- ncol(pro1_pls)                                                    # the number of projections
        h2 <- n2^(-2/9)                                                                   # bandwidth
        
        pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_num_pls) # p-value matrix
        errormat2 <- U2%*%t(U2)                                                           #residual matrix based on x2 y2
        for(l in 1:pro1_num_pls){
          x_pro2 <- x2%*%pro1_pls[,l]
          x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2       # kernel function matrix
          #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                                # Gaussian kernel
          indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
          kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                     # Epanechnikov kernel
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
      
      result <- c(pval_cauchy1_PCvm, pval_cauchy2_PCvm, pval_cauchy_PCvm,
                  pval_fisher1_PLS, pval_fisher2_PLS, pval_fisher_cauchy,
                  pval_cauchy_hybrid)
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)
  
  # toc()
  #power
  martingale_cauchy1_power <- mean(hybrid_power[,1] <= 0.05,na.rm=TRUE)
  martingale_cauchy2_power <- mean(hybrid_power[,2] <= 0.05,na.rm=TRUE)
  martingale_cauchy_power <- mean(hybrid_power[,3] <= 0.05,na.rm=TRUE)
  
  pval_fisher1_power <- mean(hybrid_power[,4] <= 0.05,na.rm=TRUE)
  pval_fisher2_power <- mean(hybrid_power[,5] <= 0.05,na.rm=TRUE)
  pval_fisher_cauchy_power <- mean(hybrid_power[,6] <= 0.05,na.rm=TRUE)
  
  hybrid_cauchy_power <- mean(hybrid_power[,7] <= 0.05,na.rm=TRUE)
  
  #return result
  return(list(martingale_cauchy1_power = martingale_cauchy1_power,
              martingale_cauchy2_power = martingale_cauchy2_power,
              martingale_cauchy_power = martingale_cauchy_power,
              pval_fisher1_power = pval_fisher1_power,
              pval_fisher2_power = pval_fisher2_power,
              pval_fisher_cauchy_power = pval_fisher_cauchy_power,
              hybrid_cauchy_power = hybrid_cauchy_power)
  )
}
Pcvm_Pls_hybrid_logit_test <- function(x,y){  
  library(MASS)
  library(glmnet)
  library(LassoSIR) 
  library(harmonicmeanp) 
  library(VariableScreening)
  library(psych) 
  
  p <- ncol(x)
  n <- nrow(x)
  #data splitting  
  index_x<-sample(1:n,floor(n/2),replace = FALSE)  
  index_x<-sort(index_x)
  
  x1 <- as.matrix(x[-index_x, , drop=FALSE])
  x2 <- as.matrix(x[index_x, , drop=FALSE])
  
  y1 <- as.matrix(y[-index_x])
  y2 <- as.matrix(y[index_x]) 
  
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  
  # =========================================================
  # Estimation: x1, y1 
  # =========================================================
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
    x1_sec <- x1[,index_beta1_non0, drop=FALSE]                                          # second estimation
    sec_model1 <- glm(y1~x1_sec, family = binomial(link = "logit"))
    beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
    beta1_hat <- beta1_lasso
    intercept1 <- unname(sec_model1$coefficients)[1] 
    pred1 = predict(sec_model1, newx = x1_sec, type="response")
    pred1 = matrix(unname(pred1), ncol = 1) 
    U1 <- y1 - pred1                                                         # residuals based on x1 y1
    # beta1_pro <- beta1_hat
    beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
    
    # Elastic Net alpha=1 
    if(sum(is.na(beta1_hat))>0){ 
      lasso_model1_1 <- cv.glmnet(x1_sec, y1, family="binomial", alpha=1) 
      beta1_lasso[index_beta1_non0] <- coef(lasso_model1_1, s = "lambda.min")[-1]  
      beta1_hat <- beta1_lasso
      intercept1 <- coef(lasso_model1_1, s = "lambda.min")[1]
      pred1 <- predict(lasso_model1_1, newx = x1_sec, s = "lambda.min", type = "response")
      U1 <- y1 - as.vector(pred1) 
      # beta1_pro <- beta1_hat
      beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
    }
  }
  
  # sir x1 U1
  screen_num1 <- floor(n1/log(n1))
  sir_Upro1 <- tryCatch({
    if(p <= screen_num1){
      sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic",
                         solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
      sir_Ubeta1<- sir_U1$beta
      sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
    }else{
      rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      # screening method based on distance correlation
      index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        # chosen index of X
      x1_Uscreen <- x1[,index_U1, drop=FALSE]
      sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic",
                              solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
      sir_Ubeta1<- sir_Ulasso1$beta
      sir_Upro1_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta1))
      sir_Upro1_t[index_U1,] <- sir_Ubeta1
      sir_Upro1_t/sqrt(colSums(sir_Upro1_t^2))
    }}, error =function(e){ matrix(NA,p,1) })
  
  # sir x1 y1
  sir_ypro1 <- tryCatch({
    if(p <= screen_num1){
      sir_ylasso1 <- LassoSIR(x1, y1, H=2, choosing.d="automatic",
                              solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
      sir_ybeta1<- sir_ylasso1$beta
      sir_ybeta1/sqrt(colSums(sir_ybeta1^2))
    }else{
      rank_y1 <- screenIID(X=x1, Y=y1, method = "MV-SIS")                      # screening method based on distance correlation
      index_y1 <- seq(1:p)[rank_y1$rank <= screen_num1]                        # chosen index of X
      x1_yscreen <- x1[,index_y1, drop=FALSE]
      sir_ylasso1 <- LassoSIR(x1_yscreen, y1, H=2, choosing.d="automatic",
                              solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
      sir_ybeta1<- sir_ylasso1$beta
      sir_ypro1_t <- matrix(0,nrow = p,ncol = ncol(sir_ybeta1))
      sir_ypro1_t[index_y1,] <- sir_ybeta1
      sir_ypro1_t/sqrt(colSums(sir_ypro1_t^2))
    }
  }, error =function(e){ matrix(NA,p,1) })
  
  # =========================================================
  # Estimation: x2, y2 
  # =========================================================
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
    x2_sec <- x2[,index_beta2_non0, drop=FALSE]
    sec_model2 <- glm(y2~x2_sec,family = binomial(link = "logit"))
    beta2_lasso[index_beta2_non0] <- unname(sec_model2$coefficients)[-1]
    beta2_hat <- beta2_lasso
    intercept2 <- unname(sec_model2$coefficients)[1]
    pred2 = predict(sec_model2, newx = x2_sec,type="response")
    pred2 = matrix(unname(pred2), ncol = 1)
    U2 <- y2-pred2
    # beta2_pro <- beta2_hat
    beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
    
    if(sum(is.na(beta2_hat))>0){ 
      lasso_model1_2 <- cv.glmnet(x2_sec, y2, family="binomial", alpha=1) 
      beta2_lasso[index_beta2_non0] <- coef(lasso_model1_2, s = "lambda.min")[-1]  
      beta2_hat <- beta2_lasso
      intercept2 <- coef(lasso_model1_2, s = "lambda.min")[1]
      pred2 <- predict(lasso_model1_2, newx = x2_sec, s = "lambda.min", type = "response")
      U2 <- y2 - as.vector(pred2) 
      # beta2_pro <- beta2_hat
      beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
    }
  }
  
  # sir x2 U2
  screen_num2 <- floor(n2/log(n2))
  sir_Upro2 <- tryCatch({
    if(p <= screen_num2){
      sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic",
                              solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
      sir_Ubeta2<- sir_Ulasso2$beta
      sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
    }else{
      rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      # screening method based on distance correlation
      index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        # chosen index of X
      x2_Uscreen <- x2[,index_U2, drop=FALSE]
      sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic",
                              solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
      sir_Ubeta2 <- sir_Ulasso2$beta
      sir_Upro2_t <- matrix(0,nrow = p,ncol = ncol(sir_Ubeta2))
      sir_Upro2_t[index_U2,] <- sir_Ubeta2
      sir_Upro2_t/sqrt(colSums(sir_Upro2_t^2))
    }}, error =function(e){ matrix(NA,p,1) })
  
  # sir x2 y2
  sir_ypro2 <- tryCatch({
    if(p <= screen_num2){
      sir_ylasso2 <- LassoSIR(x2, y2, H=2, choosing.d="automatic",
                              solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
      sir_ybeta2<- sir_ylasso2$beta
      sir_ybeta2/sqrt(colSums(sir_ybeta2^2))
    }else{
      rank_y2 <- screenIID(X=x2, Y=y2, method = "MV-SIS")                      # screening method based on distance correlation
      index_y2 <- seq(1:p)[rank_y2$rank <= screen_num2]                        # chosen index of X
      x2_yscreen <- x2[,index_y2, drop=FALSE]
      sir_ylasso2 <- LassoSIR(x2_yscreen, y2, H=2, choosing.d="automatic",
                              solution.path=FALSE, categorical=TRUE, nfolds=5, screening=TRUE)
      sir_ybeta2 <- sir_ylasso2$beta
      sir_ypro2_t <- matrix(0,nrow = p,ncol = ncol(sir_ybeta2))
      sir_ypro2_t[index_y2,] <- sir_ybeta2
      sir_ypro2_t/sqrt(colSums(sir_ypro2_t^2))
    }}, error =function(e){matrix(NA,p,1) })
  
  # =========================================================
  # Schmidt Orthogonalization 
  # =========================================================
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
  
  # ==============================================================
  # 1. PCvM Test 
  # ==============================================================
  pro2_num <- ncol(sir_Upro2_orth)                                                                                # the number of SIR-U projections
  pval_matrix1_PCvm <- matrix(nrow=1, ncol=pro2_num+1)                                                            # p-value matrix
  
  #deri_link1 <- exp(-x1 %*% beta1_hat - intercept1)/(1+exp(-x1 %*% beta1_hat - intercept1))^2      # The derivative of the link funciont \mu(x)
  deri_link1 <- plogis(x1 %*% beta1_hat + intercept1) * (1 - plogis(x1 %*% beta1_hat + intercept1 ))
  
  for(q in seq_len(pro2_num)){
    x1_pro <- x1%*%sir_Upro2_orth[,q]
    band_y1 <- cbind(U1^2, deri_link1, (x1 %*% beta1_hat)*deri_link1)
    h1 <- bandwidth_choice(x1_pro, band_y1)                                                            # choose the bandwidth
    
    Ker_inter1 <- (x1_pro %*% rep(1,n1) - rep(1,n1) %*% t(x1_pro))/h1                                  # kernel function matrix
    #kernel_1 <- exp(-0.5*Ker_inter1^2)                                                                # Gaussian kernel
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
    t_1 <- ordx1_pro[floor(0.99 * n1)]
    psi_1 <- (1/n1) * sum(U1^2 * (x1_pro <= t_1))
    PCvM_martingle1 <- (1/psi_1^2) * mean(U1^2 * (x1_pro <= t_1) * (t(martingle_sta1)^2))        # the CvM test statistic based on martingale transformation
    pval_matrix1_PCvm[,q] <- pvalue_integ_Brown(PCvM_martingle1)
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
  beta_t1 <- ordx1_betapro[floor(0.99 * n1)]
  psi_beta_1 <- (1/n1) * sum(U1^2 * (x1_betapro <= beta_t1))
  PCvM_beta_martingle1 <- (1/psi_beta_1^2) * mean(U1^2 * (x1_betapro <= beta_t1) * (t(martingle_beta_sta1)^2))  # the CvM test statistic based on martingale transformation
  pval_matrix1_PCvm[pro2_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle1)
  
  
  #construct martingale-based test statistics based on x2 y2 and projections sir_Upro1_orth based on x1 y1
  pro1_num <- ncol(sir_Upro1_orth)                                                                             # the number of projections
  pval_matrix2_PCvm <- matrix(nrow=1, ncol=pro1_num+1)                                                         # p-value matrix
  
  #deri_link2 <- exp(-x2 %*% beta2_hat - intercept2)/(1+exp(-x2 %*% beta2_hat - intercept2))^2         # The derivative of the link funciont \mu(x)
  deri_link2 <-  plogis(x2 %*% beta2_hat + intercept2 ) * (1 - plogis(x2 %*% beta2_hat + intercept2 ))
  
  for(qq in seq_len(pro1_num)){
    x2_pro <- x2%*%sir_Upro1_orth[,qq]
    band_y2 <- cbind(U2^2, deri_link2, (x2 %*% beta2_hat)*deri_link2)
    h2 <- bandwidth_choice(x2_pro, band_y2)                                                            # choose the bandwidth
    
    Ker_inter2 <- (x2_pro %*% rep(1,n2) - rep(1,n2) %*% t(x2_pro))/h2                                  # kernel function matrix
    #kernel_2 <- exp(-0.5*Ker_inter2^2)                                                                # Gaussian kernel
    
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
    t_2 <- ordx2_pro[floor(0.99 * n2)]
    psi_2 <- (1/n2) * sum(U2^2 * (x2_pro <= t_2))
    PCvM_martingle2 <- (1/psi_2^2) * mean(U2^2 * (x2_pro <= t_2) * (t(martingle_sta2)^2))       # The CvM test statistic based on martingale transformation
    pval_matrix2_PCvm[,qq] <- pvalue_integ_Brown(PCvM_martingle2)
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
  beta_t2 <- ordx2_betapro[floor(0.99 * n2)]
  psi_beta_2 <- (1/n2) * sum(U2^2 * (x2_betapro <= beta_t2))
  PCvM_beta_martingle2 <- (1/psi_beta_2^2) * mean(U2^2 * (x2_betapro <= beta_t2) * (t(martingle_beta_sta2)^2))   # The CvM test statistic based on martingale transformation
  pval_matrix2_PCvm[pro1_num+1] <- pvalue_integ_Brown(PCvM_beta_martingle2)
  
  pval_matrix_PCvm <- cbind(pval_matrix1_PCvm, pval_matrix2_PCvm)
  pval_num <- ncol(pval_matrix_PCvm)
  pval_cauchy1_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PCvm[1:(pro2_num+1)])*pi)))          # cauchy combination based on x1 y1
  pval_cauchy2_PCvm <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PCvm[1:(pro1_num+1)])*pi)))          # cauchy combination based on x2 y2
  pval_cauchy_PCvm  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PCvm[1:pval_num])*pi)))               # cauchy combination based on x y
  
  # ==============================================================
  # Hybrid and PLS
  # ==============================================================
  if(T){
    #construct hybrid test based Martingale test statistics and PLS test statistics based on x1 y1 and projections based on x2 y2
    pro2_pls_num <- ncol(pro2_orth)                                                              # the number of projections based on x2 y2
    h1 <- n1^(-2/9)                                                                              # bandwidth
    
    pval_matrix1_PLS <- matrix(nrow=1, ncol=pro2_pls_num) # p-value matrix
    errormat1 <- U1%*%t(U1)                                                                      # residual matrix based on x1 y1
    for(q in 1:pro2_pls_num){
      x_pro1 <- x1%*%pro2_orth[,q]
      x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1   # kernel function matrix
      #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                             # Gaussian kernel
      indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
      kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                 # Epanechnikov kernel
      #PLS test statistics
      Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
      pval_matrix1_PLS[,q] <- 1-pnorm(Tn1)
    }
    
    #construct hybrid test based Martingale test statistics and PLS test statistics based on x2 y2 and projections based on x1 y1
    pro1_pls_num <- ncol(pro1_orth)                                                              # the number of projections based on x1 y1
    h2 <- n2^(-2/9)                                                                              # bandwidth
    
    pval_matrix2_PLS <- matrix(nrow=1, ncol=pro1_pls_num) # p-value matrix
    errormat2 <- U2%*%t(U2)  #residual matrix based on x2 y2
    for(l in 1:pro1_pls_num){
      x_pro2 <- x2%*%pro1_orth[,l]
      x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        # kernel function matrix
      #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                                 # Gaussian kernel
      indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
      kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      # Epanechnikov kernel
      # PLS test statistics
      Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
      pval_matrix2_PLS[,l] <- 1-pnorm(Tn2)
    }
    pval_matrix_PLS <- cbind(pval_matrix1_PLS, pval_matrix2_PLS)
    pval_pls_num <- pro1_pls_num + pro2_pls_num
    
    pval_cauchy1_PLS <- 1- pcauchy(mean(tan((0.5-pval_matrix1_PLS[1:pro2_pls_num])*pi)))             
    pval_cauchy2_PLS <- 1- pcauchy(mean(tan((0.5-pval_matrix2_PLS[1:pro1_pls_num])*pi)))             
    pval_cauchy_PLS  <- 1- pcauchy(mean(tan((0.5-pval_matrix_PLS[1:pval_pls_num])*pi)))              
    
    pval_min1_PLS <- 1 - (1 - min(pval_matrix1_PLS))^pro2_pls_num
    pval_min2_PLS <- 1 - (1 - min(pval_matrix2_PLS))^pro1_pls_num
    pval_min_m <- cbind(pval_min1_PLS, pval_min2_PLS) 
    pval_min_cauchy_PLS <- 1- pcauchy(mean(tan((0.5-pval_min_m)*pi)))  
  }
  
  # pval_cauchy_hybrid
  pval_matrix_hybrid <- cbind(pval_matrix_PCvm, pval_matrix_PLS)
  pval_num_h <- ncol(pval_matrix_hybrid)
  pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
  
  # ==============================================================
  # Standalone PLS Test (Tan et al. 2025)
  # ==============================================================
  if(T){
    #construct PLS test statistic based on x1 y1 and projections based on x2 y2
    pro2_pls <- cbind(sir_Upro2, sir_ypro2)                                                       # projections based on x2 y2
    pro2_num_pls_std <- ncol(pro2_pls)                                                                # the number of projections
    h1 <- n1^(-2/9)                                                                               # bandwidth
    
    pval_matrix1_PLS_std <- matrix(nrow=1, ncol=pro2_num_pls_std) # p-value matrix
    errormat1 <- U1%*%t(U1)                                                                       # residual matrix based on x1 y1
    for(q in 1:pro2_num_pls_std){
      x_pro1 <- x1%*%pro2_pls[,q]
      x_pro1_mat <-((x_pro1)%*%matrix(1,1,n1)- matrix(1,n1,1)%*%(t(x_pro1)))/h1            # kernel function matrix
      #kermat1 <-(1/sqrt(2*pi))*exp(-(x_pro1_mat^2)/2)                                     # Gaussian kernel
      indictor1 <- ifelse(abs(x_pro1_mat) <= 1, 1, 0)
      kermat1 <- (3/4)*(1-x_pro1_mat^2)*indictor1                                          # Epanechnikov kernel
      #PLS test statistics
      Tn1 <- (sum(kermat1*errormat1)-tr(kermat1*errormat1))/sqrt(2*(sum((kermat1*errormat1)^2)-tr((kermat1*errormat1)^2)))
      pval_matrix1_PLS_std[,q] <- 1-pnorm(Tn1)
    }
    
    #construct PLS test statistics based on x2 y2 and projections based on x1 y1
    pro1_pls <- cbind(sir_Upro1, sir_ypro1)                                                       # projections based on x1 y1
    pro1_num_pls_std <- ncol(pro1_pls)                                                                # the number of projections
    h2 <- n2^(-2/9)                                                                               # bandwidth
    
    pval_matrix2_PLS_std <- matrix(nrow=1, ncol=pro1_num_pls_std) # p-value matrix
    errormat2 <- U2%*%t(U2)                                                                       #residual matrix based on x2 y2
    for(l in 1:pro1_num_pls_std){
      x_pro2 <- x2%*%pro1_pls[,l]
      x_pro2_mat <-((x_pro2)%*%matrix(1,1,n2)- matrix(1,n2,1)%*%(t(x_pro2)))/h2        # kernel function matrix
      #kermat2 <-(1/sqrt(2*pi))*exp(-(x_pro2_mat^2)/2)                                 # Gaussian kernel
      indictor2 <- ifelse(abs(x_pro2_mat) <= 1, 1, 0)
      kermat2 <- (3/4)*(1-x_pro2_mat^2)*indictor2                                      # Epanechnikov kernel
      # PLS test statistics
      Tn2 <- (sum(kermat2*errormat2)-tr(kermat2*errormat2))/sqrt(2*(sum((kermat2*errormat2)^2)-tr((kermat2*errormat2)^2)))
      pval_matrix2_PLS_std[,l] <- 1-pnorm(Tn2)
    }
    
    pval_fisher1_PLS_std <- 1- pchisq(-2*sum(log(pval_matrix1_PLS_std[1:pro2_num_pls_std])),df=2*pro2_num_pls_std)
    pval_fisher2_PLS_std <- 1- pchisq(-2*sum(log(pval_matrix2_PLS_std[1:pro1_num_pls_std])),df=2*pro1_num_pls_std)
    pval_fisher_m_std <- cbind(pval_fisher1_PLS_std, pval_fisher2_PLS_std)
    pval_fisher_cauchy_std <- 1- pcauchy(mean(tan((0.5-pval_fisher_m_std)*pi)))
  }
  
  result <- list(pval_cauchy1_PCvm = pval_cauchy1_PCvm, 
                 pval_cauchy2_PCvm = pval_cauchy2_PCvm, 
                 pval_cauchy_PCvm = pval_cauchy_PCvm,
                 pval_cauchy1_PLS = pval_cauchy1_PLS,
                 pval_cauchy2_PLS = pval_cauchy2_PLS,
                 pval_cauchy_PLS = pval_cauchy_PLS,
                 pval_min1_PLS = pval_min1_PLS,
                 pval_min2_PLS = pval_min2_PLS,
                 pval_min_cauchy_PLS = pval_min_cauchy_PLS,
                 pval_fisher1_PLS = pval_fisher1_PLS_std,
                 pval_fisher2_PLS = pval_fisher2_PLS_std,
                 pval_fisher_cauchy_PLS = pval_fisher_cauchy_std,
                 pval_cauchy_hybrid = pval_cauchy_hybrid)
  
  return(result)
}

# Example
# Pcvm_Pls_hybrid_logit(50,600,0,0)
# Pcvm_Pls_hybrid_logit(100,600,0,0)
# Pcvm_Pls_hybrid_logit(200,600,0,0) 

#rp grp tests for testing linear regression models
grp_rp_lm_parallel <- function(p,n,a,pho){
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
  if(pho == 0){
    sigma <- diag(rep(1,p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }   
  beta0 <- c(rep(1,5),rep(0,p-5))
  beta1 <- c(rep(1,10),rep(0,p-10)) # /sqrt(10)
  #parallel
  # tic()
  
  cores <- detectCores(logical=F)-2
  cl <- makeCluster(cores) 
  registerDoParallel(cl, cores=cores)
  rp_grp_power <- foreach(k = 1:s, .combine='rbind', 
                          # .packages = c('MASS','matrixcalc','RPtests','GRPtests','foreach','parallel','iterators','doParallel'),
                          .packages = c('MASS','matrixcalc','RPtests','GRPtests')) %dopar%  
    {
      set.seed(k)
      x <- mvrnorm(n, mu, sigma)
      y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    # H11
      #y <- x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n)            # H12
      # y <- x %*% beta0 + a *exp(0.25*x%*%beta1) + rnorm(n)                    # H13
      
      pval_rp <- RPtest(x, y, test="nonlin", B=49L, nperms=2, resid_type = "Lasso")
      pval_grp <- GRPtest(x, y, fam = "gaussian", nsplits = 1,RP_function = NULL,output_all=FALSE,penalize=TRUE)
      result <- c(pval_rp,pval_grp)
    }
  #end parallel
  stopImplicitCluster()
  stopCluster(cl)
  
  # toc()
  #power
  rp_power <- mean(rp_grp_power[,1]<=0.05,na.rm=TRUE)
  grp_power <- mean(rp_grp_power[,2]<=0.05,na.rm=TRUE)
  return(c(rp_power,grp_power))
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
      # z <- x %*% beta0 + a *(x[,1]*x[,2] + x[,2]*x[,3] + x[,3]*x[,4] + x[,4]*x[,5])      # H22
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
 
# This function performs repeated k-fold cross-validation to compare the predictive performance of standard GLM versus penalized GLM (glmnet).
cv_compare_models <- function(x, y, 
                              model_type = c("linear", "logistic"), 
                              nfolds = 5, 
                              nrep = 20,
                              glm_highdim_limit = 3,
                              sparse_threshold = 1000) {
  library(glmnet)
  library(caret)
  library(pROC)
  library(Matrix)
 
  # model_type <- match.arg(model_type)
  n <- nrow(x)
  p <- ncol(x)

  use_glm <- (p <= glm_highdim_limit * n)  # If it is suitable to use ordinary glm/lm, skip glm/lm
  
  
  # Initialize the summary storage
  glm_sum <- if(use_glm) matrix(0, nrow=nrep, ncol=ifelse(model_type=="linear",4,2)) else NULL
  glmnet_sum <- matrix(0, nrow=nrep, ncol=ifelse(model_type=="linear",4,2))
  
  colnames(glmnet_sum) <- if(model_type=="linear") c("MSE","NMSE","R2","MAE") else c("AUC","Accuracy")
  if(use_glm) colnames(glm_sum) <- colnames(glmnet_sum)
  
  # Sparsity is only attempted when p or n exceeds this threshold
  convert_sparse <- (p > sparse_threshold | n > sparse_threshold)
  if(convert_sparse && !inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  
  for (r in 1:nrep) {
    print(r)
    # set.seed(r)
    folds <- createFolds(y, k = nfolds, list = TRUE, returnTrain = TRUE)
    
    fold_glm <- if(use_glm) matrix(NA, nrow=nfolds, ncol=ncol(glm_sum)) else NULL
    fold_glmnet <- matrix(NA, nrow=nfolds, ncol=ncol(glmnet_sum))
    
    for (k in 1:nfolds) {
      train_idx <- folds[[k]]
      test_idx <- setdiff(1:length(y), train_idx)
      
      x_train <- x[train_idx, , drop=FALSE]
      y_train <- y[train_idx]
      x_test  <- x[test_idx, , drop=FALSE]
      y_test  <- y[test_idx]
      
      if(model_type=="linear"){
        # ----------------- Ordinary lm -----------------
        if(use_glm){
          model <- lm.fit(cbind(1, as.matrix(x_train)), y_train)
          pred_glm <- as.numeric(cbind(1, as.matrix(x_test)) %*% model$coefficients)
          mse <- mean((y_test - pred_glm)^2, na.rm=TRUE)
          nmse <- mse / var(y_test)
          r2 <- 1 - mse / var(y_test)
          mae <- mean(abs(y_test - pred_glm), na.rm=TRUE)
          if(var(y_test)==0) nmse <- r2 <- mae <- NA
          fold_glm[k, ] <- c(mse, nmse, r2, mae)
        }
        
        # ----------------- glmnet -----------------
        # x_train_sparse <- as(x_train, "dgCMatrix") 
        # cvfit <- cv.glmnet(x_train_sparse, y_train, alpha=1, nfolds=5, family="gaussian")
        cvfit <- cv.glmnet(as.matrix(x_train), y_train, alpha=1, nfolds=5, family="gaussian")
        pred_glmnet <- predict(cvfit, newx=as.matrix(x_test), s="lambda.min")
        mse <- mean((y_test - pred_glmnet)^2, na.rm=TRUE)
        nmse <- mse / var(y_test)
        r2 <- 1 - mse / var(y_test)
        mae <- mean(abs(y_test - pred_glmnet), na.rm=TRUE)
        if(var(y_test)==0) nmse <- r2 <- mae <- NA
        fold_glmnet[k, ] <- c(mse, nmse, r2, mae)
        
      } else if(model_type=="logistic"){
        
        if(length(unique(y_test))<2) next
        
        # ----------------- glm -----------------
        if(use_glm){
          model <- tryCatch(glm.fit(cbind(1, as.matrix(x_train)), y_train, family=binomial()), 
                            error=function(e) NULL)
          if(!is.null(model)){
            prob <- as.numeric(1 / (1 + exp(-(cbind(1, as.matrix(x_test)) %*% model$coefficients))))
            auc <- tryCatch(auc(y_test, prob), error=function(e) NA)
            pred <- ifelse(prob>0.5,1,0)
            acc <- mean(pred==y_test)
            fold_glm[k, ] <- c(auc, acc)
          }
        }
        
        # ----------------- glmnet -----------------
        # x_train_sparse <- as(x_train, "dgCMatrix") 
        # cvfit <- cv.glmnet(x_train_sparse, y_train, alpha=1, nfolds=5, family="binomial")
        cvfit <- cv.glmnet(as.matrix(x_train), y_train, alpha=1, nfolds=5, family="binomial")
        prob <- as.numeric(predict(cvfit, newx=as.matrix(x_test), s="lambda.min", type="response"))
        auc <- tryCatch(auc(y_test, prob), error=function(e) NA)
        pred <- ifelse(prob>0.5,1,0)
        acc <- mean(pred==y_test)
        # prob <- as.numeric(predict(cvfit, newx=as.matrix(x_train), s="lambda.min", type="response"))
        # auc <- tryCatch(auc(y_train, prob), error=function(e) NA)
        # pred <- ifelse(prob>0.5,1,0)
        # acc <- mean(pred==y_train)  
        # coef_vector <- as.matrix(coef(cvfit, s = "lambda.min"))
        # sum(coef_vector != 0) 
        fold_glmnet[k, ] <- c(auc, acc)
      }
    }
    
    # if(use_glm) results_glm[[r]] <- fold_glm
    # results_glmnet[[r]] <- fold_glmnet
    if(use_glm) glm_sum[r,] <- colMeans(fold_glm, na.rm=TRUE)
    glmnet_sum[r,] <- colMeans(fold_glmnet, na.rm=TRUE)
  }
  
  return(list(
    glm_mean = if(use_glm) colMeans(glm_sum, na.rm=TRUE) else NULL, 
    glmnet_mean = colMeans(glmnet_sum, na.rm=TRUE),
    glm_sd   = if(use_glm) apply(glm_sum, 2, sd, na.rm=TRUE) else NULL,
    glmnet_sd   = apply(glmnet_sum, 2, sd, na.rm=TRUE)
  ))
}

############################################################################### 
#########             linear: Communities+and+crime data           ############ 
###############################################################################  
data_crime <- read.csv("communities+and+crime/communities.data", header = FALSE)
data_crime <- data_crime[,-c(1:5)]
data_crime[data_crime=="?"] <- NA
data_crime <- apply(data_crime,2,as.numeric)
data_crime <- data_crime[,!is.na(apply(data_crime,2,sum))]
# na_prop <- colMeans(is.na(data_crime)) 
# data_crime <- data_crime[, na_prop <= 0.5]
dim(data_crime)
sum(is.na(data_crime))

# https://search.r-project.org/CRAN/refmans/fairml/html/communities.and.crime.html
# install.packages('fairml')
# library(fairml)
# data(communities.and.crime)
# data_crime <- communities.and.crime
# dim(data_crime)
# data_crime[,1:3] <- c()
# sum(is.na(data_crime))
# # which(is.na(data_crime))/1969
# data_crime <- na.omit(data_crime)
# y <- data_crime[,101]
# x <- data_crime[,-101] 

# dim(x)
# table(y)
y <- data_crime[,100]
x <- data_crime[,-100] 
# y_scaled<- scale(y)
# x_scaled <- as.matrix(scale(x))
# if(sum(is.na(x_scaled))>0){
#   x_scaled <- as.data.frame(lapply(x, function(x1){
#     if(is.numeric(x1)){ (x1 - mean(x1, na.rm = TRUE)) / ( 1e-8+ sd(x1, na.rm = TRUE))
#     } else { x1 } }))
# }
set.seed(123)
result_81 <- Pcvm_Pls_hybrid_linear_test(x,y) 
result_81 
# 26.09.18 
# $pval_cauchy1_Pcvm[1] 0.01185575
# $pval_cauchy2_Pcvm[1] 0.009337864
# $pval_cauchy_Pcvm[1] 0.01044729
# $pval_cauchy1_PLS[1] 3.824679e-05
# $pval_cauchy2_PLS[1] 0.1119242
# $pval_cauchy_PLS[1] 7.646853e-05
# $pval_min1_PLS[1] 3.825191e-05
# $pval_min2_PLS[1] 0.1363474
# $pval_min_cauchy_PLS[1] 7.648369e-05
# $pval_fisher1_PLS[1] 1.195892e-05
# $pval_fisher2_PLS[1] 0.001036568
# $pval_fisher_cauchy_PLS[1] 2.364504e-05
# $pval_cauchy_hybrid[1] 0.0001518262


if(F){
  # A model with added interaction items
  p <- ncol(x)
  y_scaled<- scale(y)
  x_scaled <- as.matrix(scale(x))
  if(sum(is.na(x_scaled))>0){
    x_scaled <- as.data.frame(lapply(x, function(x1){
      if(is.numeric(x1)){ (x1 - mean(x1, na.rm = TRUE)) / ( 1e-8+ sd(x1, na.rm = TRUE))
      } else { x1 } }))
  }
  lasso_model1 <- cv.glmnet(as.matrix(x_scaled), y_scaled, family="gaussian", intercept = T)
  # lasso_beta1 <- coef(lasso_model1, s = "lambda.min")[-1]
  lasso_beta1 <- coef(lasso_model1, s = "lambda.1se")[-1]
  index_beta1_non0 <- seq(1:p)[as.numeric(lasso_beta1) != 0]
  x <- data_crime[,-100]
  x_sel <- as.matrix(x[, index_beta1_non0] )
  C <- do.call(cbind, lapply(1:ncol(x_sel), function(i) x_sel[,i] * x))
  # x_2 <- cbind(x,x^2,C) 
  x_2 <- cbind(x^2,C)  
  y_scaled<- scale(y)
  x2_scaled <- as.matrix(scale(x_2))
  if(sum(is.na(x2_scaled))>0){
    x2_scaled <- as.data.frame(lapply(x_2 , function(x1){
      if(is.numeric(x1)){ (x1 - mean(x1, na.rm = TRUE)) / ( 1e-8+ sd(x1, na.rm = TRUE))
      } else { x1 } }))
  }
  x2_scaled <- cbind(x,x2_scaled)  
  dim(x2_scaled) # 1994 1584
  
  set.seed(123)
  result8_2 <- Pcvm_Pls_hybrid_linear_test(x2_scaled,y)
  result8_2 
  # 26.09.18
  # $pval_cauchy1_Pcvm[1] 0.1119547
  # $pval_cauchy2_Pcvm[1] 0.2840858
  # $pval_cauchy_Pcvm[1] 0.1640482
  # $pval_cauchy1_PLS[1] 0.1498274
  # $pval_cauchy2_PLS[1] 0.6339
  # $pval_cauchy_PLS[1] 0.2933538
  # $pval_min1_PLS[1] 0.1688748
  # $pval_min2_PLS[1] 0.7889622
  # $pval_min_cauchy_PLS[1] 0.4334903
  # $pval_fisher1_PLS[1] 0.5202078
  # $pval_fisher2_PLS[1] 0.7593538
  # $pval_fisher_cauchy_PLS[1] 0.6629945
  # $pval_cauchy_hybrid[1] 0.2132719
 
  if(F){
    res1_scaled <- cv_compare_models(x, y, model_type = "linear", nrep = 20) 
    res1_scaled 
    # $glm_mean
    # MSE       NMSE         R2        MAE 
    # 0.01886274 0.34796044 0.65203956 0.09697282  
    # $glmnet_mean
    # MSE       NMSE         R2        MAE 
    # 0.01869249 0.34472271 0.65527729 0.09525850 
    # $glm_sd
    # MSE         NMSE           R2          MAE 
    # 0.0002316463 0.0046709878 0.0046709878 0.0006713176  
    # $glmnet_sd
    # MSE         NMSE           R2          MAE 
    # 0.0001686717 0.0036706034 0.0036706034 0.0005768945
    
    res_2_scaled <- cv_compare_models(x2_scaled, y, model_type = "linear", nrep = 20) 
    res_2_scaled 
    
    # $glm_mean
    # MSE NMSE   R2  MAE 
    # NaN  NaN  NaN  NaN  
    # $glmnet_mean
    # MSE       NMSE         R2        MAE 
    # 0.01819428 0.33500777 0.66499223 0.09209830  
    # $glm_sd
    # MSE NMSE   R2  MAE 
    # NA   NA   NA   NA  
    # $glmnet_sd
    # MSE         NMSE           R2          MAE 
    # 0.0002070836 0.0040108607 0.0040108607 0.0005466594 
  }
  
  # draw
  if(F){
    if(F){
      xy <- cbind(result1$y1,result1$x1_bata,result1$U1,result1$pred1)
      xy_scaled <- cbind(result11$y1,result11$x1_bata,result11$U1,result11$pred1)
      xy_1 <- cbind(result1_1$y1,result1_1$x1_bata,result1_1$U1,result1_1$pred1)
      xy_1_scaled <- cbind(result1_11$y1,result1_11$x1_bata,result1_11$U1,result1_11$pred1)
      
      xy <- as.data.frame(xy)
      xy_scaled <- as.data.frame(xy_scaled)
      xy_1 <- as.data.frame(xy_1)
      xy_1_scaled <- as.data.frame(xy_1_scaled)
      colnames(xy) <- 
        colnames(xy_scaled) <- 
        colnames(xy_1) <- 
        colnames(xy_1_scaled) <- c("y", "X_beta", "U", "pred")
      xy$case <- "xy"
      xy_scaled$case <- "xy_scaled"
      xy_1$case <- "xy_1"
      xy_1_scaled$case <- "xy_1_scaled"
      
      library(dplyr)
      all_data <- bind_rows(xy, xy_scaled, xy_1, xy_1_scaled)
      
      df1 <- all_data %>% 
        mutate(plot_type = "X_beta vs y") %>%
        select(case, plot_type, x = X_beta, y = y)
      
      df2 <- all_data %>% 
        mutate(plot_type = "pred vs U") %>%
        select(case, plot_type, x = pred, y = U)
      
      plot_data <- bind_rows(df1, df2)
      plot_data$plot_type <- factor(plot_data$plot_type,
                                    levels = c("X_beta vs y", "pred vs U"))
      
      # Drawing: 4*2 layout
      ggplot(plot_data, aes(x = x, y = y)) +
        geom_point(alpha = 0.6) +
        facet_grid(case ~ plot_type, scales = "free") +
        theme_bw() +
        labs(x = "", y = "")
    } 
    
    if(F){
      library(ggplot2)
      df <- data.frame( xbeta = xbeta, Y = Y )
      ggplot(df, aes(x = xbeta, y = Y)) +
        geom_point(shape = 21, fill = "black", color = "black", size = 0.2, stroke = 0.5) + 
        labs(
          x = expression( hat(beta)^T * X),
          y = expression("Y")
        ) + 
        theme_classic(base_size = 16) + 
        theme(
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, color = "black"),
          axis.ticks = element_line(color = "black", size = 0.6),
          axis.line = element_line(color = "black", size = 0.6),
          panel.border = element_rect(color = "black", fill = NA, size = 0.8), 
          plot.margin = margin(10, 10, 10, 10)
        ) + 
        coord_cartesian(xlim = c(-0.01,1.01),  ylim = c(-0.01, 1.01))  # 6*4.5

     
      # residual plot 
      df <- data.frame(U = U, pred = pred) 
      ggplot(df, aes(x =pred , y = U)) +
        geom_point(shape = 21, fill = "black", color = "black", size = 0.2, stroke = 0.5) +  
        labs(x = expression(hat(Y)), 
             y = expression("Residuals")) + 
        theme_classic(base_size = 16) +   
        theme(
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, color = "black"),
          axis.ticks = element_line(color = "black", size = 0.6),
          axis.line = element_line(color = "black", size = 0.6),
          panel.border = element_rect(color = "black", fill = NA, size = 0.8),  
          plot.margin = margin(10, 10, 10, 10)
        )+  coord_cartesian(xlim = c(-0.02,1.01),  ylim = c(-0.6, 0.8))
 
    }
  }
  
}
 
####################################################################################################
##################### logit: Acute myeloid leukemia (OHSU, 2022) ####################################
####################################################################################################
# data pre-processing
if(F){
  
  library(data.table)   
  library(dplyr) 
  clinical_patient <- fread("data_clinical_patient.txt",sep = "\t",header = TRUE)
  clinical_sample <- fread("aml_ohsu_2022/data_clinical_sample.txt", sep = "\t",header = TRUE )
  mrna_zscore <- fread("data_mrna_seq_rpkm_zscores_ref_all_samples.txt", sep="\t", header=TRUE)  # mRNA Z-score 
  
  a <- make.unique(mrna_zscore$Hugo_Symbol) 
  mrna_zscore$Hugo_Symbol <- c()  
  mrna_zscore <- as.data.frame(t(mrna_zscore))
  colnames(mrna_zscore) <- a 
  mrna_zscore$Patient <- rownames(mrna_zscore)
  
  clinical_patient <- clinical_patient[5:dim(clinical_patient)[1],]
  clinical_sample <- clinical_sample[5:dim(clinical_sample)[1],] 
  # colnames(clinical_patient)
  # colnames(clinical_sample)
  length(intersect(clinical_sample$`Sample Identifier`, rownames(mrna_zscore)))
  
  table(clinical_sample$ELN2017)   #  Prognosis/Risk
  table(clinical_sample$`De Novo`)  # Is it primary
  
  clinical_sample1 <- clinical_sample[,c("Sample Identifier" ,"ELN2017", "De Novo")]
  colnames(clinical_sample1)[1] <- "Patient" 
  colnames(mrna_zscore) <- make.unique(colnames(mrna_zscore))
  data_full <- mrna_zscore %>%
    left_join(clinical_sample1, by = "Patient")
  data_full <- data_full %>%
    mutate(ELN_binary = case_when(
      ELN2017 == "Adverse" ~ 1,                                                # high ris
      ELN2017 %in% c("Favorable", "Intermediate",
                     "FavorableOrIntermediate", "IntermediateOrAdverse") ~ 0, # Non-high-risk (0)
      TRUE ~ NA_real_  # (Missing/NonAML/NonInitial) NA
    ))
  data_full_filtered <- data_full %>%
    filter(!is.na(ELN_binary)) %>%                                                       # Remove NA
    filter(!ELN2017 %in% c("NonAML", "NonInitial", "MissingKaryo", "MissingMutations"))  # Remove non-AML and missing samples
  
  colnames(data_full_filtered)[22840:22847]
  data_full_filtered$Patient <- c()
  data_full_filtered$ELN2017 <- c()
  data_full_filtered <- data_full_filtered[!data_full_filtered$`De Novo`=="unknown" |
                                             !is.na(data_full_filtered$ELN_binary),]
  # data_full_filtered <- data_full_filtered[ !is.na(data_full_filtered$ELN_binary),]
  
  data_full_filtered$`De Novo` <- ifelse(data_full_filtered$`De Novo` =="TRUE",1,0 )
  
  write.csv(data_full_filtered, file = "aml_ohsu_2022/data_AML.csv", row.names = FALSE)
} 
 
data_AML <- read.csv("aml_ohsu_2022/data_AML.csv")
dim(data_AML)
apply(data_AML[,1:40],2,var)
apply(data_AML[,1:40], 2, mean)
y <- data_AML$ELN_binary
 # table(y)
# 0   1
# 269 175 
# y <- data_AML$`De Novo`
# table(y) 
# 0   1 
# 125 319
x <- data_AML[,1:(ncol(data_AML)-2)]
x <- x[,-which(is.na(apply(data_AML,2,sum)))]
# x <- apply(x,2,as.numeric)
sum(is.na(x))   

x_matrix <- as.matrix(x)
res <- cv_compare_models(x_matrix, y, model_type = "logistic", nrep = 1) 

# seurat: dimensionality reduction
if(F){
  library(Seurat)
  y <- data_AML$ELN_binary 
  x <- data_AML[,1:(ncol(data_AML)-2)]
  x <- x[,-which(is.na(apply(data_AML,2,sum)))]
  
  data_mat <- t(x)
  seu <- CreateSeuratObject(counts = data_mat)
  seu <- FindVariableFeatures( seu, selection.method = "vst", nfeatures = 2000)
  # head(VariableFeatures(seu), 10) 
  mat_hvg <- GetAssayData(seu, layer = "counts")[hvg_genes, ]
  x_2000 <- as.data.frame(t(as.matrix(mat_hvg))) 
  
  y <- data_AML$ELN_binary
  x <- x_2000
}

set.seed(123)
result_8 <- Pcvm_Pls_hybrid_logit_test(x ,y)
result_8
# x_scale <- lapply(x, function(x1){
#   if(is.numeric(x1)){
#     (x1 - mean(x1, na.rm = TRUE)) / ( 1e-8+ sd(x1, na.rm = TRUE))
#     # (x1 - min(x1, na.rm = TRUE)) / (max(x1, na.rm = TRUE) - min(x1, na.rm = TRUE) + 1e-8)
#   } else { x1 } })  
# x_scale <- do.call(cbind,x_scale)
# result_81 <- Pcvm_Pls_hybrid_logit_test( x_scale,y)
# result_81
 
# 26.09.18 (No dimensionality reduction)
# $pval_cauchy1_PCvm[1] 0.9960244
# $pval_cauchy2_PCvm[1] 0.9969783
# $pval_cauchy_PCvm[1] 0.9964701
# $pval_cauchy1_PLS[1] 0.05521259
# $pval_cauchy2_PLS[1] 0.06294995
# $pval_cauchy_PLS[1] 0.05807016
# $pval_min1_PLS[1] 0.05485493
# $pval_min2_PLS[1] 0.05698042
# $pval_min_cauchy_PLS[1] 0.05589768
# $pval_fisher1_PLS[1] 0.003061315
# $pval_fisher2_PLS[1] 0.971497
# $pval_fisher_cauchy_PLS[1] 0.006856316
# $pval_cauchy_hybrid[1] 0.9924899


set.seed(123)
result_83 <- cv_compare_models(as.matrix(x), y, model_type = "logistic", nfolds = 5, nrep = 20 ) 
result_83

# $glm_mean
# NULL
# $glmnet_mean
# AUC  Accuracy 
# 0.9292104 0.8513649 
# $glm_sd
# NULL
# $glmnet_sd
# AUC    Accuracy 
# 0.007692602 0.011093031

# $glmnet_mean
# AUC  Accuracy 
# 0.9301585 0.8522446 
# result_84 <- cv_compare_models(x_scale, y, model_type = "logistic", nfolds = 5, nrep = 20 ) 
# result_84
# $glmnet_mean
# AUC  Accuracy 
# 0.9301585 0.8522446 

# save.image(file = "aml_ohsu_2022/r_data.RData")
# load("aml_ohsu_2022/r_data.RData")
