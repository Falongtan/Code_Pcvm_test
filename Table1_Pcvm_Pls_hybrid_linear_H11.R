# The hybrid test, the martingale-based test, and local-smoothing test for testing high dimensional linear regression models
# projections with schmidt orthogonalization
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
        
        x <- mvrnorm(n, mu, sigma)
        y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                  # H11
        # y <- x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n)         # H12
        # y <- x %*% beta0 + a *exp(0.25 * x%*%beta1) + rnorm(n)                # H13 
        
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
  beta1 <- c(rep(1,10),rep(0,p-10))# /sqrt(10)
  #parallel
  # tic()
  
  cores <- detectCores(logical=F)-2
  cl <- makeCluster(cores) 
  registerDoParallel(cl, cores=cores)
  rp_grp_power <- foreach(k = 1:s, .combine='rbind', 
                          # .packages = c('MASS','matrixcalc','RPtests','GRPtests','foreach','parallel','iterators','doParallel'),
                          .packages = c('MASS','matrixcalc','RPtests','GRPtests')) %dopar%  
    {
      
      x <- mvrnorm(n, mu, sigma)
      y <-  x %*% beta0 + a * 0.1*(x %*% beta0)^2 + rnorm(n)                    # H11
      #y <- x %*% beta0 + a * cos(0.6 * pi * x %*% beta0) + rnorm(n)              # H12
      # y <- x %*% beta0 + a *exp(0.25 * x%*%beta1) + rnorm(n)                       # H13
      
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

source("utils.R") # pvalue_integ_Brown

if(T){ 
  a <- c(0,.5, 1) 
  n <- 300                         #sample size
  p<- c(50,100,300,600,900,1200)   #dimension 
  c_h <- 1                         #bandwidths
  pho <- c(0,.4,.8)                #correlation
  r <- 9
  compar_lm_result <- matrix(0,nrow = r*length(pho)*length(a),ncol = length(p)+4)
  for (kk in 1:length(pho)){
    for (ii in 1:length(a)){
      for (jj in 1:length(p)){
        t3 <- Sys.time()
        result_pls <- Pcvm_Pls_hybrid_linear(p[jj], n, a[ii], pho[kk])
        result_grp_rp <- grp_rp_lm_parallel(p[jj], n, a[ii], pho[kk]) #  c(NA,NA) # 
        
        compar_lm_result[r*length(a)*(kk-1) + ii + 0*length(a), jj] <- result_pls$martingale_cauchy1_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 1*length(a),jj] <- result_pls$martingale_cauchy2_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 2*length(a),jj] <- result_pls$martingale_cauchy_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 3*length(a), jj] <- result_pls$pval_fisher1_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 4*length(a),jj] <- result_pls$pval_fisher2_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 5*length(a),jj] <- result_pls$pval_fisher_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 6*length(a), jj] <- result_pls$hybrid_cauchy_power
        compar_lm_result[r*length(a)*(kk-1) + ii + 7*length(a),jj] <- result_grp_rp[1]
        compar_lm_result[r*length(a)*(kk-1) + ii + 8*length(a),jj] <- result_grp_rp[2]
        t4 <- Sys.time()
        cat("p = ", p[jj],", n = ", n, ", a = ", a[ii]," , pho = ", pho[kk],
            ", martingale_cauchy1 = ", result_pls$martingale_cauchy1_power,
            ", martingale_cauchy2= ", result_pls$martingale_cauchy2_power,
            ", martingale_cauchy = ", result_pls$martingale_cauchy_power,
            ", pval_fisher1_power = ", result_pls$pval_fisher1_power,
            ", pval_fisher2_power= ", result_pls$pval_fisher2_power,
            ", pval_fisher_power = ", result_pls$pval_fisher_power,
            ", hybrid_cauchy_power = ", result_pls$hybrid_cauchy_power,
            ", rp_power = ", result_grp_rp[1],
            ", grp_power = ",result_grp_rp[2],
            ", time:",t4-t3, "\n")
        compar_lm_result[r*length(a)*(kk-1) + ii, (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 2*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 3*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 4*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 5*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 6*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 7*length(a), (length(p)+1):(length(p)+4)] <-
          compar_lm_result[r*length(a)*(kk-1) + ii + 8*length(a), (length(p)+1):(length(p)+4)] <-
          c(n, a[ii], pho[kk], c_h)
      }
    }
  }
  
  write.table(compar_lm_result, "H11.txt", sep = "  ", row.names = FALSE, col.names = FALSE)
}
