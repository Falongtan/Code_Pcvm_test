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
        set.seed(k)
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
      set.seed(k)
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
