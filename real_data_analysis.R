# The hybrid test, the martingale-based test, and local-smoothing test for testing high dimensional linear and logistical regression models

# linear + test
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
      high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9, verbose = FALSE)
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
      high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9, verbose = FALSE)
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
 
  }
  
  # pval_cauchy_hybrid 
  pval_matrix_hybrid <- cbind(pval_matrix_Pcvm, pval_matrix_PLS)
  pval_num_h <- ncol(pval_matrix_hybrid)
  pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
 
  result <- list(pval_cauchy1_Pcvm = pval_cauchy1_Pcvm, 
                 pval_cauchy2_Pcvm = pval_cauchy2_Pcvm, 
                 pval_cauchy_Pcvm = pval_cauchy_Pcvm, 
                 pval_cauchy_hybrid = pval_cauchy_hybrid) 
  
  return(result) 
}
# logit + test
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
  }
  
  # pval_cauchy_hybrid
  pval_matrix_hybrid <- cbind(pval_matrix_PCvm, pval_matrix_PLS)
  pval_num_h <- ncol(pval_matrix_hybrid)
  pval_cauchy_hybrid <- 1 - pcauchy(mean(tan((0.5-pval_matrix_hybrid[1:pval_num_h])*pi)))
 
  result <- list(pval_cauchy1_PCvm = pval_cauchy1_PCvm, 
                 pval_cauchy2_PCvm = pval_cauchy2_PCvm, 
                 pval_cauchy_PCvm = pval_cauchy_PCvm, 
                 pval_cauchy_hybrid = pval_cauchy_hybrid)
  
  return(result)
}

source("utils.R") # bandwidth_choice  pvalue_integ_Brown

# This function performs repeated k-fold cross-validation to compare the predictive performance of standard GLM versus penalized GLM (glmnet).
cv_compare_models <- function(x, y, model_type = c("linear", "logistic"), 
                              nfolds = 5, nrep = 20, glm_highdim_limit = 3, sparse_threshold = 1000) {
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
    # print(r)
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
            auc <- tryCatch(auc(y_test, prob, quiet = TRUE), error=function(e) NA)
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
        auc <- tryCatch(auc(y_test, prob, quiet = TRUE), error=function(e) NA)
        pred <- ifelse(prob>0.5,1,0)
        acc <- mean(pred==y_test) 
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
dim(data_crime)
sum(is.na(data_crime))
 
y <- data_crime[,100]
x <- data_crime[,-100] 
table(y)

set.seed(123)
t3 <- Sys.time()
result_81 <- Pcvm_Pls_hybrid_linear_test(x,y) 
t4 <- Sys.time()
cat("Dataset: Communities & Crime (Linear) | n = ", nrow(x), ", p = ", ncol(x), "\n",
    "cauchy1_pcvm = ", result_81$pval_cauchy1_Pcvm, 
    ", cauchy_pcvm = ", result_81$pval_cauchy_Pcvm, 
    ", cauchy_hybrid = ", result_81$pval_cauchy_hybrid,
    " | time:", as.numeric(difftime(t4, t3, units="secs")), "s \n\n") 

# Quadratic polynomial model
if(T){
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
  t3 <- Sys.time()
  result8_2 <- Pcvm_Pls_hybrid_linear_test(x2_scaled,y)
  t4 <- Sys.time()
  cat("Dataset: Communities & Crime (Quadratic polynomial model) | n = ", nrow(x2_scaled), ", p = ", ncol(x2_scaled), "\n",
      "cauchy1_pcvm = ", result8_2$pval_cauchy1_Pcvm,
      ", cauchy_pcvm = ", result8_2$pval_cauchy_Pcvm,
      ", cauchy_hybrid = ", result8_2$pval_cauchy_hybrid,
      " | time:", as.numeric(difftime(t4, t3, units="secs")), "s \n\n")
 
  if(F){
    res1_scaled <- cv_compare_models(x, y, model_type = "linear", nrep = 20) 
    res1_scaled 
 
    res_2_scaled <- cv_compare_models(x2_scaled, y, model_type = "linear", nrep = 20) 
    res_2_scaled  
  }
  
  # draw
  cat("\n=== The scatter plot ===\n")  
  if(T){
    library(ggplot2)
    X <- x
    Y <- y
    p <- ncol(X)
    set.seed(123)
    lasso_model <- cv.glmnet(as.matrix(X), Y, family="gaussian")
    beta_lasso <- coef(lasso_model, s = "lambda.min")[-1]
    index_beta_non0 <- seq(1:p)[as.numeric(beta_lasso) != 0]
    x_sec <- X[,index_beta_non0, drop=FALSE]                                        
    sec_model = glm(Y~x_sec-1,family = gaussian)                           
    
    pred = predict(sec_model, newx = x_sec, type="response")
    U <- Y-pred  
    xbeta <- x_sec%*%sec_model$coefficients   

    #  scatter plot
    df1 <- data.frame( xbeta = xbeta, Y = Y )
    p1 <- ggplot(df1, aes(x = xbeta, y = Y)) +
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
      ) +  coord_cartesian(xlim = c(-0.01,1.01),  ylim = c(-0.01, 1.01))  # 6*4.5
    print(p1)
    ggsave("Figure_1_scatter.pdf", plot = p1, width = 6, height = 4.5, units = "in")
    
    # residual plot  
    set.seed(123)
    lasso_model_quad <- cv.glmnet(as.matrix(x2_scaled), Y, family="gaussian")
    beta_lasso_quad <- coef(lasso_model_quad, s = "lambda.min")[-1]
    index_beta_non0_quad <- seq(1:ncol(x2_scaled))[as.numeric(beta_lasso_quad) != 0]
    
    x_sec_quad <- x2_scaled[, index_beta_non0_quad, drop=FALSE]                                        
    sec_model_quad = glm(Y~x_sec_quad-1, family = gaussian)                           
    
    pred_quad = predict(sec_model_quad, newx = x_sec_quad, type="response")
    U_quad <- Y - pred_quad  
    
    df2 <- data.frame(U = U_quad, pred = pred_quad) 
    p2 <- ggplot(df2, aes(x = pred , y = U)) +
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
      ) +  coord_cartesian(xlim = c(-0.02,1.01),  ylim = c(-0.6, 0.8))
    print(p2)
    ggsave("Figure_2_residual.pdf", plot = p2, width = 6, height = 4.5, units = "in")
  }
  cat("\n=== Plots generated and saved successfully!! ===\n")  
  
}  
 
####################################################################################################
##################### logit: Acute myeloid leukemia (OHSU, 2022) ####################################
####################################################################################################
data_AML <- read.csv("aml_ohsu_2022/data_AML.csv")
dim(data_AML)
apply(data_AML[,1:40],2,var)
apply(data_AML[,1:40], 2, mean)
y <- data_AML$ELN_binary
table(y)

x <- data_AML[,1:(ncol(data_AML)-2)]
x <- x[,-which(is.na(apply(data_AML,2,sum)))]
# x <- apply(x,2,as.numeric)
sum(is.na(x))   

# x_matrix <- as.matrix(x)
# res <- cv_compare_models(x_matrix, y, model_type = "logistic", nrep = 1) 

set.seed(123)
t3 <- Sys.time()
result_8 <- Pcvm_Pls_hybrid_logit_test(x ,y)
t4 <- Sys.time()
cat("Dataset: OHSU AML 2022 (Logistic) | n = ", nrow(x), ", p = ", ncol(x), "\n",
    "cauchy1_pcvm = ", result_8$pval_cauchy1_PCvm, 
    ", cauchy_pcvm = ", result_8$pval_cauchy_PCvm, 
    ", cauchy_hybrid = ", result_8$pval_cauchy_hybrid,
    " | time:", as.numeric(difftime(t4, t3, units="secs")), "s \n\n") 
 
set.seed(123)
t3 <- Sys.time()
result_83 <- cv_compare_models(as.matrix(x), y, model_type = "logistic", nfolds = 5, nrep = 20 ) 
t4 <- Sys.time()
cat("Dataset: OHSU AML 2022 (Logistic CV Performance) | n = ", nrow(x), ", p = ", ncol(x), "\n",
    "glmnet_AUC_mean = ", result_83$glmnet_mean["AUC"],
    ", glmnet_ACC_mean = ", result_83$glmnet_mean["Accuracy"],
    ", glmnet_AUC_sd = ", result_83$glmnet_sd["AUC"],
    ", glmnet_ACC_sd = ", result_83$glmnet_sd["Accuracy"],
    " | time:", as.numeric(difftime(t4, t3, units="secs")), "s \n\n")

# save.image(file = "aml_ohsu_2022/r_data.RData")
# load("aml_ohsu_2022/r_data.RData")
