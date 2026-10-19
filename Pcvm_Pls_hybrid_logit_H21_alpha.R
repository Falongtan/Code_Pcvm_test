# Projection test 
run_single_setting<- function(p, n, a = 0, pho = 0, s = 1000) {
  library(MASS)
  library(glmnet)
  library(LassoSIR)
  library(VariableScreening)
  library(foreach)
  library(parallel)
  library(doParallel)
  
  mu <- rep(0, p)
  if(pho == 0){
    sigma <- diag(rep(1, p))
  }else{
    v <- pho^(0:(p-1))
    sigma <- toeplitz(v)
  }
  beta0 <- c(rep(1, 5), rep(0, p - 5))
  
  cores <- detectCores(logical = FALSE) - 2
  if (is.na(cores) || cores < 1) cores <- 4
  
  random_port <- sample(10000:60000, 1)
  cl <- makeCluster(cores, port = random_port)
  registerDoParallel(cl, cores = cores) 
  hybrid_frequency <- 
    foreach(k = 1:s, #.errorhandling = "remove",
            .packages = c('MASS', 'glmnet', 'LassoSIR', 'VariableScreening')) %dopar% { 
    set.seed(k)
    x <- mvrnorm(n, mu, sigma)
    z <- x %*% beta0 + a*0.2*(x %*% beta0)^2                                             
    pr <- 1/(1 + exp(-z))
    y <- matrix(rbinom(n, 1, pr), ncol = 1)
    
    index_x <- sample(1:n, floor(n/2), replace = FALSE)
    index_x <- sort(index_x)
    
    x1 <- x[-index_x, , drop=FALSE]
    x2 <- x[index_x, , drop=FALSE]
    
    y1 <- as.matrix(y[-index_x, ])
    y2 <- as.matrix(y[index_x, ])
    
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
      x1_sec <- x1[, index_beta1_non0, drop=FALSE]                                         
      sec_model1 <- glm(y1 ~ x1_sec, family = binomial(link = "logit"))
      beta1_lasso[index_beta1_non0] <- unname(sec_model1$coefficients)[-1]
      beta1_hat <- beta1_lasso
      intercept1 <- unname(sec_model1$coefficients)[1]
      pred1 <- predict(sec_model1, newdata = as.data.frame(x1_sec), type="response")
      pred1 <- matrix(unname(pred1), ncol = 1)
      U1 <- y1 - pred1                                                       
      beta1_pro <- beta1_hat/sqrt(sum(beta1_hat^2))
      beta1_pro <- as.matrix(beta1_pro)
    }
    
    screen_num1 <- floor(n1/log(n1))
    sir_Upro1 <- tryCatch({
      if(p <= screen_num1){
        sir_U1 <- LassoSIR(x1, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta1 <- sir_U1$beta
        sir_Ubeta1/sqrt(colSums(sir_Ubeta1^2))
      }else{
        rank_U1 <- screenIID(X=x1, Y=U1, method = "DC-SIS")                      
        index_U1 <- seq(1:p)[rank_U1$rank <= screen_num1]                        
        x1_Uscreen <- x1[, index_U1, drop=FALSE]
        sir_Ulasso1 <- LassoSIR(x1_Uscreen, U1, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta1 <- sir_Ulasso1$beta
        sir_Upro1_temp <- matrix(0, nrow = p, ncol = ncol(sir_Ubeta1))
        sir_Upro1_temp[index_U1,] <- sir_Ubeta1
        sir_Upro1_temp/sqrt(colSums(sir_Upro1_temp^2))
      }
    }, error = function(e){ NA })
    
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
      x2_sec <- x2[, index_beta2_non0, drop=FALSE]
      sec_model2 <- glm(y2 ~ x2_sec, family = binomial(link = "logit"))
      beta2_lasso[index_beta2_non0] <- unname(sec_model2$coefficients)[-1]
      beta2_hat <- beta2_lasso
      intercept2 <- unname(sec_model2$coefficients)[1]
      pred2 <- predict(sec_model2, newdata = as.data.frame(x2_sec), type="response")
      pred2 <- matrix(unname(pred2), ncol = 1)
      U2 <- y2-pred2
      beta2_pro <- beta2_hat/sqrt(sum(beta2_hat^2))
      beta2_pro <- as.matrix(beta2_pro)
    }
    
    screen_num2 <- floor(n2/log(n2))
    sir_Upro2 <- tryCatch({
      if(p <= screen_num2){
        sir_Ulasso2 <- LassoSIR(x2, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta2 <- sir_Ulasso2$beta
        sir_Ubeta2/sqrt(colSums(sir_Ubeta2^2))
      }else{
        rank_U2 <- screenIID(X=x2, Y=U2, method = "DC-SIS")                      
        index_U2 <- seq(1:p)[rank_U2$rank <= screen_num2]                        
        x2_Uscreen <- x2[, index_U2, drop=FALSE]
        sir_Ulasso2 <- LassoSIR(x2_Uscreen, U2, H=5, choosing.d="automatic", solution.path=FALSE, categorical=FALSE, nfolds=5, screening=FALSE)
        sir_Ubeta2 <- sir_Ulasso2$beta
        sir_Upro2_temp <- matrix(0, nrow = p, ncol = ncol(sir_Ubeta2))
        sir_Upro2_temp[index_U2,] <- sir_Ubeta2
        sir_Upro2_temp/sqrt(colSums(sir_Upro2_temp^2))
      }
    }, error = function(e){ NA })
    
    cor1 <- numeric(0)
    if (!any(is.na(sir_Upro1))) {
      number_pro1 <- ncol(sir_Upro1)
      cor1 <- rep(0, number_pro1)
      for(kk in 1:number_pro1){
        cor1[kk] <- t(beta2_pro) %*% sir_Upro1[,kk]
      }
    }
    
    cor2 <- numeric(0)
    if (!any(is.na(sir_Upro2))) {
      number_pro2 <- ncol(sir_Upro2)
      cor2 <- rep(0, number_pro2)
      for(ll in 1:number_pro2){
        cor2[ll] <- t(beta1_pro) %*% sir_Upro2[,ll]
      }
    }
    
    cor <- as.numeric(c(abs(cor1), abs(cor2)))
    
    list(rho_alpha = cor, num_dirs = length(cor))
}
  stopImplicitCluster()
  stopCluster(cl)
  
  all_rho <- unlist(lapply(hybrid_frequency, function(res) res$rho_alpha))
  num_dirs_vec <- sapply(hybrid_frequency, function(res) res$num_dirs)
  rho_valid <- all_rho[!is.na(all_rho)]
  if(length(rho_valid) == 0) rho_valid <- c(0)
  summary_row <- data.frame( 
    n = n,
    p = p,
    pho = pho,
    mean_dirs = mean(num_dirs_vec, na.rm = TRUE),
    max_dirs = suppressWarnings(max(num_dirs_vec, na.rm = TRUE)),
    rho_mean = mean(rho_valid),
    rho_80th = unname(quantile(rho_valid, 0.80)),
    rho_85th = unname(quantile(rho_valid, 0.85)),
    rho_90th = unname(quantile(rho_valid, 0.90)),
    rho_95th = unname(quantile(rho_valid, 0.95)),
    rho_99th = unname(quantile(rho_valid, 0.99)),
    prop_over_080 = mean(rho_valid > 0.80),
    prop_over_085 = mean(rho_valid > 0.85),
    prop_over_090 = mean(rho_valid > 0.90),
    prop_over_095 = mean(rho_valid > 0.95),
    prop_over_099 = mean(rho_valid > 0.99)
  )
  if(is.infinite(summary_row$max_dirs)) summary_row$max_dirs <- 0
  return(summary_row)
}

n_seq <- c(600)  
p_seq <- c(50, 100, 300, 600, 900, 1200)
pho_seq <- c(0, 0.4, 0.8)  
s_simulations <- 1000     


final_results <- data.frame()  
save_path <- "Results_hybrid_cross.csv"

for(current_pho in pho_seq) {
  for(current_p in p_seq) {
    for(current_n in n_seq) {
      cat(sprintf("-> Runing: pho = %.1f, p = %d, n = %d, ", current_pho, current_p, current_n))
      res_row <- run_single_setting(p = current_p, n = current_n, a = 0, pho = current_pho, s = s_simulations)
      print(res_row)
      cat("\n")
      final_results <- rbind(final_results, res_row)
      if(!file.exists(save_path)) {
        write.table(res_row, save_path, sep = ",", row.names = FALSE, col.names = TRUE)
      } else {
        write.table(res_row, save_path, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      }
    }
  }
} 
rownames(final_results) <- NULL
final_results_sorted <- final_results[order(final_results$pho, final_results$n, final_results$p), ]
print(final_results_sorted)
 

 