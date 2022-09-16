# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article: Doubly-robust dynamic treatment regimen 

# Doubly-robust dynamic treatment regimen estimation via dynamic generalized linear model (dWGLM), assuming linear blip functions and logistic regression for treatment model
# install.packages("e1071")
library("e1071")
library("vioplot")
absplus <- function(x) {ifelse(x > 0, x, 0)}
dWGLM <- function(outcome.mod, blip.mod, treat.mod, tf.mod, k, data, R = 5) {
  obj <- list()
  # extract original outcome
  Y <- model.response(model.frame(outcome.mod,data))
  # work in stages starting from final stage (k)
  for (j in k:1) {
    # treatment model for weights (assuming logistic regression)
    alpha <- glm(treat.mod[[j]],binomial)
    A <- model.response(model.frame(treat.mod[[j]],data))
    # weights
    w <- abs(A - fitted(alpha))
    # regress Y on terms in treatment free model and blip model
    Hpsi <- model.matrix(blip.mod[[j]], model.frame(blip.mod[[j]], data))
    Hbeta <- model.matrix(tf.mod[[j]], model.frame(tf.mod[[j]], data))	
    Psi0 <- matrix(NA, R, dim(Hpsi)[2])
    # blip parameter estimates extracted via dimensions of Hbeta and Hpsi
    if (j == k){
      ## First step
      glm1 <- glm(Y~0+Hbeta+A:Hpsi, weights = w, family=quasibinomial(link = "logit"))
      mu <- predict(glm1, newdata=data.frame(Hbeta=Hbeta, A = 1 - A, Hpsi = Hpsi))
      w_n = w * dsigmoid(mu)
      ## Second step
      glm2 <- glm(Y~0+Hbeta+A:Hpsi, weights = w_n, family=quasibinomial(link = "logit"))
      psi <- as.numeric(glm2$coef[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])])
      # use this to identify optimal treatments
      opt <- as.numeric(Hpsi %*% psi  > 0)
      Y_final_mean <- predict(glm2, newdata=data.frame(Hbeta= model.matrix(tf.mod[[k]], model.frame(tf.mod[[k]], data))
                                                    , A = model.response(model.frame(treat.mod[[k]],data))
                                                    , Hpsi = model.matrix(blip.mod[[k]], model.frame(blip.mod[[k]], data)))) 
      # update (pseudo-) Y by adding regret
      Y_final <- expit(Y_final_mean + (opt - A)*(Hpsi %*% psi) )
    } else {
      for (rr in 1:R) {
        Y <-rbinom(n,1, Y_final)
        ## First step
        glm1 <- glm(Y~0+Hbeta+A:Hpsi, weights = w, family=quasibinomial(link = "logit"))
        mu <- predict(glm1, newdata=data.frame(Hbeta=Hbeta, A = 1 - A, Hpsi = Hpsi))
        w_n = w * dsigmoid(mu)
        ## Second step
        glm2 <- glm(Y~0+Hbeta+A:Hpsi, weights = w_n, family=quasibinomial(link = "logit"))
        psi <- as.numeric(glm2$coef[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])])
        Psi0[rr,] <- psi
      }
      psi <- as.vector(apply(Psi0, 2, mean))
      # use this to identify optimal treatments
      opt <- as.numeric(Hpsi %*% psi  > 0)
      # update (pseudo-) Y by adding regret
      Y_final  <- expit(Y_final_mean + (opt - A)*(Hpsi %*% psi) )
    }

    
    # store estimates
    obj$psi[[j]] <- psi
  }
  # return estimates
  return(obj)
}



# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:

# data setup
n <- 1000
# expit function
expit <- function(x) {1/(1+exp(-x))}

# model parameters a0, a1 (treatment model), p0, p1 (blip model)
a0 <- 0; a1 <- 1; p20 <- 0.25; p21 <- 0.5; p22 <- 0.35; th3 <- -0.5; th4 <- -0.1;
th6 <- p20; th7 <- p21; th8 <- p22; del1 <- 0.5; del2 <- 0.5

r <- 500
stage2M3 <- matrix(NA, nrow = r, ncol = 3)
colnames(stage2M3) <- c("psi0_hat","psi1_hat", "psi2_hat")

stage1M3 <- matrix(NA, nrow = r, ncol = 2)
colnames(stage1M3) <- c("psi0_hat","psi1_hat")



for (i in 1:r) {
  X1 <- rnorm(n,3,1)
  A1 <- rbinom(n,1,expit(-2.5 + 1.25*X1))
  X2 <- rnorm(n, mean = -0.5 + 0.5*X1)
  A2 <- rbinom(n,1,expit(-0.5 + 1.25*X2))
  
  O1 <- rbinom(n,1,0.5)
  O2 <- rbinom(n,1,expit(del1*O1 + del2*A1))
  f1 <- -log(abs(X1))
  f2 <- cos(pi*X2) - I(X2^3)
  #f2 <- X2^3*as.numeric(X2 > 0.25)
  
  mu.mean <- X1 + th3*A1 + th4*O1*A1 + X2 + p20*A2 + p21*I(O2*A2) + p22*I(A1*A2) + f1 + f2
  Y <- rbinom(n, 1, expit(mu.mean))
  
  # analysis
  # models to be passed to dWOLS
  outcome.mod <- Y ~ 1
  blip.mod <- list(~O1,~ O2 + A1)
  treat.mod <- list(A1~ X1, A2~ X2)
  tf.mod <- list(~ X1 + O1, ~ X1 + O1 * A1 + O2 + X2)
  mydata <- data.frame(X1,X2,O1,O2,A1,A2,Y)
  k <- 2
  
  # perform dynamic WOLS
  mod <- dWGLM(outcome.mod,blip.mod,treat.mod,tf.mod,k,data=mydata)
  stage2M3[i,] <- mod$psi[[2]]
  stage1M3[i,] <- mod$psi[[1]]
}

#library("vioplot")
#library("e1071")
par(mfrow=c(2,3))
vioplot( stage2M3[,1], names=c("dWGLM"),  
         xlab= "Stage 2 (n = 1000)", ylab = expression(paste( psi[0], " estimates")) )
abline(h = p20, col = "red")
vioplot(  stage2M3[,2], names=c("dWGLM"),  
          xlab= "Stage 2 (n = 1000)", ylab = expression(paste( psi[1], " estimates")) )
abline(h = p21, col = "red")
vioplot(  stage2M3[,3], names=c("dWGLM"),  
          xlab= "Stage 2 (n = 1000)", ylab = expression(paste( psi[2], " estimates")) )
abline(h = p22, col = "red")
p10 <- th3 + absplus(th6 + th8) - absplus(th6) + expit(del2)*(absplus(th6 + th7 + th8) - absplus(th6 + th8)) - expit(0)*(absplus(th6 + th7) - absplus(th6))
p11 <- th4 + (expit(del1 + del2) - expit(del2))*(absplus(th6 + th7 + th8) - absplus(th6 + th8)) - (expit(del1) - expit(0))*(absplus(th6 + th7) - absplus(th6))

vioplot( stage1M3[,1], names=c("dWGLM"),  
         xlab= "Stage 1 (n = 1000)", ylab = expression(paste( psi[0], " estimates")) )
abline(h = p10, col = "red")
vioplot(  stage1M3[,2], names=c("dWGLM"),  
          xlab= "Stage 1 (n = 1000)", ylab = expression(paste( psi[1], " estimates")) )
abline(h = p11, col = "red")



