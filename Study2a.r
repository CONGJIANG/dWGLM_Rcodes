# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article: Novel Robust Dynamic Treatment Regimen Estimation for Discrete Outcome: An application to smoking cessation using e-cigarettes
# dynamic weighted generalized linear model (dWGLM)
# install.packages("e1071")
library("e1071")
library("vioplot")
dWGLM <- function(outcome.mod, blip.mod, treat.mod, tf.mod, k, data, R = 10) {
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
    # blip parameter estimates extracted via dimensions of Hbeta and Hpsi
    ## First step
    glm1 <- glm(Y~0+Hbeta+A:Hpsi, weights = w, family=quasibinomial(link = "logit"))
    mu <- predict(glm1, newdata=data.frame(Hbeta=Hbeta, A = 1 - A, Hpsi = Hpsi))
    w_n = w * dsigmoid(mu)
    ## Second step
    Psi0 <- matrix(NA, R, length(glm1$coef[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])]))
    if (j == k){
      glm2 <- glm(Y~0+Hbeta+A:Hpsi, weights = w_n, family=quasibinomial(link = "logit"))
      psi <- as.numeric(glm2$coef[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])])
    } else {
      for (rr in 1:R) {
        Y <-rbinom(n,1, Y_final)
        glm2 <- glm(Y~0+Hbeta+A:Hpsi, weights = w_n, family=quasibinomial(link = "logit"))
        psi <- as.numeric(glm2$coef[(dim(Hbeta)[2]+1):(dim(Hbeta)[2]+dim(Hpsi)[2])])
        Psi0[rr,] <- psi
      }
      psi <- as.vector(apply(Psi0, 2, mean))
    }
    # use this to identify optimal treatments
    opt <- as.numeric(Hpsi %*% psi  > 0)
    
    # update (pseudo-) Y by adding regret
    Y_final <- expit(predict(glm2, newdata=data.frame(Hbeta= model.matrix(tf.mod[[k]], model.frame(tf.mod[[k]], data))
                                                      , A = model.response(model.frame(treat.mod[[k]],data))
                                                      , Hpsi = model.matrix(blip.mod[[k]], model.frame(blip.mod[[k]], data)))) 
                     + (opt - A)*(Hpsi %*% psi) )
    
    # store estimates
    obj$psi[[j]] <- psi
  }
  # return estimates
  return(obj)
}




# data setup
n <- 1000
# expit function
expit <- function(x) {1/(1+exp(-x))}

# model parameters a0, a1 (treatment model), p0, p1 (blip model)
a0 <- 0; a1 <- 1; p0 <- -2; p1 <- 1

r <- 1000
stage2M3 <- matrix(NA, nrow = r, ncol = 2)
colnames(stage2M3) <- c("psi0_hat","psi1_hat")

stage1M3 <- matrix(NA, nrow = r, ncol = 2)
colnames(stage1M3) <- c("psi0_hat","psi1_hat")

for (i in 1:r) {
  X1 <- rnorm(n,2,1)
  A1 <- rbinom(n,1,expit(-5 + a1*X1 + a1*I(X1^2)))
  
  X2 <- rnorm(n,1+ 0.5*X1, 1)
  A2 <- rbinom(n,1,expit(-2.5*X2 + a1*sin(X2) + a1*I(X2^2)))
  
  # regrets
  A1opt <- as.numeric(p0 + p1*X1 > 0)
  A2opt <- as.numeric(p0 + p1*X2 > 0)
  reg1 <- (p0 + p1*X1)*(A1opt- A1)
  reg2 <- (p0 + p1*X2)*(A2opt- A2)
  
  lgtYopt <- X1 + log(abs(X1)) + cos(pi*X1) 
  # outcome mean: optimal outcome minus regrets at each stage
  p <- expit(lgtYopt - reg1 -reg2)
  Y <- rbinom(n, 1, p)
  
  # analysis
  # models to be passed to dWGLM
  outcome.mod <- Y ~ 1
  blip.mod <- list(~X1,~X2)
  treat.mod <- list(A1~ X1 + I(X1^2), A2~ X2)
  tf.mod <- list(~ X1, ~ X1 * A1 + X2 + I(log(abs(X1))) + I(cos(pi*X1)))
  mydata <- data.frame(X1,X2,A1,A2,Y)
  k <- 2
  
  # perform dynamic WGLM
  mod <- dWGLM(outcome.mod,blip.mod,treat.mod,tf.mod,k,data=mydata)
  stage2M3[i,] <- mod$psi[[2]]
  stage1M3[i,] <- mod$psi[[1]]
}

#library("vioplot")
#library("e1071")
par(mfrow=c(2,2))
vioplot( stage2M3[,1], names=c("dWGLM"),  
         xlab= "Stage 2 (n = 1000)", ylab = expression(paste( psi[0], " estimates")) )
abline(h = p0, col = "red")
vioplot(  stage2M3[,2], names=c("dWGLM"),  
          xlab= "Stage 2 (n = 1000)", ylab = expression(paste( psi[1], " estimates")) )
abline(h = p1, col = "red")


vioplot( stage1M3[,1], names=c("dWGLM"),  
         xlab= "Stage 1 (n = 1000)", ylab = expression(paste( psi[0], " estimates")) )
abline(h = p0, col = "red")
vioplot(  stage1M3[,2], names=c("dWGLM"),  
          xlab= "Stage 1 (n = 1000)", ylab = expression(paste( psi[1], " estimates")) )
abline(h = p1, col = "red")
mtext("Blip estimates of Case 2 (study 2a)", side = 3, line = -2, outer = TRUE)







