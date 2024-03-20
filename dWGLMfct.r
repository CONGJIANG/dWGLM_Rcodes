# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article: Novel Robust Dynamic Treatment Regimen Estimation for Discrete Outcome: An application to smoking cessation using e-cigarettes
# 
#require("e1071")
#require("vioplot")
dWGLM <- function(outcome.mod, blip.mod, treat.mod, tf.mod, k, data, R = 25) {
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

