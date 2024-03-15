# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:Doubly-Robust Dynamic Treatment Regimen Estimation with Binary Outcomes
# install.packages("e1071")
install.packages("degree")
library("e1071")
library("rgenoud")
library("glmnet")
library("drgee")
# logistic regression
r <- 500
n <- 1000
res.glm00 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm00)[(ncol(res.glm00) - 3):ncol(res.glm00)] <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.glm01 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm01)[(ncol(res.glm00) - 3):ncol(res.glm00)] <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.glm02 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm02)[(ncol(res.glm00) - 3):ncol(res.glm00)] <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.AIPWE <- matrix(NA, nrow = r, ncol = 2)
colnames(res.AIPWE) <- c("psi1_hat","psi2_hat")

res.EC <- matrix(NA, nrow = r, ncol = 2)
colnames(res.EC) <- c("psi1_hat","psi2_hat")
ratio1 <- rep(NA, r); ratio2 <- rep(NA, r); ratio3 <- rep(NA, r); ratio4 <- rep(NA, r); ratio5 <- rep(NA, r)
out1 <- rep(NA, r); out2 <- rep(NA, r); out3 <- rep(NA, r); out4 <- rep(NA, r); out5 <- rep(NA, r)
# A Robust Method for Estimating Optimal Treatment Regimes
# by Baqun Zhang,* Anastasios A. Tsiatis, Eric B. Laber, and Marie Davidian
AIPWE <- function(psi) {
  X <- runif(n, 0,0.9);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3) + A*(-1 + 2*X)
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p);
  g_psi <- ifelse(cbind(1, X) %*% psi > 0, 1, 0)
  treat.mod <- glm(A ~ X + sin(X)+ (X)^2 , family = 'binomial')
  C_psi <- A* as.vector(g_psi) + (1 - A)*(1 - as.vector(g_psi))
  pi_psi <- fitted(treat.mod)*g_psi + (1 - fitted(treat.mod))*(1 - g_psi)
  mu <- glm(Y ~ X + log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), family=quasibinomial(link = "logit"))
  mu1 <- predict(mu,newdata=as.data.frame(cbind(A=1,X = X)),type="response")
  mu0 <- predict(mu,newdata=as.data.frame(cbind(A=0,X = X)),type="response")
  mu_psi <- mu1 *g_psi + mu0*(1 - g_psi)
  AIPWE <- (C_psi*Y)/pi_psi - (mu_psi)* (C_psi-pi_psi)/pi_psi
  return(mean(AIPWE))
}
X <- runif(n, 0,0.9);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
m <- X + log(abs(X)) + 1*cos(pi*X)+ 1.5*(X^3) + A*(-1 + 2*X)
p <- sigmoid(m)
p <- sigmoid(m)
Y <- rbinom(n, 1, p);
mean(Y)
Yopt <- rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                             + ifelse((cbind(1, X) %*% c(-1, 2)) > 0, (cbind(1, X) %*% c(-1, 2)), 0)))
mean(Yopt)
for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,0.9);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3) + A*(-1 + 2*X)
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p)
  simdata <- data.frame(X=X, A=A, Y=Y)
  # Step1
  treat.mod <- glm(A ~ X + sin(X)+ (X)^2 , family = 'binomial')
  w = abs(A - fitted(treat.mod))
  mu <- glm(Y ~ X + A + I(A*X), family=quasibinomial(link = "logit"))
  res.glm00[i,] <- mu$coefficients
  ratio1[i] <- sum(ifelse((cbind(1, X) %*% res.glm00[i,(ncol(res.glm00) - 1):ncol(res.glm00)]) > 0, 1, 0) == ifelse((cbind(1, X) %*% c(-1, 2)) > 0, 1, 0))/n

  out1[i] <- mean(rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                                       + ifelse((cbind(1, X) %*% res.glm00[i,(ncol(res.glm00) - 1):ncol(res.glm00)]) > 0, (cbind(1, X) %*% res.glm00[i,(ncol(res.glm00) - 1):ncol(res.glm00)]), 0))))
  res.glm01[i,] <- glm(Y ~ X+ A + I(A*X), weights = w, family=quasibinomial(link = "logit"))$coefficients
  par <- as.vector(res.glm01[i,])
  ratio2[i] <- sum(ifelse((cbind(1, X) %*% res.glm01[i,(ncol(res.glm00) - 1):ncol(res.glm00)]) > 0, 1, 0) == ifelse((cbind(1, X) %*% c(-1, 2)) > 0, 1, 0))/n
  out2[i] <- mean(rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                                       + ifelse((cbind(1, X) %*% res.glm01[i,(ncol(res.glm00) - 1):ncol(res.glm00)]) > 0, (cbind(1, X) %*% res.glm01[i,(ncol(res.glm00) - 1):ncol(res.glm00)]), 0))))
  # Step2
  mu <- cbind(1, X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dsigmoid(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res.glm02[i,] <- glm(Y ~ X + A + I(A*X), weights = w_n, family=quasibinomial(link = "logit"))$coefficients

  ratio3[i] <- sum(ifelse(cbind(1, X) %*% res.glm02[i,(ncol(res.glm00) - 1):ncol(res.glm00)] > 0, 1, 0) == ifelse(cbind(1, X) %*% c(-1, 2) > 0, 1, 0))/n
  out3[i] <- mean(rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                                       + ifelse((cbind(1, X) %*% res.glm02[i,(ncol(res.glm00) - 1):ncol(res.glm00)]) > 0, (cbind(1, X) %*% res.glm02[i,(ncol(res.glm00) - 1):ncol(res.glm00)]), 0))))
  AIPWE.res  <- genoud(AIPWE, nvars=2, max=TRUE, pop.size=200, max.generations=3,
                       wait.generations=5, gradient.check=FALSE, print=1)
  res.AIPWE[i,] <- AIPWE.res$par
  ratio4[i] <- sum(ifelse(cbind(1, X) %*% res.AIPWE[i,] > 0, 1, 0) == ifelse(cbind(1, X) %*% c(-1, 2) > 0, 1, 0))/n
  out4[i] <- mean(rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                                       + ifelse(cbind(1, X) %*% res.AIPWE[i,] > 0,cbind(1, X) %*% res.AIPWE[i,], 0))))
  dr.est<-drgee(oformula=formula(Y~X), eformula=formula(A~X + sin(X)+ (X)^2),
                iaformula=formula(~X), olink="logit",elink="logit", data=simdata,estimation.method="dr")
  EC.res <- summary(dr.est)
  res.EC[i,] <- as.vector(EC.res$coefficients[,1])
  ratio5[i] <- sum(ifelse(cbind(1, X) %*% res.EC[i,] > 0, 1, 0) == ifelse(cbind(1, X) %*% c(-1, 2) > 0, 1, 0))/n
  out5[i] <- mean(rbinom(n, 1, sigmoid(X + log(abs(X)) + 1*cos(pi*X) + 1.5*(X^3)
                                       + ifelse(cbind(1, X) %*% res.EC[i,] > 0,cbind(1, X) %*% res.EC[i,], 0))))
}


apply(na.omit(res.glm00), 2 , mean)
apply(na.omit(res.glm01), 2 , mean)
apply(na.omit(res.glm02), 2 , mean)
apply(na.omit(res.AIPWE), 2 , mean)
apply(na.omit(res.EC), 2 , mean)

apply(res.glm00, 2 , sd)
apply(res.glm01, 2 , sd)
apply(res.glm02, 2 , sd)
apply(res.AIPWE, 2 , sd)
apply(res.EC, 2 , sd)

mean(na.omit(ratio1)); mean(na.omit(ratio2)); mean(na.omit(ratio3)); mean(na.omit(ratio4)); mean(na.omit(ratio5))


