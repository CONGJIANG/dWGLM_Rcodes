# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:Novel Robust Dynamic Treatment Regimen Estimation for Discrete Outcome: An application to smoking cessation using e-cigarettes
# install.packages("e1071")
library("e1071")
options(scipen=200)
##############################################
###### Scenario 1
##############################################
###### Scenario 1
r <- 1000
res1_0 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_0) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")


res1_1 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_1) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res1_2 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_2) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2); A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi *X) + (X^3) + A*(-1 + 2*X) 
  p <- pnorm(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X, family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res1_0[i,] <- glm(Y ~ X + A + I(A*X), family=quasibinomial(link = "probit"))$coefficients
  res1_1[i,] <- glm(Y ~ X + A + I(A*X), weights = w, family=quasibinomial(link = "probit"))$coefficients
  par <- as.vector(res1_1[i,])
  # Step2
  mu <- cbind(1, X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dnorm(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res1_2[i,] <- glm(Y ~ X + A + I(A*X), weights = w_n, family=quasibinomial(link = "probit"))$coefficients
}
apply(res1_0, 2 , mean)
apply(res1_1, 2 , mean)
apply(res1_2, 2 , mean)

# logistic regression
res1_00 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_00) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res1_3 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_3) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res1_4 <- matrix(NA, nrow = r, ncol = 4)
colnames(res1_4) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi*X) + (X^3) + A*(-1 + 2*X) 
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X, family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res1_00[i,] <- glm(Y ~ X + A + I(A*X), family=quasibinomial(link = "logit"))$coefficients
  res1_3[i,] <- glm(Y ~ X + A + I(A*X), weights = w, family=quasibinomial(link = "logit"))$coefficients
  par <- as.vector(res1_3[i,])
  # Step2
  mu <- cbind(1, X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dsigmoid(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res1_4[i,] <- glm(Y ~ X + A + I(A*X), weights = w_n, family=quasibinomial(link = "logit"))$coefficients
}

apply(res1_00, 2 , mean)
apply(res1_3, 2 , mean)
apply(res1_4, 2 , mean)


library(vioplot)
par(mfrow=c(2,2))
vioplot(res1_0[,3], res1_1[,3], res1_2[,3],
        names=c( "M0","M1", "dWGLM"),xlab= "probit regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")

vioplot(res1_00[,3], res1_3[,3], res1_4[,3],
        names=c( "M0","M1","dWGLM"),xlab= "logistic regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")
vioplot(res1_0[,4], res1_1[,4], res1_2[,4],
        names=c( "M0","M1","dWGLM"),xlab= "probit regression n = 1000",  ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 

vioplot(res1_00[,4], res1_3[,4], res1_4[,4],
        names=c( "M0","M1","dWGLM"), xlab= "logistic regression n = 1000", ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 




##############################################
###### Scenario 2
##############################################
###### Scenario 2
# probit regression
#r <- 1000 # r replications

res.glm0 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm0) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")


res.glm1 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm1) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.glm2 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm2) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2); A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi *X) + (X^3) + A*(-1 + 2*X) 
  p <- pnorm(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X + sin(X) + (X)^2, family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res.glm0[i,] <- glm(Y ~ X + A + I(A*X), family=quasibinomial(link = "probit"))$coefficients
  res.glm1[i,] <- glm(Y ~ X + A + I(A*X), weights = w, family=quasibinomial(link = "probit"))$coefficients
  par <- as.vector(res.glm1[i,])
  # Step2
  mu <- cbind(1, X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dnorm(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res.glm2[i,] <- glm(Y ~ X + A + I(A*X), weights = w_n, family=quasibinomial(link = "probit"))$coefficients
}

apply(res.glm1, 2 , mean)
apply(res.glm2, 2 , mean)

# logistic regression
r <- 1000

res.glm00 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm00) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.glm3 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm3) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

res.glm4 <- matrix(NA, nrow = r, ncol = 4)
colnames(res.glm4) <- c("X_hat", "beta_hat", "psi1_hat","psi2_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi*X) + (X^3) + A*(-1 + 2*X) 
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X + sin(X)+ (X)^2 , family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res.glm00[i,] <- glm(Y ~ X + A + I(A*X), family=quasibinomial(link = "logit"))$coefficients
  res.glm3[i,] <- glm(Y ~ X + A + I(A*X), weights = w, family=quasibinomial(link = "logit"))$coefficients
  par <- as.vector(res.glm3[i,])
  # Step2
  mu <- cbind(1, X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dsigmoid(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res.glm4[i,] <- glm(Y ~ X + A + I(A*X), weights = w_n, family=quasibinomial(link = "logit"))$coefficients
}

apply(res.glm3, 2 , mean)
apply(res.glm4, 2 , mean)


library(vioplot)
par(mfrow=c(2,2))
vioplot(res.glm0[,3], res.glm1[,3], res.glm2[,3],
        names=c( "M0","M1", "dWGLM"),xlab= "probit regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")

vioplot(res.glm00[,3], res.glm3[,3], res.glm4[,3],
        names=c( "M0","M1","dWGLM"),xlab= "logistic regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")
vioplot(res.glm0[,4], res.glm1[,4], res.glm2[,4],
        names=c( "M0","M1","dWGLM"),xlab= "probit regression n = 1000",  ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 

vioplot(res.glm00[,4], res.glm3[,4], res.glm4[,4],
        names=c( "M0","M1","dWGLM"), xlab= "logistic regression n = 1000", ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 



##############################################
###### Scenario 3
##############################################
###### Scenario 3
res3_0 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_0) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")
res3_1 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_1) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")
res3_2 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_2) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2); A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi *X) + (X^3) + A*(-1 + 2*X) 
  p <- pnorm(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X, family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res3_0[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), family=quasibinomial(link = "probit"))$coefficients
  res3_1[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w, family=quasibinomial(link = "probit"))$coefficients
  par <- as.vector(res3_1[i,])
  # Step2
  mu <- cbind(1, log(abs(X)), cos(pi*X), X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dnorm(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res3_2[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w_n, family=quasibinomial(link = "probit"))$coefficients
}

apply(res3_1, 2 , mean)
apply(res3_2, 2 , mean)
# logistic regression


res3_00 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_00) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

res3_3 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_3) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

res3_4 <- matrix(NA, nrow = r, ncol = 6)
colnames(res3_4) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi*X) + (X^3) + A*(-1 + 2*X) 
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X , family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res3_00[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), family=quasibinomial(link = "logit"))$coefficients
  res3_3[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w, family=quasibinomial(link = "logit"))$coefficients
  par <- as.vector(res3_3[i,])
  # Step2
  mu <- cbind(1, log(abs(X)), cos(pi*X), X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dsigmoid(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res3_4[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w_n, family=quasibinomial(link = "logit"))$coefficients
}

apply(res3_3, 2 , mean)
apply(res3_4, 2 , mean)


library(vioplot)
par(mfrow=c(2,2))
vioplot(res3_0[,5], res3_1[,5], res3_2[,5][res3_2[,5] < 20],
        names=c( "M0","M1", "dWGLM"), xlab= "probit regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")

vioplot(res3_00[,5], res3_3[,5][res3_3[,5] < 20], res3_4[,5][res3_4[,5] < 20],
        names=c( "M0","M1","dWGLM"),xlab= "logistic regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")
vioplot(res3_0[,6], res3_1[,6][res3_1[,6] < 20], res3_2[,6][res3_2[,6] < 20],
        names=c( "M0","M1","dWGLM"),xlab= "probit regression n = 1000",  ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 

vioplot(res3_00[,6], res3_3[,6], res3_4[,6],
        names=c( "M0","M1","dWGLM"), xlab= "logistic regression n = 1000", ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 




##############################################
###### Scenario 4
##############################################
###### Scenario 4
res4_0 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_0) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")


res4_1 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_1) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

res4_2 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_2) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2); A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi *X) + (X^3) + A*(-1 + 2*X) 
  p <- pnorm(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X + sin(X) + (X)^2, family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res4_0[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), family=quasibinomial(link = "probit"))$coefficients
  res4_1[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w, family=quasibinomial(link = "probit"))$coefficients
  par <- as.vector(res4_1[i,])
  # Step2
  mu <- cbind(1, log(abs(X)), cos(pi*X), X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dnorm(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res4_2[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w_n, family=quasibinomial(link = "probit"))$coefficients
}




apply(res4_1, 2 , mean)
apply(res4_2, 2 , mean)




# logistic regression
r <- 1000

res4_00 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_00) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

res4_3 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_3) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

res4_4 <- matrix(NA, nrow = r, ncol = 6)
colnames(res4_4 ) <- c("X_hat", "beta_hat",  "beta1_hat",  "beta2_hat", "psi0_hat","psi1_hat")

for (i in 1:r) {
  # generate data:
  X <- runif(n, 0,2);  A <- rbinom(n, 1, sigmoid( -2*X + sin(X) + (X)^2) );
  m <- X + log(abs(X)) + cos(pi*X) + (X^3) + A*(-1 + 2*X) 
  p <- sigmoid(m)
  Y <- rbinom(n, 1, p);
  # Step1
  treat.mod <- glm(A ~ X + sin(X)+ (X^2) , family = 'binomial')
  w = abs(A - fitted(treat.mod))
  res4_00[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), family=quasibinomial(link = "logit"))$coefficients
  res4_3[i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w, family=quasibinomial(link = "logit"))$coefficients
  par <- as.vector(res4_3[i,])
  # Step2
  mu <- cbind(1, log(abs(X)), cos(pi*X), X, (1 - A), (1 - A)*X) %*% par
  k_mu <- dsigmoid(mu)
  w_n = abs(A - fitted(treat.mod)) * k_mu
  # Step3
  res4_4 [i,] <- glm(Y ~ log(abs(X)) + cos(pi*X) + (X^3) + A + I(A*X), weights = w_n, family=quasibinomial(link = "logit"))$coefficients
}

apply(res4_3, 2 , mean)
apply(res4_4 , 2 , mean)


library(vioplot)
par(mfrow=c(2,2))
vioplot(res4_0[,5], res4_1[,5][res4_1[,5] < 20], res4_2[,5][res4_2[,5] < 20],
        names=c( "M0","M1", "dWGLM"), xlab= "probit regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")

vioplot(res4_00[,5], res4_3[,5][res4_3[,5] < 20], res4_4 [,5][res4_4[,5] < 20],
        names=c( "M0","M1","dWGLM"),xlab= "logistic regression n = 1000", ylab = expression(paste( psi[0], " estimates")))
abline(h = -1, col = "red")
vioplot(res4_0[,6], res4_1[,6], res4_2[,6],
        names=c( "M0","M1","dWGLM"),xlab= "probit regression n = 1000",  ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 

vioplot(res4_00[,6], res4_3[,6], res4_4 [,6],
        names=c( "M0","M1","dWGLM"), xlab= "logistic regression n = 1000", ylab = expression(paste( psi[1], " estimates")))
abline(h = 2, col = "red") 




##############################################
###### Final conclusion
##############################################


##### problit
par(mfrow=c(3,2))
vioplot( res1_0[,3], res.glm0[,3], res3_0[,5], res4_0[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "ME0\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_0[,4], res.glm0[,4], res3_0[,6], res4_0[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "ME0\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")
vioplot( res1_1[,3], res.glm1[,3], res3_1[,5], res4_1[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "ME1\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_1[,4], res.glm1[,4], res3_1[,6], res4_1[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "ME1\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")
vioplot( res1_2[,3], res.glm2[,3], res3_2[,5], res4_2[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "dWGLM\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_2[,4], res.glm2[,4], res3_2[,6], res4_2[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Probit \ (n = 1000)", ylab = expression(paste( "dWGLM\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")



##### logistic
par(mfrow=c(3,2))
vioplot( res1_00[,3], res.glm00[,3], res3_00[,5], res4_00[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "ME0\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_00[,4], res.glm00[,4], res3_00[,6], res4_00[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "ME0\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")
vioplot( res1_3[,3], res.glm3[,3], res3_3[,5], res4_3[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "ME1\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_3[,4], res.glm3[,4], res3_3[,6], res4_3[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "ME1\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")
vioplot( res1_4[,3], res.glm4[,3], res3_4[,5], res4_4[,5], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "dWGLM\ ", psi[0], " estimates")) )
abline(h = -1, col = "red")
vioplot( res1_4[,4], res.glm4[,4], res3_4[,6], res4_4[,6], names=c("1", "2", "3", "4"),  
         xlab= "Scenarios \n Logistic \ (n = 1000)", ylab = expression(paste( "dWGLM\ ", psi[1], " estimates")) )
abline(h = 2, col = "red")
