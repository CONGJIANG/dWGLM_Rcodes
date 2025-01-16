###################################
###################################
###################################
# Clear the entire environment
rm(list = ls())
############################################################################
# Author: Cong Jiang (c55jiang@uwaterloo.ca)
# Article:Novel Robust Dynamic Treatment Regimen Estimation for Discrete Outcome: An application to smoking cessation using e-cigarettes
# install.packages("e1071")
#install.packages("drgee")
#install.packages("rgenoud")

library("e1071")
library("rgenoud")
library("glmnet")
library("drgee")
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}
###################################

#1. Calculate ratios
calculate_ratio <- function(X2, res) {
  # Compute predictions
  pred_res <- as.numeric((cbind(1, X2) %*% res) > 0)
  comp_vector <- as.numeric((cbind(1, X2) %*% c(4, -1)) > 0)
  # Calculate and return the ratio of matches
  match_count <- mean(pred_res == comp_vector)
  return(match_count)
}

#2. Calculate outcomes
calculate_outcome <- function(X1, X2, res) {
  # Calculate the linear prediction
  linear_pred <- X1 + exp(X2) + cos(pi * X1) + 1.5 * (X2^3)
  # Adjust linear prediction based on the sign of the linear combination with `res`
  adjustment <- ifelse(cbind(1, X2) %*% res > 0, cbind(1, X2) %*% res, 0)
  # Calculate probabilities using the sigmoid
  probs <- plogis(linear_pred + adjustment)
  # Compute and return the mean of these prob
  return(mean(probs))
}

#3. Iterative Weighted GLMs (IWGLMs)
IWGLMs <- function(X1, X2, A, Y, trt_frm, out_frm, max_iter = 1000, convergence_threshold = 1e-8) {
  # Fit the treatment model
  treatment_model <- glm(trt_frm, family = 'binomial', data = data.frame(X1 = X1, A = A))
  
  # Initialize weights using absolute difference from predicted treatment probability
  initial_weights <- abs(A - fitted(treatment_model))
  
  # Initial fit of the outcome model
  initial_fit <- glm(out_frm, data = data.frame(X1 = X1, X2 = X2, A = A, Y = Y), family = quasibinomial(link = "logit"))
  res.glm00 <- initial_fit$coef
  res.glm00mod <- initial_fit
  
  # Initialize convergence variables
  coefficients_diff <- Inf
  iteration <- 0
  
  # Iterative fitting process
  while (coefficients_diff > convergence_threshold && iteration < max_iter) {
    # Predict the outcome based on the current model
    mu <- predict(res.glm00mod, newdata = data.frame(X1 = X1, X2 = X2, A = 1 - A, Y = Y), type = "response")
    
    # Calculate the derivative of the sigmoid function at mu
    k_mu <- dsigmoid(mu)
    
    # Update weights based on the predicted outcome and treatment model
    updated_weights <- abs(A - fitted(treatment_model)) * k_mu
    
    # Fit the model with the updated weights
    updated_fit <- glm(out_frm, data = data.frame(X1 = X1, X2 = X2, A = A, Y = Y, w = updated_weights), weights = w, family = quasibinomial(link = "logit"))
    res.glm01 <- updated_fit$coef
    res.glm01mod <- updated_fit
    
    # Calculate the difference between the current and previous coefficients
    coefficients_diff <- max(abs(res.glm01 - res.glm00))
    
    # Update the previous coefficients and model for the next iteration
    res.glm00 <- res.glm01
    res.glm00mod <- res.glm01mod
    
    # Increment the iteration counter
    iteration <- iteration + 1
  }
  
  # Print convergence information
  if (iteration < max_iter) {
    cat("Converged after", iteration, "iterations.\n")
  } else {
    cat("Reached maximum iterations without convergence.\n")
  }
  
  # Return the final coefficients and model
  return(list(final_coefficients = res.glm01, final_model = res.glm01mod))
}

n = 5000
X1 <- rnorm(n, 0, 2)
X2 <- runif(n, -1, 1)
X <- cbind(X1, X2)  # Combine X1 and X2 into a matrix

A <- rbinom(n, 1, sigmoid(-2 * X1 + sin(X1) + (X2)^2))
m <- X1 + 5*exp(X1) + 3*cos(pi * X2) + 5 * (X1^3) + A * (4 - 1 * X2)
p <- sigmoid(m)
Y <- rbinom(n, 1, p)

simdata <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y)
trt_frm <- A ~ X1 + sin(X1) + X1^2
#out_frm <- Y ~ X + A + I(A * X)
out_frm <- Y ~ X1 + exp(X1) + cos(pi*X1) + X1^3+ A + I(A * X2)
(res <- IWGLMs(X1,X2, A, Y, trt_frm, out_frm))

#4. AIPWE estimation (example)
# A Robust Method for Estimating Optimal Treatment Regimes
# by Baqun Zhang,* Anastasios A. Tsiatis, Eric B. Laber, and Marie Davidian
AIPWE <- function(psi, X1, X2, A, Y, trt_frm, out_frm) {
  n <- length(Y)
  g_psi <- ifelse(cbind(1, X2) %*% psi > 0, 1, 0)
  
  treat.mod <- glm(trt_frm, family = 'binomial', data = data.frame(X1=X1, A=A))
  fitted_vals <- fitted(treat.mod)
  
  C_psi <- A * g_psi + (1 - A) * (1 - g_psi)
  pi_psi <- fitted_vals * g_psi + (1 - fitted_vals) * (1 - g_psi)
  
  mu <- glm(out_frm, family = quasibinomial(link = "logit"), data = data.frame(X1=X1,X2=X2, A=A, Y=Y))
  mu1 <- predict(mu, newdata = data.frame(A = 1, X1 = X1, X2 = X2), type = "response")
  mu0 <- predict(mu, newdata = data.frame(A = 0, X1 = X1, X2 = X2), type = "response")
  
  mu_psi <- mu1 * g_psi + mu0 * (1 - g_psi)
  AIPWE <- (C_psi * Y) / pi_psi - mu_psi * (C_psi - pi_psi) / pi_psi
  
  return(mean(AIPWE))
}

train_data <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y)
# Run genoud
genoud(AIPWE, nvars = 2, max = TRUE, pop.size = 200, max.generations = 3,
       wait.generations = 5, gradient.check = FALSE, X1 = train_data$X1, X2 = train_data$X2, A = train_data$A, Y = train_data$Y, trt_frm = trt_frm, out_frm = out_frm)

#5. Doubly Robust Estimation using GEE (drgee package needed)
new_frm <- as.formula(paste("Y ~", paste(attr(terms(out_frm), "term.labels")[-((length(attr(terms(out_frm), "term.labels"))-1):length(attr(terms(out_frm), "term.labels")))], collapse = " + ")))
(dr.est <- drgee(oformula = formula(new_frm), eformula = formula(trt_frm),
                 iaformula = formula(~X2), olink = "logit", elink = "logit",
                 data = simdata, estimation.method = "dr"))

run_simulation <- function(n, r, trt_frm, out_frm) {
  # Placeholder for storing results
  res.glm00 <- matrix(NA, nrow = r, ncol = 2)  # Adjust size as needed
  res.glm01 <- matrix(NA, nrow = r, ncol = 2)  # Adjust size as needed
  res.IWGLMs <- matrix(NA, nrow = r, ncol = 2)  # Adjust size as needed
  res.AIPWE <- matrix(NA, nrow = r, ncol = 2)  # Adjust size as needed
  res.EC <- matrix(NA, nrow = r, ncol = 2)  # Adjust size as needed
  
  ratio1 <- numeric(r); ratio2 <- numeric(r); ratio3 <- numeric(r); ratio4 <- numeric(r); ratio5 <- numeric(r)
  out1 <- numeric(r); out2 <- numeric(r); out3 <- numeric(r); out4 <- numeric(r); out5 <- numeric(r)
  
  # Main simulation loop
  for (i in 1:r) {
    # Generate data:
    X1 <- rnorm(n, 0, 2)
    X2 <- runif(n, -1, 1)
    X <- cbind(X1, X2)  # Combine X1 and X2 into a matrix
    
    A <- rbinom(n, 1, sigmoid(-2 * X1 + sin(X1) + (X2)^2))
    m <- X1 + 5*exp(X1) + 3*cos(pi * X2) + 5 * (X1^3) + A * (4 - 1 * X2)
    p <- sigmoid(m)
    Y <- rbinom(n, 1, p)
    
    simdata <- data.frame(X1 = X1, X2 = X2, A = A, Y = Y)
    
    # Split data into training and testing sets
    train_indices <- sample(1:n, size = floor(0.35 * n))
    test_indices <- setdiff(1:n, train_indices)
    
    train_data <- simdata[train_indices, ]
    test_data <- simdata[test_indices, ]
    
    # Step 1: Fit models using the training set
    mu <- glm(out_frm, family = quasibinomial(link = "logit"), data = train_data)
    vec <- mu$coefficients
    res.glm00[i, ] <- vec[(length(vec)-1):length(vec)]
    
    # Calculate ratio and out using the testing set
    test_X1 <- test_data$X1; test_X2 <- test_data$X2; test_A <- test_data$A; test_Y <- test_data$Y
    
    # Modify these functions to handle the matrix X if needed
    ratio1[i] <- calculate_ratio(test_X2, res.glm00[i, ])
    out1[i] <- calculate_outcome(test_X1, test_X2, res.glm00[i, ])

    # Step 2: Outcome model with weights
    treat.mod <- glm(trt_frm, family = 'binomial', data = train_data)
    train_data$w <- abs(train_data$A - fitted(treat.mod))
    mu1 <- glm(out_frm, family = quasibinomial(link = "logit"), weights = w, data = train_data)
    vec1 <- mu1$coefficients
    res.glm01[i, ] <- vec1[(length(vec1)-1):length(vec1)]
    ratio2[i] <- calculate_ratio(test_X2, res.glm01[i, ])
    out2[i] <- calculate_outcome(test_X1, test_X2, res.glm01[i, ])
    
    # Step 4: IWGLMs estimation
    mod3 <- IWGLMs(train_data$X1, train_data$X2, train_data$A, train_data$Y, trt_frm, out_frm)
    vec3 <- mod3$final_coefficients
    res.IWGLMs[i, ] <- vec3[(length(vec3)-1):length(vec3)]
    ratio3[i] <- calculate_ratio(test_X1, res.IWGLMs[i, ])
    out3[i] <- calculate_outcome(test_X1, test_X2, res.IWGLMs[i, ])
    
    # Step 5: AIPWE estimation (need genoud function)
    AIPWE.res <- genoud(AIPWE, nvars = 2, max = TRUE, pop.size = 200, max.generations = 3,
                        wait.generations = 5, gradient.check = FALSE, X1 = train_data$X1, X2 = train_data$X2, A = train_data$A, Y = train_data$Y, trt_frm = trt_frm, out_frm = out_frm)
    res.AIPWE[i, ] <- AIPWE.res$par
    ratio4[i] <- calculate_ratio(test_X1, res.AIPWE[i, ])
    out4[i] <- calculate_outcome(test_X1, test_X2, res.AIPWE[i, ])
    
    # Step 6: Doubly Robust Estimation using GEE (drgee package needed)
    # adjust the new outcome formula (first couple of terms)
    new_frm <- as.formula(paste("Y ~", paste(attr(terms(out_frm), "term.labels")[-((length(attr(terms(out_frm), "term.labels"))-1):length(attr(terms(out_frm), "term.labels")))], collapse = " + ")))
    
    dr.est <- drgee(oformula = formula(new_frm), eformula = formula(trt_frm),
                    iaformula = formula(~ X2), olink = "logit", elink = "logit",
                    data = train_data, estimation.method = "dr")
    EC.res <- summary(dr.est)
    res.EC[i, ] <- as.vector(EC.res$coefficients[, 1])
    ratio5[i] <- calculate_ratio(test_X1, res.EC[i, ][1:2])
    out5[i] <- calculate_outcome(test_X1, test_X2, res.EC[i, ][1:2])
  }
  
  # Return all results as a list
  return(list(res_glm00 = res.glm00, res_glm01 = res.glm01, res_IWGLMs = res.IWGLMs,
              res_AIPWE = res.AIPWE, res_EC = res.EC, ratio1 = mean(ratio1), ratio2 = mean(ratio2),
              ratio3 = mean(ratio3), ratio4 = mean(ratio4), ratio5 = mean(ratio5), out1 = mean(out1), out2 = mean(out2),
              out3 = mean(out3), out4 = mean(out4), out5 = mean(out5)))
}

trt_frm <- A ~ X1 + sin(X1) + X1^2
out_frm <- Y ~ X1 + X2 + A + I(A * X2)
#out_frm <- Y ~ X + exp(X) + cos(pi*X) + X^3 + A + I(A * X)
(res <- run_simulation(5000, 15, trt_frm = trt_frm , out_frm = out_frm))

apply(res$res_EC, 2, mean)
apply(res$res_AIPWE, 2, mean)
apply(res$res_IWGLMs, 2, mean)
apply(res$res_glm00, 2, mean)
apply(res$res_glm01, 2, mean)



# p(d^{opt} = 1) \in [0.3, 0.7]
# blip values should be bigger, i.e., V(d^{opt}) - V(1 - d^{opt}) should be bigger
# P(Y = 1) big for the purpose, treatment-free should be bigger. 



