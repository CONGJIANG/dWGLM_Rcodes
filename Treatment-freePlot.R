set.seed(1)

# Generate data for first plot
x1 <- runif(100, 0, 2)
y1 <- x1 + log(abs(x1)) + cos(pi * x1) + x1^3 + rnorm(100)

# Generate data for second plot
x2 <- runif(100, -0.5, 0.5)
y2 <- x2 + exp(x2) + cos(pi * x2) + x2^3 + rnorm(100)

# Set up the plotting area with subplots
par(mfrow=c(1,2), oma=c(0,0,2,0))

# Plot first figure
plot(y1 ~ x1, ylab = "f(x) = x + log(abs(x)) + cos(pi * x) + x^3", xlab = "x")
curve(x + log(abs(x)) + cos(pi * x) + x^3, add = TRUE, col = "red", lwd = 2)
l1 <- predict(lm(y1 ~ x1))
lines(x1, l1, col = "blue", lwd = 2)
legend(0.1, 13, legend = c("Treatment-free function", "Linear approximation"),
       col = c("red", "blue"), lty = 1, cex = 1, box.lty = 2)

# Plot second figure
plot(y2 ~ x2, ylab = "f(x) = x + exp(x) + cos(pi * x) + x^3", xlab = "x")
curve(x + exp(x) + cos(pi * x) + x^3, add = TRUE, col = "red", lwd = 2)
l2 <- predict(lm(y2 ~ x2))
lines(x2, l2, col = "blue", lwd = 2)
legend(-0.39, 4.9, legend = c("Treatment-free function", "Linear regression"),
       col = c("red", "blue"), lty = 1, cex = 1, box.lty = 2)

# Add main title in the middle of the top
mtext("Weibull distribution", line=2, side=3, outer=TRUE, cex=2)
mtext("Treatment-free function and its linear approximation", side = 3, line = - 2, outer = TRUE, cex=1.5)



