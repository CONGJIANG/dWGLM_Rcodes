############
set.seed(1)
x <- runif(100, 0,2)
y <- x + log(abs(x)) + cos(pi*x) + x^3 + rnorm(100)
plot(y ~ x, main = "Treatment-free function and its linear regression",  ylab = "f(x)")
curve(x + log(abs(x)) + cos(pi*x) + x^3, add = TRUE, col = "red",lwd=2)
l1<-predict(lm(y ~ x))
lines(x[j],l1[j],col="blue",lwd=2)
legend(0.3, 8, legend=c("Treatment-free function", "Linear regression"),
       col=c("red", "blue"), lty=1, cex=1,
       box.lty=2)
