birth_death <- function(S0 = 1, B = 0.0067, D = 0.01) {
  S <- S0
  t <- 0
  times <- c(0)
  pop <- c()
  while (t <= 4000) {
    dt <- rexp(1, ((B + D)*S))
    S <- S + 1
    t <- t + dt
    times <- c(times, t)
  }
  if (times[length(times)] > 4000) {
    times <- times[1:(length(times) - 1)]
  }
  if (t %% 5 == 0){
    pop <- c(pop, S)
  }
  return(times)
}

library(ggplot2)

B <- 0.0067
D <- B * 0.1
# D <- abs(rnorm(4000, mean = 0, sd = 0.002))
t <- seq(1,4000, by = 1)
S = 1 * exp((B-D) * t)
Pr = (1 - D/B) / (1 - (D/B) * exp(-(B - D) * t))
ggplot(NULL, aes(x = t, y = Pr)) + geom_point() 
# + scale_y_log10()


# package testing
library(DOBAD)

B <- 0.0067
D <- 0
bd1 <- birth.death.simulant(t = 1000, X0 = 20, lambda = B, mu = D, nu = 0)
N1 <- bd1@states[length(bd1@states)]
bd2 <- birth.death.simulant(t = 1000, X0 = N1, lambda = B, mu = D, nu = 0)
N2 <- bd2@states[length(bd2@states)]

# probability
beta <- (exp((B - D) * 4000) - 1)/(exp((B - D) * 4000) - (D/B))
(exp((B - D) * 4000) - 1)
(exp((B - D) * 4000) - (D/B))
(1-beta) * beta^(10^12 -1)

# predicting diversification
n <- 11000
N0 <- 1
t <- 4000
ext <- 0.9
# (log(n))/t
r <- (1/t) * log(n * (1 - ext) + ext)
beta <- (exp(r * t) - 1)/(exp(r*t) - ext)
alpha <- ext * beta
Nmean <- (N0 * exp(r *t))/(1 - alpha^N0)
N0 * exp(r * t)
