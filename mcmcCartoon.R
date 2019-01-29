#################
# Simple example of MCMC sampling
#################
set.seed(2)
library(fields)
cbPalette <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73")
colorTable <- designer.colors(500, cbPalette,
                              x = 10 * seq(1:4))

#########
# first, let's build a function that generates random numbers from a bivariate normal distribution

rho <- 0.9
x <- mvtnorm::rmvnorm(10000, sigma = matrix(c(1,sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2))

plot(x, col = 1:10000, cex = 0.1)
kde <- MASS::kde2d(x[,1], x[,2], n = 200)
image(kde, col = colorTable)

a <- x[1:35,1:2]
plot(a, ylim = c(-5,5), xlim = c(-5,5), col = "blue", cex = 2, pch = 19)
mixtools::ellipse(mu = c(0, 0),
                  sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                  alpha = .05, npoints = 250, newplot = FALSE, draw = TRUE)
mixtools::ellipse(mu = c(0, 0),
                  sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                  alpha = .01, npoints = 250, newplot = FALSE, draw = TRUE)
mixtools::ellipse(mu = c(0, 0),
                  sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2),
                  alpha = .10, npoints = 250, newplot = FALSE, draw = TRUE)

## code moved to pRoloc