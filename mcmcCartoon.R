#################
# Simple example of MCMC sampling
#################
library(fields)
cbPalette <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73")
colorTable<- designer.colors(500, cbPalette, 
                             x = 10*seq(1:4))

#########
# first, let's build a function that generates random numbers from a bivariate normal distribution

rho <- 0.9
x <- mvtnorm::rmvnorm(10000, sigma = matrix(c(1,sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2))

plot(x, col = 1:10000, cex = 0.1)
kde <- MASS::kde2d(x[,1], x[,2], n = 200)
image(kde, col = colorTable)

a <- x[1:35,1:2]
plot(a, ylim = c(-5,5), xlim = c(-5,5), col = "blue", cex = 2, pch = 19)
ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .05, npoints = 250, newplot = FALSE,
        draw = TRUE)
ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .01, npoints = 250, newplot = FALSE,
        draw = TRUE)
ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .10, npoints = 250, newplot = FALSE,
        draw = TRUE)
        

metropolisHastings <- function (n, rho) {   
  mat <- matrix(ncol = 2, nrow = n)   # matrix for storing the random samples
  x <- 4   # initial values for all parameters
  y <- 4
  cov <- matrix(c(1,sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2)
  prev <- dmvnorm(c(x, y), mu = c(0,0), sigma = cov) 
  mat[1, ] <- c(x, y)        # initialize the markov chain
  counter <- 1
  while(counter<=n) {
    propx <- x + rnorm(1, 0, 0.1)
    propy <- y + rnorm(1, 0, 0.1)
    
    newprob <- dmvnorm(c(propx, propy), sigma = cov) 
    ratio <- newprob/prev 
    
    prob.accept <- min(1, ratio)     # ap
    rand <- runif(1)
    if (rand <= prob.accept) {
      x <- propx
      y <- propy
      mat[counter, ] <- c(x, y)    # store this in the storage array 
      counter=counter+1
      prev <- newprob    # get ready for the next iteration
    }
    
  }
  return(mat)
}

met <- metropolisHastings(10000, 0.9)


plotMet <- function(r) {
  a <- x[1:35,1:2]
  plot(a, ylim = c(-5,5), xlim = c(-5,5), xlab = "Channel 1", ylab = "Channel 2", col = "blue", cex = 2, pch = 19, main = paste0("Iteration ", r))
  ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .05, npoints = 250, newplot = FALSE,
          draw = TRUE)
  ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .01, npoints = 250, newplot = FALSE,
          draw = TRUE)
  ellipse(mu = c(0, 0), sigma = matrix(c(1, sqrt(1 - rho^2),sqrt(1 - rho^2),1), ncol = 2), alpha = .10, npoints = 250, newplot = FALSE,
          draw = TRUE)
  
  points(met[100 * c(1:r) - 99,], type = "b", col = "red", cex = 2, lwd = 3, pch = 19)
}
par(mfrow = c(2,2))
plotMet(2)
plotMet(3)
plotMet(5)
plotMet(30)