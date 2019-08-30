

set.seed(1)
x <- rbeta(n = 5000, shape1 = 2, shape2 = 10)
h <- hist(x, breaks=50, plot = FALSE)
cuts <- cut(h$breaks, c(0,0.01, 0.5, 1))
plot(h, col=cuts, xlim = c(0,1), ylim = c(0, 500))

gg <- ggplot(data.frame(x), aes(x)) + geom_histogram(fill = "darkgreen", col = "black") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.border = element_rect(colour = "black", fill = NA, size = 1),
           plot.title = element_text(hjust = 0.5, size = 20),
           legend.text=element_text(size = 14)) + 
      ggtitle(label = "Histogram of simulations from the prior") + xlim(c(-0.05, 1))
gg

out <- mcmc_get_outliers(e14Tagm_converged_pooled)

propout <- out[[1]]/nrow(unknownMSnSet(E14TG2aR))
hist(propout, add = TRUE, col = "blue")

mu_post <- mean(propout)
sigma_post <- sd(propout)

z_post <- abs((mu_post - mean(x))/sigma_post)

shrink_post <- 1 - sigma_post^2/var(x)
plot(rep(shrink_post, 5000), z_post, xlim = c(0,1))

abs(mean(propout) - mean(x)/sigma_post)


getMarkers(E14TG2aR)

N <- nrow(E14TG2aR)

meanComp <- tabulate(fData(E14TG2aR)[, "markers"])[1:10] + 1
meanComp <- meanComp/sum(meanComp)

varComp <- meanComp * (1 - meanComp)/sum(tabulate(fData(E14TG2aR)[, "markers"])[1:10] + 2)

cmptable <- apply(tmp@Component, 2, tabulate)

abs(mean(cmptable[5,]/N) - meanComp[5])/sd(cmptable[5,]/N)


varCmptable <- apply(cmptable/N,1, var)
sdCmptable <- apply(cmptable/N, 1, sd)

z_post_dir <- abs((rowMeans(cmptable/N) - meanComp)/sdCmptable)
shrink_post_dir <- 1 - (varCmptable/varComp)^2


mixing_posterior_check <- function(object,
                                   params,
                                   priors,
                                   fcol = "markers") {

K <- length(priors$beta0)
N <- nrow(object)

# Compute prior mean and variance
tallyComp <- c(table(fData(object)[, fcol])[1:K] + priors$beta0)
meanComp <- tallyComp/sum(tallyComp)  
varComp <- meanComp * (1 - meanComp)/sum(tallyComp + 1)

# Compute posterior quantities
cmptable <- apply(params@Component, 2, tabulate)
varCmptable <- apply(cmptable/N, 1, var)
sdCmptable <- apply(cmptable/N, 1, sd)

post_z_mixing <- abs((rowMeans(cmptable/N) - meanComp)/sdCmptable)
post_shrink_mixing <- 1 - (varCmptable/varComp)^2

post_check_mixing <- data.frame(x = post_shrink_mixing, y = post_z_mixing, getMarkerClasses(object))
colnames(post_check_mixing) <- c("posteriorShrinkage", "posteriorZscore", "organelle")
rownames(post_check_mixing) <- getMarkerClasses(object)

cols <- getStockcol()[1:K] # get appropriate colours
gg <- ggplot(post_check_mixing, aes(x = posteriorShrinkage, y = posteriorZscore, color = organelle)) + geom_point(size = 5)
gg <- gg + theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text=element_text(size = 14)) + scale_color_manual(values = cols) + 
  ggtitle(label = "Posterior predictive check for mixing proportions") + xlim(c(0, 1))

return(gg)
}


