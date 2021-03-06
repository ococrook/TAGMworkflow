
# data("tan2009r1")
# tanres <- tagmMcmcTrain(object = tan2009r1)
# tanres <- tagmMcmcProcess(tanres)
# tan2009r1 <- tagmMcmcPredict(object = tan2009r1, params = tanres, probJoint = TRUE)
# 
# mydata <- E14TG2aR#tan2009r1
# myparams <- chains(e14Tagm_converged_pooled)[[1]]#tanres@chains@chains[[1]]
# 
# myparams2 <- chains(mcmc_pool_chains(tanres))[[1]]
# priors <- tanres@priors

nicheMeans2D <- function(object,
                         params,
                         priors,
                         dims = c(1, 2),
                         fcol = "markers",
                         aspect = 0.5){

# Make coordinates
.pca <- plot2D(object, dims = dims, plot = FALSE)
d <- 2
mcmcIter <- seq.int(by = 10, from = 1, to = params@n)
iter_len <- length(mcmcIter)
orgMeans <- array(NA, c(params@K, iter_len, ncol(object)))

  # Compute means at each iteration of MCMC algorithm
  idx <- rownames(params@Component)
  mydata_sub <- object[idx, ]
  for (i in seq.int(mcmcIter)) {
    for (j in seq.int(params@K)) {
      idxj <- (params@Component[, mcmcIter[i]] ==  j)
      
      mj <- colMeans(rbind(exprs(mydata_sub)[idxj, ],
                                   exprs(object[fData(object)[, fcol] == getMarkerClasses(object)[j], ])) )
      nj <- sum(fData(object)[, "tagm.mcmc.allocation"] == getMarkerClasses(object)[j])
      lambdaj <- params@ComponentParam@lambdak[j]
      orgMeans[j, i, ] <- (nj * mj + priors$lambda0 * priors$mu0)/lambdaj
    }
  }

  # Get data to align to 
  M <- matrix(NA, nrow = params@K, ncol = ncol(object))
  rownames(M) <- getMarkerClasses(object)
  for (j in seq.int(params@K)) {
    M[j, ] <- colMeans(exprs(object)[fData(object)[, fcol] == getMarkerClasses(object)[j],])
  }
  
  # Create coordinates 
  coords <- matrix(NA, nrow = params@K, ncol = 2)
  res0 <- prcomp(M, scale. = TRUE)
  eigs <- colnames(.pca)

  # Create inital dataset, computing average location in PCA plot
  Y0 <- matrix(NA, nrow = params@K, ncol = ncol(object))
  for (j in seq.int(params@K)) {
    Y0[j, ] <- colMeans(.pca[fData(object)[, fcol] == getMarkerClasses(object)[j], seq_len(d)])
  }
  Y0.df <- data.frame(organelle = getMarkerClasses(object), Y0)
  
  # Repeat for different monte-carlo samples
  Y.lst <- list()
  for ( i in seq.int(iter_len)) {
    res <- prcomp(orgMeans[, i, ], scale = TRUE, center = TRUE)
    res.proc <- vegan::procrustes(Y0, res$x[, dims])
    Y.df <- data.frame(organelle = getMarkerClasses(object), res.proc$Yrot) 
    Y.lst[[i]] <- Y.df
  }
  
  # create long data formats
  names(Y.lst) <- seq_len(iter_len)
  Y.lst.df <- plyr::ldply(Y.lst, .fun = function(x) x, .id = "mcmcIter")
  table(Y.lst.df$organelle)
  cols <- getStockcol()[1:params@K] # appropriate colours
  
  ## ggplot
  gg <- ggplot(
    data = Y.lst.df %>%
      mutate(organelle = factor(organelle)),
    aes(x = X1, y = X2, color = organelle)) +
    coord_fixed() + 
    geom_density2d(contour = TRUE) + 
    geom_point(alpha = 0.3) + 
    xlab(paste0(eigs[1])) + 
    ylab(paste0(eigs[2])) +
    theme(legend.position = "right", 
          text = element_text(size = 12)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    theme_minimal() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = aspect,
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.text = element_text(size = 14)) + 
    ggtitle(label = "Uncertainty in mean of organelle localisation")
  gg

  return(gg)  
}


# par(mfrow = c(1,1))
# plot2D(E14TG2aR)


## Probability variation plot
require(akima)
require(fields)


spatial2D <- function(object,
                      dims = c(1, 2),
                      cov.function = wendland.cov,
                      theta = 1,
                      derivative = 2,
                      k = 1,
                      breaks = c(0.99, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7),
                      aspect = 0.5) {

  # generate pca plot and create data from with probabilties
  .pca <- plot2D(object, dims = dims, plot = FALSE)
  probs <- data.frame(x = .pca[, 1], y = .pca[, 2], mcmc.prob = fData(object)$tagm.mcmc.joint)
  colnames(probs) <- c(c("x", "y"), getMarkerClasses(object))
  eigs <- colnames(.pca)
  
  # put data in appropriate long format
  probs.lst <- list()
  for(j in getMarkerClasses(object)) {
    probs.lst[[j]] <- probs[, c("x", "y", j)]
    colnames(probs.lst[[j]]) <- c("x", "y", "probability")
  }
  probs.lst.df <- plyr::ldply(probs.lst, .fun = function(x) x, .id = "organelle")

  # Create storage
  coords <- list()
  locations <- list()
  df <- list()
  # Create appropriate spatial grid
  for (j in getMarkerClasses(object)) {
    idxOrg <- c(probs.lst.df$organelle == j)
    coords[[j]] <- akima::interp(x = probs.lst.df$x[idxOrg],
                                 y = probs.lst.df$y[idxOrg],
                                 z = probs.lst.df$probability[idxOrg],
                                 extrap=FALSE, linear = TRUE, duplicate = TRUE) # interpolate onto appropriate grid
    coords[[j]]$z[is.na(coords[[j]]$z)] <- 0 # NaNs beyond data set to 0
    locations[[j]] <- cbind(rep(coords[[j]]$x, 40), rep(coords[[j]]$y, each = 40)) # get grid
    smoothedprobs <- fields::smooth.2d(coords[[j]]$z, x = locations[[j]],
                               cov.function = cov.function,
                               theta = theta,
                               derivative = derivative, k = k) # spatial smoothing of probabilities
    # normalisation and formatting
    zmat <- matrix(smoothedprobs$z, ncol = 1)
    zmat <- zmat/max(zmat)
    df[[j]] <- data.frame(x = rep(smoothedprobs$x, 64), y = rep(smoothedprobs$y, each = 64), z = zmat)
  }
  # format data
  df.lst <- plyr::ldply(df, .fun = function(x) x, .id = "organelle") 
  df.lst <- df.lst %>%
    mutate(organelle = factor(organelle)) 
  K <- length(getMarkerClasses(object))
  cols <- getStockcol()[1:K] # get appropriate colours
  
  
gg <- ggplot(
  data = df.lst,
  aes(x = x, y = y, z = z, color = organelle)) +
  coord_fixed() + 
  geom_contour(breaks = breaks, size = 1.2, aes(alpha = stat(level))) + 
  geom_point(alpha = 0) + 
  xlab(paste0(eigs[1])) + 
  ylab(paste0(eigs[2])) +
  scale_alpha(guide = "none") + 
  theme(legend.position = "right", 
        text = element_text(size = 12)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        aspect.ratio = aspect,
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        plot.title = element_text(hjust = 0.5, size = 20),
        legend.text=element_text(size = 14)) +
  ggtitle(label = "Spatial variation of localisation probabilities") 

  return(gg)
}



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
    ggtitle(label = "") + xlim(c(0, 1))
  
  return(gg)
}


#require(patchwork)
#gg1 + gg2 + gg3 + gg4
