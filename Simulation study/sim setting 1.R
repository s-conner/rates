### Simulation Setting 1
### PEx Rate manuscript
### Sarah Conner 12/21/22

library(MASS)

set.seed(123)
niter <- 2000
nboot <- 1000
beta0 <- 1.6
beta1 <- -.03 # continuous covariate ppFEV1
beta2 <- -1 # log RR of treatment
samp_sizes <- c(500)
res <- NULL


# --------- Analysis function --------- 

analysis <- function(dat){
  
  # --------- Observed rates --------- 
  
  obs.trt <- sum(dat$y[dat$trt==1])/sum(dat$t[dat$trt==1])
  obs.ctrl <- sum(dat$y[dat$trt==0])/sum(dat$t[dat$trt==0])
  obs.logrr <- log(obs.trt/obs.ctrl)
  
  
  # --------- Fit Models --------- 
  
  # NB Models 
  tryCatch({nb.x <- glm.nb(y ~ x + trt + offset(log(t)), data=dat); summ.nb.x <- summary(nb.x)},
           error=function(e) {skip.e1 <<- TRUE}, 
           warning=function(w) {skip.w1 <<- TRUE})
  nb <- glm.nb(y ~ trt + offset(log(t)), data=dat)
  summ.nb <- summary(nb)

  
  
  # --------- Obtain RRs and rates --------- 
  
  # Prediction profiles
  trt.cond <- data.frame(x=mean(dat$x), trt=1, t=1)
  ctrl.cond <- data.frame(x=mean(dat$x), trt=0, t=1)
  trt.marg <- data.frame(x=dat$x, trt=rep(1,n), t=dat$t)
  ctrl.marg <- data.frame(x=dat$x, trt=rep(0,n), t=dat$t)

  
  # NB conditional and marginal rates, unadjusted model
  nb.trt.cond <- predict(nb, newdata=trt.cond, type='response')
  nb.logtrtse.cond <- predict(nb, newdata=trt.cond, se=TRUE)$se.fit
  nb.ctrl.cond <- predict(nb, newdata=ctrl.cond, type='response')
  nb.logctrlse.cond <- predict(nb, newdata=ctrl.cond, se=TRUE)$se.fit
  
  nb.logrr.cond <- coef(nb)[2]
  nb.logrrse.cond <- summ.nb$coef[2,2]
  nb.logrrcil.cond <- nb.logrr.cond-1.96*nb.logrrse.cond
  nb.logrrciu.cond <- nb.logrr.cond+1.96*nb.logrrse.cond
  
  nb.trt.marg <- mean(predict(nb, newdata=trt.marg, type='response'))/mean(dat$t)
  nb.ctrl.marg <- mean(predict(nb, newdata=ctrl.marg, type='response'))/mean(dat$t)
  nb.logrr.marg <- log(nb.trt.marg/nb.ctrl.marg)
  
  
  # NB conditional and marginal rates, adjusted model - only predict if converged
  if(!skip.e1 & !skip.w1) {
    
    nb.x.trt.cond <- predict(nb.x, newdata=trt.cond, type='response') 
    nb.x.logtrtse.cond <- predict(nb.x, newdata=trt.cond, se=TRUE)$se.fit
    nb.x.ctrl.cond  <- predict(nb.x, newdata=ctrl.cond, type='response')
    nb.x.logctrlse.cond <- predict(nb.x, newdata=ctrl.cond, se=TRUE)$se.fit
    
    nb.x.logrr.cond <- coef(nb.x)[3]
    nb.x.logrrse.cond <- summ.nb.x$coef[3,2]
    nb.x.logrrcil.cond <- nb.x.logrr.cond-1.96*nb.x.logrrse.cond
    nb.x.logrrciu.cond <- nb.x.logrr.cond+1.96*nb.x.logrrse.cond
    
    nb.x.trt.marg <- mean(predict(nb.x, newdata=trt.marg, type='response'))/mean(dat$t)
    nb.x.ctrl.marg <- mean(predict(nb.x, newdata=ctrl.marg, type='response'))/mean(dat$t)
    nb.x.logrr.marg <- log(nb.x.trt.marg/nb.x.ctrl.marg)
    
  } else {
    nb.x.trt.cond <- nb.x.ctrl.cond <- nb.x.trt.marg <- nb.x.ctrl.marg <- NA
    nb.x.logtrtse.cond <- nb.x.logctrlse.cond <- nb.x.logrrse.cond <- nb.x.logrrcil.cond <- nb.x.logrrciu.cond <- NA
    nb.x.logrr.cond <- nb.x.logrr.marg <- NA
  }
  
  return(cbind(n, i,
             obs.trt, obs.ctrl, obs.logrr,
             nb.trt.cond, nb.logtrtse.cond, nb.ctrl.cond, nb.logctrlse.cond,
             nb.logrr.cond, nb.logrrse.cond, nb.logrrcil.cond, nb.logrrciu.cond,
             nb.trt.marg, nb.ctrl.marg, nb.logrr.marg,
             nb.x.trt.cond, nb.x.logtrtse.cond, nb.x.ctrl.cond, nb.x.logctrlse.cond,
             nb.x.logrr.cond, nb.x.logrrse.cond, nb.x.logrrcil.cond, nb.x.logrrciu.cond,
             nb.x.trt.marg, nb.x.ctrl.marg, nb.x.logrr.marg))
}

ana_boot <- function(dat){
  

  # --------- Fit Models --------- 
  
  # NB Models 
  tryCatch({nb.x <- glm.nb(y ~ x + trt + offset(log(t)), data=dat); summ.nb.x <- summary(nb.x)},
           error=function(e) {skip.e1 <<- TRUE}, 
           warning=function(w) {skip.w1 <<- TRUE})
  
  
  # --------- Obtain RRs and rates --------- 
  
  # Prediction profiles
  trt.marg <- data.frame(x=dat$x, trt=rep(1,n), t=dat$t)
  ctrl.marg <- data.frame(x=dat$x, trt=rep(0,n), t=dat$t)
  
  
  # NB conditional and marginal rates, adjusted model - only predict if converged
  if(!skip.e1 & !skip.w1) {
    
    nb.x.trt.marg <- mean(predict(nb.x, newdata=trt.marg, type='response'))/mean(dat$t)
    nb.x.ctrl.marg <- mean(predict(nb.x, newdata=ctrl.marg, type='response'))/mean(dat$t)
    nb.x.logrr.marg <- log(nb.x.trt.marg/nb.x.ctrl.marg)
    
  } else {
    nb.x.trt.marg <- nb.x.ctrl.marg <- nb.x.logrr.marg <- NA
  }
  
  return(cbind(n, i, nb.x.trt.marg, nb.x.ctrl.marg, nb.x.logrr.marg))
}


# --------- Loop through scenarios ---------

for(k in samp_sizes){
  
  res_k <- data.frame(matrix(NA, nrow=niter, ncol=33))

  # --------- Loop through iterations ---------
  
  for(i in 1:niter){
  
    skip.e1 <- skip.w1 <- FALSE
    if(exists('nb.x')){rm(nb.x)}
    
    
    # --------- Simulate data ---------
    
    n <- k
    x <- rnorm(n, 60, 23)
    trt <- rbinom(n,1,0.5)
    t <- runif(n, 0.8, 1.2)
    m <- as.matrix(cbind(x, trt, t))
    y <- apply(m, 1, function(a) rnbinom(1, size=0.5, mu=exp(beta0 + beta1*a['x'] + beta2*a['trt'] + log(a['t'])))) 
    dat <- data.frame(id=1:n, y=y, x=x, trt=trt, t=t)


    # --------- Analysis ---------
    
    res_i <- analysis(dat)
    
    
    # --------- Bootstrap confidence intervals for marginal RRs ---------
    
    out.boot <- matrix(NA, nboot, 5)

    for(ii in 1:nboot){
      
      #ii=490
      #set.seed(ii)
      skip.e1 <- skip.w1 <- FALSE
      if(exists('nb.x')){rm(nb.x)}
      
      # Bootstrap sample
      ids <- sample(dat$id, n, replace=TRUE)
      boot.dat <- dat[ids,]
      
      # Repeat analysis in bootstrap sample
      res.boot <- data.frame(ana_boot(boot.dat))
      
      print(c(i, k, ii, skip.w1))
      
      # Extract results
      out.boot[ii,] <- c(ii, res.boot$nb.x.trt.marg, res.boot$nb.x.ctrl.marg, res.boot$nb.x.logrr.marg, skip.w1)
    }
    
    # Boostrap sample sizes that converged
    nb.x.logtrtbootn.marg <- nboot-sum(out.boot[,5])
    nb.x.logctrlbootn.marg <- nboot-sum(out.boot[,5])
    nb.x.logrrbootn.marg <- nboot-sum(out.boot[,5])
    
    # Boostrap SEs    
    nb.x.logtrtbootse.marg <- sd(log(out.boot[,2]), na.rm=TRUE)
    nb.x.logctrlbootse.marg <- sd(log(out.boot[,3]), na.rm=TRUE)
    nb.x.logrrbootse.marg <- sd(out.boot[,4], na.rm=TRUE)
    
    # Bootstrap CIs
    nb.x.logrrbootcil.marg <- quantile(out.boot[,4], 0.025, na.rm=TRUE)
    nb.x.logrrbootciu.marg <- quantile(out.boot[,4], 0.975, na.rm=TRUE)
    
    
    # --------- Save results --------- 
    
    res_k[i,] <- cbind(res_i, nb.x.logrrbootn.marg, 
                       nb.x.logtrtbootse.marg, nb.x.logctrlbootse.marg, 
                       nb.x.logrrbootse.marg, nb.x.logrrbootcil.marg, nb.x.logrrbootciu.marg)
  }
  
  # Output k results
  colnames(res_k) <- c(colnames(res_i), 
                     'nb.x.logrrbootn.marg', 
                     'nb.x.logtrtbootse.marg', 'nb.x.logctrlbootse.marg', 
                     'nb.x.logrrbootse.marg', 'nb.x.logrrbootcil.marg', 'nb.x.logrrbootciu.marg')
  
  write.csv(res_k, paste0("sim setting 1 nb 2000k 1000boot n", k, ".csv"), row.names=FALSE)
}


