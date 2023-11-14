### Simulation Setting 1
### PEx Rate manuscript
### Sarah Conner 12/21/22

library(MASS)

set.seed(123)
niter <- 2000
nboot <- 1000

beta0 <- 0.5
beta1 <- -1.2 # binary covariate ppFEV1
beta2 <- -1.4 # log RR of treatment
beta3 <-  0.7 

#samp_sizes <- c(100,200,500,1000)
samp_sizes <- c(500)
res <- NULL


# --------- Analysis function --------- 

analysis <- function(dat){
  
  skip.e1 <- skip.w1 <- skip.e2 <- skip.w2 <- FALSE
  if(exists('nb.x')){rm(nb.x)}
  if(exists('nb.intxn')){rm(nb.intxn)}
  
  # --------- Observed rates --------- 
  
  obs.trt <- sum(dat$y[dat$trt==1])/sum(dat$t[dat$trt==1])
  obs.ctrl <- sum(dat$y[dat$trt==0])/sum(dat$t[dat$trt==0])
  obs.logrr <- log(obs.trt/obs.ctrl)
  
  
  # --------- Fit Models --------- 
  
  # NB Models 
  tryCatch({nb.intxn <- glm.nb(y ~ x*trt + offset(log(t)), data=dat); summ.nb.intxn <- summary(nb.intxn)}, 
           error=function(e) {skip.e2 <<- TRUE}, 
           warning=function(w) {skip.w2 <<- TRUE})
  tryCatch({nb.x <- glm.nb(y ~ x + trt + offset(log(t)), data=dat); summ.nb.x <- summary(nb.x)},
           error=function(e) {skip.e1 <<- TRUE}, 
           warning=function(w) {skip.w1 <<- TRUE})
  
  skip.e2 <- ifelse(abs(nb.intxn$coefficients[4])>5, TRUE, skip.e2)
  
  nb <- glm.nb(y ~ trt + offset(log(t)), data=dat)
  summ.nb <- summary(nb)
  
  
  
  # --------- Obtain RRs and rates --------- 
  
  # Prediction profiles
  trt.cond <- data.frame(x=mean(dat$x), trt=1, t=1)
  ctrl.cond <- data.frame(x=mean(dat$x), trt=0, t=1)
  trt.marg <- data.frame(x=dat$x, trt=rep(1,n), t=dat$t)
  ctrl.marg <- data.frame(x=dat$x, trt=rep(0,n), t=dat$t)
  
  trt.x1 <- data.frame(x=1, trt=1, t=1)
  ctrl.x1 <- data.frame(x=1, trt=0, t=1)  
  trt.x0 <- data.frame(x=0, trt=1, t=1)
  ctrl.x0 <- data.frame(x=0, trt=0, t=1)
  
  
  # NB rates, unadjusted model
  nb.trt.cond <- predict(nb, newdata=trt.cond, type='response')
  nb.logtrtse.cond <- predict(nb, newdata=trt.cond, se=TRUE)$se.fit
  nb.ctrl.cond <- predict(nb, newdata=ctrl.cond, type='response')
  nb.logctrlse.cond <- predict(nb, newdata=ctrl.cond, se=TRUE)$se.fit
  
  nb.logrr.cond <- coef(nb)[2]
  nb.logrrse.cond <- summ.nb$coef[2,2]
  nb.logrrcil.cond <- nb.logrr.cond-1.96*nb.logrrse.cond
  nb.logrrciu.cond <- nb.logrr.cond+1.96*nb.logrrse.cond
  
  
  # NB conditional and marginal rates, adjusted model - only predict if converged
  if(!skip.e1 & !skip.w1) {
    
    nb.x.trt.cond <- predict(nb.x, newdata=trt.cond, type='response') 
    nb.x.ctrl.cond  <- predict(nb.x, newdata=ctrl.cond, type='response') 
    nb.x.logtrtse.cond <- predict(nb.x, newdata=trt.cond, se=TRUE)$se.fit
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
  
  # NB conditional and marginal rates, interaction model - only predict if converged
  if(!skip.e2 & !skip.w2) {
    
    # Overall
    nb.intxn.trt.cond <- predict(nb.intxn, newdata=trt.cond, type='response') 
    nb.intxn.ctrl.cond  <- predict(nb.intxn, newdata=ctrl.cond, type='response') 
    nb.intxn.logtrtse.cond <- predict(nb.intxn, newdata=trt.cond, se=TRUE)$se.fit
    nb.intxn.logctrlse.cond <- predict(nb.intxn, newdata=ctrl.cond, se=TRUE)$se.fit
    
    nb.intxn.logrr.cond <- coef(nb.intxn)[3] + coef(nb.intxn)[4]*mean(dat$x)
    nb.intxn.logrrse.cond <- sqrt(summ.nb.intxn$coef[3,2]^2 + (mean(dat$x)*summ.nb.intxn$coef[4,2])^2 + 2*mean(dat$x)*summ.nb.intxn$cov.scaled[3,4])
    nb.intxn.logrrcil.cond <- nb.intxn.logrr.cond-1.96*nb.intxn.logrrse.cond
    nb.intxn.logrrciu.cond <- nb.intxn.logrr.cond+1.96*nb.intxn.logrrse.cond
    
    nb.intxn.trt.marg <- mean(predict(nb.intxn, newdata=trt.marg, type='response'))/mean(dat$t)
    nb.intxn.ctrl.marg <- mean(predict(nb.intxn, newdata=ctrl.marg, type='response'))/mean(dat$t)
    nb.intxn.logrr.marg <- log(nb.intxn.trt.marg/nb.intxn.ctrl.marg)
    
    
    # Subgroup x=1
    nb.intxn.trt.x1 <- predict(nb.intxn, newdata=trt.x1, type='response') 
    nb.intxn.ctrl.x1  <- predict(nb.intxn, newdata=ctrl.x1, type='response') 
    nb.intxn.logtrtse.x1 <- predict(nb.intxn, newdata=trt.x1, se=TRUE)$se.fit
    nb.intxn.logctrlse.x1 <- predict(nb.intxn, newdata=ctrl.x1, se=TRUE)$se.fit
    
    nb.intxn.logrr.x1 <- coef(nb.intxn)[3] + coef(nb.intxn)[4]
    nb.intxn.logrrse.x1 <- sqrt(summ.nb.intxn$coef[3,2]^2 + summ.nb.intxn$coef[4,2]^2 + 2*summ.nb.intxn$cov.scaled[3,4])
    nb.intxn.logrrcil.x1 <- nb.intxn.logrr.x1-1.96*nb.intxn.logrrse.x1
    nb.intxn.logrrciu.x1 <- nb.intxn.logrr.x1+1.96*nb.intxn.logrrse.x1
    
    
    # Subgroup x=0
    nb.intxn.trt.x0 <- predict(nb.intxn, newdata=trt.x0, type='response') 
    nb.intxn.ctrl.x0  <- predict(nb.intxn, newdata=ctrl.x0, type='response') 
    nb.intxn.logtrtse.x0 <- predict(nb.intxn, newdata=trt.x0, se=TRUE)$se.fit
    nb.intxn.logctrlse.x0 <- predict(nb.intxn, newdata=ctrl.x0, se=TRUE)$se.fit
    
    nb.intxn.logrr.x0 <- coef(nb.intxn)[3] 
    nb.intxn.logrrse.x0 <- summ.nb.intxn$coef[3,2]
    nb.intxn.logrrcil.x0 <- nb.intxn.logrr.x0-1.96*nb.intxn.logrrse.x0
    nb.intxn.logrrciu.x0 <- nb.intxn.logrr.x0+1.96*nb.intxn.logrrse.x0
    

  } else {
    nb.intxn.trt.cond <- nb.intxn.ctrl.cond <- nb.intxn.trt.marg <- nb.intxn.ctrl.marg <- NA
    nb.intxn.logtrtse.cond <- nb.intxn.logctrlse.cond <- nb.intxn.logrrse.cond <- nb.intxn.logrrcil.cond <- nb.intxn.logrrciu.cond <- NA
    nb.intxn.logrr.cond <- nb.intxn.logrr.marg <- NA
    nb.intxn.trt.x1 <- nb.intxn.ctrl.x1 <- nb.intxn.logtrtse.x1 <- nb.intxn.logctrlse.x1 <- NA
    nb.intxn.logrr.x1 <- nb.intxn.logrrse.x1 <- nb.intxn.logrrcil.x1 <- nb.intxn.logrrciu.x1 <- NA
    nb.intxn.trt.x0 <- nb.intxn.ctrl.x0 <- nb.intxn.logtrtse.x0 <- nb.intxn.logctrlse.x0 <- NA
    nb.intxn.logrr.x0 <- nb.intxn.logrrse.x0 <- nb.intxn.logrrcil.x0 <- nb.intxn.logrrciu.x0 <- NA
  }
  
  
  return(cbind(n, i,
               obs.trt, obs.ctrl, obs.logrr,
               nb.trt.cond, nb.logtrtse.cond, nb.ctrl.cond, nb.logctrlse.cond, 
               nb.logrr.cond, nb.logrrse.cond, nb.logrrcil.cond, nb.logrrciu.cond,
               nb.x.trt.cond, nb.x.logtrtse.cond, nb.x.ctrl.cond, nb.x.logctrlse.cond,
               nb.x.logrr.cond, nb.x.logrrse.cond, nb.x.logrrcil.cond, nb.x.logrrciu.cond,
               nb.x.trt.marg, nb.x.ctrl.marg, nb.x.logrr.marg,
               nb.intxn.trt.cond, nb.intxn.logrrse.cond, nb.intxn.ctrl.cond, nb.intxn.logctrlse.cond,
               nb.intxn.logrr.cond, nb.intxn.logrrse.cond, nb.intxn.logrrcil.cond, nb.intxn.logrrciu.cond,
               nb.intxn.trt.marg, nb.intxn.ctrl.marg, nb.intxn.logrr.marg,
               nb.intxn.trt.x1, nb.intxn.ctrl.x1, nb.intxn.logtrtse.x1, nb.intxn.logctrlse.x1, 
               nb.intxn.logrr.x1, nb.intxn.logrrse.x1, nb.intxn.logrrcil.x1, nb.intxn.logrrciu.x1, 
               nb.intxn.trt.x0, nb.intxn.ctrl.x0, nb.intxn.logtrtse.x0, nb.intxn.logctrlse.x0, 
               nb.intxn.logrr.x0, nb.intxn.logrrse.x0, nb.intxn.logrrcil.x0, nb.intxn.logrrciu.x0
               ))
}


ana_boot <- function(dat){
  
  skip.e1 <- skip.w1 <- skip.e2 <- skip.w2 <- FALSE
  if(exists('nb.x')){rm(nb.x)}
  if(exists('nb.intxn')){rm(nb.intxn)}
  
  
  # --------- Fit Models --------- 
  
  # NB Models 
  tryCatch({nb.intxn <- glm.nb(y ~ x*trt + offset(log(t)), data=dat); summ.nb.intxn <- summary(nb.intxn)}, 
           error=function(e) {skip.e2 <<- TRUE}, 
           warning=function(w) {skip.w2 <<- TRUE})
  tryCatch({nb.x <- glm.nb(y ~ x + trt + offset(log(t)), data=dat); summ.nb.x <- summary(nb.x)},
           error=function(e) {skip.e1 <<- TRUE}, 
           warning=function(w) {skip.w1 <<- TRUE})
  
  
  # --------- Obtain RRs and rates --------- 
  
  # Prediction profiles
  trt.marg <- data.frame(x=dat$x, trt=rep(1,n), t=dat$t)
  ctrl.marg <- data.frame(x=dat$x, trt=rep(0,n), t=dat$t)
  
  
  # NB marginal rates, adjusted model - only predict if converged
  if(!skip.e1 & !skip.w1) {
    
    nb.x.trt.marg <- mean(predict(nb.x, newdata=trt.marg, type='response'))/mean(dat$t)
    nb.x.ctrl.marg <- mean(predict(nb.x, newdata=ctrl.marg, type='response'))/mean(dat$t)
    nb.x.logrr.marg <- log(nb.x.trt.marg/nb.x.ctrl.marg)
    
  } else {
    nb.x.trt.marg <- nb.x.ctrl.marg <- nb.x.logrr.marg <- NA
  }
  
  
  # NB conditional and marginal rates, interaction model - only predict if converged
  if(!skip.e2 & !skip.w2) {
    
    nb.intxn.trt.marg <- mean(predict(nb.intxn, newdata=trt.marg, type='response'))/mean(dat$t)
    nb.intxn.ctrl.marg <- mean(predict(nb.intxn, newdata=ctrl.marg, type='response'))/mean(dat$t)
    nb.intxn.logrr.marg <- log(nb.intxn.trt.marg/nb.intxn.ctrl.marg)
    
    
  } else {
    nb.intxn.trt.marg <- nb.intxn.ctrl.marg <- nb.intxn.logrr.marg <- NA
  }
  
  return(cbind(n, i,
               nb.x.trt.marg, nb.x.ctrl.marg, nb.x.logrr.marg,
               nb.intxn.trt.marg, nb.intxn.ctrl.marg, nb.intxn.logrr.marg,
               skip.w1, skip.w2))
}


# --------- Loop through scenarios ---------


for(k in samp_sizes){
  
  res_k <- data.frame(matrix(NA, nrow=niter, ncol=63))
  
  
  # --------- Loop through iterations ---------
  
  for(i in 1:niter){
    
    
    # --------- Simulate data ---------
    
    n <- k
    x <- rbinom(n, 1, 0.5)
    trt <- rbinom(n, 1, 0.5)
    t <- runif(n, 0.8, 1.2)
    m <- as.matrix(cbind(x, trt, t))
    y <- apply(m, 1, function(a) rnbinom(1, size=0.5, mu=exp(beta0 + beta1*a['x'] + beta2*a['trt'] + beta3*a['x']*a['trt'] + log(a['t'])))) 
    dat <- data.frame(id=1:n, y=y, x=x, trt=trt, t=t)
    
    
    # --------- Analysis ---------
    
    res_i <- analysis(dat)
    nb.intxn <- glm.nb(y ~ x*trt + offset(log(t)), data=dat)
    (summ.nb.intxn <- summary(nb.intxn))

    
    # --------- Bootstrap confidence intervals for marginal RRs ---------
    
    out.boot <- matrix(NA, nboot, 9)
    
    for(ii in 1:nboot){
      
      #ii=490
      #set.seed(ii)
      
      # Bootstrap sample
      ids <- sample(dat$id, n, replace=TRUE)
      boot.dat <- dat[ids,]
      
      # Repeat analysis in bootstrap sample
      res.boot <- data.frame(ana_boot(boot.dat))
      
      print(c(i, k, ii, res.boot$skip.w1, res.boot$skip.w2))
      
      # Extract results
      out.boot[ii,] <- c(ii, 
                         res.boot$nb.x.trt.marg, res.boot$nb.intxn.trt.marg, 
                         res.boot$nb.x.ctrl.marg, res.boot$nb.intxn.ctrl.marg, 
                         res.boot$nb.x.logrr.marg, res.boot$nb.intxn.logrr.marg, 
                         res.boot$skip.w1, res.boot$skip.w2)
    }
    
    # Bootstrap sample sizes that converged
    nb.x.bootnconverge <- nboot-sum(out.boot[,8])
    nb.intxn.bootnconverge <- nboot-sum(out.boot[,9])
    
    
    # Bootstrap SEs
    nb.x.logtrtbootse.marg <- sd(log(out.boot[,2]), na.rm=TRUE)
    nb.intxn.logtrtbootse.marg <- sd(log(out.boot[,3]), na.rm=TRUE)
    
    nb.x.logctrlbootse.marg <- sd(log(out.boot[,4]), na.rm=TRUE)
    nb.intxn.logctrlbootse.marg <- sd(log(out.boot[,5]), na.rm=TRUE)
    
    nb.x.logrrbootse.marg <- sd(out.boot[,6], na.rm=TRUE)
    nb.intxn.logrrbootse.marg <- sd(out.boot[,7], na.rm=TRUE)
    
    
    # Bootstrap CIs
    nb.x.logrrbootcil.marg <- quantile(out.boot[,6], 0.025, na.rm=TRUE)
    nb.x.logrrbootciu.marg <- quantile(out.boot[,6], 0.975, na.rm=TRUE)
    
    nb.intxn.logrrbootcil.marg <- quantile(out.boot[,7], 0.025, na.rm=TRUE)
    nb.intxn.logrrbootciu.marg <- quantile(out.boot[,7], 0.975, na.rm=TRUE)
    
    
    # --------- Save results --------- 
    
    res_k[i,] <- cbind(res_i, 
                       nb.x.bootnconverge, nb.intxn.bootnconverge,
                       nb.x.logtrtbootse.marg, nb.intxn.logtrtbootse.marg,
                       nb.x.logctrlbootse.marg, nb.intxn.logctrlbootse.marg,
                       nb.x.logrrbootse.marg, nb.intxn.logrrbootse.marg,
                       nb.x.logrrbootcil.marg, nb.intxn.logrrbootcil.marg,
                       nb.x.logrrbootciu.marg, nb.intxn.logrrbootciu.marg)
  }
  
  # Output k results
  colnames(res_k) <- c(colnames(res_i), 
                       'nb.x.bootnconverge', 'nb.intxn.bootnconverge',
                       'nb.x.logtrtbootse.marg', 'nb.intxn.logtrtbootse.marg', 
                       'nb.x.logctrlbootse.marg', 'nb.intxn.logctrlbootse.marg', 
                       'nb.x.logrrbootse.marg', 'nb.intxn.logrrbootse.marg',
                       'nb.x.logrrbootcil.marg', 'nb.intxn.logrrbootcil.marg',
                       'nb.x.logrrbootciu.marg', 'nb.intxn.logrrbootciu.marg')
  
  res_k
  write.csv(res_k, paste0("sim setting 3 nb 2000k 1000boot n", k, ".csv"), row.names=FALSE)
}




