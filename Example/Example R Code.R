# Example analysis
# Marginal versus conditional rate estimation for count and recurrent event data with an estimand framework
# Sarah Conner, Yijie Zhou, and Tu Xu
# November 13 2023


library(MASS)

dat <- read.csv("C://Users//conner//OneDrive - Vertex Pharmaceuticals//Documents//Exploratory//PEx rate//Manuscript//example.csv")


# Observed RR and rates
rate_trt <- sum(dat$y[dat$trt==1])/nrow(dat[dat$trt==1,])
rate_ctrl <- sum(dat$y[dat$trt==0])/nrow(dat[dat$trt==0,])
rr <- rate_trt/rate_ctrl
cbind(rate_trt, rate_ctrl, rr)



# Unadjusted NegBin Model
nb <- glm.nb(y ~ trt + offset(log(t)), data=dat)
summary(nb)

nb.trt <- predict(nb, newdata=data.frame(trt=1, t=1), type='response') 
nb.ctrl <- predict(nb, newdata=data.frame(trt=0, t=1), type='response') 
nb.rr <- exp(coef(nb)[2])
cbind(nb.trt, nb.ctrl, nb.rr)



# Adjusted NegBin Model
nb.x <- glm.nb(y ~ x + trt + offset(log(t)), data=dat)
summary(nb.x)

trt.cond <- data.frame(x=mean(dat$x), trt=1, t=1)
ctrl.cond <- data.frame(x=mean(dat$x), trt=0, t=1)

nb.x.trt.cond <- predict(nb.x, newdata=trt.cond, type='response') 
nb.x.ctrl.cond  <- predict(nb.x, newdata=ctrl.cond, type='response')
cbind(nb.x.trt.cond, nb.x.ctrl.cond, exp(-.4))



# G-computation of adjusted model
trt.marg <- data.frame(x=dat$x, trt=rep(1,n), t=dat$t)
ctrl.marg <- data.frame(x=dat$x, trt=rep(0,n), t=dat$t)

nb.x.trt.marg <- mean(predict(nb.x, newdata=trt.marg, type='response'))/mean(dat$t)
nb.x.ctrl.marg <- mean(predict(nb.x, newdata=ctrl.marg, type='response'))/mean(dat$t)
nb.x.rr.marg <- log(nb.x.trt.marg/nb.x.ctrl.marg)
cbind(nb.x.trt.marg, nb.x.ctrl.marg, nb.x.rr.marg)
