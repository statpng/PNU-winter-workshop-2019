  ## Introdution to regularization
## Setup
### Importing the dataset and library

library(qqman)
library(glmnet)

load("[Data]PNU-winter-workshop-2019.RData")

ls()

str(workshop.data)

attach(workshop.data)


### Description of sequencing data

str(x)
str(y)
str(ref)


### Composition of glmnet

ls("package:glmnet")


### Main functions
# See below for details,
help(glmnet)
vignette("glmnet_beta", package="glmnet")

## The Lasso
### Fitting the lasso model with alpha=1

fit.lasso <- glmnet(x, y, alpha=1, nlambda=10, family="gaussian")
fit.lasso

# - nlambda: recommended to be 100 or more.
# - Df: the number of nonzero coefficients.
# - %dev: the percent deviance explained (for gaussian, this is the R-square).
# - Deviance: the difference of 2*log-likelihood between saturated model and current model.
# - Lambda: the corresponding value of lambda.

### Solution path for lasso

plot( fit.lasso, xvar="norm", label=TRUE )


### Getting coefficients for the lasso regression
#### Three functions which retrieve coefficients

dim( fit.lasso$beta )
head( fit.lasso$beta )


### Getting coefficients for the lasso regression

beta_coef <- coef( fit.lasso, s = fit.lasso$lambda[5] )
beta_pred <- predict( fit.lasso, type="coef", s = fit.lasso$lambda[5] )


### Comparison betwwen three approaches to get coefficients


tmp <- c(33, 54, 192, 350, 361, 435)

data.frame(
  beta=fit.lasso$beta[tmp, 5],
  coef=beta_coef[tmp+1,],
  pred=beta_pred[tmp+1,]
)


plot(fit.lasso, xvar="lambda", label=TRUE)
abline(v=log(116.7), cex=2, col="red")



### Selecting the tuning parameter using cross-validation

cv.lasso <- cv.glmnet(x, y, alpha=1, nlambda=100, family="gaussian", type.measure = "mse")
str(cv.lasso, max.level = 1)



### Plot for cv.glmnet

plot(cv.lasso)


# - lambda.min: the value which gives minumum mean cross-validated error (like "mse").
# - lambda.1se: the value which gives the most regularized model such that error is within one standard error of the minimum(lambda.1se).


wh.1se <- min( which( cv.lasso$cvm < cv.lasso$cvup[which.min(cv.lasso$cvm)] ))
cv.lasso$lambda[ wh.1se ]
cv.lasso$lambda.1se

fit.lasso.min <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.min)
fit.lasso.1se <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.1se)

fit.lasso.min
fit.lasso.1se



### Visualization in manhattan plot

lasso.nonzero.min <- which( fit.lasso.min$beta != 0 )
str(lasso.nonzero.min)

summary( lm(y~x[, 1]) )

lm.pvalue <- NULL
for( j in 1:ncol(x) ){
  lm.coef <- summary( lm(y~x[, j]) )$coef
  if( nrow(lm.coef)==2 ){
    lm.pvalue[j] <- lm.coef[2, 4]
  } else if ( nrow(lm.coef)==1 ){
    lm.pvalue[j] <- NA
  }
}


str( lm.pvalue )
summary( lm.pvalue )


na.lm.pvalue <- which(is.na(lm.pvalue))
table( x[, na.lm.pvalue[1]] )


df.lm.pvalue <- data.frame( ref, p=lm.pvalue)
colnames(df.lm.pvalue) <- c("SNP", "CHR", "BP", "P")
list(ref=head(ref), df.lm.pvalue=head(df.lm.pvalue))


bonferroni.cutoff <- 0.05/ncol(x)
qqman::manhattan(df.lm.pvalue, chr="CHR", bp="BP", p="P", snp="SNP", suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[lasso.nonzero.min,"SNP"], ylim=c(0,10))



### Prediction

set.seed(20190130)
idx.train <- sample(nrow(x), 0.9*nrow(x))
idx.test <- (1:nrow(x))[-idx.train]


cv.lasso.train <- cv.glmnet(x[idx.train,], y[idx.train], alpha=1, family="gaussian")

fit.lasso.train.min <- glmnet(x[idx.train,], y[idx.train], alpha=1, family="gaussian", lambda=cv.lasso.train$lambda.min) 
fit.lasso.train.1se <- glmnet(x[idx.train,], y[idx.train], alpha=1, family="gaussian", lambda=cv.lasso.train$lambda.1se) 

fit.lasso.train.min
fit.lasso.train.1se

newx <- x[idx.test, ]

newy.fit.min <- predict( fit.lasso.train.min, newx, type="response" )

head(
  data.frame(
    test.y=y[idx.test],
    fit.newy.min=as.numeric(newy.fit.min)
  )
)



## Elastic-net with alpha$\in(0, 1)$ 

# - Distinct advantages of elastic-net
#   * Grouping selection
#   * When n<<p, we can select variables up to **p**
  
### Generate binary response (y2)
y2 <- ifelse( y > median(y), 1, 0 )

table(y2)

### Fitting elastic-net model
fit.enet <- glmnet(x, y2, alpha=0.5, nlambda=10, family="binomial")
fit.enet


### Selecting the tuning parameters (alpha, lambda)

set.seed(1234)
foldid <- sample(rep(1:10, length=length(y)))



### Cross-validation with two tuning parameters

alpha.grid <- seq(0, 1, 0.1)
cv.enet <- as.list( 1:length(alpha.grid) )
for( h in 1:length(alpha.grid) ){
  cv.enet[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                            alpha=alpha.grid[h], nlambda=15, 
                            foldid=foldid)
}
# save(cv.enet, alpha.grid, file="[object]cv.enet.RData")


load("[object]cv.enet.RData")
cv.enet.cvm <- NULL
for( h in 1:length(cv.enet)){
  cv.enet.cvm <- cbind(cv.enet.cvm, cv.enet[[h]]$cvm)
}
dimnames(cv.enet.cvm) <- list( paste0("lambda.",1:14), paste0("alpha.",1:11) )

opt.tune <- which( cv.enet.cvm == min(cv.enet.cvm), arr.ind=TRUE )
opt.tune

opt.alpha <- alpha.grid[ opt.tune[2] ]
opt.lambda.min <- cv.enet[[ opt.tune[2] ]]$lambda.min
opt.lambda.1se <- cv.enet[[ opt.tune[2] ]]$lambda.1se

c(opt.alpha=opt.alpha, opt.lambda.min=opt.lambda.min, opt.lambda.1se=opt.lambda.1se)



# - we must use the same foldid so that we compare the results for each alpha.
# - You should investigate the alpha less than 0.1.


### Additional grid search for the alpha less than 0.1

alpha.grid2 <- seq(0.01, 0.1, 0.02)
cv.enet2 <- as.list( 1:length(alpha.grid2) )
for( h in 1:length(alpha.grid2) ){
  cv.enet2[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                             alpha=alpha.grid2[h], nlambda=15, 
                             foldid=foldid)
}
# save(cv.enet2, alpha.grid2, file="[object]cv.enet2.RData")



load("[object]cv.enet2.RData")
cv.enet.cvm2 <- NULL
for( h in 1:length(cv.enet2)){
  cv.enet.cvm2 <- cbind(cv.enet.cvm2, cv.enet2[[h]]$cvm)
}
dimnames(cv.enet.cvm2) <- list( paste0("lambda.",1:14), paste0("alpha.",1:5) )

opt.tune2 <- which( cv.enet.cvm2 == min(cv.enet.cvm2), arr.ind=TRUE )
opt.tune2

opt.alpha2 <- alpha.grid2[ opt.tune2[2] ]
opt.lambda.min2 <- cv.enet2[[ opt.tune2[2] ]]$lambda.min
opt.lambda.1se2 <- cv.enet2[[ opt.tune2[2] ]]$lambda.1se

c(opt.alpha2=opt.alpha2, opt.lambda.min2=opt.lambda.min2, opt.lambda.1se2=opt.lambda.1se2)

# - With repeating this procedures, the optimal alpha value is too small that the number of selected variables can be extremely large.


fit.enet.opt.min <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.min)
fit.enet.opt.1se <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.1se)
fit.enet.opt.min2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.min2)
fit.enet.opt.1se2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.1se2)

fit.enet.opt.min
fit.enet.opt.1se
fit.enet.opt.min2
fit.enet.opt.1se2

min(cv.enet.cvm)
min(cv.enet.cvm2)

# - Thus, we choose the values of opt.alpha and opt.lambda as the the optimal tuning paramters


### The significant variables of elastic-net in Mahattan plot

nonzero.enet.opt <- which( fit.enet.opt.min$beta != 0 )
bonferroni.cutoff <- 0.05/ncol(x)
qqman::manhattan(df.lm.pvalue, 
                 chr="CHR", bp="BP", p="P", snp="SNP", 
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[nonzero.enet.opt,"SNP"], 
                 ylim=c(0,10))


## Penalty factors
  
# - If we consider the variable known to be associated with response variable, 
# we involve them as covariates in the model.


### Generating arbitrary covariates

cov1 <- NULL
cov1[y2==0] <- rbinom(sum(y2==0), 1, prob=0.8)
cov1[y2==1] <- rbinom(sum(y2==1), 1, prob=0.4)
table(y2, cov1)



cov2 <- NULL
cov2[y2==0] <- rpois(sum(y2==0), lambda=40)
cov2[y2==1] <- rpois(sum(y2==1), lambda=60)
boxplot(cov2~y2)

x.cov <- cbind(cov1, cov2, x)

pf <- c(0, 0, rep(1,ncol(x)) )
table(pf)

pf.fit <- glmnet(x.cov, y2, family="binomial", penalty.factor = pf, nlambda=10 )
pf.fit



# save.image(file="part3.RData")
