## Introdution to regularization



## Setup
### Importing the dataset and library
install.packages("glmnet")

library(glmnet)

load("[Data]PNU-winter-workshop-2019.RData")

ls()

str(workshop.data)

attach(workshop.data)

### Description of sequencing data
str(x)
str(y)
str(ref)

### Elements of glmnet
ls("package:glmnet")

### Help pages for more details
help(glmnet)
vignette("glmnet_beta", package="glmnet")

## The Lasso
### Fitting the lasso model with alpha=1
fit.lasso <- glmnet(x, y, alpha=1, nlambda=10, family="gaussian")
fit.lasso

### Solution path for lasso
plot( fit.lasso, xvar="norm", label=TRUE )

### Getting coefficients for the lasso regression
dim( fit.lasso$beta )
head( fit.lasso$beta )

beta_coef <- coef( fit.lasso, s = fit.lasso$lambda[5] )
beta_pred <- predict( fit.lasso, type="coef", s = fit.lasso$lambda[5] )

### Comparison betwwen three approaches to get coefficients
tmp <- c(33, 54, 192, 350, 361, 435)

data.frame(
  fit.lasso$beta[tmp, 5],
  beta_coef[tmp+1,],
  beta_pred[tmp+1,]
)

### Coefficients in solution path for the Lasso
plot(fit.lasso, xvar="lambda", label=TRUE)
abline(v=log(116.7), cex=2, col="red")


### Selecting the tuning parameter using cross-validation
cv.lasso <- cv.glmnet(x, y, alpha=1, nlambda=100, family="gaussian", type.measure = "mse")
str(cv.lasso, max.level = 1)

### CV curve for cv.glmnet
plot(cv.lasso)

wh.1se <- min( which( cv.lasso$cvm < cv.lasso$cvup[which.min(cv.lasso$cvm)] ))
cv.lasso$lambda[ wh.1se ]
cv.lasso$lambda.1se

fit.lasso.min <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.min)
fit.lasso.1se <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.1se)

fit.lasso.min
fit.lasso.1se

lasso.nonzero.min <- which( fit.lasso.min$beta != 0 )
str(lasso.nonzero.min)

### Visualization in manhattan plot
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

#### NA is from the variables which have only one value.
na.lm.pvalue <- which(is.na(lm.pvalue))
table( x[, na.lm.pvalue[1]] )

#### Input data format for manhattan plot
df.lm.pvalue <- data.frame( ref, p=lm.pvalue)
colnames(df.lm.pvalue) <- c("SNP", "CHR", "BP", "P")

bonferroni.cutoff <- 0.05/ncol(x) # Bonferroni significant line
qqman::manhattan(df.lm.pvalue, 
                 chr="CHR", bp="BP", p="P", snp="SNP", 
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[lasso.nonzero.min,"SNP"], 
                 ylim=c(0,10))


### Prediction
#### Split training and test set
set.seed(20190130)
idx.train <- sample(nrow(x), 0.9*nrow(x))
idx.test <- (1:nrow(x))[-idx.train]

#### Training a regularization model in which a tuning paramter is selected using validation set
cv.lasso.train <- cv.glmnet(x[idx.train,], y[idx.train], family="gaussian", 
                            alpha=1, type.measure="mse")

fit.lasso.train.min <- glmnet(x[idx.train,], y[idx.train], family="gaussian", 
                              alpha=1, lambda=cv.lasso.train$lambda.min) 
fit.lasso.train.1se <- glmnet(x[idx.train,], y[idx.train], family="gaussian", 
                              alpha=1, lambda=cv.lasso.train$lambda.1se) 

fit.lasso.train.min
fit.lasso.train.1se

#### Prediction
newx <- x[idx.test, ]

newy.fit.min <- predict( fit.lasso.train.min, newx, type="response" )

head(
  data.frame(
    test.y=y[idx.test],
    fit.newy.min=as.numeric(newy.fit.min)
  )
)

##### Since regularization procedures are not prediction model,
##### we recommend "randomforest" or "boosting model" for prediction



## Elastic-net with alpha$\in(0, 1)$ 


### Generate binary response (y2)
y2 <- ifelse( y > median(y), 1, 0 )

table(y2)

### Fitting elastic-net model
fit.enet <- glmnet(x, y2, alpha=0.5, nlambda=10, family="binomial")
fit.enet


### Selecting the tuning parameters (alpha, lambda)
#### fix up a foldid}
set.seed(1234)
foldid <- sample(rep(1:10, length=length(y)))
table(foldid)

### Cross-validation with two tuning parameters
alpha.grid <- seq(0, 1, 0.1)
cv.enet <- as.list( 1:length(alpha.grid) )
for( h in 1:length(alpha.grid) ){
  cv.enet[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                            alpha=alpha.grid[h], nlambda=15, 
                            foldid=foldid)
}

cv.enet.cvm <- NULL
for( h in 1:length(cv.enet)){
  cv.enet.cvm <- cbind(cv.enet.cvm, cv.enet[[h]]$cvm)
}
dimnames(cv.enet.cvm) <- list( paste0("lambda.",1:14), paste0("alpha.",1:9) )

opt.tune <- which( cv.enet.cvm == min(cv.enet.cvm), arr.ind=TRUE )
opt.tune

opt.alpha <- alpha.grid[ opt.tune[2] ]

opt.lambda.min <- cv.enet[[ opt.tune[2] ]]$lambda.min
opt.lambda.1se <- cv.enet[[ opt.tune[2] ]]$lambda.1se

c(opt.alpha=opt.alpha, opt.lambda.min=opt.lambda.min, opt.lambda.1se=opt.lambda.1se)



### Grid search for the alpha less than 0.1
alpha.grid2 <- seq(0.01, 0.1, 0.02)
cv.enet2 <- as.list( 1:length(alpha.grid2) )
for( h in 1:length(alpha.grid2) ){
  cv.enet2[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                             alpha=alpha.grid2[h], nlambda=15, 
                             foldid=foldid)
}

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

fit.enet.opt.min <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.min)
fit.enet.opt.1se <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.1se)
fit.enet.opt.min2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.min2)
fit.enet.opt.1se2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.1se2)

fit.enet.opt.min
fit.enet.opt.min2

min(cv.enet.cvm)
min(cv.enet.cvm2)

##### We choose the values of opt.alpha and opt.lambda.min as the the optimal tuning paramters

### The significant variables of elastic-net in Mahattan plot
nonzero.enet.opt <- which( fit.enet.opt$beta != 0 )
bonferroni.cutoff <- 0.05/ncol(x)
qqman::manhattan(df.lm.pvalue, 
                 chr="CHR", bp="BP", p="P", snp="SNP", 
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[nonzero.enet.opt,"SNP"], 
                 ylim=c(0,10))


## Penalty factors

##### If we consider the variable known to be associated with response variable,
##### we involve them as covariates in a model.


### Generating arbitrary covariates
sex <- NULL
sex[y2==0] <- rbinom(sum(y2==0), 1, prob=0.8)
sex[y2==1] <- rbinom(sum(y2==1), 1, prob=0.4)
table(y2, sex)

age <- NULL
age[y2==0] <- rpois(sum(y2==0), lambda=40)
age[y2==1] <- rpois(sum(y2==1), lambda=60)
boxplot(age~y2)

x.cov <- cbind(sex, age, x)

### Penalty.factor
pf <- c(0, 0, rep(1,ncol(x)) )
table(pf)

pf.fit <- glmnet(x.cov, y2, family="binomial", penalty.factor = pf, nlambda=10 )
pf.fit

alpha.grid <- seq(0, 1, 0.1)
cv.pf <- NULL
for( h in 1:length(alpha.grid) ){
  cv.pf[[h]] <- cv.glmnet( x.cov, y2, family="binomial", type.measure="dev",
                           penalty.factor=pf, 
                           alpha=alpha.grid[h], nlambda=15 )
}

cv.pf.cvm <- NULL
for( h in 1:length(cv.pf) ){
  cv.pf.cvm <- cbind( cv.pf.cvm, cv.pf[[h]]$cvm )
}

opt.tune.pf <- which( cv.pf.cvm == min(cv.pf.cvm), arr.ind=TRUE )

opt.alpha.pf <- alpha.grid[ opt.tune.pf[2] ]
opt.lambda.min.pf <- cv.pf[[ opt.tune.pf[1] ]]$lambda.min
opt.lambda.1se.pf <- cv.pf[[ opt.tune.pf[1] ]]$lambda.1se

fit.pf.opt.min <- glmnet( x.cov, y2, family="binomial", penalty.factor=pf, 
                          alpha=opt.alpha.pf, lambda=opt.lambda.min.pf )
fit.pf.opt.1se <- glmnet( x.cov, y2, family="binomial", penalty.factor=pf, 
                          alpha=opt.alpha.pf, lambda=opt.lambda.1se.pf )

fit.pf.opt.min
fit.pf.opt.1se

##### save.image(file="part3.RData")

