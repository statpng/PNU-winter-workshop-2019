---
layout: page
title: Part 2
permalink: /part2/
---


# Part 2. Regularization for analysis of high-dimensional genomic data

<hr>
<br>

Today, we will start from taking a look at regularzation methods filtering variant calls to analyzing single
cell RNA-seq data across various contexts of high-throughput genomic
analysis.

We will assume that you have basic UNIX skills, which can be obtained
by taking a free online course
[HERE](https://www.codecademy.com/learn/learn-the-command-line), but
we we will guide from the beginning, just in case you may get lost
from the beginning.

---
<br>

### Schedule:

| Part    | Time                  | Topics                                                        |
| :-----: |   :--------------:    | :-----------------------                                      |
| II      | Regularization        | **Variable selection in high-dimensional data**               |
|         |                       | **Regularization procedures**                                 |
|         |                       | **"glmnet"**                                                  |
|         |                       | **"The Lasso"**                                               |
|         |                       | &nbsp; &nbsp; - Solution path                                 |
|         |                       | &nbsp; &nbsp; - Cross-validation                              |
|         |                       | &nbsp; &nbsp; - Variable selection                            |
|         |                       | &nbsp; &nbsp; - Comparison with univariate analysis           |
|         |                       | &nbsp; &nbsp; - Prediction                                    |
|         |                       | **Elastic-net**                                               |
|         |                       | &nbsp; &nbsp; - Cross-validation for two tuning parameters    |
|         |                       | **Covariate-adjusted model**                                  |

<br>
### Instructors:

Kipoong Kim

---



<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/highlight.js/9.5.0/styles/default.min.css">
<script src="//cdnjs.cloudflare.com/ajax/libs/highlight.js/9.5.0/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

<link rel="stylesheet" href="//cdn.jsdelivr.net/highlight.js/9.5.0/styles/default.min.css">
<script src="//cdn.jsdelivr.net/highlight.js/9.5.0/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>


<!-- dynamically load mathjax for compatibility with self-contained -->
<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

<br>

### 0. Setup
#### 0-1. Installing and calling R packages
```
install.packages("glmnet")
install.packages("qqman")
library(glmnet)
library(qqman)
```
<br>
#### 0-2. Importing the dataset
```
load("[Data]PNU-winter-workshop-2019.RData")
ls()
str(workshop.data)
attach(workshop.data)
```
<br>
#### 0-3. Descriptions of data
```
str(x)
str(y)
str(ref)
```
<br>
#### 0-4. Elements of glmnet
```
ls("package:glmnet")
```
<br>
#### 0-5. Help pages for more details
```
help(glmnet)
vignette("glmnet_beta", package="glmnet")
```
<br>
### 1. The Lasso
#### 1-1. Fitting the lasso model with alpha=1
```
fit.lasso <- glmnet(x, y, alpha=1, nlambda=10, family="gaussian")
fit.lasso
```
<br>
#### 1-2. Solution path for lasso
```
plot( fit.lasso, xvar="norm", label=TRUE )
```
<br>
#### 1-3. Getting coefficients for the lasso regression
```
dim( fit.lasso$beta )
head( fit.lasso$beta )
```
<br>
#### 1-4. Other functions to get coefficients
```
beta_coef <- coef( fit.lasso, s = fit.lasso$lambda[5] )
beta_pred <- predict( fit.lasso, type="coef", s = fit.lasso$lambda[5] )
```
<br>
#### 1-5. Comparison betwwen three approaches to retrieve coefficients
```
tmp <- c(33, 54, 192, 350, 361, 435)
data.frame(
  fit.lasso$beta[tmp, 5],
  beta_coef[tmp+1,],
  beta_pred[tmp+1,]
)
```
<br>
#### 1-6. Coefficients in solution path for the Lasso
```
plot(fit.lasso, xvar="lambda", label=TRUE)
abline(v=log(116.7), cex=2, col="red")
```
<hr>
<br>
### 2. Cross-validation
#### 2-1. Selecting the tuning parameter using cross-validation
```
cv.lasso <- cv.glmnet(x, y, alpha=1, nlambda=100, family="gaussian", type.measure = "mse")
str(cv.lasso, max.level = 1)
```
<br>
#### 2-2. CV curve for cv.glmnet
```
plot(cv.lasso)
```
<br>
#### 2-3. The optimal lambda value
```
wh.1se <- min( which( cv.lasso$cvm < cv.lasso$cvup[which.min(cv.lasso$cvm)] ))
cv.lasso$lambda[ wh.1se ]
cv.lasso$lambda.1se
```
<br>
#### 2-4. The Lasso models with the optimal lambda values ($$\lambda_{min}, \lambda_{1se}$$)
```
fit.lasso.min <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.min)
fit.lasso.1se <- glmnet(x, y, alpha=1, family="gaussian", lambda=cv.lasso$lambda.1se)
fit.lasso.min
fit.lasso.1se
```
<br>
#### 2-5. Our final Lasso model with $$\lambda_{min}$$
```
lasso.nonzero.min <- which( fit.lasso.min$beta != 0 )
str(lasso.nonzero.min)
```
<hr><br>
### 3. Visualization in Manhattan plot
#### 3-1. P-values for univariate analysis
```
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
```
<br>
#### 3-2. Treating the missing values
```
na.lm.pvalue <- which(is.na(lm.pvalue))
table( x[, na.lm.pvalue[1]] )
```
- NA is from the variables which have only one value.

<br>
#### 3-3. Formatting input data for manhattan plot
```
df.lm.pvalue <- data.frame( ref, p=lm.pvalue)
colnames(df.lm.pvalue) <- c("SNP", "CHR", "BP", "P")
```
<br>
#### 3-4. Drawing Manhattan plot with green-colored dots (variables selected by the **Lasso**)
```
bonferroni.cutoff <- 0.05/ncol(x) # Bonferroni significant line
qqman::manhattan(df.lm.pvalue,
                 chr="CHR", bp="BP", p="P", snp="SNP",
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[lasso.nonzero.min,"SNP"],
                 ylim=c(0,10))
```
<hr>
<br>

### 4. Prediction
#### 4-1. Split training and test set
```
set.seed(20190130)
idx.train <- sample(nrow(x), 0.9*nrow(x))
idx.test <- (1:nrow(x))[-idx.train]
```
<br>
#### 4-2. Training a regularization model in which a tuning paramter is selected using validation set
```
cv.lasso.train <- cv.glmnet(x[idx.train,], y[idx.train], family="gaussian",
                            alpha=1, type.measure="mse")

fit.lasso.train.min <- glmnet(x[idx.train,], y[idx.train], family="gaussian",
                              alpha=1, lambda=cv.lasso.train$lambda.min)
fit.lasso.train.1se <- glmnet(x[idx.train,], y[idx.train], family="gaussian",
                              alpha=1, lambda=cv.lasso.train$lambda.1se)
fit.lasso.train.min
fit.lasso.train.1se
```
<br>
#### 4-3. Predicting the new data
```
newx <- x[idx.test, ]
newy.fit.min <- predict( fit.lasso.train.min, newx, type="response" )

head(
  data.frame(
    test.y=y[idx.test],
    fit.newy.min=as.numeric(newy.fit.min)
  )
)
```
  - **Notification** <br>Since regularization procedures are not prediction model, we recommend "randomforest" or "boosting model" for prediction

<br>
 <hr> <hr> <br><br>

### 5. Elastic-net with alpha in [0, 1]
#### 5-1. Generating binary response "y2"
```
y2 <- ifelse( y > median(y), 1, 0 )
table(y2)
```
<br>
#### 5-2. Fitting elastic-net model
```
fit.enet <- glmnet(x, y2, alpha=0.5, nlambda=10, family="binomial")
fit.enet
```
<hr><br>
### 6. Selecting the tuning parameters ($$\alpha, \lambda$$)
#### 6-1. Fixing subsamples with "foldid"
```
set.seed(1234)
foldid <- sample(rep(1:10, length=length(y)))
table(foldid)
```
<br>
#### 6-2. Cross-validation for two tuning parameters with "dev" measure
```
alpha.grid <- seq(0, 1, 0.1)
cv.enet <- as.list( 1:length(alpha.grid) )
for( h in 1:length(alpha.grid) ){
  cv.enet[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                            alpha=alpha.grid[h], nlambda=15,
                            foldid=foldid)
}
```
<br>
#### 6-3. Combining the results of cross-validation according to $$\alpha$$.
```
cv.enet.cvm <- NULL
for( h in 1:length(cv.enet)){
  cv.enet.cvm <- cbind(cv.enet.cvm, cv.enet[[h]]$cvm)
}
dimnames(cv.enet.cvm) <- list( paste0("lambda.",1:14), paste0("alpha.",1:9) )
```
<br>
#### 6-4. Finding the optimal tuning parameters with two approaches ($$\lambda_{min}, \lambda_{1se}$$)
```
opt.tune <- which( cv.enet.cvm == min(cv.enet.cvm), arr.ind=TRUE )
opt.tune

opt.alpha <- alpha.grid[ opt.tune[2] ]
opt.lambda.min <- cv.enet[[ opt.tune[2] ]]$lambda.min
opt.lambda.1se <- cv.enet[[ opt.tune[2] ]]$lambda.1se
c(opt.alpha=opt.alpha, opt.lambda.min=opt.lambda.min, opt.lambda.1se=opt.lambda.1se)
```

<br>
#### 6-5. Further grid search for the alpha less than 0.1
```
alpha.grid2 <- seq(0.01, 0.1, 0.02)
cv.enet2 <- as.list( 1:length(alpha.grid2) )
for( h in 1:length(alpha.grid2) ){
  cv.enet2[[h]] <- cv.glmnet(x, y2, family="binomial", type.measure="dev",
                             alpha=alpha.grid2[h], nlambda=15,
                             foldid=foldid)
}
```
<br>
#### 6-6. Determining whether to view the detail ranges.
```
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
```
- Since $$\alpha=0.1$$ and $$\alpha=0.01$$ are also small value enough to have many variables, we stop here.

<br>
#### 6-7. The optimal elastic-net models with $$\alpha=0.1$$.
```
fit.enet.opt.min <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.min)
fit.enet.opt.1se <- glmnet(x, y2, family="binomial", alpha=opt.alpha, lambda=opt.lambda.1se)
fit.enet.opt.min2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.min2)
fit.enet.opt.1se2 <- glmnet(x, y2, family="binomial", alpha=opt.alpha2, lambda=opt.lambda.1se2)

fit.enet.opt.min
fit.enet.opt.min2
min(cv.enet.cvm)
min(cv.enet.cvm2)
```
  - We choose the values of opt.alpha ($$\alpha=0.1$$) and opt.lambda.min ($$\lambda_{min}$$ corresponding to $$\alpha=0.1$$) as the the optimal tuning paramters

<hr><br>

### 7. Manhattan plot with the significant variables of the optimal elastic-net model ("fit.enet.opt.min")
```
nonzero.enet.opt <- which( fit.enet.opt$beta != 0 )
bonferroni.cutoff <- 0.05/ncol(x)
qqman::manhattan(df.lm.pvalue,
                 chr="CHR", bp="BP", p="P", snp="SNP",
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[nonzero.enet.opt,"SNP"],
                 ylim=c(0,10))
```

<hr><br>

### 8. Covariate-adjusted model
  - If we consider the variable known to be associated with response variable, we involve them as covariates in a model.

<br>
#### 8-1. Generating arbitrary covariates
```
cov1 <- NULL
cov1[y2==0] <- rbinom(sum(y2==0), 1, prob=0.8)
cov1[y2==1] <- rbinom(sum(y2==1), 1, prob=0.4)
table(y2, cov1)

cov2 <- NULL
cov2[y2==0] <- rpois(sum(y2==0), lambda=40)
cov2[y2==1] <- rpois(sum(y2==1), lambda=60)
boxplot(cov2~y2)
```
<br>
#### 8-2. Penalty.factor
```
x.cov <- cbind(cov1, cov2, x)
pf <- c(0, 0, rep(1,ncol(x)) )
table(pf)

pf.fit <- glmnet(x.cov, y2, family="binomial", penalty.factor = pf, nlambda=10 )
pf.fit
```
  - We also adjust the eigen vecotrs or the n-th principle components to acommodate population structure. <br>

<br>

```
#save.image(file="part3.RData")
```
