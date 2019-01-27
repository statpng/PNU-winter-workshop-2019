---
layout: page
title: Part 3
permalink: /part3/
---


# Part 3. Selection probabilities.

<hr>
<br>

### Schedule:

| Part    | Title                   | Topics                                                     |
| :-----: |   :--------------:     | :-----------------------                                   |
| III     | Selection probabilties | **An algorithm of selection probabilities**                |
|         |                        | &nbsp; &nbsp; - Setting a grid of tuning parameters          |
|         |                        | &nbsp; &nbsp; - Applying the regularization for subsamples   |
|         |                        | **Why split data with 0.5 proportion**                     |
|         |                        | &nbsp; &nbsp; - Advantages of "subagging"                    |
|         |                        | **Algorithm summary**                                      |
|         |                        | **The stability path**                                     |
|         |                        | **Threshold to control the false positive**                |
|         |                        | **Manhattan plot with selection probabilities**            |

<br>

### Instructors:

Kipoong Kim

---

<link rel="stylesheet" href="//cdnjs.cloudflare.com/ajax/libs/highlight.js/9.6.0/styles/ir-black.min.css">
<script src="//cdnjs.cloudflare.com/ajax/libs/highlight.js/9.6.0/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

<script>
  (function () {
    var script = document.createElement("script");
    script.type = "text/javascript";
    script.src  = "https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML";
    document.getElementsByTagName("head")[0].appendChild(script);
  })();
</script>

<br>

### 0. Data setup (equivalent to part2)
```
library(glmnet)
library(qqman)

load("[Data]PNU-winter-workshop-2019.RData")

attach(workshop.data)
```
<hr><br>

### 1. An algorithm of selection probabilities

||**Algorithm** Selection probabilities with elastic-net |
|--:|:------------ |
|0: | Let us assume that a genomic data has $$n$$ samples and $$p$$ variables. |
|------------- |
|1: | For all $$\Lambda=(\alpha, \lambda)$$, $$\text{where} ~~ \alpha \in [0,1], \lambda>0$$ |
|------------------- |
|2: | **for** k=1 to K **do** |
|------------------- |
|3: | &nbsp; &nbsp; Subsample $$I_k$$ with size $$[n/2]$$ |
|----------- |
|4: | &nbsp; &nbsp; Compute $$\hat{\beta}_j^{\Lambda}(I_k)$$ with regularization model |
|----------- |
|5: | **end for** |
|----------- |
|6: | $$SP_j^\Lambda = \frac{1}{K}\#\{k\le K: \hat{\beta}_j^\Lambda(I_k) \ne 0 \} $$ |
|------------- |
|7: | $$SP_j = \underset{\Lambda}{\max}SP_j^\Lambda,~~ j=1,\cdots, p$$ |
|------------- |
|8: | **return** **SP**=$$(SP_1, \cdots, SP_p)$$ |
|------------------- |

<hr><br>
### 2. Calculation of selection probabilities
#### 2-1. Setting a grid of tuning parameters
```
alpha.grid <- seq(0.1, 0.9, 0.1)
if( FALSE ){
  lambda.mat <- array( NA, dim = c(length(alpha.grid), 10, 10),
                   dimnames=list(paste0("alpha=", alpha.grid),
                                 paste0("lambda", 1:10),
                                 paste0("REP=", 1:10)) )

  for( REP in 1:10 ){
    for( i in 1:length(alpha.grid) ){
      set.seed( 1000*REP + i )
      sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
      lambda.mat[i,,REP] <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian",
                               alpha=alpha.grid[i], nlambda=10)$lambda
    }
  }
}
# save(lambda.mat, file="[object]lambda.mat.RData")

load("[object]lambda.mat.RData")
lambda.grid <- seq( summary(lambda.mat)[2], summary(lambda.mat)[4], length.out = 100 )
range(lambda.grid)
```

<br>

#### 2-2. Applying regularization with the same subsamples for a grid of tuning paramters.
```
if( FALSE ){
  sp.array <- array(0, dim=c(ncol(x), length(alpha.grid), length(lambda.grid)) )
  K <- 50 ## recommended to be 100 or more
  for( k in 1:K ){
    if( k %% 10 == 0 ) cat("A current iteration is ", k, "\n")
    set.seed(123*k)
    sp.subset <- sample( nrow(x), floor(nrow(x)/2) )

    for( i in 1:length(alpha.grid) ){
      for( j in 1:length(lambda.grid) ){
        sp.fit <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian",
                         alpha=alpha.grid[i], lambda=lambda.grid[j])
        sp.array[,i,j] <- sp.array[,i,j] + as.numeric(sp.fit$beta!=0)/K
      }
    }
  }
}
# save(sp.array, file="[object]sp.array.RData")
```
<br>

#### 2-3. Retrieving the selection probabilities.
```
load("[object]sp.array.RData")

sp <- apply( sp.array, 1, max )

head( data.frame(sp = sort( sp, decreasing=TRUE ), variable = order( sp, decreasing=TRUE ) ) )
```
<br>

#### 2-4. Summary of an algorithm of selection probabilities
```
# library(glmnet)
# alpha.grid <- seq(0.1, 0.9, 0.1)
#
# lambda.mat <- array( NA, dim = c(length(alpha.grid), 10, 10),
#                  dimnames=list(paste0("alpha=", alpha.grid), paste0("lambda", 1:10), paste0("REP=", 1:10)) )
# for( REP in 1:10 ){
#   for( i in 1:length(alpha.grid) ){
#     set.seed( 1000*REP + i )
#     sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
#     lambda.mat[i,,REP] <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian", alpha=alpha.grid[i], nlambda=10)$lambda
#   }
# }
#
# lambda.grid <- seq( summary(lambda.mat)[2], summary(lambda.mat)[4], length.out = 100 )
# lambda.grid
#
# sp.array <- array(0, dim=c(ncol(x), length(alpha.grid), length(lambda.grid)) )
# K <- 50 ## recommended to be 100 or more
# for( k in 1:K ){
#   set.seed(123*k)
#   sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
#   
#   for( i in 1:length(alpha.grid) ){
#     for( j in 1:length(lambda.grid) ){
#       sp.fit <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian", alpha=alpha.grid[i], lambda=lambda.grid[j])
#       sp.array[,i,j] <- sp.array[,i,j] + as.numeric(sp.fit$beta!=0)/K
#     }
#   }
# }
#
# sp <- apply( sp.array, 1, max )
```
<hr><br>

### 3. Threshold($$\pi_thr$$) to control the false positive
```
sp.threshold <- function(sp, FP=5){
  if(max(sp)>1 | min(sp)<0) stop("'Sel.prob' should be in [0,1].")
  if(length( dim(sp) ) < 3) stop("'sp' should be a 3-dimensional array.")
  qhat <- sum(sp)/(dim(sp)[2]*dim(sp)[3])
  p <- dim(sp)[1]
  threshold <- qhat^2/(2*FP*p) + 0.5

  if(threshold>1) threshold <- 1
  return( threshold )
}

theoretical.threshold5 <- sp.threshold(sp.array, FP=5)
theoretical.threshold10 <- sp.threshold(sp.array, FP=10)

c(theoretical.threshold5, theoretical.threshold10)
```

<hr><br>
### 4. Manhattan plot with sel.prob $$>= \pi_thr$$
```
df.sp <- data.frame( ref, p=sp)
colnames(df.sp) <- c("SNP", "CHR", "BP", "P")

var.theoretical.threshold5 <- which( sp > theoretical.threshold5 )

qqman::manhattan(df.sp,
                 chr="CHR", bp="BP", p="P", snp="SNP", logp=FALSE,
                 suggestiveline=c(theoretical.threshold5, theoretical.threshold10),
                 highlight=df.sp[var.theoretical.threshold5, "SNP"],
                 ylim=c(0,1))
```
<hr><br>

### 5. Stability selection
#### 5-1. Computing the selection probabilities of each variable for a specific set of $$\lambda$$
```
## Stability selection
stable.lambda <- glmnet(x=x, y=y, family="gaussian", alpha=0.1, nlambda=10)$lambda

if( FALSE ){
  stable.sp.out <- array(0, dim=c(ncol(x), length(stable.lambda)) )
  K=50
  for( k in 1:K ){
    if(k%%10==0) cat("A current iteration is ", k, "\n")
    set.seed(123*k)
    sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
    for( j in 1:length(stable.lambda) ){
      sp.fit <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian",
                       alpha=0.1, lambda=stable.lambda[j])
      stable.sp.out[,j] <- stable.sp.out[,j] + as.numeric(sp.fit$beta!=0)/K
    }
  }
}
# save(stable.lambda, file="[object]stable.lambda.RData")
# save(stable.sp.out, file="[object]stable.sp.out.RData")
```
<br>

#### 5-2. The stability path
```
load("[object]stable.lambda.RData")
load("[object]stable.sp.out.RData")
var.threshold10 <- which( sp > theoretical.threshold10 )
str(var.threshold10)

matplot( log(stable.lambda), t(stable.sp.out), type="l",
         lty=ifelse(1:ncol(x) %in% var.threshold10, 1, 2),
         col=ifelse(1:ncol(x) %in% var.threshold10, 2, 1),
         xlab="log(lambda)", ylab="Selection probabilities" )
```
<br>

```{r}
# save.image( file="part4.RData" )
```

<hr><br>
