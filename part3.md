---
layout: page
title: Part 3
permalink: /part3/
---


# Part 3. Selection probabilities.

<br>

Today, we will start from subset selection to selection probabilty analyzing
high-dimensional genomic sequencing data.

We will assume that you have basic R programming skills, which can be obtained
by taking a free online course
[Datacamp](https://www.datacamp.com/courses/free-introduction-to-r) or
[Coursera](https://www.coursera.org/courses?query=r%20programming).

---
<br>

### Schedule:

| Part    | Time                   | Topics                                                     |
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

### Introdution to selection probability
Regularization procedure, 특히 lasso,는 변수선택 측면에서 매우 좋은 성능을 보여줌.
그러나 최적의 tuning parameter를 선택하는 문제에 있어서 매우 신중해야함.
tuning parameter에 따라 coefficients 값에 차이가 있는 불안정성을 보여줌.
tuning parameter에 관계없이 변수 선택에 대한 지표를 나타내는
selection probability를 제안함.

selection probability는 다음과 같은 알고리즘을 통해 계산됨.
저희는 현재 elastic-net 모형에서 lambda와 alpha 두개의 tuning parameters를
갖고 있기 때문에, 모든 alpha와 lambda 그리드에 대하여
N/2 사이즈의 subsample을 선택하여 elastic-net 모형의 coefficients를 구하는
과정을 총 K번 반복함. 이 때, K는 많을수록 좋으며 최소 100이상의 값을
가지는 것을 추천함.
K번의 반복을 통해 각가의 tuning parameter grid에서 각 변수가 선택된
횟수의 최댓값을 K로 나누어 selection probability를 계산함.
그리고 tuning parameter grid 중에서 sp의 최대값을 최종 sp로 정의함.

간단하게 말해서 selection probability는 총 K번 반복 중에서
몇번 선택되는지 비율의 tuning parameter에 대한 최댓값을 뜻함.

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
|6: | $$SP_j^\Lambda = \frac{1}{K}\#\{k<=K: \hat{\beta}_j^\Lambda(I_k) \ne 0 \} $$ |
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
```{r}
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
