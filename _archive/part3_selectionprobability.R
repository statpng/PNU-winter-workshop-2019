
## Introdution to selection probability
library(glmnet)

load("[Data]PNU-winter-workshop-2019.RData")

attach(workshop.data)

## Introdution to selection probability

library(glmnet)
alpha.grid <- seq(0.1, 0.9, 0.1)
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
# save(lambda.mat, file="[object]lambda.mat.RData")

load("[object]lambda.mat.RData")
lambda.grid <- seq( summary(lambda.mat)[2], summary(lambda.mat)[4], length.out = 100 )
range(lambda.grid)

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
# save(sp.array, file="[object]sp.array.RData")

load("[object]sp.array.RData")

sp <- apply( sp.array, 1, max )

head( data.frame(sp = sort( sp, decreasing=TRUE ), variable = order( sp, decreasing=TRUE ) ) )


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


df.sp <- data.frame( ref, p=sp)
colnames(df.sp) <- c("SNP", "CHR", "BP", "P")

var.theoretical.threshold5 <- which( sp > theoretical.threshold5 )

qqman::manhattan(df.sp, 
                 chr="CHR", bp="BP", p="P", snp="SNP", logp=FALSE, 
                 suggestiveline=c(theoretical.threshold5, theoretical.threshold10),
                 highlight=df.sp[var.theoretical.threshold5, "SNP"], 
                 ylim=c(0,1))





## Stability selection
stable.lambda <- glmnet(x=x, y=y, family="gaussian", alpha=0.1, nlambda=10)$lambda

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
# save(stable.lambda, file="[object]stable.lambda.RData")
# save(stable.sp.out, file="[object]stable.sp.out.RData")


load("[object]stable.lambda.RData")
load("[object]stable.sp.out.RData")
var.threshold10 <- which( sp > theoretical.threshold10 )
str(var.threshold10)

matplot( log(stable.lambda), t(stable.sp.out), type="l",
         lty=ifelse(1:ncol(x) %in% var.threshold10, 1, 2), 
         col=ifelse(1:ncol(x) %in% var.threshold10, 2, 1),
         xlab="log(lambda)", ylab="Selection probabilities" )



# save.image( file="part3.RData" )



