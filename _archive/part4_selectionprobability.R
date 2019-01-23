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
dim(lambda.mat)

# dimnames(lambda.mat)=list( paste0("alpha=",ALPHA), paste0("lambda",1:10), paste0("REP=", 1:10) )

for( REP in 1:10 ){
  for( i in 1:length(alpha.grid) ){
    set.seed( 1000*REP + i )
    sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
    lambda.mat[i,,REP] <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian", 
                             alpha=alpha.grid[i], nlambda=10)$lambda
  }
}

lambda.grid <- seq( summary(lambda.mat)[2], summary(lambda.mat)[4], length.out = 100 )
lambda.grid



sp.array <- array(0, dim=c(ncol(x), length(alpha.grid), length(lambda.grid)) )
K <- 50 ## recommended to be 100 or more
for( k in 1:K ){
  cat("A current iteration is ", k, "\n")
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

sp <- apply( sp.array, 1, max )


str(sp)
summary(sp)


hist(sp)




lm.pvalue <- NULL
for( j in 1:ncol(x) ){
  lm.coef <- summary( lm(y~x[, j]) )$coef
  if( nrow(lm.coef)==2 ){
    lm.pvalue[j] <- lm.coef[2, 4]
  } else if ( nrow(lm.coef)==1 ){
    lm.pvalue[j] <- NA
  }
}


df.lm.pvalue <- data.frame( ref, p=lm.pvalue)
colnames(df.lm.pvalue) <- c("SNP", "CHR", "BP", "P")
list(ref=head(ref), df.lm.pvalue=head(df.lm.pvalue))


data.frame(sp = sort( sp, decreasing=TRUE ), 
           variable = order( sp, decreasing=TRUE )
)[1:5,]


sp.threshold <- function(sp, FP=5){
  if(max(sp)>1 | min(sp)<0) stop("Sel.prob should be in [0,1].")
  if(length( dim(sp) ) < 3) stop("sp should be a 3-dimensional array.")
  qhat <- sum(sp)/(dim(sp)[2]*dim(sp)[3])
  p <- dim(sp)[1]
  threshold <- 1/FP * qhat^2/(2*p) + 0.5
  return( threshold )
}

theoretical.threshold <- sp.threshold(sp.array/K)

theoretical.threshold

bonferroni.cutoff <- 0.05/ncol(x)

var.theoretical.threshold <- which( sp > theoretical.threshold )
qqman::manhattan(df.lm.pvalue, 
                 chr="CHR", bp="BP", p="P", snp="SNP", 
                 suggestiveline=-log10(bonferroni.cutoff),
                 highlight=df.lm.pvalue[var.theoretical.threshold, "SNP"], 
                 ylim=c(0,10))





## Stability selection
stable.lambda <- glmnet(x=x, y=y, family="gaussian", alpha=0.1, nlambda=10)$lambda

# stable.enet <- glmnet(x=x, y=y, family="gaussian", alpha=0.1, lambda=stable.lambda)

stable.sp.out <- array(0, dim=c(ncol(x), length(stable.lambda)) )
K=50
for( k in 1:K ){
  cat("A current iteration is ", k, "\n")
  set.seed(123*k)
  sp.subset <- sample( nrow(x), floor(nrow(x)/2) )
  for( j in 1:length(stable.lambda) ){
    sp.fit <- glmnet(x=x[sp.subset,], y=y[sp.subset], family="gaussian", 
                     alpha=0.1, lambda=stable.lambda[j])
    stable.sp.out[,j] <- stable.sp.out[,j] + as.numeric(sp.fit$beta!=0)/K
  }
}
summary(stable.sp.out)

candidate <- c(1163,1366,1992,4287,5544,5548,5675,
               6456,7340,8494,8620,8903,9015,9892,
               10095,11313,11561,12294,13496)


# png(file="stable.enet.solutionpath.png", height=4, width=7, unit="in", res=300)
# matplot( log(stable.enet$lambda), t(stable.enet$beta), type="l", 
#          lty=ifelse(1:ncol(x) %in% candidate, 1, 2), 
#          col=ifelse(1:ncol(x) %in% candidate, 2, 1),
#          xlab="log(lambda)", ylab="Coefficients" )
# dev.off()

png(file="stable.sp.solutionpath.png", height=4, width=7, unit="in", res=300)
matplot( log(stable.lambda), t(stable.sp.out), type="l",
         lty=ifelse(1:ncol(x) %in% candidate, 1, 2), 
         col=ifelse(1:ncol(x) %in% candidate, 2, 1),
         xlab="log(lambda)", ylab="Selection probabilities" )
dev.off()


save.image( file="part4.RData" )


