load("imputed_data.RData")
set.seed(11)
nc <- ncol(imputed.data$x)
rand.sample <- sample( nc, 15000 )
Y <- imputed.data$y$Total[ !is.na( imputed.data$y$Total ) ]
X <- imputed.data$x[ !is.na( imputed.data$y$Total ), rand.sample ]
REF <- imputed.data$a[ rand.sample, ]
workshop.data <- list( x=X, y=Y, ref=REF )

save(workshop.data, imputed.data, file="[data]PNU-winter-workshop-2019.RData")
