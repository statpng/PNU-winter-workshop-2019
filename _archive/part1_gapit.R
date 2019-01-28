source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest")
biocLite("snpStats")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d")
library(multtest); library(gplots); library(LDheatmap); library(genetics)
library(ape); library(EMMREML); library(compiler) 
library(scatterplot3d)

# GAPIT can be installed by this command
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

# set working directory
setwd("C:/Users/pc/Desktop/gapitex")

# read data
library(data.table)
raw_myY <- read.table(file=".csv", sep=",", header=T)
raw_myX <- as.data.frame(fread(file=".csv", sep=",", header=F))
# load("data_for_analyze.RData")
head(raw_myY); raw_myX[1:5,1:20]
raw_myY <- raw_myY[,1:2]
dim(raw_myY); dim(raw_myX)

## pre-processing
#remove NA Y
Yna <- is.na(raw_myY[,2])
myY <- raw_myY[!Yna,]
myX <- raw_myX[,c(rep(T,11),!Yna)]
dim(myY);dim(myX)

# check proprotion of NA
genedata <-  myX[-1,-c(1:11)]
na <- apply(genedata, 1, function(x) sum(x=='NN'))
length(na)
max(na); max(na)/ncol(genedata) 

# check one genotype per SNP
one <- apply(genedata, 1, function(x) length(unique(x[x!="NN"]))==1)
sum(one)
wone <- which(one) + 1
myX <- myX[-wone,]  
dim(myX)  #161798    381

### code 1: treating heterozygosity ###
setwd("C:/Users/pc/Desktop/gapitex/result1")
result1 <- GAPIT(Y=myY, G=myX)
summary(myY)

# heterozygosity(2min)
myGD= apply(myX[-1,-(1:11)], 1,
            function(one) GAPIT.Numericalization(one, bit=2, impute="Middle", Major.allele.zero=F))
H=1-abs(myGD-1)
het.ind=apply(H,1,mean)
het.snp=apply(H,2,mean)
hist(het.ind,col="gray", main="", xlab="Heterozygosity of individuals")
hist(het.snp,col="gray", main="", xlab="Heterozygosity of markers")

sum(het.ind > 0.5)
sum(het.snp > 0.5)
w1 <- which(het.snp > 0.5) + 1
myX <- myX[-w1,]
dim(myX)  #161677    381

# PCoA
library(ggplot2)
dim(myGD)
# load()
distance.matrix <- dist(myGD)
mds.stuff <- cmdscale(distance.matrix, eig=T, x.ret=T)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.var.per[1:5]
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample=1:nrow(myY),
                       X=mds.values[,1],
                       Y=mds.values[,2])
pdf("PCoA.pdf")
gg <- ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() + 
  theme_bw() + 
  xlab(paste("MDS1-", mds.var.per[1],"%", sep="")) +
  ylab(paste("MDS2-", mds.var.per[2],"%", sep="")) +
  ggtitle("MDS plot using Euclidean distance")
print(gg)
dev.off()


### code 2: imputation and controlling MAF ###
# SNP.MAF=0.02
setwd("C:/Users/pc/Desktop/gapitex/result2")
result2 <- GAPIT(Y=myY, G=myX, SNP.impute="Middle", SNP.MAF=0.02)



### code 3: optimal PC number ###
# SNP.MAF=0.02, Model.selection=TRUE
setwd("C:/Users/pc/Desktop/gapitex/result3")
result3 <- GAPIT(Y=myY, G=myX, 
                 SNP.MAF=0.02,
                 PCA.total=100,
                 Model.selection=TRUE)
# optimal number of PC using BIC
BICdata <- read.csv("C:/Users/pc/Desktop/gapitex/result3/GAPIT.MLM.Aspartic.BIC.Model.Selection.Results.csv",
         header=T)
head(BICdata)
BIC <- BICdata[,2]
BICdata[which.max(BIC),1]
# check population structure using genomic inflation factor, lambda
pvaldata <- read.csv("C:/Users/pc/Desktop/gapitex/result3/GAPIT.MLM.Aspartic.GWAS.Results.csv", header=T)
p <- pvaldata[,4]
if (stat_type == "Z")  z = Z[, 1]
if (stat_type == "CHISQ") z = sqrt(CHI[, 1])
if (stat_type == "PVAL")  z = qnorm(p / 2)
lambda = round(median(z^2) / 0.454, 3)
lambda


### code 4 : Compressed MLM ###
setwd("C:/Users/pc/Desktop/gapitex/result4")
result4 <- GAPIT(Y=myY, G=myX, 
                        SNP.MAF=0.02,
                        group.to=370,
                        group.from=1,
                        group.by=1,
                        model="CMLM"
                                    )



### code 5 : Enriched CMLM ###
setwd("C:/Users/pc/Desktop/gapitex/result5")
result5 <- GAPIT(Y=myY, G=myX,
                        SNP.MAF=0.02,
                        kinship.cluster=c("average", "complete", "ward.D"),
                        kinship.group=c("Mean", "Max"),
                        group.from=1,
                        group.to=370,
                        group.by=1,
                        memo="ECMLM"
) 


# small p-value
pvaldata <- read.csv("C:/Users/pc/Desktop/gapitex/result5/GAPIT.CMLM.Aspartic.GWAS.Results.csv", header=T)
head(pvaldata)
o <- order(pvaldata$P.value)
opvaldata <- pvaldata[o,]
head(opvaldata)


