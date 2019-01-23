---
layout: page
---

---

### IV-1. Single cell RNA-seq analysis examples

WARNING: This step will consume >30GB of memory space, so it is
impossible that everyone run this session simultaneously. Therefore,
we will ask only a
small number of people run these steps during the practical
session, and each individual can practice these steps later on.

These steps follows the excellent guide on single cell analysis
available at
[https://github.com/broadinstitute/BipolarCell2016](https://github.com/broadinstitute/BipolarCell2016). 

---

#### a. Preparing the analysese. 

- We will first need to set up directory to store output files, and
run R
><pre>cd ~/
mkdir --p ~/out/s8
cd out/s8
R
ls -l</pre>

- First, install some packages that will be required for the analysis.
><pre>install.packages(c("ggplot2","grid","Matrix","gmodels","RANN","reshape"))</pre>
It is recommended to select Korean (Seoul) download site when
prompted.

- Next, load the prepared data and R scripts
><pre>load("~/data/10x/GSE81904_bipolar_data_Cell2016.Rdata")
source("~/data/10x/class.R")</pre>

- The dimensions of digital expression matrix is
><pre>dim(bipolar_dge) 
## [1] 24904 44994</pre>

#### b. Quality control of the batches of single cell data

- Typically, cells with a very large proportion of mitochontrial genes
are likely to contain dead cells or reads with technical artifacts. We
will filter out cells that has more than 10% of mitochondrial reads
among the total reads.
> <pre>mt.genes = grep("mt-", rownames(bipolar_dge), value = TRUE)
cells.use = colnames(bipolar_dge)[colSums(bipolar_dge[mt.genes, ])/colSums(bipolar_dge) < 0.1]
bipolar_dge = bipolar_dge[, cells.use]
dim(bipolar_dge) </pre>
  * Q. How many cells are filtered out by this stage?
  
- We can visualize the number of transcripts per cell, stratified by
  batch
><pre>dsq.bip=scDrop(count.data=bipolar_dge)
dsq.bip=initialize(dsq.bip, min.genes = 500, min.cells = 30, min.counts=60) 
print(dim(dsq.bip@data)) 
table(dsq.bip@meta$sample)
pdf('violin.pdf',width=8,height=5)
violinplot(dsq.bip,c("num.genes","num.trans"))
dev.off()</pre>

- To copy the PDF file, type the following commnad from the terminal screen.
><pre>scp kobic00@IP.AD.DR.ESS:out/s8/violin.pdf . </pre>

  
#### c. Visualizing cells with t-SNE manifolds
- The t-stochastic neighbor embedding (t-SNE) manifolds are already
produced, and they can be visualized as follows:
><pre>dsq.bip@tsne.y = tsne.y
data.plot = dsq.bip@tsne.y
data.plot$group = dsq.bip@group
pdf('tsne.pdf',width=8,height=5)
ggplot(data.plot, aes(x = tSNE1, y = tSNE2)) + geom_point(aes(colour = factor(group), size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw() 
dev.off()</pre>

- To show the gene expression for each individual gene, by
representing each individual droplet in the t-SNE manifold, one can
run the following command:
><pre>pdf('scatter.pdf',width=8,height=8)
gene.expression.scatter(dsq.bip, c("Prkca", "Glul","Scgn", "Grm6")) 
dev.off()</pre>

- To copy the output PDF files, type the following commnad from the terminal screen.
><pre>scp kobic00@IP.AD.DR.ESS:out/s8/*.pdf . </pre>

---
Feel free to ask questions to your instructor(s) if you are stuck. 
, or go back to [Day 2 Overview](../day2).

---
