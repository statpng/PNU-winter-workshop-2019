---
layout: page
---

---

### III-1. RNA-seq data from GEUVADIS

In this part, we will use RNA-seq data generated from GEUVADIS data and prepare for the eQTL
association analysis with sequencing variants from the 1000 Genomes Project. In particular,
we will extract quantification data of a specific gene, ENSG00000077044 (ensemble ID), from 
the GEUVADIS project. 

Data from the GEUVADIS project has been downloaded onto data/1000g/geuvadis. 
> <pre>
# We will use normalized and QC'ed per-gene quantification data
less -S ~/data/1000g/geuvadis/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz 
# Let's try to extract our gene of interest, ENSG00000077044, from this file
zgrep ENSG00000077044 ~/data/1000g/geuvadis/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz 
# Now, let's try to extract expression levels 
zgrep ENSG00000077044 ~/data/1000g/geuvadis/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz  |tr '\t' '\n'|head
# Extract sample IDs from header line
zcat ~/data/1000g/geuvadis/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz|head -1|tr '\t' '\n'|head </pre>

This GEUVADIS dataset includes two population groups, 373 EUR and 89 YRI. We could try to adjust for
population stratification issues using generalized mixed-models and covariates, for simplicity of the
practice, we will extract only the EUR samples for eQTL analysis.

Population and sex information of the samples are available from the 1000 Genomes Projects sample pedigree info:
> <pre>
less ~/data/1000g/panel/20130606_g1k.ped </pre>

We should extract sample IDs from EUR populations only (CEU, FIN, GBR, IBS, TSI) and then also need to assign sex
information for the anlaysis. This can be done with simple scripting (Perl or Python) but may be outside of the scope
of this workshop, so I made a pre-generated phenotype file that can be used for the eQTL analysis:
> <pre>
cp ~/data/1000g/geuvadis/analysis_results/ENSG00000077044.ped .  </pre>

Now, let's prepare genotype file for the analysis. For simplicity, we will only analyze cis-eQTL regions of the gene. 
This specific gene's genomic location is chr2:233,354,507-233,472,104 in GRCh38. We will extract +/- 1Mb from the genotype
file around the gene. We will use 1000 Genomes Project Phase 3 VCF.

> <pre>
tabix -h ~/data/1000g/panel/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	2:232354507-234472104 > cis_ENSG00000077044.vcf </pre>

Now we have genotype and phenotype files ready, so we can move on to the association analysis

---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's go to next step 
[III-2. eQTL analysis](../class-material/day2-eqtl)
, or go back to [Day 2 Overview](../day2).

---
