---
layout: page
---

---

### I-2. Variant Annotation

In this part, we will use SnpEff to annotate list of variants in a VCF file.

---

#### Preparing input (site-only) VCF

For annotation, we do not need 'genotype' fields of VCF, which could be quite large depending on
the sample size. It is more efficient to make 'site-only' VCF by cutting first 8 columns.

We will use 1000 Genomes Project's chr20 VCF file at ```data/1000g/panel/```

> <pre>
ls data/1000g/panel/
# We will use the ALL.chr20_GRCh38.genotypes.20170504.vcf.gz file 
# Because this VCF file contains ~1.8M variants, we will extract only small region
tabix -h data/1000g/panel/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz 20:10000000-15000000 \
	|cut -f 1-8 > ALL.chr20_GRCh38.part.sites.vcf
less ALL.chr20_GRCh38.part.sites.vcf 
grep -v ^# ALL.chr20_GRCh38.part.sites.vcf | wc -l</pre>

- Q. What columns did you find in the sites VCF?
- Q. How many variants are in the partial chr20 VCF?
- Q. Why ```'grep -v ^#'``` ?

#### Running SnpEff

SnpEff includes options to select different databases (species, reference build, gene database). 

> <pre>
# Checking list of databases supported by the current version of SnpEff
java -jar /BiO/apps/snpEff/snpEff.jar databases | less
# We need GRCh38 for our VCF file, so check out which versions are included
java -jar /BiO/apps/snpEff/snpEff.jar databases | grep GRCh38 </pre>

We will use Ensenble version of GRCh38 database, GRCh38.86

> <pre>
java -Xmx8g -jar /BiO/apps/snpEff/snpEff.jar GRCh38.86 ALL.chr20_GRCh38.part.sites.vcf > ALL.chr20_GRCh38.part.annotated.sites.vcf 
less -S ALL.chr20_GRCh38.part.annotated.sites.vcf</pre>

- Q. How many STOP-GAIN variants exist in this VCF?

---

Feel free to ask questions to your instructor(s) if you are stuck. 
, or go back to [Day 2 Overview](../day2).

---
