---
layout: page
---

---

### II-1. Phasing Variant Calls

In this part, we will perform a simple haplotype-aware analysis from
sequence-based genotype calls.

---

#### a. Preparing input genotype files.

- We will start with the vt/cramore variant calls from the previous
  session.
><pre>cd ~/out/s5.1
ls -l</pre>

- We will focus on PASS-filtered genotypes to produce high-quality haplotypes.
><pre>
bcftools view -f PASS -Ob \
-o PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.bcf\
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf</pre>

- We will extract genotypes only, setting genotypes with depth less
  than 10 as missing.
><pre>
cramore vcf-squeeze --minDP 10 \
--in PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.bcf \
--out PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.bcf</pre>

- Now, looking at the patterns of genotypes is easier:
><pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.bcf | cut -f 10- \
| sort | uniq -c </pre>

#### b. Phasing genotypes using EAGLE2 without reference panel.

- We will first phase the genotypes using `eagle2` software tool, without
  using a reference panel.
> <pre>eagle \
--geneticMapFile=$HOME/data/ref/genetic_map_hg38_withX.txt.gz \
--chrom=chr20 \
--outPrefix=PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_without_ref \
--vcf PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.bcf </pre>

- You can have a peek of the results using the following command:
> <pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_without_ref.vcf.gz \
| less -S</pre>

- To see the haplotypes only, try
> <pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_without_ref.vcf.gz \
| cut -f 10- | less -S</pre>

- We can compare how the phased haplotypes are concordant between the
  duplicate pairs.
> <pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_without_ref.vcf.gz \
| cut -f 10,12 | sort | uniq -c </pre>
  * Q. How concordant are the phased haplotypes between the duplicated
    pair? 
  * Q. Do you think that this version of phased haplotype is of high
    quality? Why or why not? 
	
#### c. Phasing genotypes using EAGLE2 with 1000G reference panel

- We will phase the genotypes using `eagle2` software tool, by
  leveraging the 1000 Genomes reference panel.
><pre> eagle \
--geneticMapFile=$HOME/data/ref/genetic_map_hg38_withX.txt.gz \
--chrom=chr20 \
--outPrefix=PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_with_ref\
--vcfTarget PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.bcf \
--vcfRef $HOME/data/1000g/panel/ALL.chr20_GRCh38.genotypes.20170504.bcf \
--bpStart 10000000 --bpEnd 11000000</pre>

- Agin, to see the phased haplotypes only, try
> <pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_with_ref.vcf.gz \
| cut -f 10- | less -S</pre>

- We can compare how the phased haplotypes are concordant between the
  duplicate pairs.
> <pre>
bcftools view -H \
PURtrio.vt.joint.chr20.10Mb.11Mb.filt.PASS.minDP10.phased_with_ref.vcf.gz \
| cut -f 10,12 | sort | uniq -c </pre>
  * Q. How concordant are the phased haplotypes between the duplicated
    pair? 
  * Q. Do you think that this version of phased haplotype is of high
    quality? Why or why not?
  * Q. When you have small sample sizes, would you prefer to phase
    your genotypes with or without reference panel? Why?

---
Feel free to ask questions to your instructor(s) if you are stuck. 
, or go back to [Day 2 Overview](../day2).

---
