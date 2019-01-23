---
layout: page
---

---

### I-1. Variant Filtering

In this part, we will perform a simple variant filtering by using the
TOPMed Freeze 5 variant list

---

#### a. Comparing GATK and vt/cramore callsets

- We've created two variant callsets spanning chromosome 20 10-11Mb
region. This is the same type of callset introduced yesterday, but
will 1Mb region instead of 200kb.

- The GATK callset can be found at
><pre>
ls -l ~/data/topmed/PURtrio.gatk.joint.chr20.10Mb.11Mb.vcf.gz</pre>

- The vt/cramore callset (based on TOPMed calling pipeline) can be
  found at
><pre>
ls -l ~/data/topmed/PURtrio.at.joint.chr20.10Mb.11Mb.bcf</pre>  

- Compare the basic summary statistics between the two BCF/VCF files
  using `bcftools stats' command:
><pre>
bcftools stats ~/data/topmed/PURtrio.gatk.joint.chr20.10Mb.11Mb.vcf.gz
  | less -S </pre>
><pre>
bcftools stats ~/data/topmed/PURtrio.vt.joint.chr20.10Mb.11Mb.bcf
  | less -S </pre>
  
  * Q. What are the differences between the two in terms of the total
    variant count, SNP count, indel count, and the multi-allelic
    variant count?
  * Q. How does the transition to trasverstion ratio (Ts/Tv) compare
    between the two? What is your evaluation?
  * Q. Can you explain the distribution of allele frequencies?
	
#### b. Applying external filter from TOPMed 65k callset.

- When making a variant calls with a small number of samples in a
  small region, applying VQSR or SVM filter may not be accurate due to
  small training sample size.

- In such a case, harnessing the variant filters trained by
  large-scale sequencing studies, such as TOPMed, could be helpful. We
  copied the site list from >65,000 subjects (from BRAVO variant
  browser). To see the summary of the file, try the following command:
><pre>bcftools view
  ~/data/bravo/chr20.TOPMed_freeze5_bravo_62784.vcf.gz | less -S</pre>
  * Q: How does the FILTER column look like? Can you explain what it
    means?

- We will use `vcf-apply-external-filter` software (from `apigenome`
  package) to transfer the TOPMed 65k callset filter to our
  callset. First, let's create a directory for storing output
><pre>cd ~/
mkdir --p out/s5.1
cd out/s5.1</pre>

- See how to use `vcf-apply-external-filer`, then type the command
><pre>vcf-apply-external-filter</pre>
or use the `-man option`
><pre>vcf-apply-external-filter -man</pre>

- To apply the software to the GATK callset, use the following
  command:
><pre>vcf-apply-external-filter \
-vcf ~/data/topmed/PURtrio.gatk.joint.chr20.10Mb.11Mb.vcf.gz \
-ext ~/data/bravo/chr20.TOPMed_freeze5_bravo_62784.vcf.gz \
-out PURtrio.gatk.joint.chr20.10Mb.11Mb.filt.vcf.gz \
--clear --region chr20:10000000-11000000</pre>

- A similar command may be applied for `vt/cramore` callset.
><pre>vcf-apply-external-filter \
-vcf ~/data/topmed/PURtrio.vt.joint.chr20.10Mb.11Mb.bcf \
-ext ~/data/bravo/chr20.TOPMed_freeze5_bravo_62784.vcf.gz \
-out PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf \
--clear --region chr20:10000000-11000000</pre>

- You can peek these files using `bcftools view`
><pre>
bcftools view -H PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf |less -S</pre>

#### c. Evaluating the filtered variants.

- One simple way to evaluate how good the filtered variants look like
  is to examine the Mendelian consistency of filtered genotype for
  this example. To see only PASS-filtered variants, try 
><pre>bcftools view -H PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf | grep -w PASS | less -S </pre>
To see variants failed any filter, try
><pre>bcftools view -H PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf | grep -v -w PASS | less -S </pre>

- To see the patterns of genotypes, we can reuse the complex one-liner
  from yesterday. For PASS-filtered variants, try
><pre>bcftools view -H PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf | \
grep -w PASS | cut -f 10- | perl \
-lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
| sort | uniq -c | sort -n</pre>
For the rest of variants, try
><pre>bcftools view -H PURtrio.vt.joint.chr20.10Mb.11Mb.filt.bcf | \
grep -v -w PASS | cut -f 10- | perl \
-lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
| sort | uniq -c | sort -n</pre>
  * Q. What are the differences between the two? How do you interpret
    the difference?

- Same things could be applied to the GATK callset.
  For PASS-filtered variants, try
><pre>bcftools view -H PURtrio.gatk.joint.chr20.10Mb.11Mb.filt.vcf.gz | \
grep -w PASS | cut -f 10- | perl \
-lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
| sort | uniq -c | sort -n</pre>
For the rest of variants, try
><pre>bcftools view -H PURtrio.gatk.joint.chr20.10Mb.11Mb.filt.vcf.gz | \
grep -v -w PASS | cut -f 10- | perl \
-lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
| sort | uniq -c | sort -n</pre>
  * Q. What are the differences in the results compared to
    `vt/cramore` call? Which callset appear to be better after
    applying filtering. Can you explain why?

- To produce summary statistics for only pass-filtered variants, you
  can try
><pre>bcftools view -f PASS PURtrio.gatk.joint.chr20.10Mb.11Mb.filt.vcf.gz | bcftools stats - | less -S</pre>
for GATK callset. For vt/cramore callset, try
><pre>bcftools view -f PASS PURtrio.gatk.joint.chr20.10Mb.11Mb.filt.bcf | bcftools stats - | less -S</pre>
  * Q. What differences do you see? Based on Ts/Tv statistics, which
    one appears to be better?

---
Feel free to ask questions to your instructor(s) if you are stuck. 
, or go back to [Day 2 Overview](../day2).

---
