---
layout: page
---

---

### III-1. Variant calling and joint genotuping using GATK HaplotypeCaller

In this part, we will follow the best practice of GATK
HaplotypeCaller. We used [GATK version
4.0.5.1](https://github.com/broadinstitute/gatk/releases/tag/4.0.5.1)
released on Jun 12, 2018.

---

#### a. Preambles for joint calling

- We will first create a directory to store output files.
> <pre>
$ cd ~/
$ mkdir --p out/s3.1 
$ cd out/s3.1/ </pre>

- We will use the PUR trio deeply sequenced repeatedly from two
  different centers.
> <pre>
$ ls ~/data/topmed </pre>

- Because they are CRAM files, it is important to set environment
  variable `REF_PATH`. 
> <pre>
$ export REF_PATH=~/data/ref/md5/%2s/%2s/%s 
$ env | grep REF_PATH </pre>

- To list the sample IDs, we can combine our UNIX skill sets as
  follows:
> <pre>
$ ls ~/data/topmed/ | grep cram$ | cut -d . -f 1 
NWD230091
NWD231092
NWD259170
NWD315403
NWD495157
NWD684137 </pre>
Make sure that you observe the same output

- You may use `xargs` to apply a series of command to each of these
  sample with a single command. For example,
> <pre>
$ ls ~/data/topmed/ | grep cram$ | cut -d . -f 1 | xargs -I {} ls -l ~/data/topmed/{}.recab.cram 
-rw-r--r-- 1 kobic kobic 38868223528 Jul 15 10:32 /BiO/home/kobic/data/topmed/NWD230091.recab.cram
-rw-r--r-- 1 kobic kobic 21030299096 Jul 15 11:04 /BiO/home/kobic/data/topmed/NWD231092.recab.cram
-rw-r--r-- 1 kobic kobic 22008268432 Jul 15 11:37 /BiO/home/kobic/data/topmed/NWD259170.recab.cram
-rw-r--r-- 1 kobic kobic 23927997911 Jul 15 12:14 /BiO/home/kobic/data/topmed/NWD315403.recab.cram
-rw-r--r-- 1 kobic kobic 23500707637 Jul 15 12:49 /BiO/home/kobic/data/topmed/NWD495157.recab.cram
-rw-r--r-- 1 kobic kobic 23135201299 Jul 15 13:24 /BiO/home/kobic/data/topmed/NWD684137.recab.cram </pre>
  * Q. Can you explain what the command did?
  
- Another example of using `xargs`:
> <pre>
$ ls ~/data/topmed/ | grep cram$ | cut -d . -f 1 | xargs -I {} echo "My name is {}"
My name is NWD230091
My name is NWD231092
My name is NWD259170
My name is NWD315403
My name is NWD495157
My name is NWD684137 </pre>
  * Q. Can you explain how this command works?

---

#### b. Creating per-sample GVCF files
	
- GATK HaplotypeCaller best practice recommends to create GVCF files
for each sample for joint variant calling. Due to time constrains, we
will focus the 200kb region starting from 10Mb of chromosome 20.

- Using `xargs`, we can run across all six samples with a single
command.
> <pre>
$ ls ~/data/topmed/ | grep cram$ | cut -d . -f 1 | \
xargs -I {} gatk HaplotypeCaller -R ~/data/ref/hs38DH.fa \
-I ~/data/topmed/{}.recab.cram \
-O {}.g.vcf.gz \
-ERC GVCF \
-L chr20:10000000-10200000 </pre>
Note that the command will take a while to finish.

- When everything finishes, take a look at what was generated in the
  current directory.
> <pre>
$ ls -l </pre>

- Take a look at what a GVCF look like..
> <pre>
$ zcat NWD230091.g.vcf.gz | grep -v ^## | less -S </pre>
  * Q. What is the difference between GVCF files from typical VCF files?

---

#### c. Consolidating the variant sites across multiple samples.

- The next step for joint calling is to merge the variant sites across
multiple samples. The new version of GATK HaplotypeCaller uses
_GenomicsDB_ to store collections of variant sites.

- Run the following command to merge the GVCF files to create a GenomicsDB
> <pre>$ gatk GenomicsDBImport \
-V NWD230091.g.vcf.gz \
-V NWD231092.g.vcf.gz \
-V NWD259170.g.vcf.gz \
-V NWD315403.g.vcf.gz \
-V NWD495157.g.vcf.gz \
-V NWD684137.g.vcf.gz \
--genomicsdb-workspace-path PURtrio \
--intervals chr20:10000000-10200000 </pre>

- Upon completion, take a peek at the GenomicsDB directory.
> <pre>$ ls PURtrio/ 
$ du -sh PURtrio/
</pre>
  * Do you expect the size of GenomicsDB will be large or small? How
    much, if you call whole genomes across 10,000 samples?

---

#### d. Joint genoptyping of the variant sites using GenomicsDB

- The joint genotyping can be performed using the GenomicsDB, without
requiring any other parameters, as follows:
> <pre>$ gatk GenotypeGVCFs \
-R ~/data/ref/hs38DH.fa \
-V gendb://PURtrio \
-O PURtrio.joint.vcf.gz \
-L chr20:10000000-10200000 </pre>

- Take a look at the output VCF file.
> <pre>$ zcat PURtrio.joint.vcf.gz | grep -v ^## | less -S </pre>
Use left and right arrows to navigate the columns across samples.
  * Q. Can you understand what each column means? 
  * Q. How about INFO fields? 
  * Q. How about FORMAT fields?

- To extract GT fields only across the 6 samples, we will use some
  hacky solution...
> <pre>$ zcat PURtrio.joint.vcf.gz | grep -v ^# | cut -f 10- | perl \
  -lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
  | sort | uniq -c | sort -n </pre>
  * Q. Considering these 6 samples are actually duplicate of a single
    trio, can you figure out which pairs of samples are identical? 
  * Q. Can you also figure out which one is offspring and which ones
    are parents?
  * Q. Does the genotypes mostly look Mendelian concordant or not?
  
---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's 
go to next step 
[III-2 Variant calling and joint genotuping using vt and cramore](../class-material/day1-vt-cramore.html), or 
go back to [Day 1 Overview](../day1).

---
