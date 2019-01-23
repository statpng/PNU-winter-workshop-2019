---
layout: page
---

---

### III-2. Variant calling and joint genotuping using vt and cramore

In this part, we will follow the best practice of TOPMed variant
calling pipeline. We used a freeze of [vt software
tool](https://github.com/hyunminkang/vt-topmed)
and [cramore software
tool](https://github.com/hyunminkang/apigenome) as of today.
We will also use `samtools` and `bamUtil` software tools. 

---

#### a. Creating per-sample VCF file

- We will first create another output directory, and check `REF_PATH` variable.
> <pre>$ cd ~/
$ mkdir --p out/s3.2
$ cd out/s3.2
$ env | grep REF_PATH
REF_PATH=/BiO/home/kobic/data/ref/md5/%2s/%2s/%s </pre>

- Using `xargs`, we process each sample sequentially with a single
  command:
> <pre>$ ls ~/data/topmed/ | grep cram$ | cut -d . -f 1 | xargs -I {} \
bash -c "samtools view -uh ~/data/topmed/{}.recab.cram \
  chr20:10000000-10200000 2> /dev/null | bam clipOverlap --poolSize \
  100000000 --in -.ubam --out -.ubam 2> /dev/null | vt discover2 -z -q \
  20 -b + -r ~/data/ref/hs38DH.fa -s {} -o {}.bcf" </pre>
  
- Then you will see six BCF files are created.
> <pre>$ ls -l 
total 228
-rw-rw-r-- 1 kobic kobic 36931 Jul 16 04:04 NWD230091.bcf
-rw-rw-r-- 1 kobic kobic 36294 Jul 16 04:04 NWD231092.bcf
-rw-rw-r-- 1 kobic kobic 36394 Jul 16 04:04 NWD259170.bcf
-rw-rw-r-- 1 kobic kobic 36851 Jul 16 04:04 NWD315403.bcf
-rw-rw-r-- 1 kobic kobic 37302 Jul 16 04:04 NWD495157.bcf
-rw-rw-r-- 1 kobic kobic 37430 Jul 16 04:04 NWD684137.bcf </pre>

- Take a peek at what the individual BCF looks like.
> <pre>$ bcftools view -H NWD230091.bcf | less -S </pre>

- We need to index individual BCF files
> <pre>$ ls *.bcf | xargs -I {} bcftools index {} 
$ ls</pre>

---

#### b. Consolidating the variant sites across multiple samples.

- The next step for joint calling is to merge the variant sites across
multiple samples. To do so, we need to create a file containing all
BCF files to be merged.
> <pre>$ ls *.bcf > bcflist.txt 
$ cat bcflist.txt 
NWD230091.bcf
NWD231092.bcf
NWD259170.bcf
NWD315403.bcf
NWD495157.bcf
NWD684137.bcf</pre>

- Then we can create a merged site list using the following command:
> <pre>$ cramore vcf-merge-candidate-variants --in-vcf-list bcflist.txt \
--out-vcf merged.sites.bcf --region chr20:10000000-10200000</pre>

- The merged site list can be further annotated and consolidated as
follows:
> <pre>$ vt annotate_indels -r ~/data/ref/hs38DH.fa merged.sites.bcf \
-o + 2> /dev/null | vt consolidate_variants + -o union.sites.bcf </pre>

- The union sites must be indexed.
> <pre>$ bcftools index union.sites.bcf</pre>

---

#### c. Joint genoptyping of the merged variant sites across all samples.

- Next, we need to create an index file that contains a sample ID and
  the path to the corresponding CRAM file for each sample at each
  line.
> <pre>$ ls NWD*.bcf | cut -d . -f 1 | xargs -I {} \
echo "{} /BiO/data/topmed/{}.recab.cram" > cramlist.txt 
$ cat cramlist.txt
NWD230091 /BiO/data/topmed/NWD230091.recab.cram
NWD231092 /BiO/data/topmed/NWD231092.recab.cram
NWD259170 /BiO/data/topmed/NWD259170.recab.cram
NWD315403 /BiO/data/topmed/NWD315403.recab.cram
NWD495157 /BiO/data/topmed/NWD495157.recab.cram
NWD684137 /BiO/data/topmed/NWD684137.recab.cram</pre>
_NOTE:_ For users 26-50, use `BiO2` instead of `BiO`. 

- The joint genotyping can be performed by reaching each CRAM file
  sequentially.
> <pre>$ cramore dense-genotype --in-cram-list cramlist.txt \
--in-vcf union.sites.bcf --region chr20:10000000-10200000 \
--out PURtrio.joint.bcf </pre>

- Take a look at the output VCF file.
> <pre>$ bcftools view -H PURtrio.joint.bcf  | less -S </pre>
Use left and right arrows to navigate the columns across samples.
  * Q. Can you understand what each column means? 
  * Q. How about INFO fields? 
  * Q. How about FORMAT fields?

- To extract GT fields only across the 6 samples, we will use some
  hacky solution...
> <pre>$ bcftools view -H PURtrio.joint.bcf | grep -v ^# | cut -f 10- | perl \
  -lane 'print join("\t",substr($F[0],0,3),substr($F[1],0,3),substr($F[2],0,3),substr($F[3],0,3),substr($F[4],0,3),substr($F[5],0,3),)' \
  | sort | uniq -c | sort -n </pre>
  * Q. Considering these 6 samples are actually duplicate of a single
    trio, can you figure out which pairs of samples are identical? 
  * Q. Can you also figure out which one is offspring and which ones
    are parents?
  * Q. Does the genotypes mostly look Mendelian concordant or not?
  
---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, go back to [Day 1 Overview](../day1).

---
