---
layout: page
---

---

### IV-2. SV Discovery using LUMPY

In this part, we will briefly learn how to use LUMPY Express for multi-sample SV discovery.

---

#### Running LUMPY

LUMPY has two distinct execution alternatives. LUMPY Express is a simplified wrapper for standard analyses.
LUMPY (traditional) is more customizable, for advanced users and specialized experiments.

We will use LUMPY Express, because our use is a standard use case. 

##### LUMPY Express standard use cases
- Run LUMPY Express on a single sample with pre-extracted splitters and discordants

> <pre>
    lumpyexpress \
        -B sample.bam \
        -S sample.splitters.bam \
        -D sample.discordants.bam \
        -o sample.vcf </pre>

- Run LUMPY Express jointly on multiple samples with pre-extracted splitters and discordants

> <pre>
    lumpyexpress \
        -B sample1.bam,sample2.bam,sample3.bam \
        -S sample1.splitters.bam,sample2.splitters.bam,sample3.splitters.bam \
        -D sample1.discordants.bam,sample2.discordants.bam,sample3.discordants.bam \
        -o multi_sample.vcf </pre>

- Run LUMPY Express on a tumor-normal pair
<pre>
    lumpyexpress \
        -B tumor.bam,normal.bam \
        -S tumor.splitters.bam,normal.splitters.bam \
        -D tumor.discordants.bam,normal.discordants.bam \
        -o tumor_normal.vcf </pre>

#### Running LUMPY on our samples

> <pre>
# From your sv project directory
lumpyexpress \
 -B bams/NWD230091.bam,bams/NWD231092.bam,bams/NWD259170.bam,bams/NWD315403.bam,bams/NWD495157.bam,bams/NWD684137.bam \
 -D NWD230091.discordant.bam,NWD231092.discordant.bam,NWD259170.discordant.bam,NWD315403.discordant.bam,NWD495157.discordant.bam,NWD684137.discordant.bam \
 -S NWD230091.splitter.bam,NWD231092.splitter.bam,NWD259170.splitter.bam,NWD315403.splitter.bam,NWD495157.splitter.bam,NWD684137.splitter.bam \
 -o multi.sv.vcf </pre>

- Tip: easy way to generate comma-separated filenames
> <pre>
ls bams/*.bam  |tr '\n' ',' </pre>

* Open multi.sv.vcf file. How many variants did you discover?
* Check INFO column for additional information about each variant 
* What do you recognize in the genotype field?

---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's go to next step 
[IV-3. Genotyping SV using SVTyper](../class-material/day1-svtyper)
, or go back to [Day 1 Overview](../day1).

---
