---
layout: page
---

---

### IV-1.  Preparing input files for LUMPY

In this part, we will briefly learn how to prepare for SV calling using 
LUMPY, a breakpoint detection tool for structural variations.

---

#### a. BAM files 

LUMPY uses information of split reads and discordant read pairs, stored as separate BAM files
Recent versions of LUMPY also can handle CRAM files, but it still requires split reads and 
discordant read pairs in BAM format, so we will use BAM files from TOPMed whole genome sequenced
samples. The BAM versions of a 1-Mbp region (chr22:36,000,000-37,000,000) are generated at

> <pre>
# Server 1
/BiO/data/topmed/bams
# Server 2
/BiO2/data/topmed/bams </pre>

- These BAMs are generated from original CRAM files using samtools:

> <pre>
samtools view -h -b NWD230091.recab.cram chr22:36000000-37000000  > NWD230091.bam
samtools index NWD230091.bam</pre>

- There are six (6) BAM files in the directory. 
 
Let's create a directory for your SV calling project and create a symbolic link 
to the BAM files directory under the project directory.

> <pre>
# From your home directory
mkdir sv
cd sv
# If you are using Server 1 (i.e. your home directory is under /BiO/home)
ln -s /BiO/data/topmed/bams bams 
# If you are using Server 2 (i.e. your home directory is under /BiO2/home)
ln -s /BiO2/data/topmed/bams bams </pre>


#### b. Split reads and discordant read pairs

LUMPY (Express) expects BWA-MEM aligned BAM files as input.
It automatically parses sample, library, and read group information using the @RG
tags in the BAM header.
Each BAM file is expected to contain exactly one sample.

The minimum input is a coordinate-sorted BAM file (-B), from which LUMPY Express
extracts splitters and discordants using SAMBLASTER before running LUMPY.
Optionally, users may supply coordinate-sorted splitter (-S) and discordant (-D)
BAM files which will bypass SAMBLASTER extraction for faster analysis.

If you have not yet aligned your data, you can generate these additional data 
with [SpeedSeq](https://github.com/cc2qe/speedseq), which
performs BWA-MEM alignment, marks duplicates and extracts split and discordant
read-pairs. For example, 

> <pre>
speedseq align -R "@RG\tID:id\tSM:sample\tLB:lib" \
    human_g1k_v37.fasta \
    sample.1.fq \
    sample.2.fq </pre>

Otherwise, split-read and discordant read pair data may be generated from BAM/CRAM files 
already aligned with BWA-MEM.

> <pre>
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 sample.bam |samtools sort - > sample.discordants.bam </pre>

> <pre>
# Extract the split-read alignments
samtools view -h sample.bam \
    | /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    | samtools sort -  \
	> sample.splitters.bam </pre>

- Generate discordant read pair BAM files
> <pre>
# From your 'sv' directory
samtools view -b -F 1294 bams/NWD230091.bam |samtools sort - > NWD230091.discordant.bam
samtools view -b -F 1294 bams/NWD231092.bam |samtools sort - > NWD231092.discordant.bam
samtools view -b -F 1294 bams/NWD259170.bam |samtools sort - > NWD259170.discordant.bam
samtools view -b -F 1294 bams/NWD315403.bam |samtools sort - > NWD315403.discordant.bam
samtools view -b -F 1294 bams/NWD495157.bam |samtools sort - > NWD495157.discordant.bam
samtools view -b -F 1294 bams/NWD684137.bam |samtools sort - > NWD684137.discordant.bam </pre>

- Generate split read BAM files
> <pre>
samtools view -h bams/NWD230091.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD230091.splitter.bam
samtools view -h bams/NWD231092.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD231092.splitter.bam
samtools view -h bams/NWD259170.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD259170.splitter.bam
samtools view -h bams/NWD315403.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD315403.splitter.bam
samtools view -h bams/NWD495157.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD495157.splitter.bam
samtools view -h bams/NWD684137.bam| /BiO/apps/bin/scripts/extractSplitReads_BwaMem -i stdin|samtools view -Sb - |samtools sort - > NWD684137.splitter.bam</pre>

* Q. Open one of the discordant read pair BAM files. What do you recognize?
* Q. Open one of the split read BAM files. What do you recognize?
* Q. If you have hundreds of samples, how could you easily automate running these commands?

---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's go to next step 
[IV-2. SV discovery using LUMPY](../class-material/day1-lumpy)
, or go back to [Day 1 Overview](../day1).

---
