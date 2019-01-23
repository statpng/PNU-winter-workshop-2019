---
layout: page
---

---

### II-3. Processing aligned sequence reads

In this part, we will briefly learn how to understand the format of
SAM/BAM files, and how to manipulate align sequence reads to
generate 'analysis ready' SAM/BAM file to perform variant calling.

---

#### a. Sort aligned sequence reads by genomic coordinates

- The aligned sequence reads by `bwa mem` must be sorted by
genomic coordinates prior to variant detection. It can be achieved
using `samtools sort` command.
> <pre>
samtools sort -o HG00406.sorted.bam HG00406.bwa.sam </pre>
  * Q. Try `less HG00406.sorted.bam`. Can it be read as plain text file?
    Is the format of the output file SAM or BAM?
  * Q. Try `samtools view HG00406.sorted.bam | less -S`. Are the aligned
    sequence reads sorted by genomic coordinates? 
	
- The following command can be useful to count the number of reads
aligned to each chromosome. 
> <pre> samtools view HG00406.sorted.bam | cut -f 3 | uniq -c </pre>
  * Q. What does the command do? Can you interpret the output?
  * Q. Which chromosome contains the most of aligned reads?
  
- Aligning the raw sequence reads using `bwa mem` and coordinate-wise
sorting using `samtools sort` can be done jointly in one step,
without having to store the SAM file.
> <pre>
bwa mem -t 2 -R '@RG\tID:HG00406_1KG_lowcov\tSM:HG00406\tLB:HG00406\tPL:ILLUMINA' \
~/data/ref/hs38DH.fa ~/data/1000g/fastqs/HG00406.R1.fastq.gz \
~/data/1000g/fastqs/HG00406.R2.fastq.gz | samtools sort -o HG00406.sorted.bam - </pre>
  * Q. What would be the benefit of jointly running the command this way?
  * Q. How was UNIX pipe used to perform the task jointly?
  * Q. What was the difference in the arguments compared to when the
    two steps were ran separately?
  
---

#### b. Preparing analysis-ready CRAM files.
	
- Next, we will perform the following steps to convert the aligned and
sorted sequence reads to be prepared for variant calling.
  1. Mark duplicated reads.
  2. Perform base quality score recalibration (BQSR).
  3. Collapse base quality scores into the following values: 2 (0-2),
     3 (3), 4 (4), 5 (5), 6 (6), 10 (7-12), 20 (13-22), 30 (23+)
  4. Store the output sequence reads into CRAM format for a better
     storage efficiency.

- These steps can be performed jointly using a single command using
  `bamUtil` and `samtools` as follows:
> <pre>
bam dedup_lowmem --allReadNames --binCustom --binQualS \
"0:2,3:3,4:4,5:5,6:6,7:10,13:20,23:30" --in HG00406.sorted.bam \
--recab --out -.ubam --refFile ~/data/ref/hs38DH.fa \
--dbsnp ~/data/ref/dbsnp_142.b38.vcf.gz | samtools view -h -C \
-T ~/data/ref/hs38DH.fa -o HG00406.recal.cram - </pre>
  * Q. What are the sizes of SAM file (right after `bwa mem`), and the
    sorted BAM file, and the final CRAM file? 
  * Q. How many folds of
    storage efficiencies are achieved when (compressed) BAM are used
    as opposed to uncompressed BAM, and when CRAM is used (with base
    quality binning) as opposed to compressed BAM?
  * Q. What would be the size of file if BAM was used instead of CRAM at
    the final step, assuming quality binning is still applied. How can
    this be tested?

- If `Picard` and `GATK` were used, these steps should be performed
  separately as follows:
  1. Run [`Picard
     MarkDuplicates`](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
     to mark duplicated reads from the sorted BAM file.
  2. Run [`GATK
     BaseRecalibrator`](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
     to create BQSR table.
  3. Run [`GATK
     PrintReads`](https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)
     to use the BQSR table to perform BQSR, while quantizing the
     quality scores into specific bins, using the following options:
     `--globalQScorePrior -1.0 --preserve_qscores_less_than 6
     --static_quantized_quals 10 --static_quantized_quals 20
     --static_quantized_quals 30 --disable_indel_quals`
   4. Run `samtools view -C -T ~/data/ref/hs38DH.fa -o
   [output_cram_file_name] [input_bam_file_name]` to convert the
   output BAM file into a CRAM file.
   
---

#### c. Indexing and viewing the analysis-ready CRAM files

- The analysis-ready CRAM file can be indexed using `samtools index`
  to enable random access at arbitrary genomic coordinate.
> <pre>
samtools index HG00406.recal.cram </pre>

- Viewing the CRAM formatted file requires the knowledge of reference
  genome. Otherwise, the access will be very slow because it is trying
  to access internet to query some information regarding the reference
  genome.
  
- One way to view the CRAM formated file is to include `-T [ref.fa]`
  when running `samtools view` as follows:
><pre>
samtools view -T ~/data/ref/hs38DH.fa HG00406.recal.cram chr20:0 | less -S </pre>
  
- Alternatively, one can first set an enviroment variable specifying
  local repository of MD5 sum of each contig, using the following command:
><pre>
export REF_PATH=~/data/ref/md5/%2s/%2s/%s </pre>
Then reading the CRAM file does not require `-T [ref.fa]` option
  without slowing things down. 
><pre>
samtools view HG00406.recal.cram chr20:0 | less -S </pre>
  As this approach works better with other
  software tools beyong `samtools`, setting `REF_PATH` environment
  variable is a preferred approach. To ensure `REF_PATH` enviroment
  variable is set, try to run:
><pre> env | grep REF_PATH </pre>
It should return the absolute path for `REF_PATH`, for example (varies
  by user), 
><pre> REF_PATH=/BiO/home/kobic00/data/ref/md5/%2s/%2s/%s </pre>


---

#### d. Using `samtools tview` for visualizing sequence reads

- Once the CRAM file is indexed, you can visualize the sequence reads
using `samtools tview` or `IGV` tool. Here we explain how to use
`samtools tview`.

- To invoke `samtools tview`, you need to specify the aligned
SAM/BAM/CRAM file and the reference FASTQ file as follows:
> <pre> samtools tview HG00406.recal.cram ~/data/ref/hs38DH.fa </pre>

- Then it will show the beginning of chromosome 1, which does not have
  any aligned read. To move to chromosome 20, try the following steps:
  1. Press `g` key, then you will see a prompt `Goto:` in the screen
  2. Type `chr20:100000 [Enter]`, to move chromosome 20 at 100,000bp
     position.
  3. Pressing left and right arrow will allow moving nearby genomic
     regions.
  4. Pressing `?` will show the help page. 
  5. For example, pressing `.`
     will turn on/off representing reference sequence as dot or
     comma. Pressing 'b' will color by base quality, and pressing 'n'
     will color by nucleotide. 
  6. Press `q` to end `samtools tview` and return to the shell.
  
- To compare how invoke `samtools tview`, you need to specify the aligned
SAM/BAM/CRAM file and the reference FASTQ file as follows:
> <pre> samtools tview HG00406.recal.cram ~/data/ref/hs38DH.fa </pre>
You may repeat the steps 1-6 from above and tell what the major
difference is.
  
---

Feel free to ask questions to your instructor(s) if you are stuck. 
Otherwise, let's 
<!-- go to next step 
[II-4 Quality control of aligned sequence reads](../class-material/day1-bam-quality-control.html)
, or -->
go back to [Day 1 Overview](../day1).

---
