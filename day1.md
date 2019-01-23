---
layout: page
title: Part 1 
permalink: /part1/
---

# Day 1. Analysis of High-throughput Sequence Reads

<hr>
<br>

Sequencing human genomes, transcriptomes, and epigenomes are becoming
more and more feasible and afforable solution to understand the
molecular basis of complex traits related to human
diseases.

In particular, understanding the effect of human genetic
variation within a large population is important for many applications
such as genome-wide association analysis (GWAS), expression
quantitative trait loci (eQTL) mapping, population genetic analysis,
and other related studies.

However, processing a large volume of sequence reads and jointly
analyzing them to detect and genotype variants requires sophisticated
technical skill sets.

Today, we will learn from the basics of DNA sequence reads to the
state-of-the-art software tools and pipelines to align sequence reads,
to detect genetic variants, and to genotype the variants across a
population of sequence genomes.

We will assume that you have basic UNIX skills, which can be obtained
by taking a free online course
[HERE](https://www.codecademy.com/learn/learn-the-command-line), but
we we will guide from the beginning, just in case you may get lost
from the beginning.

<br>
<hr>
<br>

### Schedule:

| Session | Time           | Topics                   |
| :-----: |:--------------:| :----------------------- |
| I       | 9:20-10:50 AM  | **Design and analysis of sequencing studies in population scale** |
|         | 10:50-11:00 AM | Coffee Break             |
| II      | 11:00-12:30 PM | **Alignment and quality control of sequenced genomes**       |
|         | 12:30-2:00 PM  | Lunch                    |
| III     | 2:00-3:30 PM   | **Calling short variants (SNPs and Indels) from sequence reads**    |
|         | 3:30-3:40 PM   | Coffee Break             |
| IV      | 3:40-5:10 PM   | **Calling structural variants (large deletions and CNVs) from sequence reads**   |

<br>

### Instructors:
Hyun Min Kang (HMK)
Goo Jun (GJ)


<hr>
<br>

## Topics:

### I) Design and analysis of sequencing studies in population scale

- This session consists of lectures only, and there will be no hands-on practice after the lecture.

—- Coffee Break [10 mins] —

### II) Alignment and quality control of sequenced genomes

This hands-on session will begin right after the completion of the
lecture part. It is separated by the following three parts.

- [II-1. Setting up the access to the compute
  server](../class-material/day1-connect-server)
- [II-2. Processing raw sequence
  reads](../class-material/day1-fastq-practice)
- [II-3. Processing aligned sequence reads
  ](../class-material/day1-bam-practice)
<!-- - [II-4. Quality control of aligned sequence
  reads](../class-material/day1-bam-quality-control) -->


—- Lunch Break [1:30 hr] —

### III) Calling short variants (SNPs and Indels) from sequence reads

- [III-1. Variant calling and joint genotyping using
  GATK HaplotypeCaller](../class-material/day1-gatk)
- [III-2. Variant calling and joint genotyping using vt+cramore](../class-material/day1-vt-cramore)  

—- Coffee Break [10 mins] —

#### IV) Calling structural variants (large deletions and CNVs) from sequence reads

- [IV-1. Preparing input files for Lumpy](../class-material/day1-prepare-svcall)
- [IV-2. SV discovery using Lumpy](../class-material/day1-lumpy)
- [IV-3. Genotyping SV using SVTyper ](../class-material/day1-svtyper)


—- End/Wrap-Up —

<br>

---
<br>

[Unix Referene Commands and Glossary](../class-material/unix-reference.html)  

[Introduction to Bash Shell Scripting](https://en.wikibooks.org/wiki/Bash_Shell_Scripting)  

[Make Install tutorial](http://www.ee.surrey.ac.uk/Teaching/Unix/unix7.html)
