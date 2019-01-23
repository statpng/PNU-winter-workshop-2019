---
layout: page
---

# Introduction to statistical genetics lecture in PNU winter workshop 2019

<hr>
<br>

**Overview:** This hands-on sessions introduce statistical methods and computational tools for analyzing high-throughput sequence data in population scale.

**Description:**

<!-- 최근 수년간의 급속한 기술발전으로 인해서 수많은 유전체 자료가 생산되고 있다.
뿐만 아니라, 작년 M&A of Illumina and PacBio, PromethION technology of ONT, and MGISEQT7 kit of BGI 등과 같은 기타 요인들은 전장유전체 해독 비용을 더욱 저렴하게 될 전망이다.
본 강의에서는 이렇게 생산된 수많은 유전체 자료를 분석하기 위한 방법들에 대해 소개하고자 한다.

첫 번째로, Genome-wide association studies(GWAS) 분야에서 가장 대표적인 univariate analysis를 수행하기 위한, 농생물 분야에서 가장 많이 사용되고 있는 GAPIT tool이다.
두 번째로, statistics와 computer science 분야에서 많이 사용되는 univariate analysis의 단점을 보완하여 더욱 좋은 성능을 갖고 있는 regularization procedure이다.

유전체 자료를 분석하기 위해서는 전처리에서부터 causal variants를 발굴하고 variant functional
유전체 자료의 예시로서 imputed wild bean dataset을 사용한다. -->

The dramatic advance of high-throughput in the last decade has produced tremendous amount genomic, epigenomic, and transcriptomic sequence data at an unprecedented scale.
Several factors such as M&A of Illumina and PacBio in last year, PromethION technology of ONT, and MGISEQT7 kit of BGI will make the price of whole genome sequencing much cheaper.
In this lecture, we introduce students to two methods for analysis of genomic sequencing data, and learn how to perform them.
The first is *GAPIT*, a tool for univariate analysis which is the most widely used for GWAS in bioinformatics.
The second is *regularization procedure*, which exceeds univariate analysis.
Using two methods, *univariate analysis* and *regulrization procedure*, we will analyze the imputed wild bean dataset and learn how to interpret the results.


<!-- **Audience:** This course is intended for researchers who are actively engaged in
genomics research and interested beginners, including laboratory scientists and clinicians with a
basic quantitative background. Ideally, participants are expected to have some basic knowledge of
human genetics (such as Mendelian inheritance), core statistical principles (such as p-values), and
basic UNIX skills (such as free contents material at
[https://www.codecademy.com/learn/learn-the-command-line](https://www.codecademy.com/learn/learn-the-command-line)) -->

<!-- **Requirements:** Participants must bring a laptop with specific [software installed]({{ site.baseurl }}/setup/). -->

**When:** January 28-30 (2019), 16:00 PM - 17:00 PM

**Where:** 313-103, Pusan National University, Pusan, Korea.

<br>
<hr>
<br>

#### [Part 0. Setup]({{ site.baseurl }}/part1/)

| Part | Time           | Topics                   |
| :-----: |   :--------------:    | :-----------------------|
| 0       | Introduction          | **1.Objectives**          |
|         |                       | **2.Methods**             |
| 0       | Setup                 | **1.Importing data**      |
|         |                       | **2.Data pre-processing** |
|         |                       | **3.Quality-control(QC)** |
|         |                       | **4.Imputation**          |



<br>
<hr>
<br>


#### [Part 1. Genome Association and Prediction Integrated Tool (GAPIT).]({{ site.baseurl }}/part1/)

| Part    |  Time                 | Topics                                                      |
| :-----: |   :--------------:    | :-----------------------                                    |
| I       | GAPIT                 | **1.Introduction of Genome-wide Association Study(GWAS)**   |
|         |                       | **2.Statistical Model of GAPIT**                            |
|         |                       | **3.Quality Control(QC) before Analysis**                   |
|         |                       | **4.Analysis using GAPIT & result**                         |
|         |                       | &nbsp; &nbsp; - code 1: treating heterozygosity             |
|         |                       | &nbsp; &nbsp; - code 2: imputation and controlling MAF      |
|         |                       | &nbsp; &nbsp; - code 3: optimal PC number                   |
|         |                       | &nbsp; &nbsp; - code 4 : Compressed MLM                     |
|         |                       | &nbsp; &nbsp; - code 5 : Enriched CMLM                      |

<br>
<hr>
<br>


#### [Part 2. Regularization procedures for variable selection.]({{ site.baseurl }}/part2/)

| Part    | Time                  | Topics                                                        |
| :-----: |   :--------------:    | :-----------------------                                      |
| II      | Regularization        | **Variable selection in high-dimensional data**               |
|         |                       | **Regularization procedures**                                 |
|         |                       | **"glmnet"**                                                  |
|         |                       | **"The Lasso"**                                               |
|         |                       | &nbsp; &nbsp; - Solution path                                 |
|         |                       | &nbsp; &nbsp; - Cross-validation                              |
|         |                       | &nbsp; &nbsp; - Variable selection                            |
|         |                       | &nbsp; &nbsp; - Comparison with univariate analysis           |
|         |                       | &nbsp; &nbsp; - Prediction                                    |
|         |                       | **Elastic-net**                                               |
|         |                       | &nbsp; &nbsp; - Cross-validation for two tuning parameters    |
|         |                       | **Covariate-adjusted model**                                  |


<br>
<hr>
<br>


#### [Part 3. Selection probabilities.]({{ site.baseurl }}/part3/)

| Part    | Time                   | Topics                                                     |
| :-----: |   :--------------:     | :-----------------------                                   |
| III     | Selection probabilties | **An algorithm of selection probabilities**                |
|         |                        | &nbsp; &nbsp; - Setting a grid of tuning parameters        |
|         |                        | &nbsp; &nbsp; - Applying the regularization for subsamples |
|         |                        | **Why split data with 0.5 proportion**                     |
|         |                        | &nbsp; &nbsp; - Advantages of "subagging"                  |
|         |                        | **Algorithm summary**                                      |
|         |                        | **The stability path**                                     |
|         |                        | **Threshold to control the false positive**                |
|         |                        | **Manhattan plot with selection probabilities**            |

<br>
<hr>
<br>


<!--## Relevant Resources
**BIOINF-575**: Programing Lab in Bioinformatics  

**Software Carpentry**: Occasional Workshops  
(Non planned for this year at UM unfortunately)
<http://software-carpentry.org> -->

<!--- Uncomment at end of course...
Add more courses when we find them.
-->
