# The promise and challenges of characterising genome-wide structural variants: A case study in a critically endangered parrot

This is the GitHub Repository supporting our submission to the Molecular Ecology Resources Special Issue *Omics Tools for Species Conservation*.  
DOI: [https://doi.org/10.1111/1755-0998.13783]( https://doi.org/10.1111/1755-0998.13783)

Here are the steps we took to trial different strategies for SV discovery using Illumina and Oxford Nanopore Technology (ONT) sequence data. The raw Illumina reads were preprocessed by [Joseph Guhlin](https://github.com/jguhlin) under the [Kākāpō125+ Sequencing Consortium](https://www.doc.govt.nz/our-work/kakapo-recovery/what-we-do/research-for-the-future/kakapo125-gene-sequencing/). Sequencing coverage for the Illumina data ranged from 9x - 38x with a mean of ~18x, and ranged from 4.2x - 15.8x with a mean of 9.7x coverage for the seven indivudals included for long-read based SV discovery.  

For this manuscript, three short-read discovery tools (Delly, Manta and SMOOVE) and two long-read SV discovery tools (CuteSV and Sniffles) were used. We also implemented the short-read genotyping tool BayesTyper to genotype outputs from CuteSV, Manta and Sniffles.  Below is a short summary of each of the steps taken in this study.  

### 1) Illumina short-read alignment with BWA
All short read data were accessed through the kākāpō 125+ sequencing consortium. These processed reads were aligned to the kākāpō reference genome ([VGP](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_004027225.2/)) using [BWA](http://bio-bwa.sourceforge.net/). The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. Other programs used in this step include: [SAMtools](https://github.com/samtools/samtools) v1.16.  

### 2) Illumina short-read alignment statistics
Programs used [Mosdepth](https://github.com/brentp/mosdepth) v0.3.3, [qualimap](http://qualimap.conesalab.org/) v2.2.2 and [SAMtools](https://github.com/samtools/samtools) v1.16. Summaries of these outputs were visualised using [MultiQC](https://github.com/ewels/MultiQC) v1.13.

### 3) ONT basecalling and read QC
Reads were basecalled using Guppy v6.3.7 using `dna_r9.4.1_450bps_sup` basecalling model. Adapters were trimmed from raw sequence reads using [Porechop](https://github.com/rrwick/Porechop) v0.2.4, lambda DNA removed using [NanoLyse](https://github.com/wdecoster/nanolyse) v1.2.0 and filtered for a minimium Q-score of 10 and length of 5kb using [NanoFilt](https://github.com/wdecoster/nanofilt) v2.8.0. Raw and filtered read quality were assessed using [NanoPlot](https://github.com/wdecoster/NanoPlot) v1.39.0.

### 4) ONT alignment and statistics
Reads were mapped using [Winnowmap](https://github.com/marbl/Winnowmap) v2.03, mapping quality was assessed using [qualimap](http://qualimap.conesalab.org/) v2.2.2, [Mosdepth](https://github.com/brentp/mosdepth) v0.3.3 and summaries visualised using [MultiQC](https://github.com/ewels/MultiQC) v1.13.  

### 5) Structural Variant discovery using the short-read discovery tool [Delly](https://github.com/dellytools/delly) v0.8.7
Called and genotyped SVs with Delly. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

### 6) Structural Variant discovery using the short-read discovery tool [Smoove](https://github.com/brentp/smoove) v0.2.8
Called and genotyped SVs with the Smoove pipeline. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

### 7) Structural Variant discovery using the short-read discovery tool [Manta](https://github.com/Illumina/manta)
Called SVs using two strategies 1) Joint and 2) Batched. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

### 8) Genotyping Manta outputs with [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper) v1.5
Genotyped SVs using BayesTyper for both Joint and Batched Manta outputs. Other programs used included [BCFtoools](http://samtools.github.io/bcftools/) and [KMC](https://github.com/refresh-bio/KMC)

### 9) Structural variant discovery using the long-read discovery tools [CuteSV](https://github.com/tjiangHIT/cuteSV) v1.0.11 and [Sniffles](https://github.com/fritzsedlazeck/Sniffles) v2.0.7
After initial SV discovery with either CuteSV or Sniffles calls were refined using [Jasmine](https://github.com/mkirsche/Jasmine)v1.1.5.

### 10) Genotyping CuteSV and Sniffles outputs with [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper) v1.5
SVs called using the long-read data were genotyped using BayesTyper.  

### 11) Identified consensus SV calls
Identified overlapping SVs with [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) v1.0.7. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

### 12) Summarising SV characteristics
Here is an R markdown outlining how SV size distributions, SV type frequency, relative frequency per individual taking into account generation and lineage were visualised. This file also outlines how DAPCs for each of the 6 methods were constructed.  
