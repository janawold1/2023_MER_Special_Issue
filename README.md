# The promise and challenges of structural variant discovery: A conservation case study in the critically endangered kākāpō (*Strigops habroptilus*)

Here are the steps we took to trial different strategies for SV discovery using Illumina and Oxford Nanopore Technology (ONT_ sequence data. The raw Illumina reads were preprocessed by [Joseph Guhlin](https://github.com/jguhlin) under the [Kākāpō125+ Sequencing Consortium](https://www.doc.govt.nz/our-work/kakapo-recovery/what-we-do/research-for-the-future/kakapo125-gene-sequencing/). Sequencing coverage for the Illumina data ranged from 9x - 38x with a mean of ~18x.

For this manuscript, three short-read discovery tools (Delly, Manta and SMOOVE) and two long-read SV discovery tools (CuteSV and Sniffles) were used.

## 1 Illumina short-read alignment with BWA
All short read data were accessed through the kākāpō 125+ sequencing consortium. These processed reads were aligned to the kākāpō reference genome ([VGP](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_004027225.2/)) using [BWA](http://bio-bwa.sourceforge.net/). The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. Other programs used in this step include: [SAMtools](https://github.com/samtools/samtools)

## 2 Illumina short-read alignment statistics
Programs used [Mosdepth](https://github.com/brentp/mosdepth), [qualimap](http://qualimap.conesalab.org/) and [SAMtools](https://github.com/samtools/samtools). Summaries of these outputs were visualised using [MultiQC](https://github.com/ewels/MultiQC).

## 3 basecalling and read QC for ONT data
Reads were basecalled using [Guppy]() under the `dna_r9.4.1_450bps_sup` basecalling model. Adapters were trimmed from raw sequence reads using [Porechop](https://github.com/rrwick/Porechop)v0.2.4, lambda DNA removed using [NanoLyse](https://github.com/wdecoster/nanolyse)v1.2.0 and filtered for a minimium Q-score of 10 and length of 5kb using [NanoFilt](https://github.com/wdecoster/nanofilt)v2.8.0. Raw and filtered read quality were assessed using [NanoPlot](https://github.com/wdecoster/NanoPlot).

## 4 ONT alignment and statistics
Reads were mapped using [Winnowmap](https://github.com/marbl/Winnowmap)v2.03, mapping quality was assessed using [qualimap](http://qualimap.conesalab.org/), [Mosdepth](https://github.com/brentp/mosdepth) and summaries visualised using [MultiQC](https://github.com/ewels/MultiQC).  

## 5 Structural Variant discovery using the short-read discovery tool [Delly](https://github.com/dellytools/delly)
Called and genotyped SVs with Delly. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 6 Structural Variant discovery using the short-read discovery tool [Smoove](https://github.com/brentp/smoove)
Called and genotyped SVs with the Smoove pipeline. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 7 Structural Variant discovery using the short-read discovery tool [Manta](https://github.com/Illumina/manta)
Called SVs using two strategies 1) Joint and 2) Batched. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 8 Structural variant discovery using the long-read discovery tools [CuteSV](https://github.com/tjiangHIT/cuteSV) and [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
After initial SV discovery, calls were refined using [Jasmine](https://github.com/mkirsche/Jasmine).

## 9 Genotyping with [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper)
Genotyped SVs using BayesTyper for both Joint and Batched Manta outputs. Other programs used included [BCFtoools](http://samtools.github.io/bcftools/) and [KMC](https://github.com/refresh-bio/KMC)

## 10 Identified consensus SV calls
Identified overlapping SVs with [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR). [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 11 Summarising SV characteristics
Plotted SV size distributions, SV type frequency, relative frequency per individual taking into account generation and lineage and assessed population structure with principal component analyses.