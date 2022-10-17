# The promise and challenges of structural variant discovery: A conservation case study in the critically endangered kākāpō (*Strigops habroptilus*)

Here are the steps I took to trial different strategies for SV discovery using Illumina and Oxford Nanopore Technology sequence data. The raw Illumina reads were preprocessed as under the kākāpō125+ sequencing consortium. Sequencing coverage for the Illumina data ranged from 9x - 38x with a mean of ~18x.

Three short-read discovery pipelines were used (Delly, Manta and SMOOVE), while two were used for long-read SV discovery (CuteSV and Sniffles).

## 1 Illumina short-read alignment with BWA
All short read data were accessed through the kākāpō 125+ sequencing consortium. These processed reads were aligned to the kākāpō reference genome ([VGP](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_004027225.2/)) using [BWA](http://bio-bwa.sourceforge.net/). The W sex chromosome (present in females) was removed prior to alignments of males as the W is highly repetitive, and may contain regions that would result in misalignment of reads. Other programs used in this step include: [SAMtools](https://github.com/samtools/samtools)

## 2 Illumina short-read alignment statistics
Programs used [Mosdepth](https://github.com/brentp/mosdepth), [qualimap](http://qualimap.conesalab.org/) and [SAMtools](https://github.com/samtools/samtools).

## 3 basecalling and read QC for Oxford Nanopore data

## 4 ONT alignment and statistics
Used winnowmap

## 5 Structural Variant discovery using the short-read discovery tool [Delly](https://github.com/dellytools/delly)
Called and genotyped SVs with Delly. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 6 Structural Variant discovery using the short-read discovery tool [Smoove](https://github.com/brentp/smoove)
Called and genotyped SVs with the Smoove pipeline. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 7 Structural Variant discovery using the short-read discovery tool [Manta](https://github.com/Illumina/manta)
Called SVs using two strategies 1) Joint and 2) Batched. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 8 Structural variant discovery using the long-read discovery tools [CuteSV]() and [Sniffles]()

## 9 Genotyping Manta, CuteSV and Sniffles outputs with [BayesTyper](https://github.com/bioinformatics-centre/BayesTyper)
Genotyped SVs using BayesTyper for both Joint and Batched Manta outputs. Other programs used included [BCFtoools](http://samtools.github.io/bcftools/) and [KMC](https://github.com/refresh-bio/KMC)

## 10 Identified consensus SV calls with [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR)
Identified overlapping SVs with SURVIVOR. [BCFtoools](http://samtools.github.io/bcftools/) was also used in this step.

## 11 Exploring SV data in R
Plotted SV size distributions, SV type frequency, relative frequency per individual taking into account generation and lineage and assessed population structure with principal component analyses.