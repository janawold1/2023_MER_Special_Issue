# Population-scale genotyping with BayesTyper
BayesTyper is a Bayesian based genotyping program. However, one draw back is that it is computationally intensive and has multiple intermediate steps to navigate.

Below I outline of how I filtered the raw Manta SV calls for genotyping with BayesTyper, create kmer counts, cluster variant kmers and finally genotype variants.  

To begin, global variables were defined as below:
```
bayes=/kakapo-data/bayestyper/
chr_ref=/kakapo-data/References/kakapo_autosomes.fa
deepV=DEEPVARIANT:/kakapo-data/bwa/manta/bayestyper/Trained_SV_scaffolds_renamed.sorted.vcf
exclude=/kakapo-data/metadata/kakapo_SVexcluded_scaffolds.bed
female=/kakapo-data/bwa/bwa_female/
male=/kakapo-data/bwa/bwa_male/
raw_data=/kakapo-data/bwa/manta/raw_variants/
ref=/kakapo-data/References/kakapo_full_ref.fa
trio=/kakapo-data/metadata/sample_trios.csv
```

## Merging sex-specific variant calls and batched calls for SV genotyping
Candidate SVs were called for the batch and joint Manta data sets.
```
bcftools merge -m all -O z -o ${raw_data}joint_total.vcf.gz \
    ${raw_data}joint_calling/female_diploidSV.vcf.gz \
    ${raw_data}joint_calling/male_diploidSV.vcf.gz

bcftools merge -m all -O z -o ${raw_data}batch_total.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch01_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch02_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch03_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch04_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch05_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch06_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_female/batch07_manta_INV_conversion_female.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch01_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch02_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch03_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch04_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch05_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch06_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch07_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch08_manta_INV_conversion_male.vcf.gz \
    ${raw_data}batch_calling/batches_male/batch09_manta_INV_conversion_male.vcf.gz
```
## SV filtering for genotyping candidates
Here, we applied thresholds to filter for SV call quality. All SVs were required to pass Manta's internal 'hard' filters (`FILTER=="PASS"`), have at least one sample pass all filters for SV calling (`FORMAT/FT=="PASS"`), did not exceed the maximum allowed depth as estimated by Manta (`FILTER!="MaxDepth"`), fail to adhere to diploid expectations (`FILTER!="Ploidy"`), all SVs <1kb in length could not have the fraction of reads with a MAPQ score of 0 exceed 0.4(`FILTER!="MaxMQ0Frac"`), and those SVs significantly larger than the paired-read length with no paired support in any sample were excluded (`FILTER!="NoPairSupport"`). Finally, all SVs were set to have a minimum length of 50bp, and all breakends were excluded.  

```
    bcftools view \
        -i '(FILTER=="PASS" & PR >= 10 & FORMAT/FT == "PASS" & FILTER!="MaxDepth" & FILTER!="Ploidy" & FILTER!="MaxMQ0Frac" & FILTER!="NoPairSupport" & SVLEN<=-50 & SVLEN > -50000 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FORMAT/FT == "PASS" & FILTER!="MaxDepth" & FILTER!="Ploidy" & FILTER!="MaxMQ0Frac" & FILTER!="NoPairSupport" & SVLEN>=50 & SVTYPE!="BND" & SVLEN < 50000)' \
        -T ^${exclude} \
        -O z -o ${out}batch_filtered/01_batch_filtered.SVs.vcf.gz \
        ${raw_data}batch_total.vcf.gz
bcftools view \
    -i 'SVLEN<50000 & SVLEN >-50000' \
    -O z -o ${mantaB}02_batch_length_filtered.SVs.vcf.gz \
    ${mantaB}01_batch_filtered.SVs.vcf.gz

bcftools view \
    -i '(FILTER=="PASS" & PR >= 10 & FORMAT/FT == "PASS" & FILTER!="MaxDepth" & FILTER!="Ploidy" & FILTER!="MaxMQ0Frac" & FILTER!="NoPairSupport" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FORMAT/FT == "PASS" & FILTER!="MaxDepth" & FILTER!="Ploidy" & FILTER!="MaxMQ0Frac" & FILTER!="NoPairSupport" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}joint_filtered/01_joint_filtered.SVs.vcf.gz \
    ${raw_data}joint_total.vcf.gz
bcftools view \
    -i 'SVLEN<50000 & SVLEN >-50000' \
    -O z -o ${mantaJ}02_joint_length_filtered.SVs.vcf.gz \
    ${mantaJ}01_joint_filtered.SVs.vcf.gz
```
It is important to note that in the process of merging multiple batches, some called SV types may not match. Breakend calls in particular are problematic as they suggest a lack of evidence to call a SV type in a given batch. To identify and remove SV calls that passed inital filtering thresholds, but may have been called as a breakend in either the joint male or female calls or the batched samples we counted these sites as per:

```
bcftools view -H ${mantaB}02_batch_length_filtered.SVs.vcf | \
    awk '{print $1, $2, $3}' | \
    grep ";" | \
    tr ";" "\t" | \
    tr ":" "\t" | \
    awk '{print $1","$2","$3","$10","$11}' > ${mantaB}batch_mistatch_calls.csv

bcftools view -H ${mantaJ}02_joint_length_filtered.SVs.vcf | \
    awk '{print $1, $2, $3}' | \
    grep ";" | \
    tr ";" "\t" | \
    tr ":" "\t" | \
    awk '{print $1","$2","$3","$10","$11}' > ${mantaJ}joint_mistatch_calls.csv
```
These files were opened in excel and non-matching SVtypes were identified. This indicated that there were 8 of these sites in the batched data and 6 in the joint data. These sites were either INDELS or Insertions/Duplications. However, there were a few sites that were paired with breakened calls (n = 4 batch; n = 1 joint). Because breakends may be indicative of unresolved complex SVs, these sites were removed prior to genotyping.

```
bgzip ${mantaB}02_batch_length_filtered.SVs.vcf
tabix ${mantaB}02_batch_length_filtered.SVs.vcf.gz
bcftools view -t ^NC_044277.2:125934488,NC_044279.2:66579444,NC_044282.2:4557370 \
    -O v -o ${mantaB}03_batch_consolidated_sites.vcf \
    ${mantaB}02_batch_length_filtered.SVs.vcf.gz

bgzip ${mantaJ}02_joint_length_filtered.SVs.vcf
tabix ${mantaJ}02_joint_length_filtered.SVs.vcf.gz
bcftools view -t ^NC_044279.2:66579444 \
    -O v -o ${mantaJ}03_joint_consolidated_sites.vcf \
    ${mantaJ}02_joint_length_filtered.SVs.vcf.gz
```

Candidate SVs were normalised prior to conversion to BayesTyper format as per:
```
bcftools norm --threads 24 -f ${chr_ref} -m -any \
    -O v -o ${mantaB}04_batch_candidates_norm.vcf ${mantaB}03_batch_consolidated.vcf

bcftools norm --threads 24 -f ${chr_ref} -m -any \
    -O v -o ${mantaJ}04_joint_candidates_norm.vcf ${mantaJ}03_joint_consolidated.vcf
```
This left 35,832 SVs in the batch call set and 32,072 in the joint call set for allele conversion and merging for BayesTyper. The `04_{batch,joint}_candidates_norm.vcf` file was used to explore SV summary characteristics as outlined below.

```
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVTYPE\tmantaB_unfiltered\n' ${raw_data}batch_total.vcf.gz > ${out}manta_summary.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVTYPE\tmantaB_SVfiltered\n' ${mantaB}04_batch_candidates_norm.vcf >> ${out}manta_summary.tsv

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVTYPE\tmantaJ_unfiltered\n' ${raw_data}joint_total.vcf.gz > ${out}manta_summary.tsv
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\t%SVTYPE\tmantaJ_SVfiltered\n' ${mantaJ}04_joint_candidates_norm.vcf >> ${out}manta_summary.tsv
```

## Convert Manta VCF to remove symbolic alleles

This converts all symbolic alleles to a full sequence. The conversion results in large VCFs, so ensure that the output is gzipped iwth the -z flag. Note, bayesTyper does NOT support bgzip files. BayesTyper benefits from including SNP data as it helps with initialising the kmers. Therefore it's recommended to conduct SNP discovery with either GATK HaplotypeCaller or FreeBayes. I used the kākāpō125+ DeepVariant SNP call set prepared by [Guhlin et al. 2022](https://www.biorxiv.org/content/10.1101/2022.10.22.513130v1.abstract).

For combining SV calls with the DeepVariant SNPs, all chromosome sizes were used to match kākāpō125+ and NCBI chromosome IDs. bcftools annotate was used to rename chromosomes to be consistent with NCBI scaffold names. SNPs were filtered down to only include regions of interest (i.e., scaffolds included for SV analysis). But after converting chromosome names, needed to sort VCF header and body of DeepVariant file to match the manta VCFs. Sorted the body with [this solution](https://www.biostars.org/p/84747/), then manually edited the header.

```
bcftools norm -c we -f ${ref} -O v -o ${out}bayestyper/Trained_SV_scaffolds_renamed_norm.vcf \
    ${out}bayestyper/Trained_SV_scaffolds_renamed.vcf

bgzip ${out}Trained_SV_scaffolds_renamed_norm.vcf; tabix -pvcf ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz

cat chr_list.txt | xargs tabix -h ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz > ${out}Trained_SV_scaffolds_renamed.sorted.vcf
```

Allele conversion and file combination was then conducted as per:

```
bayesTyperTools convertAllele \
    --variant-file ${bayes}mantaB/04_batch_candidates_norm.vcf \
    --genome-file ${chr_ref} \
    --output-prefix ${bayes}05_batch_converted
bayesTyperTools combine -v ${deepV},MANTAB:${bayes}mantaB/05_batch_converted.vcf \
    -o ${bayes}mantaB/06_batch_combined -z


bayesTyperTools convertAllele \
    --variant-file ${bayes}mantaJ/04_joint_candidates_norm.vcf \
    --genome-file ${chr_ref} \
    --output-prefix ${bayes}mantaJ/05_joint_converted
bcftools sort -T ${bayes}mantaJ/ -O v -o ${bayes}mantaJ/06_joint_converted.sorted.vcf \
    ${bayes}mantaJ/05_joint_converted.vcf
bayesTyperTools combine -v ${deepV},MANTA:${bayes}mantaJ/06_joint_converted.sorted.vcf \
    -o ${bayes}mantaJ/07_joint_combined -z
```

## Running KMC and makeBloom

It is really easy to over-resource KMC, that is to say if you provide it with too much memory/threads the program freezes. I have had the best luck with 24Gb RAM and 48 threads as per below. It also found that I needed to run KMC one individual at a time. I think threads clashed when I attempted to run multiple individuals at once since KMC uses the same temp file naming convention each time it runs.

To aid in genotyping, KMC was run using BAMs filtered for autosomal scaffolds only using the solution provided [here](https://www.biostars.org/p/302771/)

```
awk '{printf("%s\t0\t%s\n",$1,$2);}' ${ref}.fai | grep NC_ | grep -v NC_044301 | grep -v NC_044302 > /kakapo-data/reference/kakapo_autosomes.bed

for mbam in ${male}nodup_bam/*_nodup.bam
    do
    base=$(basename $mbam .bam)
    echo "Extracting autosomes for $base..."
    samtools view -@ 64 -L /kakapo-data/references/kakapo-autosomes.bed -o ${male}nodup_autosomes/${base}_autosomes.bam
    samtools index -@ 64 ${male}nodup_autosomes/${base}_autosomes.bam
    samtools stats -@ 64 ${male}nodup_autosomes/${base}_autosomes.bam > ${male}nodup_autosome_stats/${base}_autosomes.stats
done

for fbam in ${female}nodup_bam/*_nodup.bam
    do
    base=$(basename $fbam .bam)
    echo "Extracting autosomes for $base..."
    samtools view -@ 64 -L /kakapo-data/references/kakapo-autosomes.bed -o ${female}nodup_autosomes/${base}_autosomes.bam
    samtools index -@ 64 ${female}nodup_autosomes/${base}_autosomes.bam
    samtools stats -@ 64 ${female}nodup_autosomes/${base}_autosomes.bam > ${female}nodup_autosome_stats/${base}_autosomes.stats
done

multiqc ${male}nodup_autosome_stats/
multiqc ${female}nodup_autosome_stats/
```
After extracting autosomal scaffolds for all BAM files, KMC was run as per:
```
ulimit -n 2048
for mkmc in ${male}nodup_autosomes/*_nodup_autosomes.bam
    do
    indiv=$(basename ${mkmc} _nodup_autosomes.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${mkmc} ${bayes}kmc/${indiv}_KMC.res ${bayes}kmc_tmp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom -k ${bayes}kmc/${indiv}_KMC.res -p 8
done
for fkmc in ${female}nodup_autosomes/*_nodup_autosomes.bam
    do
    indiv=$(basename ${fkmc} _nodup_autosomes.bam)
    echo "Running KMC for ${indiv}..."
    kmc -k55 -ci1 -m24 -t48 -fbam ${fkmc} ${bayes}kmc/${indiv}_KMC.res ${bayes}kmc_tmp/
    echo "Running makeBloom for ${indiv}..."
    bayesTyperTools makeBloom -k ${bayes}kmc/${indiv}_KMC.res -p 8
done
```
Variant clusters were estimated for all individuals. Note, no more than 30 indivs can be genotyped at once. The <samples>.tsv file should contain one sample per row with columns:
<sample_id>, <sex> and <kmc_output_prefix>
No header for the <samples>.tsv file.

```
batch=/kakapo-data/bayestyper/mantaB/
joint=/kakapo-data/bayestyper/mantaJ/

for samps in ${out}sample_batches/sample_batch*.tsv
    do
    base=$(basename ${samps} .tsv)
    for svs in $batch $joint
        do
        echo "Running cluster for $base for $svs..."
        bayesTyper cluster \
            --variant-file ${svs}05_*_combined.vcf.gz \
            --samples-file ${samps} --genome-file ${chr_ref} \
            --output-prefix ${svs}clusters/${base} --threads 24
        bayesTyper genotype \
            --variant-clusters-file ${svs}clusters/${base}_unit_1/variant_clusters.bin \
            --cluster-data-dir ${svs}clusters/${base}_cluster_data/ \
            --samples-file ${samps} --genome-file ${chr_ref} \
            --output-prefix ${svs}genotypes/${base}_genotypes --threads 24 -z
    done
done
```

## Merging SV genotype batches
To expedite the process of merging genotype batches, sites called in either the batch or joint Manta call sets were extracted, compressed with bgzip and indexed with tabix.
```
for vcf in ${out}bayestyper/joint_filtered/sample_batch{1..6}/*_genotypes.vcf
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Compressing $vcf with bgzip..."
    bgzip ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf.gz"
    tabix ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o ${joint}06_mantaJ_genotypes.vcf.gz \
    ${joint}sample_batch{1..6}_mantaJ_genotypes.vcf.gz

for vcf in ${batch}*_genotypes.vcf.gz
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Converting $vcf to bgzip compression..."
    gunzip $vcf
    bgzip ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf"
    tabix ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o ${batch}06_manta_genotypes.vcf.gz \
    ${batch}genotypes/sample_batch{1..6}_genotypes.vcf.gz
```

## Finding SV type

BayesTyper removes symbolic alleles from the genotype files. To make comparisons among SV tools similar, we attempted to relate SV types called by Manta.

As recommended [here](https://github.com/bioinformatics-centre/BayesTyper/issues/41), we tried to use the convertSeqToAlleleId script included in the source installation of BayesTyper (not available through Bioconda install).  

The format of this command is:  
<convertSeqToAlleleID> <input_VCF> <output_prefix> <minimum_SV_length>

Because SVs were already quality filtered prior to their inclusion in BayesTyper genotyping and to capture as many SVs in the output, we chose not to implement a minimum size threshold.  
```
convertSeqToAlleleID 04_raw_cuteSV_genotypes.vcf.gz 05_cuteSV_converted 0
convertSeqToAlleleID 04_raw_sniffles_genotypes.vcf.gz 05_sniffles_converted 0
```
The breakdown of SV types output by `convertSeqToAlleleID` command are as per:
|   SV Type   | Manta - Batch | Manta - Joint |
|:-----------:|:-------------:|:-------------:|
|  Deletion   |      969      |      900      |
| Duplication |      289      |      268      |
|  Insertion  |       0       |       0       |
|  Inversion  |     32,281    |     30,181    |
|    Total    |     33,539    |     31,349    |

Comparing with the output of SVs in the normalised file for genotyping, we see that the counts of SV types were not consistent with the expectations from the normalised file before running `bayesTyperTools combine`.  

|   SV Type   | Manta - Batch | Manta - Joint |
|:-----------:|:-------------:|:-------------:|
|  Deletion   |     1,981     |     1,419     |
| Duplication |      882      |      503      |
|  Insertion  |      752      |      221      |
|  Inversion  |     34,173    |     31,338    |
|    Total    |     37,791    |     33,481    |

However, the total number of SVs in the Genotype file prior to conversion with `convertSeqToAlleleID` for each call set was 34,134 and 31,747 for the batch and joint call sets repectively. While the conversion doesn't work for insertions, even adding these doesn't match with the expected total outputs. Taking the Batched call set as an example:  

33,539 total SVs converted + 752 expected insertions = 34,291 Total SVs  

This is 157 more SVs than expected from counts of the number of SVs in the combined candidates file (34,291 Total converted SVs - 34,134 Total not converted SVs).  

This leaves uncertainty around whether the conversion by `convertSeqToAlleleID` changed the type of SVs called by Manta (e.g., Insertions coverted to duplications). To explore this further, attempted to reconcile genotype variants with their original SV class called by Manta.  

### Identifying the locations of genotyped variants

Here we used `bcftools query` to extract the chromosome, start and end positions of genotyped variants.  

```
bcftools query \
    -f '%CHROM\t%POS\t%END\n' \
    ${batch}08_batch_filtered_trios.vcf > ${batch}batch_genos
bcftools query \
    -f '%CHROM\t%POS\t%END\n' \
    ${joint}08_joint_filtered_trios.vcf > ${joint}joint_genos
```
We also ensured that extracted regions were unique, and did not represent multiple SVs in the same region.

```
cat ${batch}batch_genos | sort | uniq | wc -l # Total unique records
cat ${batch}batch_genos | wc -l # Total records

cat ${joint}joint_genos | sort | uniq | wc -l # Total unique records
cat ${joint}joint_genos | wc -l # Total records
```
In both instances, the number of regions remained consistent.  

### Extracting overlaping SVs

Here we attempted to extract SVs from the normalised manta call sets prior to removal of symbolic alleles with `bayesTyperTools convert`.  

```
bcftools query \
    -T ${batch}batch_genos \
    -f '%CHROM\t%POS\t%END\t%SVTYPE\n' \
    ${batch}batch_filtered/02_batch_filtered_norm.vcf.gz | sort | uniq -c | sort -n | tail

bcftools query \
    -T ${joint}joint_genos \
    -f '%CHROM\t%POS\t%END\t%SVTYPE\n' \
    ${joint}joint_filtered/02_joint_filtered_norm.vcf.gz | sort | uniq -c | sort -n | tail
```

Found that 76 SVs don't overlap in the batch data and 126 SVs don't overlap in the joint data.  

### Annotate VCF with SV type  
First prepare the annotation file by identifying the resolvable SVs and their types:
```
bcftools query -f '%CHROM\t%POS\n' batch_filtered/07_batch_filtered_genotypes.vcf > batch_geno_sites
bcftools query -T batch_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' batch_filtered/02_batch_filtered_norm.vcf.gz > batch_conversion

bcftools query -f '%CHROM\t%POS\n' joint_filtered/07_joint_filtered_genotypes.vcf > joint_geno_sites
bcftools query -T joint_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' joint_filtered/02_joint_filtered_norm.vcf.gz > joint_conversion
```

Then edit the `batch_conversion` file so the first line in this file contains `#CHROM POS SVTYPE  SVLEN` with nano.  

And create a file that captures the information to be appended into the header of the VCF.

Contents of annots.hdr for example:
```
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant called by Manta">
##INFO=<ID=SVLEN,Number=1,Type=String,Description="Length of structural variant called by Manta">
```
Now time to annotate the file to continue with the genotype based analyses.
```
bgzip batch_conversion
bgzip joint_conversion

tabix -s1 -b2 -e2 batch_conversion.gz
tabix -s1 -b2 -e2 joint_conversion.gz

bcftools annotate -a batch_conversion.gz -h annots.hdr \
    -c CHROM,POS,SVTYPE,SVLEN -O v -o 09_batch_annotated.vcf \
    ${out}bayestyper/batch_filtered/07_batch_filtered_genotypes.vcf

bcftools annotate -a joint_conversion.gz -h annots.hdr \
    -c CHROM,POS,SVTYPE,SVLEN -O v -o 09_joint_annotated.vcf \
    ${out}bayestyper/joint_filtered/07_joint_filtered_genotypes.vcf
```

## Genotype quality filtering

SVs were simultaneously filtered for missingness and to ensure variant sites were included. As with the Manta SV call sets, for all sites where at least one individual carried an alternate allele (`N_PASS(GT=="alt") >= 1`), we allowed an individual missingness of 20% (`N_PASS(GT=="mis") < 35`) and required that at least 80% of all samples passed all four of BayesTyper's 'hard' filtering thresholds (`N_PASS(SAF!=0) < 35`).
```
 bcftools view --threads 24 \
    -i '((N_PASS(GT=="mis") < 35) & (N_PASS(GT=="alt") >= 1) & (N_PASS(SAF!=0) < 35))' \
    -O z -o ${cute}05_cuteSV_SV_genotypes.vcf.gz \
    ${cute}04_raw_cuteSV_DeepV_genotypes.vcf.gz

 bcftools view --threads 24 \
    -i '((N_PASS(GT=="mis") < 35) & (N_PASS(GT=="alt") >=1 ) & (N_PASS(SAF!=0) < 35))' \
    -O z -o ${sniff}05_sniffles_SV_genotypes.vcf.gz \
    ${sniff}04_raw_sniffles_DeepV_genotypes.vcf.gz
```
We are now ready to test Mendelian trios and prepare summary files.

## Mendelian Inheritance Tests
Tests of Mendelian Inheritance were conducted as per: 

```
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf \
    ${out}bayestyper/joint_filtered/08_joint_annotated.vcf
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf \
    ${out}bayestyper/batch_filtered/08_batch_annotated.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf



bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf
```
## Summarising the number of SVs carried by individuals

```
while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var"\tmantaB"}' >> ${out}bayestyper/summary/manta_batch_generations.tsv
    bcftools view -s ${indiv} ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var"\tmantaJ"}' >> ${out}bayestyper/summary/manta_joint_generations.tsv
done < /kakapo-data/metadata/generations.tsv
```


## Lineage Comparisons
```
mkdir -p ${out}bayestyper/lineage_{batch,joint}_comparisons

bcftools view -s M /kakapo-data/bwa/manta/raw_variants/batch_total.vcf.gz | \
    bcftools view -i 'GT=="alt"' -O z -o ${out}bayestyper/lineage_batch_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz
bcftools view -s M /kakapo-data/bwa/manta/raw_variants/joint_total.vcf.gz | \
    bcftools view -i 'GT=="alt"' -O z -o ${out}bayestyper/lineage_joint_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz

bcftools view -s ^M,G,F,R,N,P,O,S /kakapo-data/bwa/manta/raw_variants/batch_total.vcf.gz |\
    bcftools view -i 'GT=="alt"' -O z -o ${out}bayestyper/lineage_batch_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz
bcftools view -s ^M,G,F,R,N,P,O,S /kakapo-data/bwa/manta/raw_variants/joint_total.vcf.gz |\
    bcftools view -i 'GT=="alt"' -O z -o ${out}bayestyper/lineage_joint_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz

bcftools index ${out}bayestyper/lineage_batch_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz
bcftools index ${out}bayestyper/lineage_batch_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz

bcftools index ${out}bayestyper/lineage_joint_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz
bcftools index ${out}bayestyper/lineage_joint_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz

bcftools isec ${out}bayestyper/lineage_batch_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz \
    ${out}bayestyper/lineage_batch_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz \
    -p ${out}bayestyper/lineage_batch_comparisons/unfiltered/

bcftools isec ${out}bayestyper/lineage_joint_comparisons/unfiltered/RH_unfiltered_variants.vcf.gz \
    ${out}bayestyper/lineage_joint_comparisons/unfiltered/SI_unfiltered_variants.vcf.gz \
    -p ${out}bayestyper/lineage_joint_comparisons/unfiltered/

bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_batch_comparisons/unfiltered/0000.vcf > ${out}bayestyper/lineage_batch_comparisons/unfiltered/Fiordland_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_batch_comparisons/unfiltered/0001.vcf > ${out}bayestyper/lineage_batch_comparisons/unfiltered/Rakiura_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_batch_comparisons/unfiltered/0002.vcf > ${out}bayestyper/lineage_batch_comparisons/unfiltered/Shared_unfiltered_sites.txt

bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_joint_comparisons/unfiltered/0000.vcf > ${out}bayestyper/lineage_joint_comparisons/unfiltered/Fiordland_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_joint_comparisons/unfiltered/0001.vcf > ${out}bayestyper/lineage_joint_comparisons/unfiltered/Rakiura_unfiltered_private_sites.txt
bcftools query -f '%CHROM\t%POS\n' \
    ${out}bayestyper/lineage_joint_comparisons/unfiltered/0002.vcf > ${out}bayestyper/lineage_joint_comparisons/unfiltered/Shared_unfiltered_sites.txt

while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_batch_comparisons/unfiltered/Fiordland_unfiltered_private_sites.txt raw_variants/batch_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tFiordland_unfiltered_lineage\n' >> ${out}mantaB_lineage_counts.tsv
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_batch_comparisons/unfiltered/Rakiura_unfiltered_private_sites.txt raw_variants/batch_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tRakiura_unfiltered_lineage\n' >> ${out}mantaB_lineage_counts.tsv
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_batch_comparisons/unfiltered/Shared_unfiltered_sites.txt raw_variants/batch_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tShared_unfiltered_lineage\n' >> ${out}mantaB_lineage_counts.tsv
done < /kakapo-data/metadata/generations.tsv

while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_joint_comparisons/unfiltered/Fiordland_unfiltered_private_sites.txt ${out}raw_variants/joint_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tFiordland_unfiltered_lineage\n' >> ${out}mantaJ_lineage_counts.tsv
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_joint_comparisons/unfiltered/Rakiura_unfiltered_private_sites.txt ${out}raw_variants/joint_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tRakiura_unfiltered_lineage\n' >> ${out}mantaJ_lineage_counts.tsv
    bcftools view -s ${indiv} -R ${out}bayestyper/lineage_joint_comparisons/unfiltered/Shared_unfiltered_sites.txt ${out}raw_variants/joint_total.vcf.gz | bcftools query -i 'GT=="alt"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVTYPE\tShared_unfiltered_lineage\n' >> ${out}mantaJ_lineage_counts.tsv
done < /kakapo-data/metadata/generations.tsv
```