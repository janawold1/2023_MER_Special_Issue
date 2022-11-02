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

```
bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}joint_filtered/01_joint_filtered.SVs.vcf.gz \
    ${raw_data}joint_total.vcf.gz

bcftools view -i '(FILTER=="PASS" & FMT/PR >= 10 & FORMAT/FT == "PASS" & SVLEN<=-50 & SVTYPE=="DEL") | (FILTER=="PASS" & PR >= 10 & FMT/FT == "PASS" & SVLEN>=50 & SVTYPE!="BND")' \
    -T ^${exclude} \
    -O z -o ${out}batch_filtered/01_batch_filtered.SVs.vcf.gz \
    ${raw_data}batch_total.vcf.gz
```

## Convert Manta VCF to remove symbolic alleles
This converts all symbolic alleles to a full sequence. The conversion results in large VCFs, so ensure that the output is gzipped iwth the -z flag. Note, bayesTyper does NOT support bgzip files. BayesTyper benefits from including SNP data as it helps with initialising the kmers. Therefore it's recommended to conduct SNP discovery with either GATK HaplotypeCaller or FreeBayes. I used the kākāpō125+ DeepVariant SNP call set prepared by (Guhlin et al. 2022)[https://www.biorxiv.org/content/10.1101/2022.10.22.513130v1.abstract].

For combining SV calls with the DeepVariant SNPs, all chromosome sizes were used to match kākāpō125+ and NCBI chromosome IDs. bcftools annotate was used to rename chromosomes to be consistent with NCBI scaffold names. SNPs were filtered down to only include regions of interest (i.e., scaffolds included for SV analysis). But after converting chromosome names, needed to sort VCF header and body of DeepVariant file to match the manta VCFs. Sorted the body with [this solution](https://www.biostars.org/p/84747/), then manually edited the header.

```
bcftools norm -c we -f ${ref} -O v -o ${out}bayestyper/Trained_SV_scaffolds_renamed_norm.vcf \
    ${out}bayestyper/Trained_SV_scaffolds_renamed.vcf

bgzip ${out}Trained_SV_scaffolds_renamed_norm.vcf; tabix -pvcf ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz

cat chr_list.txt | xargs tabix -h ${out}Trained_SV_scaffolds_renamed_norm.vcf.gz > ${out}Trained_SV_scaffolds_renamed.sorted.vcf
```
Allele conversion and file combination was then conducted as per:

```
for manta in ${out}{batch,joint}_filtered/
    do
    file=$(basename ${manta})
    echo "Normalising 01_$file.SVs.vcf"
    bcftools norm --threads 24 -f ${chr_ref} -m -any \
        -O v -o ${manta}02_${file}_bayestyper_candidates_norm.vcf ${manta}01_$file.SVs.vcf
    echo "Running convertAllele for 02_${file}_norm.vcf..."
    bayesTyperTools convertAllele --variant-file ${manta}02_bayestyper_candidates_norm.vcf --genome-file ${chr_ref} \
        --output-prefix ${manta}03_${file}_converted
    echo "Finally combining DeepVariant and Manta alleles into 04_${file}_mantaDV_candidates.vcf.gz"
    bayesTyperTools combine -v ${deepV},MANTA:${manta}03_${file}_converted.vcf.gz \
        -o ${manta}04_${file}_combined_variants -z
done
```

## Running KMC and makeBloom
It is really easy to over-resource KMC, that is to say if you provide it with too much memory/threads the program freezes. I have had the best luck with 24Gb RAM and 48 threads as per below. It also found that I needed to run KMC one individual at a time. I think threads clashed when I attempted to run multiple individuals at once since KMC uses the same temp file naming convention each time it runs.

To aid in genotyping, KMC was run using BAMs filtered for autosomal scaffolds only using the solution provided (here)[https://www.biostars.org/p/302771/].

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
            --variant-file ${svs}04_*_combined_variants.vcf.gz \
            --samples-file ${samps} --genome-file ${chr_ref} \
            --output-prefix ${svs}clusters/${base} --threads 24
        bayesTyper genotype \
            --variant-clusters-file ${svs}clusters/${base}_unit_1/variant_clusters.bin \
            --cluster-data-dir ${svs}clusters/${base}_cluster_data/ \
            --samples-file ${samps} --genome-file ${chr_ref} \
            --output-prefix ${svs}genotypes/${base}_genotypes --threads 24
    done
done
```

## Merging SV genotype batches

```
for vcf in ${out}bayestyper/joint_filtered/sample_batch{1..6}/*_genotypes.vcf.gz
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Converting $vcf to bgzip compression..."
    gunzip $vcf
    bgzip ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf"
    tabix ${out}bayestyper/joint_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o joint_filtered/06_manta_genotypes.vcf.gz \
    joint_filtered/sample_batch1/sample_batch1_genotypes.vcf.gz \
    joint_filtered/sample_batch2/sample_batch2_genotypes.vcf.gz \
    joint_filtered/sample_batch3/sample_batch3_genotypes.vcf.gz \
    joint_filtered/sample_batch4/sample_batch4_genotypes.vcf.gz \
    joint_filtered/sample_batch5/sample_batch5_genotypes.vcf.gz \
    joint_filtered/sample_batch6/sample_batch6_genotypes.vcf.gz

 bcftools view --threads 24 -i '(ACO="MANTA") & ((N_PASS(GT=="mis") < 17) & (N_PASS(GT="alt")>=1))' \
    -O v -o joint_filtered/07_joint_filtered_genotypes.vcf \
    joint_filtered/06_joint_genotypes.vcf.gz

for vcf in ${out}bayestyper/batch_filtered/sample_batch{1..6}/*_genotypes.vcf.gz
    do
    batch=$(basename $vcf _genotypes.vcf.gz)
    echo "Converting $vcf to bgzip compression..."
    gunzip $vcf
    bgzip ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf
    echo "Now indexing $vcf"
    tabix ${out}bayestyper/batch_filtered/${batch}/${batch}_genotypes.vcf.gz
done

bcftools merge -m id -O z -o batch_filtered/06_manta_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch1/sample_batch1_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch2/sample_batch2_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch3/sample_batch3_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch4/sample_batch4_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch5/sample_batch5_genotypes.vcf.gz \
    ${out}bayestyper/batch_filtered/sample_batch6/sample_batch6_genotypes.vcf.gz

 bcftools view --threads 24 -i '(ACO="MANTA") & ((N_PASS(GT=="mis") < 17) & (N_PASS(GT="alt")>=1))' \
    -O v -o ${out}bayestyper/batch_filtered/07_batch_filtered_genotypes.vcf \
    ${out}bayestyper/batch_filtered/06_batch_genotypes.vcf.gz

```

## Finding SV type

BayesTyper removes symbolic alleles from the genotype files. To make comparisons among SV tools similar, the SV type called by Manta was resoved as per below.
1) First identify the locations of genotyped variants:
```
bcftools query -f '%CHROM\t%POS\n' ${out}bayestyper/batch_filtered/08_batch_filtered_trios.vcf > ${out}bayestyper/batch_genos
bcftools query -f '%CHROM\t%POS\n' ${out}bayestyper/joint_filtered/08_joint_filtered_trios.vcf > ${out}bayestyper/joint_genos
```

2) Count overlaping SVtypes:
```
bcftools query -T ${out}bayestyper/batch_genos -f '%SVTYPE\n' \
    ${out}bayestyper/batch_filtered/02_batch_filtered_norm.vcf.gz | sort | uniq -c

bcftools query -T ${out}bayestyper/joint_genos -f '%SVTYPE\n' \
    ${out}bayestyper/joint_filtered/02_joint_filtered_norm.vcf.gz | sort | uniq -c
```

Found that 76 SVs don't overlap in the batch data and 126 SVs don't overlap in the joint data. 

3) Find non-overlapping sites (i.e., those present in genotyped output, but absent in call set):
```

```

4) Annotate VCF for individual counts:
First prepare the annotation file by identifying the resolvable SVs and their types:
```
bcftools query -f '%CHROM\t%POS\n' batch_filtered/07_batch_filtered_genotypes.vcf > batch_geno_sites
bcftools query -T batch_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' batch_filtered/02_batch_filtered_norm.vcf.gz > batch_conversion

bcftools query -f '%CHROM\t%POS\n' joint_filtered/07_joint_filtered_genotypes.vcf > joint_geno_sites
bcftools query -T joint_geno_sites -f '%CHROM\t%POS\t%SVTYPE\t%SVLEN\n' joint_filtered/02_joint_filtered_norm.vcf.gz > joint_conversion
```

Then edit the ```batch_conversion``` file so the first line in this file contains ```#CHROM POS SVTYPE  SVLEN``` with nano. 

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

## Mendelian Inheritance Tests
Tests of Mendelian Inheritance were conducted as per: 

```
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf \
    ${out}bayestyper/joint_filtered/08_joint_annotated.vcf
bcftools +mendelian -m a -T ${trio} -O v -o ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf \
    ${out}bayestyper/batch_filtered/08_batch_annotated.vcf

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_joint_manta_genofilter\n' \
    ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv



bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.05' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.05_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.1' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.1_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv

bcftools query -i '(MERR / N_PASS(GT!="mis")) <=0.2' \
    -f '%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\t0.2_fail_batch_manta_genofilter\n' \
    ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf >> ${out}bayestyper/summary/manta_mendel.tsv
```
## Summarising the number of SVs carried by individuals

```
while read -r line
    do
    indiv=$(echo $line | awk '{print $1}')
    gen=$(echo $line | awk '{print $2}')
    echo "Counting SVs for $indiv..."
    bcftools view -s ${indiv} ${out}bayestyper/batch_filtered/09_batch_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var}' >> ${out}bayestyper/summary/manta_batch_generations.tsv
    bcftools view -s ${indiv} ${out}bayestyper/joint_filtered/09_joint_annotated_trios.vcf | bcftools query -i 'GT!= "RR" & GT!="mis"' -f '[%SAMPLE]\t%CHROM\t%POS\t%END\t%SVLEN\t%SVTYPE\n' | awk -v var="$gen" '{print $0"\t"var}' >> ${out}bayestyper/summary/manta_joint_generations.tsv
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