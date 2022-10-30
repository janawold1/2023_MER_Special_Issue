# Genotyping with Paragraph
## Creating manifest files
Paragraph recommends running population-scale genotyping in the single-sample mode by:
    1) Creating a manifest for each single sample.
    2) Run Paragraph for each manifest with the `-M` option for each sample to its depth.
    3) Merge all genotypes.vcf.gz to create final VCF with `bcftools merge`.

Sample manifests have four required columns outlined as per below:

|      id      |                path             |  depth   |  read length  |
|:------------:|:-------------------------------:|:--------:|:-------------:|
|   SAMPLE_1   |    /path/to/bam/SAMPLE_1.bam    |    20    |      125      |
|     ...      |               ...               |    ...   |      ...      |
|   SAMPLE_N   |    /path/to/bam/SAMPLE_N.bam    |    20    |      125      |

To simplify genotyping, bam files were subsetted to autosomal chromosomes.
```
awk '{printf("%s\t0\t%s\n",$1,$2);}' ${ref}.fai | grep NC_ | grep -v NC_044301 | grep -v NC_044302 > /kakapo-data/reference/kakapo_autosomes.bed

male=/kakapo-data/bwa/bwa_male/
female=/kakapo-data/bwa/bwa_female/

for bam in ${male}nodup_bam/batch*/*_nodup.bam
    do
    base=$(basename $bam .bam)
    echo $base 
    samtools view -@ 64 -L /kakapo-data/references/kakapo_autosomes.bed -o ${male}nodup_autosomes/${base}_autosomes.bam $bam
    samtools index -@ 64 ${male}nodup_autosomes/${base}_autosomes.bam
    samtools stats -@ 64 ${male}nodup_autosomes/${base}_autosomes.bam > ${male}nodup_autosome_stat/${base}_autosomes.stats
done
```
And sample manifests were created as per:
```
for stat in nodup_autosome_stat/*.stats
    do
    data=$(cat $stat | grep ^RL | cut -f 2-)
    samp=$(basename $stat _nodup_autosomes.stats)
    while read -r line
        do
        read=$(echo $line | awk '{print $1}')
        count=$(echo $line | awk '{print $2}')
        awk -v var="$count" -v var2="$read" 'BEGIN{for(c=0;c<var;c++) print var2}' >> nodup_autosome_stat/${samp}_readlength.txt
    done < <(printf "$data\n")
    avg=$(cat nodup_autosome_stat/${samp}_readlength.txt | awk '{sum += $1};END {print sum/NR}')
    depth=$(samtools depth -a nodup_autosomes/${samp}_nodup_autosomes.bam | awk '{sum+=$3}; END {print sum/NR}')
    printf "id\tpath\tdepth\tread length\n$samp\t/kakapo-data/bwa/bwa_female/nodup_autosomes/${samp}_nodup_autosomes.bam\t$depth\t$avg\n" > /kakapo-data/paragraph/manifests/${samp}.tsv
done
```
## Running Paragraph
Once file manifests were generated for all samples, Paragraph was run for SVs called by CuteSV and Sniffles. However, Paragraph cannot genotype SVs that occur within the average short-read length of either end of the chromosome. The largest mean short-read length was 150bp, so all SVs called within the first or last 150 of each chromosome were filtered.
```
awk '{print $1"\t"$2"\t"$2+150}' /kakapo-data/references/kakapo_autosomes.bed > rmv.tsv
awk '{print $1"\t"$3-150"\t"$3}' /kakapo-data/references/kakapo_autosomes.bed >> rmv.tsv

```
This removed 2 deletions, 5 duplications and 1 insertion from the CuteSV callset, whereas 5 insertions fell within this range and an additional 17 insertions in the Sniffles callset were excluded due to missing `INFO/SEQ` and a large inversion that was consistently problematic for genotyping (NC_044293.2:7407591-7441145). These SVs were found using:
```
bcftools query -f '%POS\t%CHROM\t%SVTYPE\n' paragraph/cuteSV/07_paragraph_candidates.vcf.gz | head -n 6 | awk '{print $2"\t"$1}' > paragraph/cuteSV/rmv.tsv

bcftools view -T ^paragraph/cuteSV/rmv.tsv -O v -o paragraph/cuteSV/08_paragraph_filteredcandidates.vcf paragraph/cuteSV/07_paragraph_candidates.vcf.gz
```
Finally, samples were genotyped as per below. As [recommended](https://github.com/Illumina/paragraph#TestExample), the maximum number of supporting reads `-M` was set to 20x the depth.
```
cute=/kakapo-data/ONT/winnowmap/cuteSV/08_paragraph_filteredcandidates.vcf
snif=/kakapo-data/ONT/winnowmap/sniffles/08_paragraph_filteredcandidates.vcf
ref=/kakapo-data/references/kakapo_autosomes.fa
out=/kakapo-data/paragraph/

for samp in /kakapo-data/paragraph/manifests/*.tsv
    do
    indiv=$(basename ${samp} .tsv)
    M=$(cat $samp | awk '{print $3*20}')
    echo "Running Paragraph using SVs called with CuteSV for ${indiv}..."
    multigrmpy.py -i $cute \
        -m $samp \
        -r $ref \
        -o ${out}cuteSV/genotypes/${indiv} \
        -t 16 \
        -M ${M} &
    echo "Running Paragraph using SVs called with Sniffles for ${indiv}..."
    multigrmpy.py -i $snif \
        -m $samp \
        -r $ref \
        -o ${out}sniffles/genotypes/${indiv} \
        -t 16 \
        -M ${M}
    wait
done
```
# SV summaries
SV characteristics were captured as below for plotting in R. See 12_SV_summaries.Rd for more details.

```
printf "chrom\tpos\tsvlen\tsvtype\tdata\n" > ${cute}cuteSV_summary.tsv
printf "chrom\tpos\tsvlen\tsvtype\tdata\n" > ${snif}sniffles_summary.tsv

bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_unfiltered\n' ${cute}01_cuteSV_rawSVs.vcf >> ${cute}cuteSV_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_SVfiltered\n' ${cute}07_paragraph_candidates.vcf >> ${cute}cuteSV_summary.tsv
#bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_genofiltered\n' ${cute}01_cuteSV_rawSVs.vcf >> ${cute}cuteSV_summary.tsv

bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_unfiltered\n' ${snif}01_sniffles_rawSVs.vcf >> ${snif}sniffles_summary.tsv
bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_SVfiltered\n' ${snif}07_paragraph_candidates.vcf >> ${snif}sniffles_summary.tsv
#bcftools query -f '%CHROM\t%POS\t%SVLEN\t%SVTYPE\tcuteSV_genofiltered\n' ${cute}01_cuteSV_rawSVs.vcf >> ${snif}sniffles_summary.tsv

while read -r line
    do
    scaf=$(echo ${line} | awk '{print $1}')
    chro=$(echo ${line} | awk '{print $2}')
    echo Converting $scaf to $chro
    sed -i "s/$scaf/$chro/g" sniffles_summary.tsv
    sed -i "s/$scaf/$chro/g" cuteSV_summary.tsv
done < convert_chr.txt
```