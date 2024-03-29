# Long-read SV discovery and genotyping
Here are some brief notes of how I conducted long-read SV discovery with cuteSV and Sniffles. All SV calls by CuteSV and Sniffles were individually refined with Jasmine prior to comparisons.

## Initial SV calls with cuteSV and Sniffles 
Sniffles was run using sensitive parameters an all supporting reads were reported as recommended below:

```
align=/kakapo-data/ONT/winnowmap/alignment/bam/
snif=/kakapo-data/ONT/winnowmap/sniffles/
cute=/kakapo-data/ONT/winnowmap/cuteSV/
mref=/kakapo-data/references/kakapo_no_Wchromosome.fna
fref=/kakapo-data/references/kakapo_full_ref.fa

for male in bird_ID
    do
    echo "Running sniffles for ${male}..."
    sniffles --input ${align}${male}_3kb.bam \
        --reference ${mref} \
        --snf ${snif}snf/${male}.snf \
        --non-germline \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64 &
    sniffles --input ${align}${male}_3kb.bam \
        --reference ${mref} \
        --vcf ${snif}vcf_for_jasmine/${male}.vcf \
        --non-germline \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64
    bcftools view -T /kakapo-data/metadata/autosomes.txt \
        -O v -o ${snif}vcf_for_jasmine/${male}_autosomes.vcf \
        ${snif}vcf_for_jasmine/${male}.vcf
done
wait
for female in bird_ID
    do
    echo "Running sniffles for ${female}..."
    sniffles --input ${align}${female}_3kb.bam \
        --reference ${fref} \
        --snf ${snif}snf/${female}.snf \
        --non-germline \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64 &
    sniffles --input ${align}${female}_3kb.bam \
        --reference ${fref} \
        --vcf ${snif}vcf_for_jasmine/${female}.vcf \
        --non-germline \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64
    bcftools view -T /kakapo-data/metadata/autosomes.txt \
        -O v -o ${snif}vcf_for_jasmine/${female}_autosomes.vcf \
        ${snif}vcf_for_jasmine/${female}.vcf
done
wait
sniffles --input ${snif}snf/*.snf --vcf ${snif}01_sniffles_SVs.vcf
bcftools view -T autosomes.txt -O v -o ${snif}02_sniffles_autosomeSVs.vcf ${snif}01_sniffles_SVs.vcf
```

## Running initial cuteSV calls
Important to note, cuteSV cannot use gzipped reference genomes. 

```
for bam in ${align}*_3kb.bam
    do
    base=$(basename $bam _3kb.bam)
    echo "Running cuteSV for ${base}..."
    if [[ $base == @(Bill|Blades|Gulliver|Rangi) ]]
        then
        cuteSV ${bam} ${mref} \
            ${cute}vcf_for_jasmine/${base}.vcf $cute \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
        bcftools view -T /kakapo-data/metadata/autosomes.txt \
            -O v -o ${snif}vcf_for_jasmine/${base}_autosomes.vcf \
            ${cute}vcf_for_jasmine/${base}.vcf
        bcftools norm --check-ref s -f $ref \
            -O z -o ${cute}rawSVs/${base}_autosomes_norm.vcf.gz \
            ${cute}vcf_for_jasmine/${base}_autosomes.vcf
        else
        cuteSV ${bam} ${fref} \
            ${cute}vcf_for_jasmine/${base}.vcf $cute \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
        bcftools view -T /kakapo-data/metadata/autosomes.txt \
            -O v -o ${cute}vcf_for_jasmine/${base}_autosomes.vcf \
            ${cute}vcf_for_jasmine/${base}.vcf
        bcftools norm --check-ref s -f $ref \
            -O z -o ${cute}rawSVs/${base}_autosomes_norm.vcf.gz \
            ${cute}vcf_for_jasmine/${base}_autosomes.vcf
    fi
done
```
Because Jasmine refines SV calls on the individual level, the total number of raw SVs called by CuteSV were estimated by merging all the normalised initial individual calls and counting the number of SVs.
```
bcftools merge -O v -o ${cute}01_cute_SV_rawSVs.vcf ${cute}rawSVs/*norm.vcf.gz
```
## Running Jasmine
After SV discovery calling, calls for each individual were concatenated into a single VCF for breakpoint refinement and SV normalisation. Below is an example of how jasmine was run for Sniffles, with this being repeated for CuteSV calls. Only SV calls on autosomal scaffolds were used to reduce the affects of poorly resolved regions (e.g., unplaced scaffolds, sex chromosomes).  
```
ref=/kakapo-data/References/kakapo_full_ref.fa
bam=/kakapo-data/ONT/winnowmap/alignment/bam/

echo "Refining merged SV calls"

for vcf in ${snif}vcf_for_jasmine/*_autosomes.vcf
    do
    indiv=$(basename $vcf _autosomes.vcf)
    echo "Running Jasmine and IRIS for $indiv..."
    jasmine --dup_to_ins file_list=${vcf} out_file=${snif}vcf_for_jasmine/${indiv}_dupToIns.vcf \
        --comma_filelist genome_file=${ref}
    iris threads=24 vcf_in=${snif}vcf_for_jasmine/${indiv}_dupToIns.vcf \
        genome_in=${ref} reads_in=${bam}${indiv}_3kb.bam vcf_out=${snif}vcf_for_jasmine/${indiv}_irisRefined.vcf
    jasmine --pre_normalize file_list=${snif}vcf_for_jasmine/${indiv}_irisRefined.vcf --comma_filelist \
        out_file=${snif}vcf_for_jasmine/${indiv}_norm.vcf
    jasmine file_list=${snif}vcf_for_jasmine/${indiv}_norm.vcf --comma_filelist \
        --mark_specific spec_reads=10 spec_len=50 \
        out_file=${snif}vcf_for_jasmine/${indiv}_mark_spec.vcf 
    jasmine file_list=${snif}vcf_for_jasmine/${indiv}_mark_spec.vcf --comma_filelist max_dist=200 \
        out_file=${snif}vcf_for_jasmine/${indiv}_dupRmv.vcf --nonlinear_dist
done

ls ${snif}vcf_for_jasmine/*_dupRmv.vcf > ${snif}jasmine_merge.txt
jasmine file_list=${snif}jasmine_merge.txt out_file=${snif}02_jasmine_raw.vcf
jasmine --dup_to_ins file_list=${snif}02_jasmine_raw.vcf --comma_filelist out_file=${snif}03_jasmine_refined.vcf genome_file=${ref}
```
Unfortunately, the long-read data assessed here varied widely from 5x - 12x coverage for autosomal scaffolds. A minimum depth of 20x would have been ideal.

The low sequence depth of our samples makes it somewhat challenging to curate high quality calls with the jasmine pipeline. This is due to the recommended minimum number of reads a variant needs to pass (`spec_reads`) is about 25% of the average coverage, which would place some samples < 2x coverage to call a high-quality variant. To conservatively call SVs, the default minimum `spec_reads` of 10 was used for all 7 samples included here.

Finally, low-confidence / imprecise call were removed. Notably, could not filter for `IS_SPECIFIC` in Sniffles calls as no SVs remained for comparisons.
```
cat ${snif}03_jasmine_refine.vcf | grep -v 'IMPRECISE;' > ${snif}04_filtered.vcf
cat ${cute}03_jasmine_refine.vcf | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' > ${cute}04_filtered.vcf

bcftools sort -O v -o ${snif}05_filtered_sorted.vcf ${snif}04_filtered.vcf
bcftools sort -O v -o ${cute}05_filtered_sorted.vcf ${cute}04_filtered.vcf
```

## Filtering for genotyping with BayesTyper and Paragraph
Genotyping SVs was challenging due to formatting SV calls to be compatible with Paragraph. First, Transversions and unplaced scaffolds were excluded while the VCF was normalised.
```
bcftools view -e 'SVTYPE="TRA"' -O v -o ${snif}06_filtered_noTRAV.vcf ${snif}05_filtered_sorted.vcf
bcftools view -e 'SVTYPE="TRA"' -O v -o ${cute}06_filtered_noTRAV.vcf ${cute}05_filtered_sorted.vcf

bcftools norm --check-ref s -D -f ${ref} -O v -o ${snif}07_paragraph_candidates.vcf ${snif}06_filtered_noTRAV.vcf
bcftools norm --check-ref s -D -f ${ref} -O v -o ${cute}07_paragraph_candidates.vcf ${cute}06_filtered_noTRAV.vcf
```
This VCF was used to estimate the number of SVs that passed filtering thresholds. However, not all of these SVs are compatible with the genotyping tools used here. More details are provided in the relevant Markdown files.  
