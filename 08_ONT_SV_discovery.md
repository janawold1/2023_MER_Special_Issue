# Long-read SV discovery and genotyping
Here are some brief notes of how I conducted long-read SV discovery with Sniffles and refined calles with Jasmine. 

## Initial SV calls with cuteSV and Sniffles 
To prepare called SVs for refinement using the Jasmine pipeline, Sniffles was run using sensitive parameters an all supporting reads were reported as recommended below:

```
for male in Bill Blades Gulliver Rangi Sass Smoko
    do
    echo "Running sniffles for ${male}..."
    sniffles --input ${data}winnowmap/alignment/bam/${male}.bam \
        --reference ${mref} \
        --vcf ${data}winnowmap/sniffles/${male}.vcf \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64
done

for female in Bella Kuia Margaret-Maree Sue
    do
    echo "Running sniffles for ${female}..."
    sniffles --input ${data}winnowmap/alignment/bam/${female}.bam \
        --reference ${fref} \
        --vcf ${data}winnowmap/sniffles/${female}.vcf \
        --minsvlen 50 \
        --minsupport 2 \
        --threads 64
done
```

## Running initial cuteSV calls

```
for bam in ${data}winnowmap/alignment/bam/*.bam
    do
    echo "Running cuteSV for ${base}..."
    if [[ $base == @(B|C|F|I|J) ]]
        then
        cuteSV ${bam} ${mref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
        else
        cuteSV ${bam} ${fref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            cuteSV ${bam} ${fref} \
            ${data}winnowmap/cuteSV/${base}.vcf cuteSV/ \
            --max_cluster_bias_INS 100 \
            --diff_ratio_merging_INS 0.3 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24 \
            --max_cluster_bias_DEL 100 \
            --diff_ratio_merging_DEL 0.3 \
            --threads 24
    fi &
done
```
## Running Jasmine
After SV discovery calling, calls for the same individual using both Sniffles and cuteSV were concatenated into a single VCF for breakpoint refinement and SV normalisation.
```
ref=/kakapo-data/References/kakapo_full_ref.fa

echo "Refining merged SV calls"
cd ${data}winnowmap/merged_SVs/
for vcf in *.vcf
    do
    indiv=$(basename $vcf .vcf)
    echo "Running Jasmine and IRIS for $indiv..."
    jasmine --dup_to_ins --preprocess_only file_list=${vcf} --comma_filelist genome_file=${ref}
    iris threads=24 vcf_in=output/${indiv}_dupToIns.vcf genome_in=${ref} reads_in=${data}winnowmap/alignment/bam/${indiv}.bam vcf_out=output/${indiv}_dupToIns_irisRefined.vcf
    echo "Normalising SV types for $indiv..."
    jasmine --preprocess_only --pre_normalize file_list=output/${indiv}_dupToIns_irisRefined.vcf --comma_filelist
done
```
After the normalisation of the files, high-confidence calls were marked. Unfortunately, the long-read data assessed here varied widely from 4x - 12x coverage. A minimum depth of 20x would have been ideal.

The low sequence depth of our samples makes is somewhat challenging to curate high quality calls with the jasmine pipeline. This is due to the recommended minimum number of reads a variant needs to pass (```spec_reads```) is about 25% of the average coverage, which would place some samples < 2x coverage to call a high-quality variant. In light of this, a minimum ```spec_reads``` of 2 was used for all 8 samples with < 10x coverage and 3x was used for the remaining 2 samples with >=10x coverage.

```
for indiv in D C2 F G H I J L
    do
    echo "Calling high-quality calls for ${indiv}..."
    jasmine file_list=output/${indiv}_dupToIns_irisRefined_normalizeTypes.vcf --comma_filelist --preprocess_only --mark_specific --spec_reads=2 spec_len=50
done

for indiv in B C
    do
    echo "Calling high-quality calls for ${indiv}..."
    jasmine file_list=output/${indiv}_dupToIns_irisRefined_normalizeTypes.vcf --comma_filelist --preprocess_only --mark_specific --spec_reads=3 spec_len=50
done

for vcf in output/*_dupToIns_irisRefined_normalizeTypes_markedSpec.vcf
    do
    indiv=$(basename $vcf _dupToIns_irisRefined_normalizeTypes_markedSpec.vcf )
    echo "Removing duplicate calls for ${indiv}..."
    jasmine file_list=$vcf --comma_filelist max_dist=200 --allow_intrasample out_file=${indiv}_jasmine.vcf --nonlinear_dist
done
```
SV calls were merged across all samples, insertions were converted back to duplications.
```
ls *_jasmine.vcf > vcf_list.txt

jasmine file_list=vcf_list.txt out_file=kakapo_mergedSVs.vcf

jasmine --dup_to_ins --postprocess_only out_file=kakapo_mergedSVs.vcf
```
The output of the conversion step of insertions back to duplications for the SV call set was: ```Number of insertions converted back to duplications: 135 out of 25322 total variants```

Finally, low-confidence / imprecise call were removed.
```
cat kakapo_mergedSVs.vcf | grep -v 'IMPRECISE;' > kakapo_mergedSVs_precise.vcf
```

## Genotyping with Paragraph
Genotyping SVs was challenging due to formatting SV calls to be compatible with Paragraph. First, Transversions and unplaced scaffolds were excluded while the VCF was normalised.
```
bcftools view -e 'SVTYPE="TRAV"' -O v -o sniffles_merged_autosomes_noTRA.vcf.gz sniffles_merged_autosomes.vcf
bcftools norm -T ../metadata/kakapo_chromosome_scaffolds.bed --threads 24 --check-ref s -D -f /kakapo-data/References/kakapo_full_ref.fa -O v -o candidates.vcf sniffles_merged_autosomes_noTRA.vcf.gz
```

Once this was complete, all calls within 125 bp of chromosomal scaffold start and ends were exluded (n = 2).
```
bcftools view -T ^bad_sites.tsv -O v -o candidates_noBadSites.vcf candidates.vcf
```
