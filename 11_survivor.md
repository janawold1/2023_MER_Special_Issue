# Finding overlapping SV calls with Survivor
Survivor was used to identify SVs that overlapped in the data filtered for SV call quality among all six SV discovery tools/strategies. Survivor uses a config file that provides the path to VCFs for comparison.
```
#!/bin/bash -e
data=/kakapo-data/survivor/

for file in ${data}input/*.txt
    do
    num=$(wc -l ${file} | awk '{print $1}')
    base=$(basename ${file} .txt)
    for i in {0,50,500,1000}
        do
        echo "Running SURVIVOR for ${base} with $i interval..."
        SURVIVOR merge $file ${i} ${num} 1 1 0 50 ${data}outputs/${base}_${i}_distance.vcf
    done &
done
```

```
for vcf in ${data}outputs/*_distance.vcf
    do
    base=$(basename $vcf _distance.vcf)
    echo "Counting SVs for $base..."
    total=$(bcftools query -f '%SVTYPE\n' ${vcf} | wc -l)
    del=$(bcftools query -i 'SVTYPE=="DEL"' -f '%SVTYPE\n' ${vcf} | wc -l)
    dup=$(bcftools query -i 'SVTYPE=="DUP"' -f '%SVTYPE\n' ${vcf} | wc -l)
    ins=$(bcftools query -i 'SVTYPE=="INS"' -f '%SVTYPE\n' ${vcf} | wc -l)
    inv=$(bcftools query -i 'SVTYPE=="INV"' -f '%SVTYPE\n' ${vcf} | wc -l)
    printf "${base}\t${total}\t${del}\t${dup}\t${ins}\t${inv}\n" >> ${data}overlap_counts.csv
done
```

# Overlaps in genotype outputs
```
for file in ${data}geno_input/*.txt
    do
    num=$(wc -l ${file} | awk '{print $1}')
    base=$(basename ${file} .txt)
    for i in {0,50,500,1000}
        do
        echo "Running SURVIVOR for ${base} with $i interval..."
        SURVIVOR merge $file ${i} ${num} 0 1 0 50 ${data}outputs/${base}_${i}bp_distance.vcf
    done &
done

for vcf in ${data}geno_output/*_distance.vcf
    do
    base=$(basename $vcf _distance.vcf)
    echo "Counting SVs for $base..."
    total=$(bcftools query -f '%SVTYPE\n' ${vcf} | wc -l)
    del=$(bcftools query -i 'SVTYPE=="DEL"' -f '%SVTYPE\n' ${vcf} | wc -l)
    dup=$(bcftools query -i 'SVTYPE=="DUP"' -f '%SVTYPE\n' ${vcf} | wc -l)
    ins=$(bcftools query -i 'SVTYPE=="INS"' -f '%SVTYPE\n' ${vcf} | wc -l)
    inv=$(bcftools query -i 'SVTYPE=="INV"' -f '%SVTYPE\n' ${vcf} | wc -l)
    printf "${base}\t${total}\t${del}\t${dup}\t${ins}\t${inv}\n" >> ${data}geno_overlap_counts.tsv
done
```

## Genotype success
To assess the number of SVs in the all v all comparison of SVs genotyped for call quality, we created a bed file
```
bcftools query -f '%CHROM\t%POS\t%END\n' ${data}outputs/allvall_50_distance.vcf > ${data}allvall_50bp_distance.bed
```

Then for each of the SVs compared the number of SVs by type to the number of expected. 
```
for vcf in ${data}geno_vcfs/*.vcf
    do
    echo $vcf
    bcftools query -T ${data}allvall_50bp_distance.bed -f '%SVTYPE\n' $vcf | sort | uniq -c
done
```