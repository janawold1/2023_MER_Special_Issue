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

# Not cosidering SV type in overlap comparisons.
```
for file in ${data}input/*.txt
    do
    num=$(wc -l ${file} | awk '{print $1}')
    base=$(basename ${file} .txt)
    for i in {0,50,500,1000}
        do
        echo "Running SURVIVOR for ${base} with $i interval..."
        SURVIVOR merge $file ${i} ${num} 0 1 0 50 ${data}outputs/${base}_notype_${i}_distance.vcf
    done &
done

for vcf in ${data}outputs/*_notype_{0,50,500,1000}_distance.vcf
    do
    base=$(basename $vcf _distance.vcf)
    echo "Counting SVs for $base..."
    total=$(bcftools query -f '%SVTYPE\n' ${vcf} | wc -l)
    del=$(bcftools query -i 'SVTYPE=="DEL"' -f '%SVTYPE\n' ${vcf} | wc -l)
    dup=$(bcftools query -i 'SVTYPE=="DUP"' -f '%SVTYPE\n' ${vcf} | wc -l)
    ins=$(bcftools query -i 'SVTYPE=="INS"' -f '%SVTYPE\n' ${vcf} | wc -l)
    inv=$(bcftools query -i 'SVTYPE=="INV"' -f '%SVTYPE\n' ${vcf} | wc -l)
    echo "${base},${total},$del,$dup,$ins,$inv" >> ${data}overlap_counts_min50bp_notype.csv
done
```