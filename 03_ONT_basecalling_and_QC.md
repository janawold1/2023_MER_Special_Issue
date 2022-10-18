# ONT basecalling, read quality and trimming
This markdown file is a summary of scripts used to basecall and filter ONT data for alignment.

Raw MinION fast5 files were basecalled using guppy v6.0.1 under the super high accuracy model. This particular step was conducted on a separate VM capable of GPU basecalling.
```
config=dna_r9.4.1_450bps_sup.cfg
guppy_basecaller -i ${data} \
    -s ${output} \
    --config ${config} \
    --verbose_logs \
    --compress_fastq \
    --device cuda:all:100%\
    --recursive \
    --calib_detect \
    --detect_mid_strand_adapter
```
A summary of raw sequencing outputs was assessed using [NanoPlot](https://github.com/wdecoster/NanoPlot).
```

```
After basecalling, reads were then moved to the VM intended for all downstream analysis. Global variables were defined as below:
```
data=/kakapo-data/ONT/
lambda=${data}DNA_CS.fasta
source ~/anaconda3/etc/profile.d/conda.sh
```
## Read trimming
Adapters, including those occuring midstrand, were trimmed using [Porechop](https://github.com/rrwick/Porechop). Lambda DNA was removed using [NanoLyse](https://github.com/wdecoster/nanolyse) and [NanoFilt](https://github.com/wdecoster/nanofilt) were used to ensure a minimum q score of 10 and filter reads to four different minimum read length thresholds (3, 4, 5 and 10kb) prior to alignment.
```
conda activate nanofilt
for fq in ${data}fastq/*.fastq.gz
    do
    base=$(basename ${fq} .fastq.gz)
    echo "Trimming reads for ${base}..."
    porechop -i ${fq} -o ${data}porechop/${base}_porechop.fastq.gz --discard_middle
    for length in {3,4,5,10}
        do
        gunzip -c ${data}porechop/${base}_porechop.fastq.gz | \
        NanoLyse -r ${lambda} | \
        NanoFilt -q 10 -l ${length}000 > ${data}trimmed/${base}_q10_${length}kbtrim.fastq
        gzip ${data}trimmed/${base}_q10_${length}kbtrim.fastq
    done &
done
conda deactivate nanofilt
```
## Read statistics
Summary statistics for each of these filtering trials were assessed using [NanoPlot](https://github.com/wdecoster/NanoPlot).
```
conda activate nanoplot
for indiv in A B C D E F G H I J K L
    do
    for length in {3,4,5,10}
        do
        echo "Running NanoPlot for ${indiv}..."
         NanoPlot -t 24 -o ${data}NanoPlot/$indiv \
            -f pdf --N50 \
            --title "${indiv} q10 ${length}kb trim" \
            -p ${indiv}_q10_${length}kbtrim \
            --fastq trimmed/${indiv}_q10_${length}kbtrim.fastq.gz
    done &
done
conda deactivate nanoplot
```