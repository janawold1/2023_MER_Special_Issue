#!/bin/bash/
#####################################################################################
#Script for running the MANTA software package. Intended for structural
#variant calling. MANTA is aprogram for the detection of structural variants
#from paired-end sequence data. Because running all individuals together is taxing on 
#MANTA, this program was run in 7 batches for females and 9 batches for males.
#####################################################################################
##Setting fixed variables
fref=~/Kakapo/References/kakapo_full_ref.fa
mref=/~Kakapo/References/kakapo_noW.fa
work=~/Kakapo/alignments_female/bam/sorted/
out=~/Kakapo/manta/
MANTA_INSTALL_PATH=~/software/manta-1.6.0.centos6_x86_64/
#To create inputs
#ls -lh alignments_female/bam/sorted/ | awk '{print $9}' | grep -v .bai | grep bam | sort | sed -e 's%^%--bam ~/Kakapo/alignments_female/bam/sorted/%' | awk '{print $0, "\"}'
#ls -lh alignments_male/bam/sorted/ | awk '{print $9}' | grep -v .bai | grep bam | sort | sed -e 's%^%--bam ~/Kakapo/alignments_male/bam/sorted/%' | awk '{print $0, "\"}'

##First need to create the configuration file
${MANTA_INSTALL_PATH}/bin/configManta.py \
--bam ~/Kakapo/alignments_female/bam/sorted/N.sorted.bam \
--referenceFasta ${fref} \
--runDir ${out}female/

${MANTA_INSTALL_PATH}/bin/configManta.py \
--bam ~/Kakapo/alignments_female/bam/sorted/N.sorted.bam \
--referenceFasta ${mref} \
--runDir ${out}male/

echo "Running MANTA workflow for females..."
~/Kakapo/manta/female/runWorkflow.py

echo "Running MANTA workflow for males..."
~/Kakapo/manta/male/runWorkflow.py