## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International$
## Author: Francisco J. Romero-Campero
## Email: fran@us.es

#$ -S /bin/bash
#$ -N index
#$ -V
#$ -cwd
#$ -j yes
#$ -o index

## Dowloading reference genome and annotation
cd /home/<grupo>/<exp>/genome 
wget -O genome.fa.gz <URL_GENOME>
gunzip genome.fa.gz

cd /home/<grupo>/<exp>/annotation
wget -O annotation.gtf.gz <URL_ANNOTATION>
gunzip annotation.gtf.gz

## Building reference genome index
cd /home/<grupo>/<exp>/genome
extract_splice_sites.py ../annotation/annotation.gtf > splices_sites.ss
extract_exons.py ../annotation/annotation.gtf > exons.exon
hisat2-build --ss splices_sites.ss --exon exons.exon genome.fa index
