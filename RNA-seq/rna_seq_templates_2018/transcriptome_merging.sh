## This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0$
## Author: Francisco J. Romero-Campero
## Email: fran@us.es

#$ -S /bin/bash
#$ -N merge
#$ -V
#$ -cwd
#$ -j yes
#$ -o merge


## Accessing results folder
cd /home/<grupo>/<exp>/results

## Merging sample transcriptomes
stringtie --merge -G ../annotation/annotation.gtf -o stringtie_merged.gtf merge_list.txt

## Comparing our assembly with the reference
gffcompare -r ../annotation/annotation.gtf -G -o comparison stringtie_merged.gtf
