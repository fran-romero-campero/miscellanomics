## Script para determinar los genes dianas de un factor de transcripción
## a partir del fichero narrowPeak generado por MaCS2.

## Autor: Francisco J. Romero-Campero - fran@us.es
## Fecha: Octubre 2019

## Instalar chipseeker y paquete de anotación de Arabidopsis thaliana
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Athaliana.BioMart.plantsmart28")

library(ChIPseeker)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28


## Leer fichero de picos
prr5.peaks <- readPeakFile(peakfile = "prr5_peaks.narrowPeak",header=FALSE)

## Definir la región que se considera promotor entorno al TSS
promoter <- getPromoters(TxDb=txdb, 
                         upstream=1000, 
                         downstream=1000)

## Anotación de los picos
prr5.peakAnno <- annotatePeak(peak = prr5.peaks, 
                             tssRegion=c(-1000, 1000),
                             TxDb=txdb)

plotAnnoPie(prr5.peakAnno)
plotDistToTSS(prr5.peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Convertir la anotación a data frame
prr5.annotation <- as.data.frame(prr5.peakAnno)
head(prr5.annotation)

target.genes <- prr5.annotation$geneId[prr5.annotation$annotation == "Promoter"]

write(x = target.genes,file = "prr5_target_genes.txt")

