###############################################################
## Biología Molecular de Sistemas                            ##
## Grado en Bioquímica                                       ##         
## Universidad de Sevilla                                    ##
##                                                           ##
## Unidad 2. Estudios Transcriptómicos Masivos Basados en    ##
##           RNA-seq                                         ##
##                                                           ##
## Análisis de los resultados de HISAT2 y STRINGTIE          ##
## usando los paquetes de R ballgown y limma.                ##  
###############################################################

## El paquete de bioconductor ballgown proporciona las funciones necesarias para 
## realizar un análisis de expresión génica diferencial y visualizar los resultados
## a partir de procesamiento de los datos brutos de secuenciación realizados con 
## hisat2 y stringtie. 

## Para ejecutar con éxito este script es necesario descargar la carpeta samples
## completa a tu ordenador, mover este script a la carpeta samples y fijar el 
## Working Directory To Source File Location. 

## Instalación y carga de los paquetes necesarios. Sólo es necesario instalar los
## paquetes la primera vez que se ejecuta este script en un ordenador el resto de las
## veces bastará cargar los paquetes simplemente. 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ballgown", version = "3.8")
BiocManager::install("genefilter", version = "3.8")

library(ballgown)
library(genefilter)


## Para cargar los datos es necesario crear previamente un fichero tabular
## que contenga como primera columna los nombres de las carpetas donde se guarda
## cada muesra típicamente sample1, sample2, etc ... El resto de columnas
## haran referencia al genotipo, tratamiento y demás caracteriśticas de cada muestra. 
pheno.data <- read.csv("pheno_data.csv")
pheno.data

## La función ballgown se usa para cargar o leer los datos. Es necesario especificar
## el directorio donde se encuentra las muestras. En nuestro caso especificamos .
## para indicar que se encuentran en el actual directorio. 
bg.data <- ballgown(dataDir = ".", samplePattern = "sample", pData=pheno.data)
bg.data
sampleNames(bg.data)

## La función gexpr extrae los niveles de expresión génicos en cada muestra
## medidos como FPKM (Fragments Per Kilobase of exon and Million of mapped reads)
gene.expression <- gexpr(bg.data)
head(gene.expression)
dim(gene.expression)

## Nombramos las columnas con los nombres de nuestras muestras. 
colnames(gene.expression) <- c("col0_1","col0_2","abc_1","abc_2")

## Previsualizamos la similitud entre las réplicas
plot(log2(gene.expression[,1]+1),log2(gene.expression[,2]+1),pch=19,cex=0.7,xlab="col0_1",ylab=substitute(italic("col0_2")),cex.lab=1.25)
plot(log2(gene.expression[,3]+1),log2(gene.expression[,4]+1),pch=19,cex=0.7,xlab="abc_1",ylab=substitute(italic("abc_2")),cex.lab=1.25)

## Construimos un boxplot para comprobar que las distribuciones globales de las
## muestras son similares y comparables.
boxplot(log2(gene.expression + 1),col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5)

## En muchas ocasiones no será necesario realizar ninguna normalización. En cualquier caso
## las instrucciones que siguen realizan una normalización de cuartil superior. Esta es la técnica de
## normalización más comúnmente usada.

## Cálculo de los cuartiles superiores para cada columna.
upper.quantiles <- vector(mode="numeric",length=ncol(gene.expression))

for(i in 1:ncol(gene.expression))
{
  upper.quantiles[i] <- quantile(gene.expression[,i],probs=0.75)
}

## División de cada columna por su cuartil superior. 
for(i in 1:ncol(gene.expression))
{
  gene.expression[,i] <- gene.expression[,i] / upper.quantiles[i]
}

## Transformación logarítmica y visualización de los datos normlaizados.
log.gene.expression <- log2(gene.expression+1)
boxplot(log.gene.expression,col=rainbow(ncol(gene.expression)),ylab="log2(FPKM + 1)",cex.lab=1.5,outline=F)

## Alternativamente la normalización se puede realizar con el paquete normalyzer. 
# ## En este punto ballgown no realiza ninguna normalización de los datos. 
# ## Utilizamos el paquete de R Normalyzer para esta tarea. Para ello es neceario generar
# ## un fichero con un formato específico.
# 
# normalyzer.data <- data.frame(rownames(gene.expression),gene.expression)
# dim(normalyzer.data)
# head(normalyzer.data)
# 
# colnames(normalyzer.data) <- NULL
# rownames(normalyzer.data) <- NULL
# 
# normalyzer.table <- rbind(c(0,rep(1:2,each=2)),                                  
#                     rbind(c("Gene",colnames(gene.expression)),normalyzer.data))
# 
# head(normalyzer.table)
# 
# write.table(normalyzer.table,file="normalyzer_table.tab",col.names = F, quote = F, row.names = F, sep="\t")
# 
# library(Normalyzer)
# library(grid)
# normalyzer(datafile = "normalyzer_table.tab", getjob = "data_normalization")
# 
# normalized.data <- read.table(file="data_normalization/Quantile-normalized.txt", header=T)
# head(normalized.data)
# 
# normalized.data[is.na(normalized.data)] <- 0
# 
# ## Testeamos la normalización
# boxplot(normalized.data[,2:5],col=rep(c("red","blue"),each=2))

## Calculamos la matrix de expresión media. 
col0 <- (log.gene.expression[,"col0_1"] + log.gene.expression[,"col0_2"])/2
abc <- (log.gene.expression[,"abc_1"] + log.gene.expression[,"abc_2"])/2

mean.expression <- matrix(c(col0,abc),ncol=2)
colnames(mean.expression) <- c("col0","abc")
rownames(mean.expression) <- names(col0)
head(mean.expression)

## Previsualizamos el efecto de la mutación en un scatterplot.
plot(col0,abc,pch=19,cex=0.7,xlab="Col0",ylab=substitute(italic("abc")),cex.lab=1.25)

##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

## Especificamos el diseño experimental

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2)))
colnames(experimental.design) <- c("col0","abc")

##A continuación, ajustamos la estimación de los niveles de expresión de cada
##gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
##fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit <- lmFit(log.gene.expression, experimental.design)

##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:

contrast.matrix <- makeContrasts(abc-col0,levels=c("col0","abc"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

nrow(log.gene.expression)

col0.abc <- topTable(contrast.results, number=7507,coef=1,sort.by="logFC")
head(col0.abc)

fold.change <- col0.abc$logFC
genes.ids <- rownames(col0.abc)

activated.genes <- genes.ids[fold.change > 1]
repressed.genes <- genes.ids[fold.change < - 1]

length(activated.genes)
length(repressed.genes)

## Código para desarrollar una función gráfico de barras

gen <- activated.genes[1]

original.data <- 2^log.gene.expression

expr.1 <- unlist(c(original.data[gen, 1:2]))
expr.2 <- unlist(c(original.data[gen, 3:4]))

mean.1 <- mean(expr.1)
mean.2 <- mean(expr.2)

sd.1 <- sd(expr.1)
sd.2 <- sd(expr.2)

means <- c(mean.1, mean.2)
sds <- c(sd.1, sd.2)

arrow.top <- means + sds
arrow.bottom <- means - sds

xpos <- barplot(means,ylim=c(0,1.5*max(arrow.top)),col=rainbow(2))
arrows(xpos, arrow.top, xpos, arrow.bottom,code = 3,angle=90,length=0.2)
