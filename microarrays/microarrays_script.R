###############################################################
## Biología Molecular de Sistemas                            ##
## Grado en Bioquímica                                       ##         
## Universidad de Sevilla                                    ##
##                                                           ##
## Unidad 1. Estudios Transcriptómicos Masivos: Análisis de  ##
##           Microarrays                                     ##
##                                                           ##
###############################################################

## El paquete affy proporciona las funciones básicas para manipular
## datos de microarray de affymetrix
library(affy)

## Los datos de microarray que vamos a analizar corresponden al estudio 
## Long et al. (2010) The bHLH transcription factor POPEYE regulates response 
## to iron deficiency in Arabidopsis roots. Plant Cell 22(7):2219-36.
## Estos datos se encuentran disponibles de forma pública en la base de datos
## GEO identificados con el número de acceso GSE21582.
## Fijamos el espacio de trabajo a la carpeta que contenga los datos de 
## microarrays y utilizamos la función ReadAffy para cargar los datos de 
## de microarrays brutos contenidos en la correspondiente carpeta
help(ReadAffy)
microarray.raw.data <- ReadAffy(verbose=TRUE)

## Si evaluamos la correspondiente variable obtenemos información sobre el 
## tamaño del microarray, el diseño de la placa, el número de muestras y 
## el número de genes.
microarray.raw.data

## Obtenemos el diseño de placa con la función cdfName
cdfName(microarray.raw.data)

## Los paquetes simpleaffy y affyPLM contienen funciones para el análisis de 
## la calidad de microarrays
library(simpleaffy)
library(affyPLM)

## La función image nos permite visualizar la luminiscencia de la placa
## de microarray
for(i in 1:8)
{
  image(microarray.raw.data[,i],col=rainbow(100))
  readline(prompt = "Pause. Press <Enter> to continue...")
}

## La función qc aplicada sobre los datos crudos nos permite realizar un 
## análisis de la calidad. Un resumen gráfico del análisis de la calidad 
## se puede representar con la función plot
quality.analysis <- qc(microarray.raw.data)
plot(quality.analysis)

## Antes de realizar ningún preprocesamiento de los datos comprabamos si
## son comparables utilizando como descriptores globales de la distribución
## de los niveles de expresión boxplot (cajas y bigotes) e histogramas.
boxplot(microarray.raw.data,col=rainbow(8),las=2,ylab="Fluorescence")
hist(microarray.raw.data,col=rainbow(8))

## Para lograr que todos los microarrays sean comparables utilizamos el 
## Robust Multiarray Analysis que realiza corrección de la fluorescencia
## del fondo, normalización, sumación y estimación de los niveles
## de expresión en log2. 
microarray.processed.data <- rma(microarray.raw.data)

## Una vez procesados los datos comprobamos si son comparables.
boxplot(microarray.processed.data,col=rainbow(8),las=2,ylab="Fluorescence")
hist(microarray.processed.data,col=rainbow(8))

## La estimación de los niveles de expresión realizada por RMA puede 
## extraerse como una matriz usando la función exprs
expression.level <- exprs(microarray.processed.data)
head(expression.level)

## Obtenemos los nombres de las sondas y nombramos apropiadamente las
## columnas con el nombre de las muestras
probe.names <- rownames(expression.level)
probe.names
sampleID <- c("WT_with_Fe_1","WT_with_Fe_2","WT_no_Fe_1","WT_no_Fe_2","pye_with_Fe_1","pye_with_Fe_2","pye_no_Fe_1","pye_no_Fe_2")
colnames(expression.level) <- sampleID
head(expression.level)

## Calculamos los valores medios de expresión para cada genotipo/condición 
## sumando las correspondientes columnas y dividiendo por el número de réplicas
wt.with.fe <- (expression.level[,"WT_with_Fe_1"] + expression.level[,"WT_with_Fe_2"])/2
wt.no.fe <- (expression.level[,"WT_no_Fe_1"] + expression.level[,"WT_no_Fe_2"])/2
pye.with.fe <- (expression.level[,"pye_with_Fe_1"] + expression.level[,"pye_with_Fe_2"])/2
pye.no.fe <- (expression.level[,"pye_no_Fe_1"] + expression.level[,"pye_no_Fe_2"])/2

## Creamos una matriz que contenga por columna la expresión media para cada 
## condición o genotipo. Nombramos las filas con el nombre de las sondas (genes) 
## y la columnas con la condición o genotipo. 
mean.expression <- matrix(c(wt.with.fe,wt.no.fe,pye.with.fe,pye.no.fe),ncol=4)
conditions.id <- c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")
rownames(mean.expression) <- names(wt.with.fe)
colnames(mean.expression) <- conditions.id
head(mean.expression)

## Scatterplots o gráficos de dispersión para la comparción de distintos 
## genotipos/condiciones. Este tipo de gráficos nos permite obtener una visión 
## global de la comparación entre genotipos/condiciones.
plot(wt.with.fe,wt.no.fe,xlab="WT with Fe",ylab="WT no Fe")
plot(pye.with.fe,pye.no.fe,xlab="pye with Fe",ylab="pye no Fe")

plot(wt.with.fe,pye.with.fe,xlab="WT with Fe",ylab="pye with Fe")
plot(wt.no.fe,pye.no.fe,xlab="WT no Fe",ylab="pye no Fe")

## El paquete limma proporciona las funciones necesarias para 
## determinar los genes expresados de forma diferencial (DEGs).
library(limma)

## Generamos una matriz que represente el diseño experimental, es decir,
## marcamos las muestras que constituyen réplicas biológicas de un
## mismo genotipo o condición con el mismo número.
experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3,4,4)))
colnames(experimental.design) <- c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")

## Ajustamos la estimación de los niveles de expresión de cada gen a un 
## modelo lineal teniendo en cuenta el diseño experimental (réplicas)
linear.fit <- lmFit(expression.level, experimental.design)

## Para determinar los genes expresados de forma diferencial debemos 
## especificar los contrastes (comparaciones entre condiciones) que se 
## van a considerar. En nuestro ejemplo buscamos determinar los genes 
## que ven su expresión afectada por la deficiencia de hierro en el 
## genotipo WT y pye. Además queremos determinar el efecto del gen 
## POPEYE. Por lo tanto, realizaremos todas las siguientes comparaciones:
##
##      WT no Fe / WT with Fe
##      pye no Fe / pye with Fe
##      WT with Fe / pye with Fe
##      WT no Fe / pye no Fe
##      Pericycle with Fe / Pericycle no Fe

## Para especificar los constrastes a realizar utilizamos la función 
## makeContrasts que recibe como entrada los contrastes a realizar 
## separados por comas y especificados con los nombres de las dos 
## condiciones correspondientes separadas por un guión -. También 
## recibe el argumento levels, un vector con el nombre de las condiciones:
contrast.matrix <- makeContrasts(WT_no_Fe-WT_with_Fe,pye_no_Fe-pye_with_Fe,WT_with_Fe-pye_with_Fe,WT_no_Fe-pye_no_Fe,levels=c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe"))

## Calculamos el fold-change y los p-valores correspondientes para cada
## gen en cada uno de los constrastes especificados utilizando las 
## funciones constrasts.fit y eBayes
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

## Con la función topTable podemos obtener un marco de datos con 
## información sobre la expresión diferencial de los genes. Esta función
## recibe como entrada la salida de eBayes, el número de genes a mostrar 
## (argumento number) y el identificador de la comparación a tener en 
## cuenta (argumento coef). El parametro sort.by nos permite ordenar las filas del
## marco de datos según el fold-change o el p-valor. 
WT.with.no.Fe <- topTable(contrast.results, number=22810,coef=1,sort.by="logFC")
head(WT.with.no.Fe)

## Para la selección de DEGs extraemos el fold change, p valores e identificadores de 
## cada gen.
fold.change.WT.with.no.Fe <- WT.with.no.Fe[["logFC"]]
p.value.WT.with.no.Fe <- WT.with.no.Fe[["adj.P.Val"]]
genes.ids <- rownames(WT.with.no.Fe)

# dependiendo de la versión de R quizás debas utilizar esta instrucción para 
## obtener los identificadores de los genes
## genes.ids <- WT.with.no.Fe[["ID"]]  

## Para la selección de DEGs seguimos el criterio del fold-change fijando 
## un umbral de 1 que corresponde a genes que se expresan más del doble o menos 
## de la mitad. 
activated.genes.WT.with.no.Fe.1 <- genes.ids[fold.change.WT.with.no.Fe > 1]
repressed.genes.WT.with.no.Fe.1 <- genes.ids[((fold.change.WT.with.no.Fe < - 1))]

length(activated.genes.WT.with.no.Fe.1)
length(repressed.genes.WT.with.no.Fe.1)

## La librería annaffy nos proporciona funciones para generar listas de
## DEGs en formato html o txt con la anotación disponible sobre los distintos 
## genes. 
library(annaffy)

## La funicón affTableAnn nos permite generar una tabla con DEGs y su anotación. 
## Las funciones saveHTML y saveText escriben la tabla anterior en formato HTML y txt.
activated.genes.WT.with.no.Fe.1.table <- aafTableAnn(activated.genes.WT.with.no.Fe.1, "ath1121501.db", aaf.handler())
saveHTML(activated.genes.WT.with.no.Fe.1.table, file="activated_genes_WT_with_no_Fe_1.html")
saveText(activated.genes.WT.with.no.Fe.1.table, file="activated_genes_WT_with_no_Fe_1.txt")

repressed.genes.WT.with.no.Fe.1.table <- aafTableAnn(repressed.genes.WT.with.no.Fe.1, "ath1121501.db", aaf.handler())
saveHTML(repressed.genes.WT.with.no.Fe.1.table, file="repressed_genes_WT_with_no_Fe_1.html")
saveText(repressed.genes.WT.with.no.Fe.1.table, file="repressed_genes_WT_with_no_Fe_1.txt")

## Para seleccionar los DEGs en el resto de comparaciones se realizarían de forma
## repetitiva las mismas instrucciones pero con datos diferentes por lo tanto es 
## apropiado definir una siguiente función que recibe como entrada el resultado de la
## función eBayes (contrast.results), el número de genes total (gene.number), el  
## identificador numérico de la comparación a considerar (comparison.number) y el umbral
## de corte en el fold change (fold.change.threshold).

## Uno de los gráficos más básicos que muestran la comparación entre dos 
## transcriptomas junto con los genes seleccionados como activados y reprimidos 
## son los gráficos de dispersión o scatterplots. Estos gráficos se usan
## principalmente cuando se determinan los DEGs según el criterio del fold-change
plot(wt.with.fe,wt.no.fe,pch=19,cex=0.5,col="grey",xlab="WT with Fe",ylab="WT no Fe")
points(wt.with.fe[activated.genes.WT.with.no.Fe.1],wt.no.fe[activated.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="red")
points(wt.with.fe[repressed.genes.WT.with.no.Fe.1],wt.no.fe[repressed.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="blue")
text(wt.with.fe["254550_at"]+0.3,wt.no.fe["254550_at"]+0.3,"IRT1", col="black", cex=0.7)
text(wt.with.fe["252427_at"]+0.3,wt.no.fe["252427_at"]+0.3,"PYE", col="black", cex=0.7)
text(wt.with.fe["257062_at"]+0.3,wt.no.fe["257062_at"]+0.3,"BTS", col="black", cex=0.7)
text(wt.with.fe["251109_at"]+0.3,wt.no.fe["251109_at"]+0.3,"FER1", col="black", cex=0.7)

plot(mean.expression[,"pye_with_Fe"],mean.expression[,"pye_no_Fe"],pch=19,cex=0.5,col="grey",xlab="pye with Fe",ylab="pye no Fe")
points(mean.expression[diff.genes.2[["activated.genes"]],"pye_with_Fe"],mean.expression[diff.genes.2[["activated.genes"]],"pye_no_Fe"],pch=19,cex=0.5,col="red")
points(mean.expression[diff.genes.2[["repressed.genes"]],"pye_with_Fe"],mean.expression[diff.genes.2[["repressed.genes"]],"pye_no_Fe"],pch=19,cex=0.5,col="blue")

plot(mean.expression[,"WT_with_Fe"],mean.expression[,"pye_with_Fe"],pch=19,cex=0.5,col="grey",xlab="WT with Fe",ylab="pye with Fe")
points(mean.expression[diff.genes.3[["activated.genes"]],"WT_with_Fe"],mean.expression[diff.genes.3[["activated.genes"]],"pye_with_Fe"],pch=19,cex=0.5,col="blue")
points(mean.expression[diff.genes.3[["repressed.genes"]],"WT_with_Fe"],mean.expression[diff.genes.3[["repressed.genes"]],"pye_with_Fe"],pch=19,cex=0.5,col="red")

  
## Selección de DEGs utilizando el criterio de inferencia estadística. Fijamos como
## umbral de significancia un p-valor menor a 0.01. Los genes activados son aquellos
## cuyo fold change es mayor que 0 y su p-valor menor 0.01. De forma semejante 
## determinamos los genes reprimidos.
activated.genes.WT.with.no.Fe.2 <- genes.ids[((fold.change.WT.with.no.Fe > 0) & p.value.WT.with.no.Fe < 0.01)]
repressed.genes.WT.with.no.Fe.2 <- genes.ids[((fold.change.WT.with.no.Fe < 0) & p.value.WT.with.no.Fe < 0.01)]

length(activated.genes.WT.with.no.Fe.2)
length(repressed.genes.WT.with.no.Fe.2)

## Definir una función para la selección de DEGs según un criterio de inferencia 
## estadística. 

diff.genes.2 <- DEGs.p.value(contrast.results,gene.number=22819,comparison.number=2,p.value.threshold=0.01)
length(diff.genes.2[["activated.genes"]])
length(diff.genes.2[["repressed.genes"]])

diff.genes.3 <- DEGs.p.value(contrast.results,gene.number=22819,comparison.number=3,p.value.threshold=0.01)
length(diff.genes.3[["activated.genes"]])
length(diff.genes.3[["repressed.genes"]])

diff.genes.4 <- DEGs.p.value(contrast.results,gene.number=22819,comparison.number=4,p.value.threshold=0.01)
length(diff.genes.4[["activated.genes"]])
length(diff.genes.4[["repressed.genes"]])

## DEG basados en inferencia estadística y fold-change
activated.genes.WT.with.no.Fe.3 <- genes.ids[((fold.change.WT.with.no.Fe > 1) & p.value.WT.with.no.Fe < 0.01)]
repressed.genes.WT.with.no.Fe.3 <- genes.ids[((fold.change.WT.with.no.Fe < -1) & p.value.WT.with.no.Fe < 0.01)]

length(activated.genes.WT.with.no.Fe.3)
length(repressed.genes.WT.with.no.Fe.3)

## Definir una función que combina los criterios de selección de DEGs basados en
## en fold change y en la inferencia estadística. 
diff.genes.2 <- DEGs.p.value.fold.change(contrast.results,gene.number=22810,comparison.number=2,fold.change.threshold=1,p.value.threshold=0.01)
length(diff.genes.2[["activated.genes"]])
length(diff.genes.2[["repressed.genes"]])

diff.genes.3 <- DEGs.p.value.fold.change(contrast.results,gene.number=22810,comparison.number=3,fold.change.threshold=1,p.value.threshold=0.01)
length(diff.genes.3[["activated.genes"]])
length(diff.genes.3[["repressed.genes"]])

diff.genes.4 <- DEGs.p.value.fold.change(contrast.results,gene.number=22810,comparison.number=4,fold.change.threshold=1,p.value.threshold=0.01)
length(diff.genes.4[["activated.genes"]])
length(diff.genes.4[["repressed.genes"]])

## Los volcano plots nos permiten representar en un único gráfico el fold change
## y la significancia estadística. Estos gráficos se utilizan para representar los
## DEGs cuando se usa el método que combina fold-change e inferencia estadística
WT.with.no.Fe <- topTable(contrast.results, number=22810,coef=1,sort.by="logFC")
head(WT.with.no.Fe)

fold.change.WT.with.no.Fe <- WT.with.no.Fe[["logFC"]]
names(fold.change.WT.with.no.Fe) <- rownames(WT.with.no.Fe)
#names(fold.change.WT.with.no.Fe) <- WT.with.no.Fe[[1]]
log.p.value.WT.with.no.Fe <- -log10(WT.with.no.Fe[["adj.P.Val"]])
names(log.p.value.WT.with.no.Fe)  <- rownames(WT.with.no.Fe)
#names(log.p.value.WT.with.no.Fe)  <- WT.with.no.Fe[[1]]

plot(fold.change.WT.with.no.Fe,log.p.value.WT.with.no.Fe,pch=19,cex=0.5,col="grey",ylab="-log10(p value)",xlab="log2 fold change")
points(fold.change.WT.with.no.Fe[activated.genes.WT.with.no.Fe.3],log.p.value.WT.with.no.Fe[activated.genes.WT.with.no.Fe.3],pch=19,cex=0.5,col="red")
points(fold.change.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.3],log.p.value.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.3],pch=19,cex=0.5,col="blue")

## Los mapas de calor o heatmap realizados sobre datos normalizados nos permiten 
## visualizar patrones de expresión si tenemos el suficiente número de condiciones/genotipos.
diff.genes.1 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=1,fold.change.threshold=1)
diff.genes.2 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=2,fold.change.threshold=1)
diff.genes.3 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=3,fold.change.threshold=1)
diff.genes.4 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=4,fold.change.threshold=1)

## Recopilamos todos los DEGs
complete.DEGs <- c(diff.genes.1[["activated.genes"]],diff.genes.1[["repressed.genes"]],diff.genes.2[["activated.genes"]],diff.genes.2[["repressed.genes"]],diff.genes.3[["activated.genes"]],diff.genes.3[["repressed.genes"]],diff.genes.4[["activated.genes"]],diff.genes.4[["repressed.genes"]])
## Eliminamos repeticiones en el conjunto de todos los DEGs.
complete.DEGs <- unique(complete.DEGs)
length(complete.DEGs)

## Extraemos los niveles de expresión de todos los DEGs en las condiciones/genotipos
## de interés. 
DEG.expression <- mean.expression[complete.DEGs,c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")]

## Normalizamos los datos teniendo en cuenta que scale normaliza por columnas y nosotros
## queremos hacerlo por filas. Calculamos con t la traspuesta de nuestra matriz. 
normalized.DEG.expression <- t(scale(t(DEG.expression)))

## Con la función heatmap.2 de gplots construímos el correspondiente heatmap.
library(gplots)
help(heatmap.2)
heatmap.2(normalized.DEG.expression,Colv=FALSE,dendrogram="row",labRow=c(""),density.info="none",trace="none",col=heat.colors(100)[100:1],margins = c(8,8),cexCol=1.2)


