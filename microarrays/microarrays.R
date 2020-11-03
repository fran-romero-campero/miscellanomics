# Unidad 1. Estudios Transcriptómicos Masivos: Análisis de Microarrays  
# Bilogía Molecular de Sistemas, 3º Grado en Bioquímica
# Francisco J. Romero-Campero - email: fran@us.es
# YouTube: https://www.youtube.com/channel/UCRBDDVQHHisLcZtLPlYvmow
# Twitter: https://twitter.com/fran_rom_cam
# IG: @greennetworks https://www.instagram.com/greennetworks/
# Dpto Ciencias de la Computación e Inteligencia Artificial
# Instituto de Bioquímica Vegetal y Fotosíntesis
# Universidad de Sevilla 
  
# Lectura de los datos brutos.
  
#El paquete **affy** proporciona las funciones básicas para manipular 
##datos de microarrays de **affymetrix**.

library(affy)

##Los datos de microarray que vamos a analizar corresponden al estudio 
##Long et al. (2010) The bHLH transcription factor POPEYE regulates response 
##to iron deficiency in Arabidopsis roots. Plant Cell 22(7):2219-36. 
##Estos datos se encuentran disponibles de forma pública en la base de datos 
##Gene Expression Omnibus (GEO) identificados con el número de acceso GSE21582.

##Fijamos el espacio de trabajo a la carpeta que contenga los datos de 
##microarrays. Siempre es aconsejable guardar en la misma carpeta el script 
##de análisis y los datos analizados. Si se sigue este consejo basta fijar 
##el espacio de trabajo a la localización donde se encuentra el script o 
##source file. Para ellos usar el menu **Session** y la opción **Set Working Directory**
##seguida de **To Source File Location**. 

##Utilizamos la función **ReadAffy** para cargar o leer los datos brutos de 
##microarrays  (fluorescencia sin procesar) contenidos en la correspondiente carpeta. 
##Los datos brutos deben estar contenidos en ficheros con **formato CEL**. 
##Se recomienda fijar el argumento **verbose** a TRUE para que muestra información 
##por pantalla sobre la lectura de los datos.


microarray.raw.data <- ReadAffy(verbose=TRUE)
microarray.raw.data

##Si evaluamos la variable que almacena los datos brutos obtenemos información 
##sobre el tamaño del microarray, el diseño de la placa, el número de muestras 
##y el número de transcritos (o sondas) que contiene la placa correspondiente. 

microarray.raw.data

##Si sólo estamos interesados en obtener el diseño de placa podemos usar la 
##función **cdfName**.


cdfName(microarray.raw.data)

## Control de la Calidad.

##Los paquetes **simpleaffy** y **affyPLM** contienen funciones para el análisis 
##de la calidad de microarrays.

library(simpleaffy)
library(affyPLM)

##La función **image** nos permite visualizar una reproducción de la foto 
##tomada por el scanner de la fluorescencia de la placa de microarray. Usando 
##esta visualización podemos determinar si se han producido daños físicos 
##durante el manejo de la placa de microarrays. 

image(microarray.raw.data[,1],col=rainbow(100))
image(microarray.raw.data[,2],col=rainbow(100))

##La función **qc** aplicada sobre los datos crudos nos permite realizar un 
##análisis de la calidad. Un resumen gráfico del análisis de la calidad se 
##puede representar con la función plot.

quality.analysis <- qc(microarray.raw.data)
plot(quality.analysis)

##Cada fila del gráfico QC representa una muestra. El primer valor númerico 
#se llama *porcentaje de deteccion* y corresponde al porcentaje de sondas para 
##el cual se ha detectado fluorescencia. Se espera que estos valores sean 
##similares en todas las muestras. El según valor númerico hace referencia a 
##la *fluorescencia del fondo* medida en posiciones del microarray donde no se 
##encuentra ninguna sonda de nucleótidos. Esta fluorescencia de fondo mide la 
##hibridación inespecífica de cDNA a la superficie del microarray. Típicamente 
##la fluorescencia del fondo aparece marcada como problemática ya que presenta 
##valores dispares en las distintas muestras. Sin embargo, este valor puede 
##considerarse como un "blanco" que es fácilmente corregible durante el 
##preprocesamiento de los datos.  Los símbolos de círculo y triángulo 
##representan la proporición de fluorescencia entre los extremos 3' y 5' de las
##sondas control para la actina y el GAPDH. Debido a los procesos de degradación 
##del RNA y de síntesis del cDNA a partir del RNA siempre se espera que el 
##extremo 3' de los transcritos esté sobrerepresentado con respecto al extremo 
##5'. Sin embargo, si se han dado procesos masivos de degradación del RNA y de 
##incompletitud durante la síntesis del cDNA el extremo 3' estará muy altamente
##sobrerepresando con respectoa l extremo 5' lo que indicará problemas en la
##calidad. 

## Preprocesamiento y Estimación de los Niveles de Expresión génicos. 

##En todo expermiento existen dos fuentes de variabilidad. Por una parte tendremos 
##la variabilidad biológica que estamos interesados en estudiar. Por otra parte 
##tendremos a variabilidad experimental producida por factores técnicos. La 
##variablidad experimental se suele denominar ruido y no es deseada. 
##Típicamente en todo análisis de datos se realiza un preprocesamiento de los 
##datos brutos que busca eliminar el ruido o variabilidad experimental mientras 
##se mantiene la variabilidad biológica.

##Antes de realizar ningún preprocesamiento de los datos comprabamos si son 
##comparables utilizando como descriptores globales de la distribución de los 
##niveles de expresión boxplot (cajas y bigotes) e histogramas.

boxplot(microarray.raw.data,col=rainbow(14),las=2,ylab="Fluorescence")
hist(microarray.raw.data,col=rainbow(14))

##Para lograr que todos los microarrays sean comparables utilizamos el 
##**Robust Multiarray Analysis** que realiza corrección de la fluorescencia 
##del fondo, normalización, sumación y estimación de los niveles de expresión 
##en **log2**. Este algoritmo realiza el correspondiente preprocesamiento de 
##los datos brutos de microarrays. 

microarray.processed.data <- rma(microarray.raw.data)

##Una vez preprocesados los datos comprobamos si son comparables.

boxplot(microarray.processed.data,col=rainbow(8),las=2,ylab="Fluorescence A.U.")
hist(microarray.processed.data,col=rainbow(8))

##La estimación de los niveles de expresión realizada por RMA puede extraerse 
##como una matriz usando la función **exprs**. Las columnas de esta matriz 
##corresponden a los ficheros de datos brutos .CEL. Por lo tanto, las columnas 
##representan muestras y las filas sondas o transcritos.

expression.level <- exprs(microarray.processed.data)
head(expression.level)
dim(expression.level)

##Por defecto, el nombre de las columnas coincide con el nombre del 
##correspondiente fichero CEL. Estos nombres no son informativos y conducen 
##fácilmente a error. Por lo tanto, renombramos apropiadamente las columnas 
##con el nombre de las muestras. 

sampleID <- c("WT_with_Fe_1","WT_with_Fe_2","WT_no_Fe_1","WT_no_Fe_2",
              "pye_with_Fe_1","pye_with_Fe_2","pye_no_Fe_1","pye_no_Fe_2")
colnames(expression.level) <- sampleID
head(expression.level)

##Calculamos los valores medios de expresión para cada genotipo/condición
##sumando las correspondientes columnas y dividiendo por el número de 
##réplicas.

wt.with.fe <- (expression.level[,"WT_with_Fe_1"] + 
                 expression.level[,"WT_with_Fe_2"])/2

wt.no.fe <- (expression.level[,"WT_no_Fe_1"] + 
               expression.level[,"WT_no_Fe_2"])/2

pye.with.fe <- (expression.level[,"pye_with_Fe_1"] + 
                  expression.level[,"pye_with_Fe_2"])/2

pye.no.fe <- (expression.level[,"pye_no_Fe_1"] + 
                expression.level[,"pye_no_Fe_2"])/2

##Creamos una matriz que contenga por columna la expresión media para cada 
##condición o genotipo. Nombramos las filas con el nombre de las sondas 
##(transcritos) y la columnas con la condición o genotipo. 

mean.expression <- matrix(c(wt.with.fe,wt.no.fe,
                            pye.with.fe,pye.no.fe),ncol=4)
conditions.id <- c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")
rownames(mean.expression) <- names(wt.with.fe)
colnames(mean.expression) <- conditions.id
head(mean.expression)

##Para realizar una previsualización comparativa entre los transcriptomas de 
##dos genotipos/condiciones diferentes se suelen usar gráficos de dispersión 
##o scatterplots. Este tipo de gráficos nos permite obtener una visión global
##de la comparación entre genotipos/condiciones. Cada punto representa un 
##transcrito, la coordenada x representa el nivel expresión en el primer 
##genotipo/condición mientras que la coordenada y representa el nivel de 
##expresión en el segundo genotipo/condición. Por lo tanto, genes cuyos 
##transcritos estén representados por puntos localizados en la diagonal no se 
##expresan de forma diferencial. Genes cuyos transcritos estén representados por 
##puntos localizados en la parte superior a la diagonal estarán activados en 
##la segunda condición/genotipo con respecto a la primera. Genes cuyos 
##transcritos estén representados por puntos localizados en la parte inferior 
##a la diagonal estarán reprimidos en la segunda condición/genotipo con 
##respecto a la primera.

##Gráficos de dispersión o scatterplots que muestren una gran dispersión 
##indica efectos notables en el transcriptoma. Mientras que poca dispersión 
##indica efectos leves. Por lo tanto, en los casos donde se observe mucha 
##dispersión se suelen fijar criterios restrictivos para la determinación de 
##genes expresados de forma diferencial.

plot(wt.with.fe,wt.no.fe,xlab="WT with Fe",ylab="WT no Fe",pch=19,cex=0.5)
  plot(pye.with.fe,pye.no.fe,xlab="pye with Fe",ylab="pye no Fe",pch=19,cex=0.5)
plot(wt.with.fe,pye.with.fe,xlab="WT with Fe",ylab="pye with Fe",pch=19,cex=0.5)
plot(wt.no.fe,pye.no.fe,xlab="WT no Fe",ylab="pye no Fe",pch=19,cex=0.5)

## Selección de Genes Expresados de Forma Diferencial. 

##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 

library(limma)

##Generamos una matriz que represente el diseño experimental. 
##Cada condición/genotipo se identifica con un número entero. En nuestro caso 
##tenemos seis condiciones/genotipos: tipo silvestre con hierro (1), tipo 
##silvestre sin hierro (2), popeye con hierro (3) y popeye sin hierro (4). 
##Seguidamente, se marca cada muestra con el número que identifica cada 
##condición/genotipo. De esta forma, cada identificador se repite el número 
##de réplicas de la condición/genotipo correspondiente (no necesariamente de 
##forma consecutiva). Usamos la función **model.matrix** con la sintáxis 
##siguiente para especificar el diseño experimental. Las columnas de esta 
##matrix se nombran con los nombres de cada condición, WT_with_Fe, WT_no_Fe, 
##pye_with_Fe y pye_no_Fe. Muy importante, a partir de este punto sólo 
##podremos referirnos a cada condición/genotipo usando exactamente los nombres 
##anteriores. 

experimental.design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3,4,4)))
colnames(experimental.design) <- c("WT_with_Fe","WT_no_Fe",
                                   "pye_with_Fe","pye_no_Fe")

##A continuación, ajustamos la estimación de los niveles de expresión de cada
##gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
##fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit <- lmFit(expression.level, experimental.design)

##Para determinar los genes expresados de forma diferencial debemos especificar 
##los contrastes (comparaciones entre condiciones) que se van a considerar.
##Típicamente se suelen considerar comparaciones entre condiciones donde sólo 
##se ha dado un cambio genotípico o de condición para poder estrablecer de 
##forma inequívoca una relación causal. En nuestro ejemplo buscamos determinar
##los genes que ven su expresión afectada por la deficiencia de hierro en el
##genotipo WT y pye. Además queremos determinar el efecto del gen POPEYE. Por
##lo tanto, realizaremos todas las siguientes comparaciones:
  
##  * WT no Fe / WT with Fe
## * pye no Fe / pye with Fe
## * WT with Fe / pye with Fe
## * WT no Fe / pye no Fe

##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:
  
contrast.matrix <- makeContrasts(WT_no_Fe-WT_with_Fe,
                                 pye_no_Fe-pye_with_Fe,
                                 pye_with_Fe-WT_with_Fe,
                                 pye_no_Fe-WT_no_Fe,
                                 levels=c("WT_with_Fe","WT_no_Fe",
                                          "pye_with_Fe","pye_no_Fe"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

##Con la función *topTable* podemos obtener un marco de datos con información
##sobre la expresión diferencial de los genes. Esta función recibe como entrada
##la salida de eBayes, el número de genes a mostrar (argumento number) y el 
##identificador de la comparación a tener en cuenta (argumento coef). 
##El parametro sort.by nos permite ordenar las filas del marco de datos según 
##el fold-change o el p-valor. 
nrow(expression.level)
WT.with.no.Fe <- topTable(contrast.results, number=22810,coef=1,sort.by="logFC")
head(WT.with.no.Fe)

pye.with.no.Fe <- topTable(contrast.results, number=22810,coef=2,sort.by="logFC")
head(pye.with.no.Fe)


##Cuando se comparan los transcriptomas de dos genotipos diferentes o de un 
##mismo genotipo bajo distintas condiciones existen diversos métodos para 
##determinar genes expresados de forma diferencial o differentially expressed
##genes (DEGs) en inglés:
  
### Método basado en el fold-change (factor de proporcionalidad): 
  
##Se fija un umbral para el fold-change típicamente 2, 4 u 8 que en log2 
##corresponde a 1, 2 ó 3. Los DEGs son aquellos que incrementan (o decrementan) 
##su expresión por encima de dicho umbral (por debajo de menos dicho umbral).
##Este método es biológicamente interpretable de forma directa y no requiere un
##alto número de réplicas biológicas. Se aplica especialmente a estudios con 
##organismo modelos donde no son necesarias muchas réplicas.

##Para la selección de DEGs extraemos el fold change e identificadores de cada
##sonda.

fold.change.WT.with.no.Fe <- WT.with.no.Fe$logFC
genes.ids.WT.with.no.Fe <- rownames(WT.with.no.Fe)

##Fijamos un umbral de 2 que corresponde a genes que se expresan más del doble
##o menos de la mitad. Ya que nuestros datos están transformados por log2 en
##nuestro código usamos log2(2) = 1 

activated.genes.WT.with.no.Fe.1 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe > 1]
repressed.genes.WT.with.no.Fe.1 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe < - 1]

length(activated.genes.WT.with.no.Fe.1)
length(repressed.genes.WT.with.no.Fe.1)



fold.change.pye.with.no.Fe <- pye.with.no.Fe$logFC
genes.ids.pye.with.no.Fe <- rownames(pye.with.no.Fe)

activated.genes.pye.with.no.Fe.1 <- genes.ids.pye.with.no.Fe[fold.change.pye.with.no.Fe > 1]
repressed.genes.pye.with.no.Fe.1 <- genes.ids.pye.with.no.Fe[fold.change.pye.with.no.Fe < - 1]

length(activated.genes.pye.with.no.Fe.1)
length(repressed.genes.pye.with.no.Fe.1)




##Los gráficos que se usan principalmente para representar los DEGs según el
##criterio del fold-change son los gráficos de dispersión o scatterplots. 

plot(wt.with.fe,wt.no.fe,pch=19,cex=0.5,col="grey",xlab="WT with Fe",ylab="WT no Fe")

points(wt.with.fe[activated.genes.WT.with.no.Fe.1],
       wt.no.fe[activated.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="red")

points(wt.with.fe[repressed.genes.WT.with.no.Fe.1],
       wt.no.fe[repressed.genes.WT.with.no.Fe.1],pch=19,cex=0.5,col="blue")

text(wt.with.fe["254550_at"]+0.3,wt.no.fe["254550_at"]+0.3,"IRT1", col="black", cex=0.7)
text(wt.with.fe["252427_at"]+0.3,wt.no.fe["252427_at"]+0.3,"PYE", col="black", cex=0.7)
text(wt.with.fe["257062_at"]+0.3,wt.no.fe["257062_at"]+0.3,"BTS", col="black", cex=0.7)
text(wt.with.fe["251109_at"]+0.3,wt.no.fe["251109_at"]+0.3,"FER1", col="black", cex=0.7)


##Anotación de los genes expresados de forma diferencial

##Podemos extraer información sobre los genes expresados de forma diferencial.
##La librería o paquete *annaffy* nos proporciona funciones para generar listas 
##de DEGs en formato html o txt con la anotación disponible sobre los distintos 
##genes. 

library(annaffy)

##La funicón *affTableAnn* nos permite generar una tabla con DEGs y su 
##anotación. Las funciones *saveHTML* y *saveText* escriben la tabla anterior 
##en formato HTML y txt.

activated.genes.WT.with.no.Fe.1.table <- aafTableAnn(activated.genes.WT.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(activated.genes.WT.with.no.Fe.1.table, 
         file="activated_genes_WT_with_no_Fe_1.html")
saveText(activated.genes.WT.with.no.Fe.1.table, 
         file="activated_genes_WT_with_no_Fe_1.txt")

repressed.genes.WT.with.no.Fe.1.table <- aafTableAnn(repressed.genes.WT.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_WT_with_no_Fe_1.html")
saveText(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_WT_with_no_Fe_1.txt")



activated.genes.pye.with.no.Fe.1.table <- aafTableAnn(activated.genes.pye.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(activated.genes.pye.with.no.Fe.1.table, 
         file="activated_genes_pye_with_no_Fe_1.html")
saveText(activated.genes.pye.with.no.Fe.1.table, 
         file="activated_genes_pye_with_no_Fe_1.txt")

repressed.genes.pye.with.no.Fe.1.table <- aafTableAnn(repressed.genes.pye.with.no.Fe.1, 
                                                     "ath1121501.db", aaf.handler())
saveHTML(repressed.genes.pye.with.no.Fe.1.table, 
         file="repressed_genes_pye_with_no_Fe_1.html")
saveText(repressed.genes.WT.with.no.Fe.1.table, 
         file="repressed_genes_pye_with_no_Fe_1.txt")




## Método que combina inferencia estdística con fold-change:

##Para aplicar este método es necesario tener un alto número de réplicas
##biológicas. Para cada gen y para cada pareja de genotipos/condiciones a 
##comparar se formula un contraste de hipótesis sobre igualdad de medias. 
##Normalmente este contraste de hipótesis utiliza un estadístico similar a la 
##t-student. Se fija un nivel de significancia y se calcula el correspondiente 
##p-valor (y p-valor corregido para el testeo múltiple o q-valor). Si dicho 
##p-valor (o q-valor) es menor que el nivel de significancia y además el
##correspondiente fold-change cumple el umbral se asume que el correspondiente
##gen se expresa de forma diferencial en los genotipos/condiciones estudiadas.

##Para la selección de DEGs según este criterio extraemos el fold change, 
##el p-valor ajustado y los identificadores de cada sonda.

fold.change.WT.with.no.Fe <- WT.with.no.Fe$logFC
p.value.WT.with.no.Fe <- WT.with.no.Fe$adj.P.Val
genes.ids.WT.with.no.Fe <- rownames(WT.with.no.Fe)

##Fijamos como umbral de significancia un p-valor menor a 0.05 y como umbral 
##para el fold change 2 (que corresponde a 1 en log2).

activated.genes.WT.with.no.Fe.2 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe > 1 & 
                                                             p.value.WT.with.no.Fe < 0.05]

repressed.genes.WT.with.no.Fe.2 <- genes.ids.WT.with.no.Fe[fold.change.WT.with.no.Fe < -1 & 
                                                             p.value.WT.with.no.Fe < 0.05]

length(activated.genes.WT.with.no.Fe.2)
length(repressed.genes.WT.with.no.Fe.2)

##Los volcano plots nos permiten representar en un único gráfico el fold change
##y la significancia estadística. Estos gráficos se utilizan para representar 
##los DEGs cuando se usa el método que combina fold-change e inferencia
##estadística.

names(fold.change.WT.with.no.Fe) <- genes.ids.WT.with.no.Fe
log.p.value.WT.with.no.Fe <- -log10(p.value.WT.with.no.Fe)
names(log.p.value.WT.with.no.Fe)  <- genes.ids.WT.with.no.Fe

plot(fold.change.WT.with.no.Fe,log.p.value.WT.with.no.Fe,
     pch=19,cex=0.5,col="grey",ylab="-log10(p value)",xlab="log2 fold change",xlim=c(-6,6))

points(fold.change.WT.with.no.Fe[activated.genes.WT.with.no.Fe.2],
       log.p.value.WT.with.no.Fe[activated.genes.WT.with.no.Fe.2],
       pch=19,cex=0.5,col="red")

points(fold.change.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.2],
       log.p.value.WT.with.no.Fe[repressed.genes.WT.with.no.Fe.2],
       pch=19,cex=0.5,col="blue")

text(fold.change.WT.with.no.Fe["254550_at"],
     log.p.value.WT.with.no.Fe["254550_at"]+0.3,"IRT1", col="black")

text(fold.change.WT.with.no.Fe["251109_at"]+0.3,
     log.p.value.WT.with.no.Fe["251109_at"],"FER1", col="black")

## Mapas de Calor 

##Los mapas de calor nos permiten visualizar los niveles de expresión de los
##DEGs en todas las condiciones analizadas utilizando un código de colores.
##Esta visualización global nos permite identificar genes con comportamientos 
##similares en las condicones analizadas. 

##Para realizar un mapa de calor el primer paso consiste en recopilar todos
##los genes expresados de forma diferencial en un vector y eliminar
##repeticiones.

## Para seleccionar los DEGs en el resto de comparaciones se realizarían de forma
## repetitiva las mismas instrucciones pero con datos diferentes por lo tanto es 
## apropiado definir la siguiente función que recibe como entrada el resultado de la
## función eBayes (contrast.results), el número de genes total (gene.number), el  
## identificador numérico de la comparación a considerar (comparison.number) y el umbral
## de corte en el fold change (fold.change.threshold).
DEGs.fold.change <- function(contrast.results,
                             gene.number,
                             comparison.number,
                             fold.change.threshold)
{
  ## Extraer con topTable la información de DEGs
  transcriptome.comparison <- topTable(contrast.results, number=gene.number,coef=comparison.number,sort.by="logFC")
  ## Extraer el fold change de cada gen y los nombres de las sondas (genes).
  current.fold.changes <- transcriptome.comparison[["logFC"]]
  gene.ids <- rownames(transcriptome.comparison)
  #gene.ids <- transcriptome.comparison[[1]]
  
  ## Extraer el nombre de los genes cuyo fold change excede el umbral prefijado (genes
  ## activados) o que es menor que menos el umbral prefijado (genes reprimidos).
  current.activated.genes <- gene.ids[current.fold.changes > fold.change.threshold]
  current.repressed.genes <- gene.ids[current.fold.changes < - fold.change.threshold]
  
  ## Devolver una lista cuyas componentes se identifican con los nombres "activated.genes"
  ## y "repressed.genes" y que contienen los correspondientes vectores. 
  return(list(activated.genes = current.activated.genes, repressed.genes = current.repressed.genes))
}

diff.genes.1 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=1,fold.change.threshold=1) 
diff.genes.2 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=2,fold.change.threshold=1) 
diff.genes.3 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=3,fold.change.threshold=1) 
diff.genes.4 <- DEGs.fold.change(contrast.results,gene.number=22810,comparison.number=4,fold.change.threshold=1) 

complete.DEGs <- c(diff.genes.1$activated.genes,
                   diff.genes.1$repressed.genes,
                   diff.genes.2$activated.genes,
                   diff.genes.2$repressed.genes,
                   diff.genes.3$activated.genes,
                   diff.genes.3$repressed.genes,
                   diff.genes.4$activated.genes,
                   diff.genes.4$repressed.genes)

complete.DEGs <- unique(complete.DEGs)
length(complete.DEGs)

##A continuacón, extraemos los niveles de expresión de todos los DEGs en las
##condiciones/genotipos de interés y normalizamos los datos teniendo en cuenta 
##que la función scale normaliza por columnas y nosotros queremos hacerlo por
##filas. Calculamos con la función t la traspuesta de nuestra matriz. 

DEG.expression <- mean.expression[complete.DEGs,
                                  c("WT_with_Fe","WT_no_Fe","pye_with_Fe","pye_no_Fe")]

normalized.DEG.expression <- t(scale(t(DEG.expression)))

##Con la función heatmap.2 del paquete gplots construímos el correspondiente 
##heatmap.

library(gplots)
help(heatmap.2)
heatmap.2(normalized.DEG.expression,Colv=FALSE,dendrogram="row",
          labRow=c(""),density.info="none",trace="none",
          col=redgreen(100),margins = c(8,8),cexCol=1.2)

