#  Pipeline Transcriptomic analysis using Salmon/HISAT2--DESeq2/EdgeR_limmavoom

By: Esteban Velásquez Agudelo

This document is a guide designed specifically for linux architectures to make a differential transcriptomics analysis with the most recomended programs to date according to some literature in 2020 [[1](#7.-referencias)]. These softwares and packages have had recent updates and changes, like DESeq2 with lfcshrink [[2](#7.-referencias)]. It is advised to read the last reviews with the latest updates in posterior dates. This guide is done using ubuntu 18.04 and with 6 samples:  2 Naturally resistant, 2 Artificially resistant and 2 susceptible strands. Apart from this guide one could use cufflinks, but not with tophat because is deprecated, and an useful guide to use it can be found here[^9](#6.-footnotes).


# Methodology
 A transcriptomic analysis in general is done in 5 steps, illustrated in Figure 1. in general these steps are:


1. Quality analysis of the data and posterior trimming.
2. Data alignment to a reference genome or transcriptome [^1](#6.-footnotes).
3. Diferential trancription analysis from gene count or transcript quantification.
4. Ontological and pathway analysis of these differential transripted genes

<p align="center">
  <img width="460" height="300" src=TRANSCRIPTOMICANALYSIS.svg>
  <p align="center"> Figure 1. Methodology of a transcriptome analysis <p>
</p>

# Index
<p align="center">

<!-- TOC depthFrom:1 depthTo:6 updateOnSave:false-->
[Begginers](#begginers)<br> [1. Quality data analysis and filtering](#1-quality-data-analysis-and-filtering)<br> 
- [1.1 Quality analysis](#11-quality-analysis)<br>
- [1.2 Filtering and adapter trimming using trimmomatic](#12-filtering-and-adaptor-trimming-using-trimmomatic)<br>
  
[2. Filtered read alignment](#2filtered-read-alignment)<br>
- [2.1 Read alignment using HISAT2](#21-read-alignment-using-hisat2)<br>
- [2.2 Salmon alignment](#22-salmon-alignment)

[3. Diferrential transcriptional analysis](#3-diferrential-transcriptional-analysis)
- [3.1 Installing and preparing r studio](#31-installing-and-preparing-r-studio)
- [3.2 DESeq2 analysis from hisat2 using featurecounts](#32-deseq2-analysis-from-hisat2-using-featurecounts)
- [3.3 DESeq2 analysis from salmon](#33-deseq2-analysis-from-salmon)
- [3.4 Analysis using limma from featurecounts](#34-analysis-using-limma-from-featurecounts)
- [3.5 Analysis using limma from salmon](#35-analysis-using-limma-from-salmon)


[4. Posterior analysis](#4-posterior-analysis)

- [4.1 Obtaining the archive of GO terms and the fasta files of the differentailly transcripted genes(DTE)](#41-obtaining-the-archive-of-go-terms-and-the-fatsa-files-of-the-diferentially-transcripted-genesdte)
- [4.2 R graphs of the GO terms](#42-r-graphs-of-the-go-terms)
- [4.3 Set analysis of DTG using the fisher test ](#43-set-analysis-of-dtg-using-the-fisher-test)

- [4.4 Pathway analysis using KEEG-KAAS](#44-pathway-analysis-using-keeg-kaas))

[5. Optionals](#5-optionals)

- [5.1 Counting with HT-Seq](#51-counting-with-ht-seq)

- [5.2 Generating a gtf from a gff](#52-generating-a-gtf-from-a-gff)

- [5.3 Changing the formats of some reads with biopython](#53-changing-the-formats-of-some-reads-with-biopython)
- [5.4 Changing file permisions](#54-to-change-permisiion-in-files)
- [5.5 Transcriptome de novo reconstruction](#55-transcriptome-de-novo-reconstruction)



[6. Footnotes](#6.-footnotes). 
<!-- /TOC -->
</p>

# Begginers

The first thing is to learn how to create environments inside of linux to be more efficient in computer CPU usage and in RAM usage. Environment creation is done with Conda, which is a package manager, but it also can be done without it. Conda is used here beacuse its easier to understand and use.

* To install conda use:
```bash 
bash Anaconda-latest-Linux-x86_64.sh
```

* with this command environments are created:
```bash 
conda create -n environment_name -c channel %program_to_download
```
* With these two comands the environment activates and deactivates
```bash 
conda activate environmente_name 
conda deactivate
```

Apart from this, it is necesary to know how to move inside a linux terminal with basic commands that move through folders, copying them and erasing them.

# 1. Quality data analysis and filtering

## 1.1 Quality analysis: 
 It is recomended to use FASTQC to all of the data you are going to analyse; and see the graph inside the html output of FASTQ to see which actions must be done. Some of the things one must verify are the percentage of bases in each position of the read, adapters content, quality in each position and the length of the reads.

In order to create an environment and run FASTQC for one sample you should use these commands, in which P102.fastq are the reads that are going to be analysed :

```bash 
conda create -n fastqc fastqc
conda activate fastqc
fastqc -o P102 -t 4 P102.fastq
```
<p align="center">
  <img width="460" height="300" src=fastqc.png>
  <p align="center"> Figure 2. Quality of each position in the FASTQC output<p>
</p>

In the case of the example, one can see that at the end of the reads the quality goes very low, so it can be filtered, but we also have to take into account  waht we want when we align the reads, and explore the parameters in order to maximize that quality. in this specific case we see that the first 15 reads are equal for all the reads, so that is an adapter to be removed. we are also going to maximize unique read alignment by changing parameters in trimmomatic.

##  1.2 filtering and adaptor trimming using trimmomatic

One should download the binary version and uncompress it to the folder in which the data is, it can be downloaded here :http://www.usadellab.org/cms/?page=trimmomatic

The commands to filter the data and get the adapter trimmed; in which the first one is for Single reads and the second is for Paired reads are (the option to get the adapters trimmed is  HEADCROP and the best option to vary for data filtering in this dataset is the last number in SLIDINGWINDOW) :

```bash 
java -jar trimmomatic-0.39.jar SE -phred33 -threads 8 P101T.fastq P101Trim.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:15

java -jar trimmomatic-0.39.jar PE -phred33 -threads 8 P101_F.fastq P101_R.fastq%(acá es necesario poner primero el forward y despues el reverse)% P101_forward_Paired P101_forward_Unpairedtrim.fq P101_reverse_Pairedtrim.fq  P101_forward_Unpairedtrim.fq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINlEN:36 HEADCROP:15
```
#  2.Filtered read alignment

## 2.1 Read alignment using HISAT2

 HISAT2 is the best genome aligner with STAR, although it uses a lot  of memmory and time. To use it one must first create an environment with HISAT2 installed, then one must create an index to align to and then do the alignment with the first commmand for single reads and the second one for paired reads.

 ```bash
 conda create -n hisat hisat2

 hisat2-build TriTrypDB-46_TcruziSylvioX10-1_Genome.fasta Tcruzi

 hisat2 -q -p 8 -x Tcruzi -U P101T.fastq -S P101T.sam

 hisat2 -q -p 8 -x Tcruzi -1 P101T_F.fastq -2 P101T_R -U P101T_Unpaired.fastq -S P101TP.sam

```
### 2.1.1 Changing the format from sam to bam

To do this we should use the package samtools and the commands are the following:

```bash
conda create -n samtools samtools

samtools sort -@ 8 -o P101T12.bam P101Trim12.sam
 
samtools index ./P101/P101T12.bam

```

__Output__ of HISAT2:
```
134810 reads; of these:<p>
  134810 (100.00%) were unpaired; of these:<p>
    58794 (43.61%) aligned 0 times<p>
    44963 (33.35%) aligned exactly 1 time<p>
    31053 (23.03%) aligned >1 times<p>
56.39% overall alignment rate<p>
```
The most important quantity is the aligned exactly one time at least for the <em>T. cruzi </em>, in which the genome is pretty repetitive. One must change the SLIDINGWINDOW to maximize the mayority of alignments of exactly 1 time in the reads in this dataset.


## 2.2 Salmon Alignment

Salmon is the most preferable tool to do alignment and quantification due to the short time of usage and the alignment being done in respect to the transcriptome[[3](#7.-referencias)]. It has 2 probable downsides which are the bad analysis for the small RNA and the lowly expressed RNA[[4](#7.-referencias)]. To use the last version of Salmon, which in this case is the best, one must download the binary and uncompress it in the folder containing the filtrered reads here: https://github.com/COMBINE-lab/salmon/releases/tag/v1.3.0 and follow the next steps:

* uncompress the folder:
```bash
tar xzvf salmon-1.2.0_linux_x86_64.tar.gz
```
* Creating a decoy with the transcriptome and the genome: (https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)[^2](#6.-footnotes)
```bash
grep "^>" <TriTrypDB-46_TcruziSylvioX10-1_Genome.fasta | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
cat TriTrypDB-48_TcruziSylvioX10-1_AnnotatedTranscripts.fasta TriTrypDB-46_TcruziSylvioX10-1_Genome.fasta > gentrome.fasta

```
* Creating an index
```bash
./salmon-late_linux_x86_64/bin/salmon index -t gentrome.fasta -d decoys.txt -p 8 -i salmon_index --gencode
```
* Read quantification for single reads and paired read respectively, in this case theA is used for automatic library identification :
```bash
./salmon-latest_linux_x86_64/bin/salmon quant -i salmon_index -l A -r ./P101/fastqs/P101T12.fastq --validateMappings -o transcripts_quantP101

./salmon-latest_linux_x86_64/bin/salmon quant -i transcripts_index -l A -1 P101_1.fq -2 P101_2.fq --validateMappings -o transcripts_quantP101

```
After these steps one goes directly to R to make the analysis with DESSeq2 and/or EdgeR limma-voom.

# 3. Diferrential transcriptional analysis

## 3.1 Installing and preparing R studio

All the steps followed here are done in R  4.0 in R studio, and to install it one must follow these commands in the linux terminal:
```bash
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'

sudo apt update

sudo apt install r-base
```
* Now to install r studio one has to download  the r studio version of the distribution one is using, which for this example is ubuntu 18.04; it can be downloaded here :https://rstudio.com/products/rstudio/download/#download
  
```bash
sudo apt install gdebi-core

sudo gdebi  rstudio-1.3.1093-amd64.deb
```

The last step is to install some dependencies in the linux terminal that r needs for certain packages that one will use:

```bash
sudo apt-get install libcurl4-gnutls-dev
sudo apt-get install libxml2-dev
sudo apt-get install libssl-dev
```

It is recommended to open rstudio from terminal with sudo permisions, and to install packages inside of r one must write in the r terminal:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rsubread")
```

## 3.2 DESeq2 analysis from HISAT2 using featurecounts

In R the first thing that has to be done is to count the mapped reads to the genome, which can be done eassily with featurecounts (this can also be done with HT-Seq, but its much slower, a guide of this is presented in optional steps). To install it one must install the package Rsubread with the comand mentioned earlier.

### 3.2.1 counting mapped reads
 The following command is the one to count mapped reads, and we must be careful with the featureType option, which is the name of the gene that we are going to save. it can be inspected in the gtf and it must be in the following format TcYSYL_0000001. The last command is to verify that we did it nicely:
```r
library(Rsubread)

fc <- featureCounts(files=c("./P101/P101T12.bam", "./P102/P102T12.bam", "./RCl3/RCl3T12.bam", "./RCL4/RCl4T12.bam","./RL11/RL11T12.bam","./RL12/RL12T12.bam"),nthreads=6,annot.ext="TriTrypDB-46_TcruziSylvioX10-1.gff",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="ID")

head(fc$counts)
```
__Output__:
```
                      P101T12.bam P102T12.bam RCl3T12.bam RCl4T12.bam RL11T12.bam RL12T12.bam
exon_TcSYL_0117800-E1           1           0           7           1           3           3
exon_TcSYL_0014440-E1           0           0           0           0           1           0
exon_TcSYL_0077180-E1           1           0           3           0           1           2
exon_TcSYL_0058970-E1           0           1           0           0           0           1
exon_TcSYL_0058900-E1           0           0           0           0           1           0
exon_TcSYL_0045600-E1           0           0           0           0           1           0
```
### 3.2.2 Importing featurecounts to DESeq2
* This specific guide is based on another guide, which is the best one according to the DESEq2 creator [^8](#6.-footnotes). 

#### 3.2.2.1 Creating a vector that summarizes the experimental design:

 One should verify that the order of the condition is the same as the fc  with the last command, and the commands to use are[^3](#6.-footnotes):
```r
condition<-c(rep("Artificial",4),"Natural","Natural")
coldata <- matrix(c(rep("Artificial",4),"Natural","Natural"), dimnames = list(colnames(fc$counts), 'condition'))


coldata
```
__Output__:
```
            condition   
P101T12.bam "Artificial"
P102T12.bam "Artificial"
RCl3T12.bam "Artificial"
RCl4T12.bam "Artificial"
RL11T12.bam "Natural"   
RL12T12.bam "Natural"   
```
#### 3.2.2.2  To import the data one must use the following command:
To save the information of the genes the command of the second paragraph is used and the last command is used to verify that everything went smooth:
```r
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData=fc$counts, colData=coldata(fc$counts), design=~condition)


featureData <- data.frame(gene=rownames(fc))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

colData(dds)
```
__Output__:
```r 
             condition
              <factor>
P101T12.bam Artificial
P102T12.bam Artificial
RCl3T12.bam Artificial
RCl4T12.bam Artificial
RL11T12.bam Natural   
RL12T12.bam Natural  
```

### 3.2.3 Generating the differential analysis

the steps to be followed are:

#### 3.2.3.1. Agregate the factors well
```r
colData(dds)$condition <- factor(colData(dds)$condition, levels=c("Artificial","Natural"))

dds$condition <- relevel(dds$condition, ref = "Artificial")
```
#### 3.2.3.2. Previsualizing the data

The colors are the number of experiments
```r
boxplot(log10(counts(dds)+1),col=c("red","green","blue","yellow",'brown','gray'),las=2)
title(main="Unnormalized data and library size DEseq2", ylab="Log-counts") 
```
<p align="center">
  <img width="460" height="300" src=unnormDESeq2.png>
  <p align="center"> Figure 3. Data previsualization using DESeq2<p>
</p>

#### 3.2.3.3. Transforming the data to PCA visualization:
When to many differences are expected the option of blind must be false.
```r
rld <- rlogTransformation(dds, blind = TRUE)
plotPCA(rld, intgroup = c("condition"))
```
<p align="center">
  <img width="460" height="300" src=PCADESeq2.png>
  <p align="center"> Figura 4. PCA DESeq2 <p>
</p>

#### 3.2.3.4. Biological variance estimation 
This is done to see the options of the different fitTypes used, be it local or parametric, this ca be done with fitType='local' with a comma in estimateDispersion [^7](#6.-footnotes).
```r
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)
```
<p align="center">
  <img width="460" height="300" src=dispDESeq2.png>
  <p align="center"> Figure 5. Dispersion DESeq2 <p>
</p>

#### 3.2.3.5. Determinend the diferential transcription,
* One can change the fitType but one must be carefull, the default and generally the best one is parametric , but it's better to analyze it with the dispersion plot and this can be changed with fitType='local'. One can ask like this guy in the bioconductor supportto see which option is the best :https://support.bioconductor.org/p/81094/
```r
dds<-estimateSizeFactors(dds)
dds <- estimateDispersions(dds,fitType='local')
dds <-nbinomWaldTest(dds)
res <- results(dds, contrast = c('condition', 'Natural', 'Artificial'),alpha=0.05)
summary(res)
```
__Output__:
```
out of 16240 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1, 0.0062%
LFC < 0 (down)     : 39, 0.24%
outliers [1]       : 2, 0.012%
low counts [2]     : 13657, 84%
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

#### 3.2.3.6. Removing the size effect with the last updated package:
One can change the hypothesis used here by changing the lfc threshold[^4](#6.-footnotes)
```r
coefi<-resultsNames(dds)
lfcutoff<-0.5
res.lfc2 <- lfcShrink(dds, coef=coefi[2], type="apeglm",res=res) 
res<-res.lfc2
summary(res.lfc2)
```
__Output__:
```
out of 16240 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1, 0.0062%
LFC < 0 (down)     : 39, 0.24%
outliers [1]       : 2, 0.012%
low counts [2]     : 13657, 84%
(mean count < 4)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

#### 3.2.3.7. Visualizing the results with an MAplot:
For publication one can use a single graph
```r
par(mfrow=c(1,2))
lims <- c(-8,8)
xlab <- "mean of normalized counts"
plotMA(res.lfc2, ylim=c(-2,2), main="apeglm, LFC NaturalvsArtificial", alpha=0.05)
hline()
plotMA(dds,ylim=c(-2,2),main="DESeq2")
hline()
```

<p align="center">
  <img width="460" height="300" src=MADESeq2.png>
  <p align="center"> Figure 6. MA plot <p>
</p>

#### 3.2.3.8. Making an interactive graph for better analysis:
```r
library("Glimma")
status <- as.numeric(res$padj < .1)
anno <- data.frame(GeneID=rownames(res), symbol=res$padj)
glMDPlot(res, status=status, counts=counts(dds,normalized=TRUE), groups=dds$condition, transform=FALSE, samples=colnames(dds), anno, folder="glimma", launch=TRUE)
```
#### 3.2.3.9. Creating a list of differentially positive and negative trancripted genes
One obtains a list of positive an a list of negative making the nomenclature the same as TriTryp like TcYSyl_0000001.
```r
df<-data.frame(res)
DOWN<-which(df$padj<0.05  & df$log2FoldChange<0 )
UP<-which(df$padj<0.05  & df$log2FoldChange>0 )
DEDOWN<-rownames(df[DOWN,])
DEUP<-rownames(df[UP,])


Nomgenes <- function(DEDOWN) {
  for(i in 1:length(DEDOWN)){
    DEDOWN[i]<-paste("TcSYL_",sub(sprintf(".*%s(.*)%s.*","_", "-"), "\\1", DEDOWN[i]))
  }
  return(DEDOWN)
}
trimgenes <- function(DOWNshrink) {
  for (i in 1:length(DOWNshrink)){
    DOWNshrink[i]<-gsub(" " ,"",DOWNshrink[i])
  }
  return(DOWNshrink)
}

DEDOWN<-trimgenes(Nomgenes(DEDOWN))
DEUP<-trimgenes(Nomgenes(DEUP))


write(paste(DEDOWN, collapse = ','), 'DESeq2.DEGdownlocal.results.txt') 
write(paste(DEUP, collapse = ','), 'DESeq2.DEGuplocal.results.txt')
```
## 3.3 DESeq2 analysis from Salmon 

### 3.3.1 Creating a matrix of gene-transcript
It's better to have the .gtf in the same folder as where the analysis is going to be done, the packages must be installed beforehand
```r
library(tximport)
library("readr")
library("tximportData")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("TriTrypDB-46_TcruziSylvioX10-1.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
```

### 3.3.2 Creating a matrix of the experimental design 
To create it one must write the next commands, it is also to neccesary to verify that the last output is true with the last command:
```r
vv <- c("transcripts_quantP101", "transcripts_quantP102", "transcripts_quantRCl3", "transcripts_quantRCl4","transcripts_quantRL11","transcripts_quantRL12",rep("Artificial",4),"Natural","Natural",1,2,3,4,5,6)
samples <- matrix(vv,nrow= 6, ncol = 3, byrow = FALSE)
colnames(samples) <- c("sample","Resistance","Experiment")
samples<-data.frame(samples)
samples$condition <- factor(c(rep("Artificial",4),rep('Natural',2)))
files <- file.path(".", samples$sample, "quant.sf")
names(files) <- paste0(samples$sample)
all(file.exists(files)

```
### 3.3.2 Importing the salmon files
```r
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
names(txi)
library(DESeq2)
condition <- factor(c(rep("Artificial",4),rep('Natural',2)))
dds <- DESeqDataSetFromTximport(txi, samples, ~condition)
```
now we should do the steps from [3.2.3.1](#3.2.3.1.-agregar-los-factores-bien)

## 3.4 Analysis using limma from featurecounts

### 3.4.1 Counting using featurecounts:

```r
library(Rsubread)
fc <- featureCounts(files=c("./P101/P101T12.bam", "./P102/P102T12.bam", "./RCl3/RCl3T12.bam", "./RCL4/RCl4T12.bam","./RL11/RL11T12.bam","./RL12/RL12T12.bam"),nthreads=8,annot.ext="TriTrypDB-46_TcruziSylvioX10-1.gff",isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="ID")
head(fc$counts)
```
### 3.4.2 Creating a DGEobject from featurecounts
```r
library(edgeR)
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
```
### 3.4.3 Data previsualization
```r
lcpm <- cpm(x, log=TRUE)
boxplot(lcpm, las=2,main="",col=c("red","green","blue","brown","yellow","gray"))
title(main="Unnormalized data and library size", ylab="Log-cpm")
```
<p align="center">
  <img width="460" height="300" src=boxplotlimma.png>
  <p align="center"> Figure 7. Previsualization of data with Limma<p>
</p>

### 3.4.4 Data filtering
```r
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]
```

### 3.4.5 Normalizing the data with  TMM

```r
x2 <- calcNormFactors(x, method = "TMM")
lcpm2 <- cpm(x2, log=TRUE)
boxplot(lcpm2, las=2,main="",col=c("red","green","blue","brown","yellow","gray"))
title(main="Normalised data", ylab="Log-cpm")

```

<p align="center">
  <img width="460" height="300" src=boxplotnormlimma.png>
  <p align="center"> Figure 8. Visualizing normalized data with limma<p>
</p>


### 3.4.6 Creating a designed matrix

```r
vv <- c("P101", "P102", "RCl3", "RCl4","RL11","RL12",rep("Artificial",4),"Natural","Natural",1,2,3,4,5,6)
targets <- matrix(vv,nrow= 6, ncol = 3, byrow = FALSE)
colnames(targets) <- c("Filename","Resistance","experiment")
celltype <- factor(targets[,2])
design <- model.matrix(~0+celltype)
colnames(design) <- levels(celltype)
```
### 3.4.7 Normalization based on voom

```r
y1<- voom(x2,design,plot=TRUE)
```
### 3.4.8 Clusterization with PCA

```r
opar <- par(no.readonly = TRUE)
par(xpd = TRUE, mar = par()$mar + c(0, 0, 0, 5))
plotMDS(y1,xlim=c(-2.5,3),labels = c("P101", "P102", "RCl3", "RCl4","RL11","RL12"),col=c(rep("green",4),"red","red"))
legend(par("usr")[2], par("usr")[4], c("Artificial","Natural"), pch = 16, col = c("green","red"), bty = "n")
```
<p align="center">
  <img width="460" height="300" src=PCAfilter.png>
  <p align="center"> Figure 9. PCA with Limma<p>
</p>

### 3.4.9 Aplying the linear model for the differential transcription:
<br>
In this case the option Trend and robust are used. Robust filters the genes taht have a variance very different fron the mean, and trend is used when the difference of mean of the size of each library is not greater than 3 fold. 

```r
fit1 <- lmFit(y1,design)
contr <- makeContrasts(ArtificialvsNatural=Artificial-Natural,levels=design)
fit1.contr <- eBayes(contrasts.fit(fit1,contr),trend = TRUE,robust = TRUE) 
dt1 <- decideTests(fit1.contr)
summary(dt1)
```
### 3.4.10 Comparing the models:

Now the dispersions are compared from each model iwth the command [[3.4.9](#-349-Aplying-the-linear-model-for-the-differential-rasncription)] taking the robust and trend to false independently and then comparing the different SA plots to se which one is better to use:.
```r
plotSA(fit1.contr)
```
<p align="center">
  <img width="460" height="300" src=SAplot.png>
  <p align="center"> Figure 10. Dispersion with Limma trend and robust<p>
</p>

### 3.4.11 Visualizing the data with a MA plot and an interactive MAplot using limma 
 The Ma is for publications and the glimma is used to analize carefully the data and later analysis.
```r
tfit1 <- treat(fit1, lfc=1)
plotMD(tfit1, column=1, status=dt1[,1], main=colnames(tfit1)[1], xlim=c(0,15),ylim=c(0,12))
library(Glimma)
glMDPlot(tfit1, coef=1, status=dt, main=colnames(tfit1)[1],side.main="ENTREZID", counts=lcpm, groups=celltype, launch=FALSE)
```
<p align="center">
  <img width="460" height="300" src=MAplotlimma.png>
  <p align="center"> Figure 11. MAPlot of limma with robust and trend<p>
</p>

### 3.4.12 Creating output files
Saving the list of the diferentially transcripted genes for GO and pathways analysis. This aproximation is a little brute force so maybe you can upgrade it.

#### 3.4.12.1 Creating a table of diferentially transcripted genes

```r
voom.tt <- topTable(fit1.contr, adjust.method = "BH",number= Inf)
voom.limma.results <- data.frame(gene.id = rownames(voom.tt), p.value = voom.tt$P.Value, q.value = voom.tt$adj.P.Val, log.fc = voom.tt$logFC, test.stat = voom.tt$t, row.names = NULL)
```
#### 3.4.12.2 selecting the genes down and up regulated
```r
De.limma.down<- which(dt1[,1]==-1)
De.limma.up<-which(dt1[,1]==1)
de.limma.down=c()
for(i in De.limma.down){
  de.limma.down<-c(de.limma.down, dimnames(dt1[i,0])[[1]])
}
de.limma.down
de.limma.up=c()
for(i in De.limma.up){
  de.limma.up<-c(de.limma.up, dimnames(dt1[i,0])[[1]])
}
de.limma.up
position<-voom.limma.results$q.value<0.05
DEgenesdown<-c()
```
#### 3.4.12.3 changing the format to TriTryp
 The names of the genes are searched in the tritryp format specified above, without spaces
```r
for(i in 1:length(position)){
  if(position[i]==TRUE & voom.limma.results$gene.id[i] %in% de.limma.down){
    DEgenesdown <- c(DEgenesdown, trimws(paste("TcSYL_",sub(sprintf(".*%s(.*)%s.*","_", "-"), "\\1", voom.limma.results$gene.id[i]))))
  }
}
DEgenedown=c()
for (x in DEgenesdown){
  DEgenedown<-c(DEgenedown, gsub(" ", "", x, fixed = TRUE))
}
DEgenesup<-c()
for(i in 1:length(position)){
  if(position[i]==TRUE & voom.limma.results$gene.id[i] %in% de.limma.up){
    DEgenesup <- c(DEgenesup, trimws(paste("TcSYL_",sub(sprintf(".*%s(.*)%s.*","_", "-"), "\\1", voom.limma.results$gene.id[i]))))
  }
}
DEgeneup=c()
for (x in DEgenesup){
  DEgeneup<-c(DEgeneup, gsub(" ", "", x, fixed = TRUE))
}
```
#### 3.4.12.4 veryfication and writing
The data are visually inspected and then saved
```r
DEgeneup
DEgenedown 
write(paste(DEgenedown, collapse = ','), 'limma.DEGdown.results.txt')
write(paste(DEgeneup, collapse = ','), 'limma.DEGup.results.txt')
```

## 3.5 Analysis using limma from salmon

A matix of pairs gene-transcripts is created the same way as in the DESeq2 in [3.3.1](#331-Creating-a-matrix-of-gene-transcript)

### 3.5.2 Importing the data to create a DGE object:

```r
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "lengthScaledTPM")
library(limma)
library(edgeR)
y <- DGEList(txi$counts)
```
### 3.5.3 Filtering the data ad creating a condition then normalizing
```r
isexpr <- rowSums(cpm(y) > 5) >= 4
y <- y[isexpr,]
condition <- factor(c(rep("Artificial",4),rep('Natural',2)))
y <- calcNormFactors(y,method = "TMM")
lcpm2 <- cpm(y, log=TRUE)
boxplot(lcpm2, las=2,main="",col=c("red","green","blue","brown","yellow","gray"))
title(main="Normalised data", ylab="Log-cpm")
y <- voom(y, design)

```
<p align="center">
  <img width="460" height="300" src=normalizedcounts.png>
  <p align="center"> Figure 12. Boxplot of normalized using salmon and limma<p>
</p>

The next steps are the same as in the [3.4.8](#348-Clusterization-with-PCA) with a little difference at the end ([3.4.12](#3.4.12-crear-archivos-de-salida)) because the names are from transcripts and not genes:
```r
DEgenesdown <- c(DEgenesdown,substr(voom.limma.results$gene.id[i],1,nchar(voom.limma.results$gene.id[i])-3))
```

# 4. Posterior analysis

## 4.1  Obtaining the archive of GO terms and the fatsa files of the Diferentially transcripted genes(DTE).

First one must enter to the site :https://tritrypdb.org/tritrypdb/app/ then click in the box  __Search for__  after that click in the __Genes__, afterwards __Annotation, curaton and identifiers__ and lastly __Gene ID__ . In this window one should select __upload a text file__  now you name your search and click in __Get Answer__. Susequently one must click in Download[^5](#6.-footnotes). <br>
After that one selects depending of the analysis: 

### 4.1.1 Go analysis:

 One must select : __Tab- or comma-delimited (openable in Excel) - choose a pre-configured table__. then check  __function prediction y GO Terms__. following that one slecets __Comma-delimited (.csv) file__ and finally  __Get__.

### 4.1.2 Pathway analysis:

One must create a fasta file of the up and down regulated genes independently in the following way: check in __FASTA - sequence retrieval, configurable__. Then slect __Protein__, after that __Text File__ and lastly __Get Sequences__.

## 4.2 R graphs of the GO terms
Once we have all the .csv we proceed to make some graphs summarizing the results in r in this way:

### 4.2.1 reading the data of a csv file in the same folder as the script

```r
DOWNDEG<-read.csv("./salmonDeseqresults/SalmonDESeq2down.csv")
```

### 4.2.2 Pie grapgh of the porcentage of each Go term

#### 4.2.2.1 for the down regulated:
```r
BP<-sum(DOWNDEG$Ontology=="Biological Process")/length(DOWNDEG$Ontology)
MF<-sum(DOWNDEG$Ontology=="Molecular Function")/length(DOWNDEG$Ontology)
CC<-sum(DOWNDEG$Ontology=="Cellular Component")/length(DOWNDEG$Ontology)

pie(c(BP,MF,CC),labels=c(paste(round(100*BP,2),"%"),paste(round(100*MF,1),"%"),paste(round(100*CC,1),"%")),explode=0.1, main="",col=c("blue", "lightblue", "cyan"))
legend("topright", c("Biological Process","Molecular Function","Cellular Component"), cex = 0.56,fill = c("blue", "lightblue", "cyan"),box.lty=0,inset = 0.062)

```
<p align="center">
  <img width="460" height="300" src=PieDown.png>
  <p align="center"> Figure 13. Pie chart of down regulated genes in salmon by DESeq2<p>
</p>

#### 4.2.2.1 For the down regulated:
```r
UPDEG<-read.csv("./salmonDeseqresults/SalmonDESeq2up.csv")
BP<-sum(UPDEG$Ontology=="Biological Process")/length(UPDEG$Ontology)
MF<-sum(UPDEG$Ontology=="Molecular Function")/length(UPDEG$Ontology)
CC<-sum(UPDEG$Ontology=="Cellular Component")/length(UPDEG$Ontology)
pie(c(BP,MF,CC),labels=c(paste(round(100*BP,1),"%"),paste(round(100*MF,1),"%"),paste(round(100*CC,1),"%")),explode=0.1, main="",col=c("blue", "lightblue", "cyan"))
legend("topleft", c("Biological Process","Molecular Function","Cellular Component"), cex = 0.58,fill = c("blue", "lightblue", "cyan"),box.lty=0,inset =  0.04)

```
<p align="center">
  <img width="460" height="300" src=PieUp.png>
  <p align="center"> Figure 14. Pie chart of up regulated genes in salmon by DESeq2<p>
</p>


### 4.2.3 Graph of principal 15 Go terms being down and up regulated:
First the 15 most cited terms are filtered
#### 4.2.3.1 For Down regulated:
```r
terms<-unique(DOWNDEG$GO.Term.Name)
counts<-c(rep(0,length(terms)))
for( i in 1:length(terms)){
  for(j in DOWNDEG$GO.Term.Name){
    if(terms[i]==j){
      counts[i]=counts[i]+1
    }
  }
}
df<-data.frame(counts,terms)
df1 <- df[order(-counts),] 
df1> df[with(df, order(-counts, terms)), ]
df2<-df1[1:15,]
```
#### 4.2.3.1 for the up regulated
One counts the different processes and then the quantity:
```r
termsup<-unique(UPDEG$GO.Term.Name)
countsup<-c(rep(0,length(termsup)))
for( i in 1:length(termsup)){
  for(j in UPDEG$GO.Term.Name){
    if(termsup[i]==j){
      countsup[i]=countsup[i]+1
    }
  }
}
countsup
dfup<-data.frame(countsup,termsup)
df1up <- dfup[order(-countsup),] 
df1up> dfup[with(df, order(-countsup, termsup)), ]
df2up<-df1up[1:15,]
```
### 4.2.4 Unifying the data frames in a single graph:
```r
names(df2up)<-names(df2)
df2up<-na.omit(df2up)
df2<-na.omit(df2)
df2all<-rbind(df2, df2up)
```
### 4.2.5 making the graph:
[^10](#6.-footnotes)
```r
library(ggplot2)
library(dplyr)
library(forcats)
colors<-c(rep("#f68060",length(df2[,1])),rep("#60d6f6",length(df2up[,1])))
rownames(df2all) <- c()
df2all["Tgene"]<-c(rep(1,length(df2[,1])),rep(2,length(df2up[,1])))
df2all["Type"]<-c(rep("DOWN",length(df2[,1])),rep("UP",length(df2up[,1])))
df2all
df2all %>%
  ggplot( aes(x=reorder(terms,Tgene), y=counts,fill=Type)) +
  geom_bar(stat="identity", alpha=.8, width=.5) +
  coord_flip() +
  xlab("") +
  theme_bw()+
  ylab("Number of genes")+
  scale_fill_manual(values = c("#f68060", "#60d6f6"))
```
<p align="center">
  <img width="460" height="300" src=Goanalysis.png>
  <p align="center"> Figure 15. Graph of the 15 most important GO terms<p>
</p>


## 4.3 set analysis of DTG using the Fisher test
 This analysis is usefull when there are a lot of DTG >20, generally the ones from limma to make more biologically significant set of genes that represent a GO term.
### 4.3.1 Enrichment analysis using the fisher test in TriTryp 

* One enters https://tritrypdb.org/tritrypdb/app/ and the slecets my strategies since one already made some searches or follow [4.1](#41-obtaining-the-archive-of-go-terms-and-the-fatsa-files-of-the-diferentially-transcripted-genesdte). Being inthe results tab we click in __Analize results__ an then select __Gene ontology Enrichment__.<br>

* after the page did the search we select __cellular component__, then we click in __Download__, after that we do the same but with  __Molecular Function__ and __Biological Process__ .[^6](#6.-footnotes)

 ### 4.3.2 reading the data from Tritryp
```r
library(data.table)
CCDOWN<-read.delim("./salmonlimmaresults/GoEnrichmentCCDOWNlimma.tab",header = T, sep = "\t")
MFDOWN<-read.delim("./salmonlimmaresults/GoEnrichmentMFDOWNlimma.tab",header = T, sep = "\t")
MFUP<-read.delim("./salmonlimmaresults/GoEnrichmentMFUPlimma.tab",header = T, sep = "\t")
BPDOWN<-read.delim("./salmonlimmaresults/GoEnrichmentBPDOWNlimma.tab",header = T, sep = "\t")
BPUP<-read.delim("./salmonlimmaresults/GoEnrichmentBPUPlimma.tab",header = T, sep = "\t")
```
### 4.3.2 Creating tables with the lower P-value
```r
library(gridExtra)
df<-data.frame(CCDOWN$Name[1:10],MFDOWN$Name[1:10],MFUP$Name[1:10],BPDOWN$Name[1:10],BPUP$Name[1:10])
df1<-data.frame(CCDOWN$Name[1:10],MFDOWN$Name[1:10])
names(df1)<-c('CCUP','MFUP')
grid.table(df1)
df2<-data.frame(BPDOWN$Name[1:10])
names(df2)<-c('BPUP')
grid.table(df2)
df3<-data.frame(MFUP$Name[1:10],BPUP$Name[1:10])
names(df3)<-c('MFDOWN','BPDOWN')
grid.table(df3)
```

<p align="center">
  <img width="460" height="300" src=enrichment.png>
  <p align="center"> Figure 16. Enriched genes table<p>
</p>

## 4.4 Pathway analysis using KEEG-KAAS

Enter the following page https://www.genome.jp/kegg/kaas/ and then select the __Partial Genome__. Select also __GHOSTX__ ans down check __File Upload__ Enter the fasta files downloaded in [4.1.2](#412-pathway-analysis). Afterwards one must give a unique name to the query and enter a valid mail to verify later. For the organisms, one should erase all and for <em>T.cruzi</em>  one should click __Organism List__ and then select kinetoplastids in the pop-up. Finally click in __Compute__  and wait until you get the mail notification. When the first notification comes, one should click it to activate the compute, and a compute can only be done one at a time.

### 4.4.1 Analysis of the KEEG results

What one should analyze are first routes. One should also analyze after selecting __Execute__, the brite and module should also be analyzed.

# 5. Optionals:

## 5.1 Counting with HT-seq
HT-seq is a tool for counting mapped reads to genes in bam files with respect to a reference genome, where the package should be downloaded. The first command is for single reads the second is for paired reads (HT-seq doesn't work with alignments that quantify single and paired reads at the same time so one shoud separate these):

```bash
conda create -n htseq HTSeq

htseq-count -s no -f bam -n 4 ./ERR3501950/accepted_hits.bam.SE.bam TriTrypDB-46_TcruziSylvioX10-1.gtf > SE_50_count.txt

htseq-count -s no -r name -f bam -n 4 ./ERR3501950/accepted_hits.bam.PE.bam TriTrypDB-46_TcruziSylvioX10-1.gtf > PE_50_counts.txt
```
 * To sepate a bam archive in a bam of single reads and a bam of paired reads do:

```bash
conda activate samtools
samtools view -bf 1 foo.bam > foo.paired-end.bam
samtools view -bF 1 foo.bam > foo.single-end.bam
```

## 5.2 Generating a gtf from a gff
```bash
conda create -n agat -c bioconda agat

agat_convert_sp_gff2gtf.pl --gff TriTrypDB-46_TcruziSylvioX10-1.gff -o 1_TriTrypDB-46_TcruziSylvioX10-1.gtf
```
## 5.3 Changing the formats of some reads with BioPython
one installs biopython and then:
```python
from Bio import SeqIO
count = SeqIO.convert("G5ABAL401.RCl3.sff", "sff", "RCl3.fastq", "fastq")
```

## 5.4 to change permisiion in files

To make everything available inside a folder use:

```bash
sudo chmod -R a+rwX file/
```
## 5.5 Transcriptome de novo reconstruction

* https://www.rna.uni-jena.de/compared-the-best-de-novo-transcriptome-assembly-tools/

Spades is very good for assembly so i would use that option
# 6. Footnotes

[^1] : In some cases the we do a De novo transcriptome ensemble in which there are other steps.<br>
[^2] : To obtain the transcriptome and the  gtf from <em>T. cruzi</em> we select in data, then Download data, current release yand search frothe species we're going to analyze in https://tritrypdb.org/tritrypdb/app <br>
[^3] for more information about this see: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-matrix-input <br>
[^4] for more information about this see : https://www.biostars.org/p/293212/ <br>
[^5] A visual tutorial about this is avalilable in: https://www.youtube.com/watch?v=npgkkychkrI&list=LLUISwaoqeJZF8tvTzj2aX9w&index=1 <br>
[^6] A visual tutorial about this is avalilable in: https://www.youtube.com/watch?v=qGLkAEzQpqU&list=LLUISwaoqeJZF8tvTzj2aX9w&index=3&t=0s<br>
[^7]for more information about this discussion see: https://support.bioconductor.org/p/81094/<br>
[^8] For a compressed guide of DESeq2 there is: https://bioconductor.github.io/BiocWorkshops/rna-seq-data-analysis-with-deseq2.html<br>
[^9] Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., ... & Pachter, L. (2012). Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks. Nature protocols, 7(3), 562-578.
[^10]for colors https://www.colorhexa.com/f68060

# 7. References
[1] Williams, C. R., Baccarella, A., Parrish, J. Z., & Kim, C. C. (2017). Empirical assessment of analysis workflows for differential expression analysis of human samples using RNA-Seq. BMC bioinformatics, 18(1), 38.<br>
[2] Zhu, A., Ibrahim, J. G., & Love, M. I. (2019). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics, 35(12), 2084-2092.<br>
[3]https://cgatoxford.wordpress.com/2016/08/17/why-you-should-stop-using-featurecounts-htseq-or-cufflinks2-and-start-using-kallisto-salmon-or-sailfish/<br>
[4] Yi, L., Liu, L., Melsted, P., & Pachter, L. (2018). A direct comparison of genome alignment and transcriptome pseudoalignment. BioRxiv, 444620.<br>
-Anders, S., Pyl, P. T., & Huber, W. (2015). HTSeq—a Python framework to work with high-throughput sequencing data. Bioinformatics, 31(2), 166-169.<br>
-Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.<br>
-Aslett, M., Aurrecoechea, C., Berriman, M., Brestelli, J., Brunk, B. P., Carrington, M., ... & Gardner, M. J. (2010). TriTrypDB: a functional genomic resource for the Trypanosomatidae. Nucleic acids research, 38(suppl_1), D457-D462.<br>
-Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.<br>
Brian Bushnel. (2016). HISAT2 or Tophat2. Retrieved from https://www.biostars.org/p/230724/<br>
-Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & De Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.<br>
-Cock, P. J., Fields, C. J., Goto, N., Heuer, M. L., & Rice, P. M. (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic acids research, 38(6), 1767-1771.<br>
-Cruz-Saavedra, L., Muñoz, M., Patiño, L. H., Vallejo, G. A., Guhl, F., & Ramírez, J. D. (2020). Slight temperature changes cause rapid transcriptomic responses in Trypanosoma cruzi metacyclic trypomastigotes. Parasites & vectors, 13, 1-16.<br>
-Garber, M., Grabherr, M. G., Guttman, M., & Trapnell, C. (2011). Computational methods for transcriptome annotation and quantification using RNA-seq. Nature methods, 8(6), 469-477.<br>
-Kim, D., Langmead, B., & Salzberg, S. L. (2015). HISAT: a fast spliced aligner with low memory requirements. Nature methods, 12(4), 357-360.<br>
-Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome biology, 15(2), R29.<br>
-Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.<br>
-Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550.<br>
-Ma’ayan Avii (2014). The Fisher Exact test, Network Analysis in Systems Biology . Coursera BC : Mount Sinai Center for Bioinformatics. <br>
-Moriya, Y., Itoh, M., Okuda, S., Yoshizawa, A. C., & Kanehisa, M. (2007). KAAS: an automatic genome annotation and pathway reconstruction server. Nucleic acids research, 35(suppl_2), W182-W185.<br>
-Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), 417-419.<br>
-Team, R. C. (2013). R: A language and environment for statistical computing.<br>
-Trapnell, C., Hendrickson, D. G., Sauvageau, M., Goff, L., Rinn, J. L., & Pachter, L. (2013). Differential analysis of gene regulation at transcript resolution with RNA-seq. Nature biotechnology, 31(1), 46-53.


