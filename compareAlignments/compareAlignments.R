source("https://bioconductor.org/biocLite.R")
biocLite('tximportData')
biocLite('tximport')
biocLite('EnsDb.Hsapiens.v86')

library(dplyr)
library(readr)
library(tximportData)
library(tximport)
library(EnsDb.Hsapiens.v86)
library(gtools)
library(rnaseqcomp)

#sample_metadata.csv is a metadata table containing the sample names, batch, ID, etc.
samples.meta <- read.csv('sample_metadata.csv')

#tx2geneFromGTF returns a DataFrame of S4Vectors
#This is needed for tximport to convert transcripts to matching genes
tx2geneFromGTF <- function(gtf.file){

  ensDbFromGtf(gtf=gtf.file)
  gs <- gsub(pattern = '.gtf.gz',replacement = '', gtf.file)
  sqlite <- cat(gs, '.sqlite',sep = '')
  edb <- EnsDb(sqlite)

  Tx.ensemble <- transcripts(edb, columns = c("tx_id", "gene_id", "gene_name"), return.type = "DataFrame")
  tx2gene <- Tx.ensemble[,c(1,2)]
  return(tx2gene)

}

#Generates a table of consolidated, gene-level TPM output from an alignment tool: 1 gene name column, 1 gene length column, and 1 TPM column per sample.
#Example Case: Kallisto is used to generate transcript-level quantification from 30 samples, yielding 30 output files.  This function will take the transcript-level TPM
#values from each sample, convert them to gene-level (using tximport), and consolidate them into a single object.

generateTPMTable <- function(             
  sample.dir,                                 #Directory containing directories of output data from alignment. Ex: "./kallisto-output/sample1/abundance.tsv, "./kallisto-output/sample1/abundance.tsv"
  align.tool=NULL,                            #Alignment tool used. Ex: rsem, kallisto, salmon
  tx2gene.map,                                #tx2gene object mapping transcripts to genes 
  align.output.filename = NULL                #Specification of output file name (ex: 'abundance.tsv'), if align.tool is NULL
  )
  
  {
  file.list <- list.dirs(path = sample.dir, full.names=TRUE, recursive = FALSE)
  file.list <- mixedsort(file.list)
  
  #This should probably be a dictionary, but I'm keeping it this way for now in case of later alignment tools having more complex directory structures.
  
  if(align.tool == "rsem"){
    file.list <- file.path(file.list, "Quant.genes.results")
  }
  if(align.tool == "kallisto"){
    file.list <- file.path(file.list, "abundance.tsv")    
  }
  if(align.tool == "salmon"){
    file.list <- file.path(file.list, "quant.sf")
  }
  if(!is.null(align.output.filename)){
    file.list <- file.path(file.list, align.out.filename)
  }
  
  tx.output <- tximport(file.list, type = align.tool, tx2gene = tx2gene.map, reader = read_tsv, ignoreTxVersion = TRUE)
  counts <- tx.output$counts
  counts <- as.data.frame(counts)
  
  name.list <- list.dirs(path=sample.dir, full.names=FALSE, recursive=FALSE)
  name.list <- mixedsort(name.list)
  names(counts) <- paste0(align.tool,"_", name.list)
  
  counts <- cbind(gene_length = tx.output$length[,1], counts)
  counts <- cbind(gene_name = rownames(counts), counts)
  return(counts)
}

rsem.counts <- generateTPMTable('../AlignmentOutput/rsem_output/', align.tool = 'rsem', tx2gene)
kallisto.counts <- generateTPMTable('../AlignmentOutput/kallisto_output/', align.tool = 'kallisto', tx2gene)
salmon.counts <- generateTPMTable('../AlignmentOutput/salmon_output/', align.tool = 'salmon', tx2gene)

#Joining count tables to limit dataset to shared genes
rsem.kallisto.counts <- inner_join(rsem.counts, kallisto.counts, by = "gene_name")
rsem.kallisto.salmon.counts <- inner_join(rsem.kallisto.counts, salmon.counts, by = "gene_name")

#Retrieving list of housekeeping genes from simdata dataset
#Subset from: https://www.ncbi.nlm.nih.gov/pubmed/23810203
#Housekeeping genes will be used to normalize comparisons between alignment tools
data(simdata)
df.simdata <- as.data.frame(simdata)

housek <- subset(df.simdata, meta.house == TRUE)
housek <- housek[,c("meta.gene", "meta.house")]
housek <- housek %>% mutate(meta.gene = gsub("\\.[0-9]+","", meta.gene))
names(housek) <- c("gene_name", "house")
house.rsem.kallisto.salmon.counts <- unique(left_join(rsem.kallisto.salmon.counts, housek))

#true.fc <- df.simdata
#true.fc <- true.fc[,c("meta.gene", "meta.positive", "meta.fcsign", "meta.fcstatus")]
#true.fc <- true.fc %>% mutate(meta.gene = gsub("\\.[0-9]+","", meta.gene))

true.fc <- simdata$meta
true.fc <- true.fc[,c("gene", "positive", "fcsign", "fcstatus")]
true.fc <- true.fc %>% mutate(gene = gsub("\\.[0-9]+","", gene))
names(true.fc) <- c("tx", "gene_name", "positive", "fcsign", "fcstatus")

true.house.align.counts <- unique(left_join(true.fc, house.rsem.kallisto.salmon.counts))

#Input into signalCalibrate requires a list of matrices of equal size, 1 for each tool
rsem.mat <- data.matrix(true.house.align.counts[,c(3:32)], rownames.force=TRUE)
rownames(rsem.mat) <- true.house.align.counts$gene_name

kallisto.mat <- data.matrix(true.house.align.counts[,c(34:63)], rownames.force=TRUE)
rownames(kallisto.mat) <- true.house.align.counts$gene_name

salmon.mat <- data.matrix(true.house.align.counts[,c(65:94)], rownames.force=TRUE)
rownames(salmon.mat) <- true.house.align.counts$gene_name

mat.list <- list(rsem.mat, kallisto.mat, salmon.mat)
names(mat.list) <- c("rsem", "kallisto", "salmon")

cellInfo <- factor(samples.meta$cell)
batInfo <- factor(samples.meta$batch)

#evaluationFeature = features being evaluated (in this case all)
#calibrationFeature = features used for normalization between samples (in this case housekeeping genes)

evaluationFeature <- rep(TRUE, nrow(true.house.align.counts))
calibrationFeature <- true.house.align.counts$house
unitReference = 1

dat <- signalCalibrate(mat.list, cellInfo, batInfo, evaluationFeature, calibrationFeature, unitReference, calibrationFeature2 = calibrationFeature)

plotSD(dat, ylim=c(0,0.5))

plotNE(dat,xlim=c(0.5,1))

plotROC(dat,simdata$meta$positive,simdata$meta$fcsign,ylim=c(0,0.8))

#6/5/2017:  For ROC, needed the 'true' expression status of genes (for FPR, etc).  Confused how this was derived for the demo set,
#6/5/2017: so I used the demo set, which had limited overlap (~3000).  The plots appear to work, but the differences are very minute,
#6/5/2017: probably because of very high sample count and low feature count.  
#6/5/2017: To-do: Puzzle out how they derived 'true' expression status and repeat it for a larger set of our genes.

