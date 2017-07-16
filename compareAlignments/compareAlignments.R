source("https://bioconductor.org/biocLite.R")
biocLite('tximportData')

library(readr)
library(dplyr)
library(tximportData)
library(tximport)
library(EnsDb.Hsapiens.v86)
library(gtools)
library(rnaseqcomp)

#sample_metadata.csv is a metadata table containing the sample names, batch, ID, etc.
samples.meta <- read.csv('compareAlignments/sample_metadata.csv')

samples.samp <- paste0('sample', samples.meta$samp)
#sampList <- c(samples.meta$samp)

#tx2geneFromGTF returns a DataFrame of S4Vectors
#This is needed for tximport to convert transcripts to matching genes
tx2geneFromGTF <- function(gtf){

  ensDbFromGtf(gtf)
  edb <- EnsDb("Homo_sapiens.GRCh38.88.sqlite")

  Tx.ensemble <- transcripts(edb, columns = c("tx_id", "gene_id", "gene_name"), return.type = "DataFrame")
  tx2gene <- Tx.ensemble[,c(1,2)]
  return(tx2gene)


}

tx2gene <- tx2geneFromGTF('Homo_sapiens.GRCh38.88.gtf.gz')

#consolidateSampleCounts generates a consolidated counts data frame from all samples
consolidateSampleCounts <- function(             
  sample.dir,                                 #Directory containing directories of output data from alignment.
  align.tool=NULL,                            #Alignment tool used (ex: rsem, kallisto, salmon) 
  tx2gene.map,                                #tx2gene object mapping transcripts to genes 
  align.output.filename = NULL                #Specification of output file name (ex: 'abundance.tsv'), if align.tool is NULL
  )
  
  {
  file.list <- list.dirs(path = sample.dir, full.names=TRUE, recursive = FALSE)
  file.list <- mixedsort(file.list)
  
  if(align.tool == "rsem"){
    file.list <- file.path(file.list, "Quant.genes.results")
    transcript.level <- TRUE
  }
  if(align.tool == "kallisto"){
    file.list <- file.path(file.list, "abundance.tsv")
    transcript.level <- TRUE
  }
  if(align.tool == "salmon"){
    file.list <- file.path(file.list, "quant.sf")
    transcript.level <- TRUE
  }
  if(align.tool == "stringtie"){
    file.list <- file.path(file.list, "gene_abundances.tsv")
    transcript.level <- FALSE
  }
  if(!is.null(align.output.filename)){
    file.list <- file.path(file.list, align.out.filename)
  }
  
  if(transcript.level == TRUE){
    tx.output <- tximport(file.list, type = align.tool, tx2gene = tx2gene.map, reader = read_tsv, ignoreTxVersion = TRUE)
    counts <- tx.output$counts
    counts <- as.data.frame(counts)
  }
 
  name.list <- list.dirs(path=sample.dir, full.names=FALSE, recursive=FALSE)
  name.list <- mixedsort(name.list)
  names(counts) <- paste0(align.tool,"_", name.list)
  
  counts <- cbind(gene_length = tx.output$length[,1], counts)
  colnames(counts)[1] <- paste0(align.tool, ".", "gene_length")
  
  counts <- cbind(gene_name = rownames(counts), counts)
  return(counts)
}

tx2gene <- tx2geneFromGTF('Homo_sapiens.GRCh38.88.gtf.gz')

rsem.counts <- consolidateSampleCounts('AlignmentOutput/rsem_output/', align.tool = 'rsem', tx2gene)
kallisto.counts <- consolidateSampleCounts('AlignmentOutput/kallisto_output/', align.tool = 'kallisto', tx2gene)
salmon.counts <- consolidateSampleCounts('AlignmentOutput/salmon_output/', align.tool = 'salmon', tx2gene)

#tximport does not currently support stringtie.  It's tricky to add because Stringtie does not directly output counts or effective length,
#so they must be calculated.  The stringtie authors have written a python script that I used to calculate counts and consolidate all samples
#to one table (which is exactly what I would want tximport to do).

stringtie.counts <- read.csv('AlignmentOutput/hisat2_stringtie_output/gene_count_matrix.csv')
colnames(stringtie.counts)[1] <- "gene_name"

#Putting columns in ascending order
for(i in 1:30){
  coluq[i] <- paste0("sample",i)
}
stringtie.counts <- stringtie.counts[c("gene_name", coluq)]

#Joining count tables to limit dataset to shared genes
joined.counts <- join_all(list(rsem.counts, kallisto.counts, salmon.counts, stringtie.counts), by="gene_name", type="inner")

#Retrieving list of housekeeping genes from simdata dataset
#Subset from: https://www.ncbi.nlm.nih.gov/pubmed/23810203
#Housekeeping genes will be used to normalize comparisons between alignment tools
data(simdata)
df.simdata <- as.data.frame(simdata)

housek <- subset(df.simdata, meta.house == TRUE)
housek <- housek[,c("meta.gene", "meta.house")]
housek <- housek %>% mutate(meta.gene = gsub("\\.[0-9]+","", meta.gene))
names(housek) <- c("gene_name", "house")
house.joined.counts <- unique(left_join(joined.counts, housek))

#Input into signalCalibrate requires a list of matrices of equal size, 1 for each tool
rsem.mat <- data.matrix(house.joined.counts[,c(3:32)], rownames.force=TRUE)
rownames(rsem.mat) <- house.joined.counts$gene_name

kallisto.mat <- data.matrix(house.joined.counts[,c(34:63)], rownames.force=TRUE)
rownames(kallisto.mat) <- house.joined.counts$gene_name

salmon.mat <- data.matrix(house.joined.counts[,c(65:94)], rownames.force=TRUE)
rownames(salmon.mat) <- house.joined.counts[,1]

stringtie.mat <- data.matrix(house.joined.counts[,c(95:124)], rownames.force=TRUE)
rownames(stringtie.mat) <- house.joined.counts[,1]

mat.list <- list(rsem.mat, kallisto.mat, salmon.mat, stringtie.mat)
names(mat.list) <- c("rsem", "kallisto", "salmon", "stringtie")

cellInfo <- factor(samples.meta$cell)
batInfo <- factor(samples.meta$batch)

#evaluationFeature = features being evaluated (in this case all)
#calibrationFeature = features used for normalization between samples (in this case housekeeping genes)

#house.length.rsem.kallisto.salmon.counts <- left_join(rsem.kallisto.salmon.counts, housek, by = "gene_name")

evaluationFeature <- rep(TRUE, nrow(joined.counts))
calibrationFeature <- house.joined.counts$house
unitReference = 1

dat <- signalCalibrate(mat.list, cellInfo, batInfo, evaluationFeature, calibrationFeature, unitReference, calibrationFeature2 = calibrationFeature)

plotSD(dat, ylim=c(0,1.4))




