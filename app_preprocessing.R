library(shiny)
library(DESeq2)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(plotly)
library(biomaRt)

# required functions
GetAnnotationFromBiomart <- function(species=NULL) {
  if (species == "hg38") {
    dataset="hsapiens_gene_ensembl"
  } else if (species == "mm10") {
    dataset="mmusculus_gene_ensembl"
  } else if (species == "rn6") {
    dataset="rnorvegicus_gene_ensembl"
  }
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=dataset, host="aug2017.archive.ensembl.org")
  annotation_biomart <- getBM(attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position",
                                             "end_position", "gene_biotype", "description"),
                              filters = "ensembl_gene_id", values = rownames(dds), mart = mart)
  rownames(annotation_biomart) <- annotation_biomart$ensembl_gene_id
  annotation_biomart <- subset(annotation_biomart, select=c(-`ensembl_gene_id`))
  return(annotation_biomart)
}

AnnotateDifflist <- function(annotations=NULL, difflist=NULL) {
  annotations <- annotations[rownames(difflist), ]
  print(all(rownames(annotations) == rownames(difflist)))
  annot_difflist <- cbind(difflist, annotations)
  annot_difflist <- subset(annot_difflist, select=c(7,1,2,3,4,5,6,8,9,10,11,12))
  colnames(annot_difflist) <- c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                "padj", "chromosome", "gene_start",
                                "gene_end", "gene_type", "description")
  return(annot_difflist)
}

ExtractDifferentialListsMultiFactor <- function(var1=NULL, var2=NULL) {
  if (is.null(var2)) {
    difflist <- results(dds, name=var1)
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    return(difflist)
    
  } else {
    difflist <- results(dds, contrast = list(c(var1, var2)))
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    return(difflist)
  }
}


ObtainNormalizedCounts <- function() {
  data_normalized = counts(dds, normalized = T)
  data_normalized <- as.data.frame(data_normalized)
  return(data_normalized)
}


counts <- read.table("htseq_counts_matrix.txt", sep="\t", header=T, row.names=1)
coldata <- read.table("metadata_wo_outlier.txt", sep="\t", header=T, row.names=1)

coldata$Genotype <- str_replace(coldata$Genotype, "RGS20_", "")
counts <- counts[rowSums(counts >= 1) >= 2, ]

all(rownames(coldata) %in% colnames(counts))
counts = counts[, rownames(coldata)]
all(rownames(coldata) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ Genotype + Treatment + Genotype:Treatment)

dds$Treatment <- relevel(dds$Treatment, 'Naive')
dds$Genotype <- relevel(dds$Genotype, 'WT')

dds <- DESeq(dds)

annotation <- GetAnnotationFromBiomart("mm10")
species="mm10"

interaction <- ExtractDifferentialListsMultiFactor("GenotypeKO.TreatmentCFA")
norm_counts <- ObtainNormalizedCounts()

coldata$Genotype <- relevel(factor(coldata$Genotype), ref="WT")
coldata$Treatment <- relevel(factor(coldata$Treatment), ref="Naive")

CFA_vs_Naive_in_WT <- ExtractDifferentialListsMultiFactor("Treatment_CFA_vs_Naive")
CFA_vs_Naive_in_KO <- ExtractDifferentialListsMultiFactor("Treatment_CFA_vs_Naive", "GenotypeKO.TreatmentCFA")

save(interaction, norm_counts, coldata, CFA_vs_Naive_in_WT, CFA_vs_Naive_in_KO, file="app_Input.RData")

