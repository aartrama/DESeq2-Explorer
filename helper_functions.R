GetAnnotationFromBiomart <- function(species) {
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

AnnotateDifflist <- function(annotations, difflist) {
  annotations <- annotations[rownames(difflist), ]
  print(all(rownames(annotations) == rownames(difflist)))
  annot_difflist <- cbind(difflist, annotations)
  annot_difflist <- subset(annot_difflist, select=c(7,1,2,3,4,5,6,8,9,10,11,12))
  colnames(annot_difflist) <- c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", 
                                "padj", "chromosome", "gene_start",
                                "gene_end", "gene_type", "description")
  return(annot_difflist)
}

ExtractDifferentialListsPairwise <- function(filename, groupname, group2, group1) {
  difflist <- results(dds, contrast=c(groupname, group2, group1))
  difflist <- AnnotateDifflist(annotation, difflist) 
  difflist <- difflist[order(difflist$pvalue), ]
  difflist$gene_id <- rownames(difflist)
  write.table(difflist, file = filename, sep = "\t", quote = F, row.names=F)
}

ExtractDifferentialListsPairwiseExcel <- function(filename, groupname, group2, group1) {
  difflist <- results(dds, contrast=c(groupname, group2, group1))
  difflist <- AnnotateDifflist(annotation, difflist) 
  difflist <- difflist[order(difflist$pvalue), ]
  difflist$gene_id <- rownames(difflist)
  write.xlsx(difflist, filename)
}

ExtractDifferentialListsMultiFactor <- function(filename, var1, var2=NULL) {
  if (is.null(var2)) {
    difflist <- results(dds, name=var1)
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    write.table(difflist, file = filename, sep = "\t", quote = F, row.names=F)
    return(difflist)
    
  } else {
    difflist <- results(dds, contrast = list(c(var1, var2)))
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    write.table(difflist, file = filename, sep = "\t", quote = F, row.names=F)
    return(difflist)
  }
}

ExtractDifferentialListsMultiFactorExcel <- function(filename, var1, var2=NULL) {
  if (is.null(var2)) {
    difflist <- results(dds, name=var1)
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    write.xlsx(difflist, filename)
    
  } else {
    difflist <- results(dds, contrast = list(c(var1, var2)))
    difflist <- AnnotateDifflist(annotation, difflist) 
    difflist <- difflist[order(difflist$pvalue), ]
    difflist$gene_id <- rownames(difflist)
    write.xlsx(difflist, filename)
  }
}

GeneratePCA <- function(filename, group1, group2=NULL) {
  if (is.null(group2)) {
    pcaData <- plotPCA(vsd.fast, intgroup=c(group1), returnData=TRUE)
  } else {
    pcaData <- plotPCA(vsd.fast, intgroup=c(group1, group2), returnData=TRUE)
  }
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  if (is.null(group2)) {
    plot1 <- ggplot(pcaData, aes(PC1, PC2, color=!!ensym(group1)))
  } else {
    plot1 <- ggplot(pcaData, aes(PC1, PC2, color=!!ensym(group1), shape=!!ensym(group2))) 
  }
  pdf = pdf(filename)
  plot1 <- plot1 + geom_point(size=3) + geom_text(aes(label = name), color="black", size = 2, 
                                                  position = position_nudge(y = 0.2))+
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  print(plot1)
  dev.off()
  
}

ObtainNormalizedCounts <- function() {
  data_normalized = counts(dds, normalized = T)
  return(data.frame(data_normalized))
}

ObtainRPKMs <- function(species) {
  annotation <- GetAnnotationFromBiomart(species)
  # RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
  annotation$gene_length <- annotation$end_position - annotation$start_position + 1 
  # this is how FeatureCounts calculates gene lengths
  annotation <- subset(annotation, select=c("gene_length", "external_gene_name"))
  mcols(dds)$basepairs <- annotation$gene_length
  fpkm_values <- fpkm(dds, robust=F) # robust = F: use colsums for lib size norm 
  # (as opposed to size factors calculated by DESEQ2)
  fpkm_values <- as.data.frame(fpkm_values) 
  Gene.name <- annotation$external_gene_name
  fpkm_values <- cbind(Gene.name=annotation$external_gene_name, fpkm_values)
  write.table(fpkm_values, file = "RPKM_values.txt", sep = "\t", quote = F, col.names = NA)
}

ObtainPvalueLogfcDESeq2 <- function(difflist) {
  compname <- strsplit(difflist, ".txt")
  data1 <- read.delim(difflist, sep = "\t", header =T, stringsAsFactors = F)
  rownames(data1) <- paste0(data1$gene_name, ",", data1$gene_id)
  data1.subset <- cbind(data1$log2FoldChange, data1$pvalue, data1$padj)
  rownames(data1.subset) <- rownames(data1)
  colnames(data1.subset) <- c("log2FoldChange", "pvalue", "padj")
  colnames(data1.subset) <- c(paste0(compname, ".log2FoldChange"), paste0(compname, ".pvalue"), paste0(compname, ".padj"))
  data1.subset
}

RunMultiFactor <- function(WT_list, KO_list, interaction_list, outfile_tsv, pval_or_padj, final_outfile_xlsx) {

    use_python("/Users/aarthi/miniconda3/bin/python3", required=T)
    py_config()
    source_python("helper_functions.py")

    list_of_files = c(WT_list, KO_list, interaction_list)

    result <- lapply(list_of_files, ObtainPvalueLogfcDESeq2)

    intersect_genes <- Reduce(intersect, list(rownames(result[[1]]),
                                              rownames(result[[2]]),
                                              rownames(result[[3]])))

    result[[1]] <- result[[1]][intersect_genes,]
    result[[2]] <- result[[2]][intersect_genes,]
    result[[3]] <- result[[3]][intersect_genes,]

    table(row.names(result[[1]])==row.names(result[[2]]))
    table(row.names(result[[2]])==row.names(result[[3]]))

    result_last <- do.call(cbind, result)

    write.table(result_last, file = outfile_tsv, sep = "\t", quote = F, col.names = NA)

    parse_diff_list_table(outfile_tsv, pval_or_padj) # "pvalue" or "padj"
    integrate_up_down_csv_files_into_excel(final_outfile_xlsx)

}