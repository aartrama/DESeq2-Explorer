library(shiny)
library(DESeq2)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(plotly)
library(biomaRt)

# required functions
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

ExtractDifferentialListInteraction <- function(termname) {
  difflist <- results(dds, name=termname)
  difflist <- AnnotateDifflist(annotation, difflist)
  difflist <- difflist[order(difflist$pvalue), ]
  difflist$gene_id <- rownames(difflist)
  difflist <- difflist[complete.cases(difflist), ]
  return(difflist)
}

ObtainNormalizedCounts <- function() {
  data_normalized = counts(dds, normalized = T)
  return(data.frame(data_normalized))

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
interaction <- ExtractDifferentialListInteraction("GenotypeKO.TreatmentCFA")
norm_counts <- ObtainNormalizedCounts()
coldata$Genotype <- relevel(factor(coldata$Genotype), ref="WT")
coldata$Treatment <- relevel(factor(coldata$Treatment), ref="Naive")


ui <- fluidPage(
  
  fluidRow(column(4,
                  selectizeInput('chosen_gene', 'Choose gene:', 
                                 choices=NULL,
                                 multiple=TRUE, 
                                 selected=c("Ttr", "Rgs22")))
  ),
  
  fluidRow(
    plotOutput('ggplot_img', height = "500px"),
    tableOutput('table_output')
    
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(session, 'chosen_gene', choices=interaction$gene_name, server=TRUE)
  factor_names <- colnames(coldata)
  genes_of_interest <- reactive(input$chosen_gene)
  interaction_data <- reactive(interaction[interaction$gene_name %in% genes_of_interest(), ])
  geneIDs_of_interest <- reactive(row.names(interaction_data()))
  
  selectedData <- reactive({ 
    validate(need(genes_of_interest(), "Please select a gene"))
    
    gene_ID <- geneIDs_of_interest() 
    exp_test_gene <- norm_counts[row.names(norm_counts) %in% gene_ID, ]
    exp_test_gene <- t(as.data.frame(exp_test_gene))
    
    if (all(row.names(coldata) == row.names(exp_test_gene))) {
      
      exp_test_gene <- melt(cbind(exp_test_gene, coldata))
      
    } else {
      
      stop("Samples in coldata do not match samples in Normalized counts")
      
    }
    exp_test_gene
    
  })
  
  
  gene_ids <- reactive({
    specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
    gene.name <- paste0(genes_of_interest(), " (p-adj=", specify_decimal(interaction_data()$padj, 3), ")")
    names(gene.name) <- geneIDs_of_interest()
    gene.name
  })
  
  plotHeight <- reactive({
    number_of_genes <- length(gene_ids())
    return(number_of_genes * 300)
  })
  
  output$ggplot_img <- renderPlot(
    
    ggplot(selectedData(), aes_string(x=factor_names[2], y="value", color=factor_names[2])) +
      geom_point(size=3) + 
      theme_linedraw(base_size = 16) + 
      ylab("Normalized Expression") + 
      facet_grid( reformulate( factor_names[1], "variable"), 
                  labeller = labeller(variable = gene_ids()), 
                  scales="fixed") + 
      expand_limits(y = 0) +
      theme(legend.position="none") +
      stat_summary(
        fun = median,
        geom = "line",
        color="black",
        aes(group=factor_names[1]), position = position_dodge(width=0.9)
      ),
    height = plotHeight, res= 96 )
  
  
}


shinyApp(ui, server)
