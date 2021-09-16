library(shiny)
library(stringr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(S4Vectors)
library(parallel)
library(BiocGenerics)
library(stats4)

###### manually change this
load("app_Input.RData")
######


ui <- fluidPage(

  titlePanel(project_title),
  
  sidebarLayout(
    
    sidebarPanel(
      selectizeInput('selected_gene', 'Choose gene:', choices=NULL, multiple=TRUE),
      textOutput('padj_explanation'),
      br(),
      tableOutput("pval_legend"),
      width=3),
    
    mainPanel(
      plotOutput('ggplot_figure')
    )
    
  )
)


server <- function(input, output, session) {
  updateSelectizeInput(session, 'selected_gene', choices=interaction$gene_name, server=TRUE)
  
  factor_names <- colnames(coldata)
  condition_names <- levels(coldata[, factor_names[1]])
  
  genes_of_interest <- reactive(input$selected_gene)
  interaction_data <- reactive(interaction[interaction$gene_name %in% genes_of_interest(), ])
  geneIDs_of_interest <- reactive(row.names(interaction_data()))
  
  selectedData <- reactive({ 
    validate(need(genes_of_interest(), "Please select a gene"))
    
    gene_ID <- geneIDs_of_interest() 
    exp_test_gene <- norm_counts[row.names(norm_counts) %in% gene_ID, ]
    exp_test_gene <- t(as.data.frame(exp_test_gene))
    
    if (all(row.names(coldata) == row.names(exp_test_gene))) {
      
      exp_test_gene <- reshape2::melt(cbind(exp_test_gene, coldata))
      
    } else {
      
      stop("Samples in coldata do not match samples in Normalized counts")
      
    }
    subset_difflist1 <- difflist1[unique(exp_test_gene$variable), ]
    subset_difflist2 <- difflist2[unique(exp_test_gene$variable), ]
    exp_test_gene[, pvalue_or_padj] = rep(NA, nrow(exp_test_gene))
    
    for(g in unique(exp_test_gene$variable)) {
      exp_test_gene[exp_test_gene[, factor_names[1]] == condition_names[1] & exp_test_gene$variable == g, ][, pvalue_or_padj] <- makeStar(subset_difflist1[g, ][, pvalue_or_padj])
      exp_test_gene[exp_test_gene[, factor_names[1]] == condition_names[2] & exp_test_gene$variable == g, ][, pvalue_or_padj] <- makeStar(subset_difflist2[g, ][, pvalue_or_padj])
    }
    
    exp_test_gene
    
  })
  
  output$pval_legend <- renderTable({
    validate(need(genes_of_interest(), " "))
    df <- data.frame(star=c('***', '**', '*', 'N.S.'),
                     range=c('0 - 0.001', '0.001 - 0.01', '0.01 - 0.05', '0.05 - 1.0'))
    colnames(df) <- c("Symbol", glue::glue("Range of {pvalue_or_padj}"))
    df
  })
  
  indiv_pvalue_df <- reactive({
    df <- data.frame(
      Genotype=selectedData()[,factor_names[1]],
      Treatment=selectedData()[,factor_names[2]],
      variable=selectedData()[,'variable'],
      label=selectedData()[, pvalue_or_padj])
    
    return(df)
    
  })
  
  makeStar <- function(pval) {
    star <- c('***', '**', '*', 'N.S.')
    breaks <- c(0, 0.001, 0.01, 0.05, 1)
    interval <- findInterval(pval, breaks)
    star[interval]
  }
  
  
  gene_ids <- reactive({
    gene.name <- paste0(interaction_data()$gene_name, " (", signif(interaction_data()[, pvalue_or_padj], 3), ")")
    names(gene.name) <- geneIDs_of_interest()
    gene.name
  })
  
  plotHeight <- reactive({
    number_of_genes <- length(gene_ids())
    return(number_of_genes * 300)
  })

  output$padj_explanation <- reactive({
    validate(need(genes_of_interest(), " "))
    glue::glue("Note - Number next to the gene represents the {pvalue_or_padj} of the gene's difference between {condition_names[1]} and {condition_names[2]}")
  })
  
  output$ggplot_figure <- renderPlot(
    
    ggplot(selectedData(), aes_string(x=factor_names[2], y="value", 
                                      color=factor_names[2])) + geom_point() +
    geom_point(size=3) +
      theme_linedraw(base_size = 16) + 
      facet_grid( reformulate( factor_names[1], "variable"), 
                  labeller=labeller(variable=gene_ids()),
                  scales="free") + 
      ylab("Normalized Expression") +
      stat_summary(
        fun = median,
        geom = "line",
        color="black",
        aes(group=factor_names[1]), 
        position = position_dodge(width=0.9)) + 
      geom_text(
        data=indiv_pvalue_df() ,
        mapping = aes(x=1.5, y=Inf, label=label), 
        hjust   = 0.5, vjust= 2,
        color = 'black',
        size=4.5
      ),
    height = plotHeight, res=96 )
  
} 


shinyApp(ui, server) 
