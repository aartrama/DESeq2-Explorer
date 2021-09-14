library(shiny)
library(ggplot2)
library(reshape2)

load("app_Input.RData")

ui <- fluidPage(
  
  titlePanel("CFA/Naive RGS20 WT/KO"),
  
  sidebarLayout(
    
    sidebarPanel(
      selectizeInput('selected_gene', 'Choose gene:', choices=NULL, multiple=TRUE),
      textOutput('padj_explanation'),
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
  
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  
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
    exp_test_gene
    print(exp_test_gene)
    
  })
  
  
  gene_ids <- reactive({
    gene.name <- paste0(interaction_data()$gene_name, " (P-adj=", specify_decimal(interaction_data()$padj, 3), ")")
    names(gene.name) <- geneIDs_of_interest()
    gene.name
  })
  
  plotHeight <- reactive({
    number_of_genes <- length(gene_ids())
    return(number_of_genes * 300)
  })
  
  output$padj_explanation <- reactive({
    validate(need(genes_of_interest(), " "))
    glue::glue("* P-adj value represents the signficance of the gene's difference between {condition_names[1]} and {condition_names[2]}")
  })
  
  output$ggplot_figure <- renderPlot({
    
    p <- ggplot(selectedData(), aes(Treatment, value)) + geom_point() 
    p <- p + 
      facet_grid( variable ~ Genotype) + 
      geom_text(
        data=data.frame( # make df reactive
          Genotype=c("WT", "KO","WT", "KO"),
          Treatment=c("CFA", "Naive","CFA", "Naive"),
          variable=c("ENSMUSG00000061808", "ENSMUSG00000061808", "ENSMUSG00000031661", "ENSMUSG00000031661"),
          label=c("2342342", "120", "adfd", "sfdsf")
        ),
        mapping = aes(x=Inf,y=Inf,label=label),
        hjust   = 1.05,
        vjust   = 1.5
      )
    print(p)
    height = plotHeight })
  
} 


shinyApp(ui, server) 
