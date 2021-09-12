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
      
      exp_test_gene <- melt(cbind(exp_test_gene, coldata))
      
    } else {
      
      stop("Samples in coldata do not match samples in Normalized counts")
      
    }
    exp_test_gene
    
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

  output$ggplot_figure <- renderPlot(
    
    ggplot(selectedData(), aes_string(x=factor_names[2], y="value", color=factor_names[2])) +
      geom_point(size=3) + 
      theme_linedraw(base_size = 16) + 
      ylab("Normalized Expression") + 
      facet_grid( reformulate( factor_names[1], "variable"), 
                  labeller = labeller(variable = gene_ids()), 
                  scales="free") + 
      expand_limits(y = 0) +
      stat_summary(
        fun = median,
        geom = "line",
        color="black",
        aes(group=factor_names[1]), position = position_dodge(width=0.9)
      ),
    height = plotHeight, res= 96 )
  
  
}


shinyApp(ui, server) 
