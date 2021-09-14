library(shiny)
library(ggplot2)
library(reshape2)
library(ggprism)

load("app_Input.RData")

# p.values <- c(9.5e-15, 0.02)
# Signif <- symnum(p.values, corr = FALSE, na = FALSE, 
#                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
#                  symbols = c("***", "**", "*", ".", " "))

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
  treatment_names <- levels(coldata[, factor_names[2]])
  
  genes_of_interest <- reactive(input$selected_gene)
  interaction_data <- reactive(interaction[interaction$gene_name %in% genes_of_interest(), ])
  # WT_data <- reactive(CFA_vs_Naive_in_WT[CFA_vs_Naive_in_WT$gene_name %in% genes_of_interest(), ])
  # KO_data <- reactive(CFA_vs_Naive_in_KO[CFA_vs_Naive_in_KO$gene_name %in% genes_of_interest(), ])
  geneIDs_of_interest <- reactive(row.names(interaction_data()))
  
  final_df <- reactive({
    df <- data.frame(
    Treatment = c("CFA"), 
    Pval = c("0.01", "0.01"))
    print(df)
    df
  })
  
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
    print(exp_test_gene)
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
  
  # get_label_positions <- reactive({
  #   genes <- geneIDs_of_interest()
  #   genotype_levels <- condition_names
  #   treatment_levels <- treatment_names
  # 
  #   gene_names_df <- c()
  #   condition_names_df <- c()
  #   y_vals <- c()
  #   
  #   for(items in genes) {
  #     for(items1 in genotype_levels) {
  #       dat_f <- selectedData()
  #       tempdf <- dat_f[dat_f$Genotype == items1 & dat_f$variable == items, ]
  #       print(tempdf)
  #       df1 <- tempdf[tempdf$Treatment == 'Naive', ]
  #       df2 <- tempdf[tempdf$Treatment == 'CFA', ]
  #       gene_names_df <- append(gene_names_df, items)
  #       condition_names_df <- append(condition_names_df, items1)
  #       y_vals <- append(y_vals, ((median(df1$value) + median(df2$value))/2) + 2)
  #     }
  #   }
  #   
  #   final_df <- data.frame(gene_id=gene_names_df,
  #                          Genotype=condition_names_df,
  #                          median_vals=y_vals,
  #                          labels_for_plot=rep("text", length(y_vals)))
  #   
  #   return(final_df)
  #   print(final_df)
  #   
  # })
  
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
        aes(group=factor_names[1]), 
        position = position_dodge(width=0.9)
      ) 
)}


shinyApp(ui, server) 





