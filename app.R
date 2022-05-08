## Author: Ryan Yordanoff
## BU BF591
## Final Project

library(shiny)
library(tidyverse)
library(DT)
library(colourpicker)
options(shiny.maxRequestSize=30*1024^2)

#++++++++++++++++++++++++++++++++UI++++++++++++++++++++++++++++++++++++++++++++
# Define UI
ui <- fluidPage(
  
  mainPanel(
    tabsetPanel(
      #---------------------------Tab 1: Samples------------------------------------
      # Tab 1: Samples
      tabPanel("Samples",
               
               # Sidebar layout
               sidebarLayout(
                 sidebarPanel(
                   
                   # Input: Select file
                   fileInput(inputId='fileuploadt1',
                             label= 'Load differential expression results:',
                             placeholder = 'DE_results.csv'),
                   
                 ),
                 
                 
                 mainPanel(
                   # Tabs
                   tabsetPanel(
                     # Tab 1: Summary
                     tabPanel("Summary",
                        tableOutput('table_t1')),
                     # Tab 2: Table
                     tabPanel("Table",
                        tableOutput('table_t1_2')),
                     # Tab 3: Plots
                     tabPanel("Plots",
                        plotOutput('hist_plots')),
                   )
                 ),
                 
                 
               )
      ),
      
      #---------------------------Tab 2: Counts------------------------------------
      # Tab 2: Counts
      tabPanel("Counts",
               
               # Sidebar layout
               sidebarLayout(
                 sidebarPanel(
                   
                   # Input: Select file
                   fileInput(inputId='',
                             label= 'Load differential expression results:',
                             placeholder = 'results.csv'),
                   
                   # Radio button 1: 
                   radioButtons(inputId = '',
                                label = 'Choose x-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'log2FoldChange'),
                   
                   # Radio button 2: 
                   radioButtons(inputId = '',
                                label = 'Choose y-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'padj'),
                   
                   # Slider: 
                   sliderInput(inputId = '',
                               label = 'Select the magnitude of the p-adjusted coloring:',
                               min = -300,
                               max = 0,
                               value = -10),
                   
                   # Plot Button
                   actionButton(inputId = "",
                                label = "Plot",
                                icon = icon('chart', class="fa fa-area-chart") )
                   
                 ),
                 
                 
                 
                 mainPanel(
                   # Tabs
                   tabsetPanel(
                     # Tab 1: Samples
                     tabPanel("1",
                     ),
                     # Tab 2: Counts
                     tabPanel("2",
                     ),
                   )
                 ),
                 
                 
               )
               
               
      ),
      
      #---------------------------Tab 3: Differential Expression-------------------
      # Tab 3: Differential Expression
      tabPanel("DE",
               
               
               # Sidebar layout
               sidebarLayout(
                 sidebarPanel(
                   
                   # Input: Select file
                   fileInput(inputId='fileuploadt3',
                             label= 'Load differential expression results:',
                             placeholder = 'results.csv'),
                   
                   # Radio button 1: x-axis
                   radioButtons(inputId = 'xaxis',
                                label = 'Choose x-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'log2FoldChange'),
                   
                   # Radio button 2: y-axis
                   radioButtons(inputId = 'yaxis',
                                label = 'Choose y-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'padj'),
                   
                   # Color Selection button 1: Base Color
                   colourInput(inputId = 'basecolor',
                               label = 'Base Point Color:',
                               value = "#C6C6DE",
                               showColour = "both",
                               palette = 'square',
                               closeOnClick = 'TRUE'),
                   
                   # Color Selection button 2: Above Threshold Color
                   colourInput(inputId = 'threshcolor',
                               label = 'Above Threshold Color:',
                               value = "#5A50A1",
                               showColour = "both",
                               palette = 'square',
                               closeOnClick = 'TRUE'),
                   
                   # Slider: p-adjusted threshold
                   sliderInput(inputId = 'sliderp',
                               label = 'Select the magnitude of the p-adjusted coloring:',
                               min = -300,
                               max = 0,
                               value = -10),
                   
                   # Plot Button
                   actionButton(inputId = "plotButton",
                                label = "Plot",
                                icon = icon('chart', class="fa fa-area-chart") )
                   
                 ),
                 
                 
                 
                 mainPanel(
                   # Tabs
                   tabsetPanel(
                     # Tab 1: Samples
                     tabPanel("Plot",
                              plotOutput('volcano')),
                     # Tab 2: Counts
                     tabPanel("Data Table",
                              tableOutput('table')),
                   )
                 ),
                 
                 
               )
      ),
      
      #---------------------------Tab 4: Gene Set Enrichment Analysis-------------------
      # Tab 4: Gene Set Enrichment Analysis
      tabPanel("GSEA",
               
               
               # Sidebar layout
               sidebarLayout(
                 sidebarPanel(
                   
                   # Input: Select file
                   fileInput(inputId='',
                             label= 'Load differential expression results:',
                             placeholder = 'results.csv'),
                   
                   # Radio button 1: 
                   radioButtons(inputId = '',
                                label = 'Choose x-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'log2FoldChange'),
                   
                   # Radio button 2: 
                   radioButtons(inputId = '',
                                label = 'Choose y-axis variable:',
                                choices = c('baseMean','log2FoldChange','lfcSE',
                                            'stat', 'pvalue','padj'),
                                selected = 'padj'),
                   
                   # Slider: 
                   sliderInput(inputId = '',
                               label = 'Select the magnitude of the p-adjusted coloring:',
                               min = -300,
                               max = 0,
                               value = -10),
                   
                   # Plot Button
                   actionButton(inputId = "",
                                label = "Plot",
                                icon = icon('chart', class="fa fa-area-chart") )
                   
                 ),
                 
                 
                 
                 mainPanel(
                   # Tabs
                   tabsetPanel(
                     # Tab 1: Normalized Enrichment Score
                     tabPanel("NES",
                              
                        sidebarPanel(
                          # Slider: Slider to adjust number of top pathways to plot by adjusted p-value
                          sliderInput(inputId = '',
                                      label = 'Select the magnitude of the p-adjusted coloring:',
                                      min = -300,
                                      max = 0,
                                      value = -10),
                          
                        ),
                        mainPanel(
                          
                        ),
                     ),
                     # Tab 2: Results Table
                     tabPanel("Results Table",
                              
                        sidebarPanel(
                          # Radio button 1: 
                          radioButtons(inputId = '',
                                       label = 'Choose x-axis variable:',
                                       choices = c('baseMean','log2FoldChange','lfcSE',
                                                   'stat', 'pvalue','padj'),
                                       selected = 'log2FoldChange'),
                          
                          # Radio button 2: 
                          radioButtons(inputId = '',
                                       label = 'Choose y-axis variable:',
                                       choices = c('baseMean','log2FoldChange','lfcSE',
                                                   'stat', 'pvalue','padj'),
                                       selected = 'padj'),
                          
                          # Slider: Slider to adjust number of top pathways to plot by adjusted p-value
                          sliderInput(inputId = '',
                                      label = 'Select the magnitude of the p-adjusted coloring:',
                                      min = -300,
                                      max = 0,
                                      value = -10),
                          
                          # Download Button: Export current filtered and displayed table results
                          downloadButton(''),
                          
                        ),
                        mainPanel(
                          
                        ),
                     ),
                     # Tab 3: NES Scatter Plot
                     tabPanel("NES Scatter Plot",
                              
                      sidebarPanel(
                        # Slider: Slider to adjust number of top pathways to plot by adjusted p-value
                        sliderInput(inputId = '',
                                    label = 'Select the magnitude of the p-adjusted coloring:',
                                    min = -300,
                                    max = 0,
                                    value = -10),
                        
                      ),
                      mainPanel(
                        
                      ),
                     ),
                   )
                 ),
               )
      ),
    )
  )
)



#++++++++++++++++++++++++++++++++SERVER+++++++++++++++++++++++++++++++++++++++++
# Define server logic
server <- function(input, output, session) {
  
  
  
  #---------------------------Tab 1: Samples-------------------------------------
  load_data_t1 <- reactive({
    DE_results_summary <- read.csv(input$fileuploadt1$datapath)
    return(DE_results_summary)
  })
  
  
  
  
  summary_table <- function(data) {
    sample_info <- data
    df <- data.frame(sapply(sample_info, class)) %>%
      as_tibble(rownames='Column Name') %>%
      rename(Type = 2)
    
    df2 <- data.frame(sapply(sample_info, function(x) str_glue('{round(mean(x),2)} (+/- {round(sd(x), 2)}))'))) %>%
      as_tibble(rownames='Column Name') %>%
      rename(`Mean (sd)` = 2)
    
    df3 <- data.frame(sapply(sample_info, function(x) str_glue('{length(unique(x))} Distinct Values')))%>%
      as_tibble(rownames='Column Name') %>%
      rename(`Distinct Values` = 2)
    
    output_df <- cbind(df, (df2[2]), (df3[2])) %>%
      mutate(`Mean (sd)` = na_if(`Mean (sd)`, "NA (+/- NA))"), 
             `Mean (sd)` = coalesce(`Mean (sd)`, `Distinct Values`)) %>%
      select(-`Distinct Values`) %>%
      rename(`Mean (sd) or Distinct Values` = `Mean (sd)`)
    
    
    return(output_df)
  }
  
  multi_hist <- function(data) {
    return(data %>% select_if(is.numeric) %>%
             hist.data.frame())
  }
  
  
  
  
  
  output$table_t1 <- renderTable({req(input$fileuploadt1)
    summary_table(load_data_t1())
  })
  
  output$table_t1_2 <- renderTable({req(input$fileuploadt1)
    load_data_t1()
  })
  
  output$hist_plots <- renderPlot({req(input$fileuploadt1)
    multi_hist(load_data_t1())
  })
  
  
  
  
  
  #---------------------------Tab 2: Counts--------------------------------------
  
  #---------------------------Tab 3: Differential Expression----------------------  
  #' load_Data
  #'
  #' @details Ok
  load_data_t3 <- reactive({
    DE_results <- read.csv(input$fileuploadt3$datapath) %>%
      rename(gene = 1)
    return(DE_results)
  })
  
  
  #' Volcano plot
  #'
  #' @details I
  #'
  #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      
      p <- dataf %>%
        ggplot(aes()) +
        geom_point(mapping = aes(x=!!sym(x_name), y=-log10(!!sym(y_name)),color=padj<(1*10^(slider)))) +
        theme_light() +
        labs(x=x_name,
             y=str_glue('-log10({y_name})'), color=str_glue('{y_name} < 1 x 10^{slider}')) +
        scale_color_manual(values=c(color1, color2)) +
        theme(legend.position="bottom")
      
      return(p)
    }
  
  #' Draw and filter table
  #' @details Same
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_table <- function(dataf, slider) {
    filtered <- filter(dataf, padj < (1*10^slider)) %>%
      mutate(pvalue = formatC(.$pvalue, digits=2, format='e')) %>%
      mutate(padj = formatC(.$padj, digits=2, format='e'))
    return(filtered)
  }
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  output$volcano <- renderPlot({req(input$fileuploadt3)
    volcano_plot(load_data_t3(), input$xaxis, input$yaxis, input$sliderp, input$basecolor, input$threshcolor)
  })
  
  # Same here, just return the table as you want to see it in the web page
  output$table <- renderTable({req(input$fileuploadt3)
    draw_table(load_data_t3(), input$sliderp)
  })
  
  
  
#---------------------------Tab 4: Differential Expression----------------------
  
  
  
  
#---------------------------------------------------------------------------------  
}






# Run the application
shinyApp(ui = ui, server = server)