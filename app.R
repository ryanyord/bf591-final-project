## Author: Ryan Yordanoff
## BU BF591
## Final Project

library(shiny)
library(tidyverse)
library(DT)
library(colourpicker)

#++++++++++++++++++++++++++++++++UI++++++++++++++++++++++++++++++++++++++++++++
# Define UI
ui <- fluidPage(
  
  mainPanel(
    tabsetPanel(
#---------------------------Tab 1: Samples------------------------------------
      # Tab 1: Samples
      tabPanel("Samples",
      ),
      
#---------------------------Tab 2: Counts------------------------------------
      # Tab 2: Counts
      tabPanel("Counts",
      ),
      
#---------------------------Tab 3: Differential Expression-------------------
      # Tab 3: Differential Expression
      tabPanel("DE",
               
               
               # Sidebar layout
               sidebarLayout(
                 sidebarPanel(
                   
                   # Input: Select file
                   fileInput(inputId='fileupload',
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
      ),
    )
  )
)



#++++++++++++++++++++++++++++++++SERVER+++++++++++++++++++++++++++++++++++++++++
# Define server logic
server <- function(input, output, session) {
  
  #' load_Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is 
  #' still a "function", but it will take no arguments. The `reactive({})` bit 
  #' says "if any of my inputs (as in, input$...) are changed, run me again". 
  #' This is useful when a user clicks a new button or loads a new file. In 
  #' our case, look for the uploaded file's datapath argument and load it with 
  #' read.csv. Return this data frame in the normal return() style.
  load_data <- reactive({
    DE_results <- read.csv(input$fileupload$datapath) %>%
      rename(gene = 1)
    return(DE_results)
  })
  
#---------------------------Tab 1: Samples-------------------------------------
  
#---------------------------Tab 2: Counts--------------------------------------
  
#---------------------------Tab 3: Differential Expression----------------------  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  #' @details I bet you're tired of these plots by now. Me too, don't worry.
  #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
  #' Write a normal volcano plot using geom_point, and integrate all the above 
  #' values into it as shown in the example app. The testing script will treat 
  #' this as a normal function.
  #' 
  #' !!sym() may be required to access column names in ggplot aes().
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
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will 
  #' evaluate it normally. Not only does this function filter the data frame to 
  #' rows that are above the slider magnitude, it should also change the format 
  #' of the p-value columns to display more digits. This is so that it looks 
  #' better when displayed on the web page. I would suggest the function 
  #' `formatC()`
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
  output$volcano <- renderPlot({req(input$fileupload)
    volcano_plot(load_data(), input$xaxis, input$yaxis, input$sliderp, input$basecolor, input$threshcolor)
  })
  
  # Same here, just return the table as you want to see it in the web page
  output$table <- renderTable({req(input$fileupload)
    draw_table(load_data(), input$sliderp)
  })
  
  
  
#---------------------------Tab 4: Differential Expression----------------------
}






# Run the application
shinyApp(ui = ui, server = server)