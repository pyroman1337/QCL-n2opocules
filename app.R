###########################################################################
### Soil GreenHouseGas Flux Visualisation and calculation tool          ###
###########################################################################
## Author: Roman Hüppi
## Date : August 2018
## Version: 0.2

# libraries ---------------------------------------------------------------
library(shiny)
library(ggplot2)
library(plotly)
library(data.table)
library(DT)

numericInput3<-function (inputId, label, value = "",...) 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "number", value = value,...))
}
# Define UI for application that draws a histogram
ui <- fluidPage(
  img(src = "ETH_logo.jpg", height = 70, width = 200, align = "right"),
  
   # Application title
  titlePanel("Mixing model simulation for nitrous oxide isotopologue measurements"),
  
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        h4("Specify simulation parameters"),
        sliderInput(inputId = "carlo",
                    label = "Number of MonteCarlo simulations",
                    min = 5000,
                    max = 50000,
                    value = 6000),
        sliderInput(inputId = "Dppb",
                    label = "Simulated concentration range [ppb]",
                    min = 300,
                    max = 12000,
                    value = c(333.5,3000)),

        # Horizontal line ----
        tags$hr(),
        h4("Experimental input:"),
        # sliderInput("conc", "Start and End concentration [ppb]",
        #             min = 325.0, max = 4500.0, value = c(333.5, 4500)),

        numericInput3(inputId="alpha.start", label=HTML("&alpha; start: &emsp;"), value = 15.5, width = 200),
        numericInput3(inputId="alpha.end"  , label=HTML("&emsp;&alpha; end:"), value = -24.35, width = 200),
        
        numericInput3(inputId="beta.start", label=HTML("&beta; start: &emsp;"), value = -2.5),
        numericInput3(inputId="beta.end"  , label=HTML("&emsp;&beta; end:"), value = -22.94),
        
        numericInput3(inputId="d18O.start", label=HTML("&delta;<sup>18</sup>O start:"), value = 44, width = 100),
        numericInput3(inputId="d18O.end"  , label=HTML("&delta;<sup>18</sup>O end:"), value = -31.79),
        # Stephens increase by 700 ppb  (LGR) ; single 5 min average measurements

        # Horizontal line ----
        tags$hr(),
        h4("Instrument precision (standard deviation [\u2030]):"),
        # numericInput3(inputId="conc_start_sd", label="initial bulk SD", value = 0.3),
        # numericInput3(inputId="conc_end_sd"  , label="bulk end SD", value = 0.9),
        
        numericInput3(inputId="alpha_start_sd", label=HTML("&alpha; SD start:&emsp;"), value = 10),
        numericInput3(inputId="alpha_end_sd"  , label=HTML("&emsp;&alpha; SD end:"), value = 4),
        
        numericInput3(inputId="beta_start_sd", label=HTML("&beta; SD start:&emsp;"), value = 10),
        numericInput3(inputId="beta_end_sd"  , label=HTML("&emsp;&beta; SD end:"), value = 6),
        
        numericInput3(inputId="d18O_start_sd", label=HTML("&delta;<sup>18</sup>O SD start:&nbsp;"), value = 12),
        numericInput3(inputId="d18O_end_sd"  , label=HTML("&delta;<sup>18</sup>O SD end:"), value = 10),
        
        numericInput3(inputId="conc_start_sd", label=HTML("[N<sub>2</sub>O] SD start:"), value = 0.3),
        numericInput3(inputId="conc_end_sd"  , label=HTML("[N<sub>2</sub>O] SD end:"), value = 0.9),
        
        # Horizontal line ----
        tags$hr(),
        sliderInput(inputId = "target",
                    label = "Enter your targeted precision level [permille SD]",
                    min = 0.1,
                    max = 10,
                    value = 1),
        actionButton("sim.go", "Start Simulation")
        
      ),
        
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(
          
          tabPanel(p(icon("line-chart"),  "MonteMattiPlot"),
  
                   
           plotlyOutput("precisionSD"), br(),br(),br(),br(),br(),br(),
           verbatimTextOutput("value")
           
          ), # end of "Visualize the Data" tab panel
          tabPanel(p(icon("table"), "Dataset"),
                   # dataTableOutput(outputId="dTable")
                   # verbatimTextOutput("value"),
                   downloadButton("downloadData", "Download"),
                   
                   DT::dataTableOutput("plot.table")
                   
          ) #, end of "Dataset" tab panel
        ) # end tab set panel
      ) # end mainPanel
   ),  # end sidebarLayout
  img(src = "Empa_logo.png", height = 75, align = "right")
  
)    # end fluidPage

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # get.input <- reactive({
  #   
  # })
  simdata <- reactive({
  input$sim.go
  isolate({
  # observeEvent(input$sim.go, {
    
  # set.seed(52)
  # sim.input <- data.table(Dppb = seq(0,250,10), Dsd = SD.init  )# defines simulated increase in concentration (start, end, steps)
  Dppb.max <- input$Dppb[2]
    # Dppb.max <- 5000
  Dppb  <-  seq(0,Dppb.max,Dppb.max/100)
  
  sim  <- data.table(id = rep(Dppb, each = input$carlo) ) # defines number of monte carlo replicates  , Dppb = Dppb,  C = rnorm(2e3, mean = C00, sd = sdGC),
  

  # Stephens increase by 700 ppb  (LGR) ; single 5 min average measurements
  sim[,':=' (alpha_start = rnorm(length(id), mean = input$alpha.start,sd = input$alpha_start_sd),
             alpha_end   = rnorm(length(id), mean = input$alpha.end , sd = input$alpha_end_sd),    # theoretische stabw von Messgeräten
             beta_start  = rnorm(length(id), mean = input$beta.start, sd = input$beta_start_sd),
             beta_end    = rnorm(length(id), mean = input$beta.end  , sd = input$beta_end_sd),
             d18O_start  = rnorm(length(id), mean = input$d18O.start, sd = input$d18O_start_sd),
             d18O_end    = rnorm(length(id), mean = input$d18O.end  , sd = input$d18O_end_sd),
             conc_start  = rnorm(length(id), mean = input$Dppb[1],    sd = input$conc_start_sd),
             conc_end    = rnorm(length(id), mean = input$Dppb[1]+id, sd = input$conc_end_sd)
  )]
  
# 
#   sim[,':=' (alpha_end = rnorm(length(id), mean = -24   , sd = 4),   alpha_start = rnorm(length(id), mean = 15.5, sd = 10), # theoretische stabw von Messgeräten
#              beta_end  = rnorm(length(id), mean = -22.94   , sd = 6),   beta_start  = rnorm(length(id), mean = -2.5, sd = 10),
#              conc_start  = rnorm(length(id), mean = 333,    sd = 2),
#              conc_end    = rnorm(length(id), mean = Dppb[1]+id, sd = 8),
#              # conc_end  = rnorm(length(id), mean = 333.5+id , sd = 0.9 ), conc_start  = rnorm(length(id), mean =  333.5, sd = 0.3),
#              d18O_end  = rnorm(length(id), mean = 31.79    , sd = 10  ),d18O_start  = rnorm(length(id), mean = 44,  sd = 12))]


   
  sim[,':=' (SP_start =  alpha_start - beta_start, SP_end   =  alpha_end - beta_end,
             bulk_start = (alpha_start + beta_start)/2, bulk_end = (alpha_end + beta_end)/2) ]
  
  # sim.mean <- sim[,.(SP_source    = mean((SP_end   * conc_end - SP_start   * conc_start)/(conc_end - conc_start)),
  #                    bulk_source  = mean((bulk_end * conc_end - bulk_start * conc_start)/(conc_end - conc_start)),
  #                    alpha_source = mean((alpha_end* conc_end - alpha_start* conc_start)/(conc_end - conc_start)),
  #                    beta_source  = mean((beta_end * conc_end - beta_start * conc_start)/(conc_end - conc_start)),
  #                    d18O_source  = mean((d18O_end * conc_end - d18O_start * conc_start)/(conc_end - conc_start)) 
  #                    ),id]

  sim.sd   <- sim[,.(SP_source    = sd((SP_end   * conc_end - SP_start   * conc_start)/(conc_end - conc_start)),
                     bulk_source  = sd((bulk_end * conc_end - bulk_start * conc_start)/(conc_end - conc_start)),
                     alpha_source = sd((alpha_end* conc_end - alpha_start* conc_start)/(conc_end - conc_start)),
                     beta_source  = sd((beta_end * conc_end - beta_start * conc_start)/(conc_end - conc_start)),
                     d18O_source  = sd((d18O_end * conc_end - d18O_start * conc_start)/(conc_end - conc_start))
                     ),id]

  
  sim.melt  <- melt(sim.sd, id.vars = c("id"), 
                         measure.vars = c("SP_source","bulk_source","alpha_source","beta_source","d18O_source"), 
                    value.name = "precision", variable.name = "variable")
  
  # sim.melt$target <- sim.sd[alpha_source < input$target, min(id)]
  
  
  target.bulk  <- sim.sd[bulk_source  < input$target, min(id)]
  target.alpha <- sim.sd[alpha_source < input$target, min(id)]
  target.beta  <- sim.sd[beta_source  < input$target, min(id)]
  target.SP    <- sim.sd[SP_source    < input$target, min(id)]  
  target.18O   <- sim.sd[d18O_source  < input$target, min(id)]
  }) # end isolate
  
  list(sim.melt = sim.melt, target.bulk = target.bulk, target.alpha = target.alpha, target.SP = target.SP, target.18O = target.18O, target.beta = target.beta)
    # }) # end observeEvent
  })
  
  output$value <- renderText({
    # input$sim.go
    paste("Concentration increase in ppb when targed precision is reached: \n - alpha:",simdata()$target.alpha,
          "\n - beta:",simdata()$target.beta, "\n - bulk:",simdata()$target.bulk, 
          "\n - SP:\t",simdata()$target.SP, "\n - 18O:\t",simdata()$target.18O)   # ,simdata()$target[1])
    
  })
  
  output$precisionSD <- renderPlotly({
    # input$sim.go
      # generate bins based on input$bins from ui.R
     ggplotly(
         ggplot(data = simdata()$sim.melt, aes_string(x ="id", y = "precision", colour = "variable") ) + 
           geom_point(size = 0.95, show.legend = T) +
           # ggplot(data = ghg.file.plot, aes(x = date.strp, y = mean, colour = group.var)) + 
           # geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2) +
           geom_line() +  ylim(-10,50) +   
           geom_hline(yintercept = input$target, colour = "darkgrey", linetype = 2) + 
           # guides(fill = guide_legend(label.hjust = 3, title = "chose species to show on graph:")) +
           # scale_linetype_manual("targed threshold", values = 2) +
           # guides(fill = guide_legend(override.aes = list(linetype = 2))) +
         
           # geom_smooth(method = "loess") +
           theme_classic() + #+  theme_dark() +
           # theme(legend.position = "bottom") +  # this does not work in ggplotly
           ylab(paste("SD δ<sup>15</sup>N source [\u2030]")) + xlab("ΔN<sub>2</sub>O [ppb]")  #expression(  δ
         , height = 500, dynamicTicks = T)  # )
   })  # end renderPlotly
  
  output$plot.table <- DT::renderDataTable(datatable(dcast(simdata()$sim.melt, id ~ variable), options = list(pageLength = 20,lengthChange=T)) 
                                           %>% formatSignif(c(2:9),3))
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    # content = function(file) {
    filename = function() {
      paste("MonteMattiSim.csv", sep="")
    },
    content = function(file) {
      fwrite(dcast(simdata()$sim.melt, id ~ variable), file, row.names = FALSE, sep = ",")
    }
    
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

