#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



source("HiPerMAb.R")

####################################### Libraries #################################################
if("Hmisc" %in% rownames(installed.packages())){
  library(Hmisc)} else{
    install.packages("Hmisc")
    library(Hmisc)}

if("pROC" %in% rownames(installed.packages())){
  library(pROC)} else{
    install.packages("pROC")
    library(pROC)}

if("discretization" %in% rownames(installed.packages())){
  library(discretization)} else{
    install.packages("discretization")
    library(discretization)}

if("plotly" %in% rownames(installed.packages())){
  library(plotly)} else{
    install.packages("plotly")
    library(plotly)}

 if("Biocomb" %in% rownames(installed.packages())){
   library(Biocomb)} else{
     install.packages("Biocomb")
     library(Biocomb)}

if("data.table" %in% rownames(installed.packages())){
  library(data.table)} else{
    install.packages("data.table")
    library(data.table)}


if("DT" %in% rownames(installed.packages())){
  library(DT)} else{
    install.packages("DT")
    library(DT)}

####################################### User Interface ############################################

ui <- fluidPage(
  fileInput(inputId = "data", label = "Choose a file.", accept = c("text/csv", "text/comma-separated-values, text/plain", ".csv")), 
  
  # Application title
  titlePanel("High Performance Measurement Abundance (HiPerMAb)"),
  tabPanel("HiPerMAb",
       tabsetPanel(
       tabPanel("Performance",
        selectInput(inputId = "pfM.method", label = "performance method",choices =list("entropy","mAUC","AAC","HUM","misClassRate")),
        #selectInput(inputId = "random.simulation", label = "simulate random data",choices = c("Monte Carlo","Permutation")),
        selectInput(inputId = "imput.method", label = "impute missing values",choices = c("median","random")),
                
        #selectInput(inputId = "is.positive", label = "is.positive",choices = c( "FALSE", "TRUE")),
        selectInput(inputId = "corrected.method", label = "corrected p-value",choices = c( "FWER", "FDR"))
      ,
        splitLayout(
          textInput("text", "positive class", value = "NULL")
          #,
          #numericInput(inputId = "no.simulations", label = "no.simulations",value = 1000, min = 100, max = 1000000, step = 100),
          #numericInput(inputId = "Con.Interval", label = "Conf.Interval",value = 0.95, min = 0, max = 1, step = 0.05)
        ),
        br(),
        
                tabPanel("Performance table",
                         
                         actionButton( "performanceTable", "Performance table"),
                         dataTableOutput("performance.table")
                         
                )),

   tabPanel("HiPerMAb table and curves",
            selectInput(inputId = "pfM.method", label = "performance method",choices =list("entropy","mAUC","AAC","HUM","misClassRate")),
            selectInput(inputId = "random.simulation", label = "simulate random data",choices = c("Monte Carlo","Permutation")),
            selectInput(inputId = "imput.method", label = "impute missing values",choices = c("median","random")),
           
            selectInput(inputId = "is.positive", label = "is positive",choices = c( "FALSE", "TRUE")),
            
            splitLayout(
              textInput("text", "positive class", value = "NULL"),
              numericInput(inputId = "no.simulations", label = "no. simulations",value = 1000, min = 100, max = 1000000, step = 100),
              numericInput(inputId = "Con.Interval", label = "Confidence Interval",value = 0.95, min = 0, max = 1, step = 0.05)
            ),
            br(),  
            
  tabPanel("HiPerMAb table",
           
           actionButton( "HiPerMAbTable", "HiPerMAb table"),
           dataTableOutput("HiPerMAb.table")
  ),
  br(), 
  tabPanel("HiPerMAb curves",
           
           actionButton( "HiPerMAbCurve", "HiPerMAb curves"),
           plotlyOutput("HiPerMAb.Curve")
  )
),
   tabPanel("HiPerMAb power",
            splitLayout("no. Biomarker Candidates",
                numericInput(inputId = "bc.no1", label = "from ",value = 100, min = 0, max = 1000000, step = 50),
                numericInput(inputId = "bc.no2", label = "to ",value = 1000, min = 50, max = 1000000, step = 50)
                
            ), 
            splitLayout("p-value of true Biomarker Candidates",
              numericInput(inputId = "pv.no1", label = "from ",value = 0.001, min = 0.001, max = 0.05, step = 0.001),
              numericInput(inputId = "pv.no2", label = "to ",value = 0.05, min = 0.001, max = 0.05, step = 0.001)
            ) ,
            splitLayout("Confidence Interval (1-alpha)",
                        numericInput(inputId = "alpha", label = "alpha",value = 0.05, min = 0.001, max = 0.05, step = 0.002),
                        numericInput(inputId = "power", label = "power",value = 0.8, min = 0, max = 1, step = 0.05)
            ),
            # splitLayout("3D plot",
            #             numericInput(inputId = "theta", label = "theta",value = 30),
            #             numericInput(inputId = "phi", label = "phi",value = 50)
            # ),
            
            #selectInput(inputId = "grid.plot", label = "grid.plot",choices = c( "TRUE","FALSE")),
            tabPanel("show plot",
                     
                     actionButton( "Dplot", "show  plot"),
                     plotlyOutput("D.plot")
            )
            
   )
 )))

####################################### Server ####################################################
server <- function(input, output) {
  
  
  
  perfo.table<- eventReactive(input$performanceTable, {
    data.input <- input$data
    if(!is.null(data.input)){
      isolate(performance(read.csv2(data.input$datapath), input$random.simulation, input$imput.method, 
                          input$pfM.method,input$no.simulations,input$pos.class,input$Con.Interval, 
                          input$is.positive,input$corrected.method))
      
    }
    
  })
  output$performance.table<- DT::renderDataTable({
    perfo.table()
  })
  
  curve<- eventReactive(input$HiPerMAbCurve, {
    data.input <- input$data
    if(!is.null(data.input)){
      isolate( hipermab(read.csv2(data.input$datapath),input$random.simulation, input$imput.method, input$pfM.method,input$no.simulations,
                        input$pos.class,input$Con.Interval, input$is.positive)$ByPlot)
    }
  })
  
  output$HiPerMAb.Curve <- renderPlotly({
    
    curve()
    
  })
  
  hipermab.table<- eventReactive(input$HiPerMAbTable, {
    data.input <- input$data
    if(!is.null(data.input)){
      isolate(hipermab(read.csv2(data.input$datapath),input$random.simulation, input$imput.method, input$pfM.method,input$no.simulations,
                       input$pos.class,input$Con.Interval, input$is.positive)$ByNumbers)
      
    }  
  })
  output$HiPerMAb.table<- DT::renderDataTable({
    hipermab.table()
  })
  
  Dplot<- eventReactive(input$Dplot, {
    
    isolate(plot.no.required.true.biomarkers(seq(from=input$bc.no1,to=input$bc.no2,by=50),
   seq(from=input$pv.no1,to=input$pv.no2,by=0.002),alpha= input$alpha,theta=input$theta,phi=input$phi,grid.plot=input$grid.plot))
    
            })
  
  output$D.plot <- renderPlotly({
    
    Dplot()
    
  })
  
  
}


####################################### Run the application ####################################### 
shinyApp(ui = ui, server = server)
