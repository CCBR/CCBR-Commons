library(shiny)

ui <- fluidPage(
  titlePanel("CCBR Microarray analysis workflow", windowTitle="CCBR Microarray analysis workflow"),
  h5("(For human or mouse data)"),
sidebarLayout(
  sidebarPanel(width=12,
  fluidRow(align="Top",
    
    column(3,
           textInput("ProjectID", label=h6("Project ID:"), value="CCBR-", width="150px"),
           checkboxGroupInput("Species", label=h6("Select species"), choices=c("Human", "Mouse"),selected=TRUE)
          ),
    column(3,
           numericInput("NoSamples", label=h6("Total number of samples"), value="2",width="150px"),
           numericInput("NoContrasts", label=h6("Number of contrasts"),value="1", width="150px")
           ),
    column(3,
           numericInput("NoGroups", label=h6("Number of groups"), value="2", width="150px"),
           h6("Group names"),
           tags$style(type="text/css", "textarea {width:100%}"),
           tags$textarea(id = "GroupNames", label=h6("Group names"),placeholder = "Type here", rows = 4, "")
            ),
    column(3,
           fileInput("CELfiles", label=h6("Choose CEL files"), multiple= T, accept = c('.cel','.CEL')),
           selectInput("Platform", label = h6("Select platform"), choices = list("Choice 1" = 1, "Choice 2" = 2, "Choice 3" = 3), selected = 1)
    )
  ),
  submitButton(text="Click to assign samples to groups and create contrasts")
),
mainPanel(
  tabsetPanel(type = "tabs", 
              tabPanel("ArrayQualityMetricsReport", plotOutput("plot")), 
              tabPanel("Plots", verbatimTextOutput("summary")), 
              tabPanel("Tables", tableOutput("table")),
              tabPanel("Inputs", 
                       textOutput("projectid"), 
                       textOutput("Species"), 
                       textOutput("NoSamples"), 
                       textOutput("NoGroups"), 
                       textOutput("NoContrasts"))
  )
)
),
div(style="display:inline-block",submitButton("Generate PDF report"), style="float:center")
)

server <- function(input, output) {
  output$projectid=renderText({paste("Project ID: ",input$ProjectID)})
  output$Species=renderText({paste("Species: ",input$Species)})
  output$CELfiles=renderText(input$CELfiles)
  output$NoSamples=renderText({paste("No. of samples: ",input$NoSamples)})
  output$NoContrasts=renderText({paste("No. of contrasts: ",input$NoContrasts)})
  output$NoGroups=renderText({paste("No. of groups: ",input$NoGroups)})
}



shinyApp(ui = ui, server = server)
