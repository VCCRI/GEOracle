
appCSS <- "
#loading-content {
position: absolute;
background: #000000;
opacity: 0.9;
z-index: 100;
left: 0;
right: 0;
height: 100%;
text-align: center;
color: #FFFFFF;
}
"

########################################################################
# User Interface
########################################################################

shinyUI(fluidPage(
  
  useShinyjs(),
  inlineCSS(appCSS),
  
  # Loading message
  div(
    id = "loading-content",
    h2("GEOracle is loading...")
  ),
  
  tags$head(
    tags$style(HTML("
.shiny-progress .progress {
  position: absolute;
                    width: 100%;
                    top: 0px;
                    height: 20px;
                    margin: 0px;
                    }
                    
  .shiny-progress .bar {
  opacity: 0.6;
                    transition-duration: 250ms;
height: 20px;
                    }
                    
.shiny-progress .progress-text {
    position: absolute;
    right: 10px;
    height: 50px;
    width: 600px;
    background-color: #eef8ff;
      margin: 0px;
    padding: 2px 3px;
    opacity: 0.85;
  }

.shiny-progress .progress-text .progress-message {
  padding: 0px 3px;
  font-weight: bold;
  font-size: 200%;
}

.shiny-progress .progress-text .progress-detail {
  padding: 0px 3px;
  font-size: 200%;
}


  "))
    ,
    tags$head(tags$style("#curtab{color: white;
                                 font-size: 5px;
                                 font-style: italic;
                                 }"
    )
    )
    
    ,
    tags$head(tags$style("#curtabmessage{color: white;
                                 font-size: 5px;
                                 font-style: italic;
                                 }"
    )
    )
    
    
  ),
  
  
  
  sidebarPanel(
    
    
    conditionalPanel("(output.curtab == 'Verify' || output.curtab == 'Add')",
                     tags$div(title="GEOracle", img(src='http://georacle.victorchang.edu.au/georacle.png', align = "left", width="90%"))
                     ,
                     fileInput('uploadFile', "Upload list of GSE IDs")
                     ,
                     
                     a(href="http://georacle.victorchang.edu.au/TGFB_Case_Study_GSE_IDs.txt", "Sample input file (Right click / save as)", target = "blank")
                     ,
                     
                     bsTooltip(id = "uploadFile", title = "<p>e.g. GSE123 GSE345 GSE567 ...<\\p>", 
                               placement = "right", trigger = "hover")
                     
                     ,
                     textOutput('numGSEs')
                     ,
                     conditionalPanel("typeof output.numGSEs === 'string'",
                                      h4('Set filters', align = "center")
                                      ,
                                      selectInput('filter', 'Strictness', c("Default", "Custom"), selected = "Default")
                                      ,
                                      ### custom filter panel
                                      conditionalPanel("input.filter == 'Custom'",
                                                       
                                                       htmlOutput("speciesUI")
                                                       ,
                                                       div(style="display:inline-block",
                                                           #sliderInput("clusterrange", "Cluster size ange:", min = 1, max = 15, value = c(2,5), width=150))
                                                           numericInput("minClusterSize", "Min. Cluster Size", value = 2, min = 1, step = 1, width = 100))
                                                       ,
                                                       div(style="display:inline-block",
                                                           numericInput("maxClusterSize", "Max. Cluster Size", value = 5, min = 2, step = 1, width= 100))
                                                       ,
                                                       div(style="display:inline-block",
                                                           checkboxInput("allbetween", "All clusters correct size", value = 1))
                                                       ,
                                                       div(style="display:inline-block",
                                                           checkboxInput("simpleOnly", "Simple (2 clusters) only", value = 0))
                                                       ,
                                                       div(style="display:inline-block",
                                                           numericInput("maxcomps", "Max. Comparisons", value = 10, min = 1, step = 1, width= 100))
                                                       ,
                                                       div(style="display:inline-block",
                                                           checkboxInput("singlePlatformOnly", "1 platform only", value = 1))
                                                       ,
                                                       div(style="display:inline-block",
                                                           checkboxInput("onechannel", "1 channel array", value = 1))
                                                       ,
                                                       div(style="display:inline-block",
                                                           checkboxInput("predpert", "Predicted perturbation", value = 1))
                                      )
                                      ,
                                      actionButton("go", " COMPUTE", icon("wrench"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                     )
                     ,
                     conditionalPanel("input.go > 0",
                                      
                                      hr()
                                      ,
                                      h4('Processed GSEs', align = "center")
                                      ,
                                      h5("Click the GSE IDs below to verify results before proceeding\n")
                                      ,
                                      hr()
                                      ,
                                      DT::dataTableOutput('processedGSEs')
                                      ,
                                      tags$script(HTML(
                                        "Shiny.addCustomMessageHandler('pager',function(page) {
      $('#'+'processedGSEs').find('table').DataTable().page(page).draw(false);
    })"
                                      )),
                                      hr()
                                      ,
                                      conditionalPanel("input.processedGSEs_row_last_clicked > 0",
                                                       h5("Once you have renamed and verified all comparisons click here to create a GRN\n")
                                                       ,
                                                       
                                                       actionButton("finished", "NEXT STEP", icon("play"), style="color: #fff; background-color: #29C28A; border-color: #2e6da4")
                                      )
                     )
    )
    ,
    conditionalPanel("output.curtab == 'GRN'",
                     tags$div(title="GEOracle", img(src='http://georacle.victorchang.edu.au/georacle.png', align = "left", width="90%"))
                     ,
                     hr()
                     ,
                     actionButton("goback", "GO BACK", icon("step-backward"), style="color: #fff; background-color: rgb(156, 5, 5); border-color: #2e6da4")
                     ,
                     hr()
                     ,
                     textInput("folder", "Output Folder Name", value = gsub("[: -]","_",Sys.time()))
                     ,
                     actionButton("process", "Calculate D.E.", icon("object-align-bottom"), style="color: #fff; background-color: rgb(18, 88, 24); border-color: rgb(141, 252, 0)")
    )
    ,
    
    uiOutput("curtabs")
    
    
    
  ),
  
  mainPanel(
    conditionalPanel("output.curtab == 'Verify'",
                     uiOutput('matchedPairs')
    )
    ,
    conditionalPanel("input.process > 0",
                     uiOutput('makeGRN')
    )
    ,
    conditionalPanel("output.curtab == 'Add'",
                     uiOutput('addComp')
    )
  )
))
