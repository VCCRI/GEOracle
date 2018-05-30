
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


jscode <- "shinyjs.closeWindow = function() { window.close(); }"

########################################################################
# User Interface
########################################################################

shinyUI(fluidPage(
  
  useShinyjs(),
  inlineCSS(appCSS),
  extendShinyjs(text = jscode, functions = c("closeWindow")),
  
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
                     #tags$div(title="GEOracle", img(src='http://georacle.victorchang.edu.au/georacle.png', align = "left", width="90%"))
                     imageOutput("GEoracleLogo", height = "auto")
                     ,
                     hr()
                     ,
                     uiOutput('analysis_select')
                     ,
                     
                     tabsetPanel(id = "inTabset",
                                 tabPanel("Search GEO",
                                          #renderText("You can search the GEO database for datasets to process here."),
                                          textInput('searchTerm', label = "Enter a search term:"),
                                          
                                          #bsTooltip(id = "searchTerm", title = "You can use & or | for multiple keywords.",
                                          #placement = "bottom", trigger = "hover"),
                                          
                                          div(id="species_select_search_div",
                                              selectInput('selectSpecies', 'Species:', choices = c("A","C"))
                                          )
                                          ,
                                          
                                          div(id="microarray_checkbox_div",
                                              checkboxInput("GSE_platform_filter", label ="Show only microarray GSE", value=TRUE)#,
                                          )
                                          ,
                                          
                                          actionButton('find', "Search GEO", 
                                                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                          hr(),
                                          
                                          uiOutput("searchStatus"),
                                          dataTableOutput(outputId="GEOtable"),
                                          hr(),
                                          
                                          div(id="keep_button_div", style="display:inline-block",
                                              uiOutput('keep.button'))
                                          ,
                                          div(id="remove_button_div", style="display:inline-block",
                                              uiOutput('remove.button'))
                                          ,
                                          hr()
                                          ,
                                          div(id="download_GSE_button_div",
                                              uiOutput('download.button'),
                                              uiOutput('batch.download.button')
                                          )
                                          ,
                                          div(id="analyse_GSE_button_div",
                                              uiOutput('direct.analyse.button')
                                          )
                                          
                                          
                                          
                                          
                                 )
                                 ,
                                 tabPanel("Process GSE list",
                                          hr()
                                          ,
                                          
                                          fluidRow(
                                            column(6,
                                                   
                                                   textAreaInput('GSEtext', "Direct GSE input", placeholder = "GSE1234 GSE5678 GSE9123"), actionButton("submitGSETtext", "Submit")
                                            )
                                            ,
                                            column(6,
                                                   div(id="upload_file_div",
                                                       fileInput('uploadFile', "Upload list of GSE IDs")
                                                   )
                                                   
                                                   ,
                                                   bsTooltip(id = "uploadFile", title = "e.g. GSE123 GSE345 GSE567 ...",
                                                             placement = "bottom", trigger = "hover")
                                                   ,
                                                   
                                                   a(href="http://georacle.victorchang.edu.au/TGFB_Case_Study_GSE_IDs.txt", "Sample input file (Right click / save as)", target = "blank")
                                            )
                                          )
                                          ,
                                          hr()
                                          ,
                                          textOutput('numGSEs')
                                          ,
                                          conditionalPanel("typeof output.numGSEs === 'string'",
                                                           h4('Set filters', align = "center")
                                                           ,
                                                           uiOutput('filter_select')
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
                                                                                checkboxInput("allbetween", "All clusters correct size", value = 0))
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
                                                                            ,
                                                                            div(style="display:inline-block",
                                                                                checkboxInput("gene_symbols_only", "Pre-annotated genes only", value = 0))
                                                                            ,
                                                                            checkboxInput("see_all_GSE", "Let every GSE through", value = F)
                                                           )
                                                           ,
                                                           actionButton("go", " COMPUTE", icon("wrench"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                          )
                                 )
                                 
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
                     imageOutput("GEoracleLogo2", height = "auto")
                     ,
                     hr()
                     ,
                     actionButton("goback", "GO BACK", icon("step-backward"), style="color: #fff; background-color: rgb(156, 5, 5); border-color: #2e6da4")
                     ,
                     hr()
                     ,
                     textInput("folder", "Output Folder Name", value = gsub("[: -]","_",Sys.time()))
                     ,
                     hr()
                     ,
                     h2("Set differential expression significance thresholds:")
                     ,
                     selectInput("logfc", "Absolute log2 FC", choices = c(0.5,1,2), selected = 1)
                     ,
                     selectInput("pval", "B.H. adjusted P-value", choices = c(0.1,0.05,0.01,0.001), selected = 0.05)
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
