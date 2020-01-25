#
# SHINY APP CREATED BY *****ANDRES CAMILO MENDEZ****
# THIS APP COMES WITH ABSOLUTLEY NO WARRANTY
# FEEL FREE TO DISTRIBUTE
#
# User interfaz for betaBS classifier
#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#


suppressMessages(if(!require(yaml)){install.packages('yaml'); library(yaml)} else {library(yaml)})
suppressMessages(if(!require(shiny)){install.packages("shiny");library(shiny)}else{library(shiny)})
suppressMessages(if(!require(shinydashboard)){install.packages("shinydashboard");library(shinydashboard)}else{library(shinydashboard)})
suppressMessages(if(!require(shinyBS)){install.packages("shinyBS");library(shinyBS)}else{library(shinyBS)})
suppressMessages(if(!require(shinyjs)){install.packages('shinyjs'); library(shinyjs)} else {library(shinyjs)})

#suppressMessages(if(!require(ggplot2)){install.packages('ggplot2'); library(ggplot2)} else {library(ggplot2)})
#suppressMessages(if(!require(plotly)){install.packages('plotly'); library(plotly)} else {library(plotly)})
source("www/helpers.R", local = TRUE)
source("r_functions.R", local = TRUE)


# Define UI for application that draws a histogram
header <- dashboardHeader(
  title = HTML("<h4><b> Dengue App </b></h4>"),titleWidth = 230
  
)

sidebar<-dashboardSidebar( width = 230,
                           
                           sidebarMenu( id = "menu",
                                        #menuItem("Beans Geenepools",tabName = "geneBeans", icon = icon("dashboard"), startExpanded = TRUE, ),
                                        menuItem("Calculate probs",tabName="item1", icon = icon( "fa fa-angle-double-right"), selected = TRUE),
                                        menuItem("Add new records", tabName = "item2" ,icon = icon("far fa-angle-double-right")  )
                                        
                                        #End menu item beans genepool
                                        
                           )
                           
                           
)

body <- dashboardBody(
  
  tabItems(
    
    tabItem(tabName = "item1",
            
            fluidRow(
              HTML("<h5><b>1- Upload database </b></h5>")
            ,useShinyjs()
             ,fileInput("data_in", 
                         "",
                        multiple = FALSE, 
                        buttonLabel = "Buscar",
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
              )
            ,tags$hr()
            ,HTML('<h5> <b>2- Preview database </b> </h5>')
            ,dataTableOutput("data_in_prev")
            ,HTML('<h5> <b>3- Database Logfile  </b> </h5>')
            ,htmlOutput("logfile", container = tags$ul, style = "width:80%; height:10px" )
            ,tags$hr()
            ,HTML('<h5> <b>4- Calculate probs taking into account Hemogram vars?  </b></h5>')
            ,radioButtons("quest", label = NULL, choices = c("No" = 1, "yes" = 2), inline = TRUE)
            ,conditionalPanel("input.quest == 2"
            ,numericInput("Leucos_pc", label = HTML("<h5><b>5- Enter the cutpoint for leukocytes </b></h5>"), value = 4200, min=0, max = 18200)
            ,numericInput("plaq_pc", label = HTML("<h5><b>6- Enter the cutpoint for platelets </b></h5>"), value = 165000, min=0, max = 474000)
            ,selectInput("hmog",label = HTML("<h5><b>7- select the way of build the hemogram variable (H) </b></h5>"), choices = c(
              "H+ = (+,+) " = 1,
              "H+ = (+,+) ; (+,-)" = 2,
              "H+ = (+,+) ; (-,+)" = 3,
              "H+ = (+,+) ; (+,-) ; (-,+)" = 4
            ), width = "350px")#end selectInput
            )#end conditional panel
            ,HTML('<h5><b> 8- Calculate probabilities </b> </h5>')
            ,radioButtons("prType", label = NULL, choices = c("Discrete" = 1, "Continuos" = 2), inline = TRUE)
            ,conditionalPanel("input.prType == 1"
                              ,sliderInput("prev", label = "Prevalence of dengue disease:", min = 0, max = 1, value = 0.1, step = 0.1, round = TRUE)
                              )#end conditional panel
            ,withBusyIndicatorUI(
             bsButton("cal_probs", label = "Calculate", style = "primary")
            )
            ,HTML('<h5><b>9- Preview  the probabilities </b> </h5>')
            ,dataTableOutput("arbol_prev")
            ,HTML("<h5><b>10- Download probabilites  </b></h5>")
            ,withBusyIndicatorUI(
              downloadButton("save_probs", label = " Download results")
            )
            ,HTML('<h5><b> 11- Summary table </b></h5>')
            ,tableOutput("summ")
            ,HTML("<h5><b>12- Download summary table  </b></h5>")
            ,withBusyIndicatorUI(
              downloadButton("save_summ", label = " Download summary")
            )
            ,HTML("<h5><b>13- ROC curve plot  </b></h5>")
            ,plotOutput("roc_curve", width = "800px", height = "400px")
            ,HTML("<h5><b>(OPTIONAL) Download ROC curve data frame  </b></h5>")
            ,withBusyIndicatorUI(
              downloadButton("save_roc", label = " Download ROC curve data")
            )

            
            )#end fluidRow
            
            
            
            
            )#end tabItem1
    
    ,tabItem(tabName = "item2"
             
             ,fluidRow(
               HTML("<h3><b> Add new individual (just for continuos probabilities) </b></h3>")
               ,HTML("<h5><b>1- Upload probabilities database (arbol) to be updated </b></h5>")
               ,useShinyjs()
               ,fileInput("arbol_in", 
                         "",
                         multiple = FALSE, 
                         buttonLabel = "Buscar",
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")
               )
               ,HTML('<h5> <b>2- Preview probabilities database </b> </h5>')
               ,dataTableOutput("arbol_in_prev")
               ,HTML('<h5> <b>3- Upload file with the new records </b> </h5>')
               ,fileInput("new_records", 
                          "",
                          multiple = FALSE, 
                          buttonLabel = "Buscar",
                          accept = c("text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
               )
               ,dataTableOutput("newr_in_prev")
               ,tags$hr()
               ,HTML('<h5> <b>New individuals database logfile </b> </h5>')
               ,htmlOutput("logfile2", container = tags$ul, style = "width:80%; height:10px" )
               ,tags$hr()
               ,tags$hr()
               ,HTML('<h5> <b>4- Should the App take into account hemogram variables?  </b></h5>')
               ,radioButtons("quest2", label = NULL, choices = c("No" = 1, "yes" = 2), inline = TRUE)
               ,conditionalPanel("input.quest2 == 2"
                                 ,numericInput("Leucos_pc2", label = HTML("<h5><b>5- Enter the cutpoint for leukocytes </b></h5>"), value = 4200, min=0, max = 18200)
                                 ,numericInput("plaq_pc2", label = HTML("<h5><b>6- Enter the cutpoint for platelets </b></h5>"), value = 165000, min=0, max = 474000)
                                 ,selectInput("hmog2",label = HTML("<h5><b>7- select the way of build the hemogram variable (H) </b></h5>"), choices = c(
                                   "H+ = (+,+) " = 1,
                                   "H+ = (+,+) ; (+,-)" = 2,
                                   "H+ = (+,+) ; (-,+)" = 3,
                                   "H+ = (+,+) ; (+,-) ; (-,+)" = 4
                                 ), width = "350px")#end selectInput
               )
               ,HTML('<h5> <b>8- Assign probabilities to new records </b> </h5>')
               ,withBusyIndicatorUI(
                 bsButton("assign_probs", label = "Assign probabilities ", style = "primary")
               )
               ,HTML('<h5> <b>9- Preview Assigned probabilities to new records </b> </h5>')
               ,dataTableOutput("assigned_in_prev")
               ,HTML('<h5> <b>10- Download new records with assigned probabilities  </b> </h5>')
               ,withBusyIndicatorUI(
                 downloadButton("save_assignation", label = " Download probs assigned")
               )
               ,HTML('<h5> <b>11- Do you want to calculate the performance measure?  </b></h5>')
               ,radioButtons("quest3", label = NULL, choices = c("No" = 1, "yes" = 2), inline = TRUE)
               ,conditionalPanel("input.quest3 == 2"
                                 ,numericInput("bayes_pc", label = HTML("<h5><b>5- Enter the cutpoint for bayes probability </b></h5>"), value = 0.5, min=0, max = 1, step = 0.1)
                                 ,conditionalPanel("input.quest2 == 2"
                                   ,numericInput("bayes.hemo.pos_pc", label = HTML("<h5><b>6- Enter the cutpoint for bayes.hemo.pos  </b></h5>"), value = 0.5, min=0, max = 1, step = 0.1)
                                   ,numericInput("bayes.hemo.neg_pc", label = HTML("<h5><b>6- Enter the cutpoint for bayes.hemo.neg  </b></h5>"), value = 0.5, min=0, max = 1, step = 0.1)
                                 )#end nested conditional panel
                                 ,withBusyIndicatorUI(
                                   bsButton("calc_summ2", label = " Calculate measures", style = "primary")
                                 )
                                 ,HTML('<h5> <b> Performance Measures  </b></h5>')
                                 ,tableOutput("summ2")
               )#end conditional panel
               ,HTML('<h5> <b> 12- Do you want to update probabilities? </b></h5>')
               ,radioButtons("quest4", label = NULL, choices = c("No" = 1, "yes" = 2), inline = TRUE)
               ,conditionalPanel("input.quest4 == 2"
                                 ,withBusyIndicatorUI(
                                   bsButton("update_prs", label = " Update probabilities", style = "primary")
                                 )
                                 ,HTML('<h5> <b> 13- Preview of updated probabilities  </b></h5>')
                                 ,dataTableOutput("updated_probs_prev")
                                 ,HTML('<h5> <b>  Download updated probabilities  </b></h5>')
                                 ,withBusyIndicatorUI(
                                   downloadButton("save_updated_probs", label = " Download updated probs")
                                 )
             
                                 ,HTML('<h5><b> 14- Summary table </b></h5>')
                                 ,tableOutput("summ3")
                                 ,HTML("<h5><b> Download summary table  </b></h5>")
                                 ,withBusyIndicatorUI(
                                   downloadButton("save_summ3", label = " Download summary")
                                 )
                                 ,HTML("<h5><b>15- ROC curve plot  </b></h5>")
                                 ,withBusyIndicatorUI(
                                   bsButton("calc_croc2", label = " Plot ROC curve", style = "primary")
                                 )
                                 ,plotOutput("roc_curve2", width = "800px", height = "400px")
                                 ,HTML("<h5><b>(OPTIONAL) Download ROC curve data  </b></h5>")
                                 ,withBusyIndicatorUI(
                                   downloadButton("save_roc2", label = " Download ROC curve data")
                                 )
                                 )#end conditonal panel
             
               
             )#end fluidRow
             
     
             )#end tabItem2
  
    
  )#end tabitemsss
 
)#end DAshboardBody

dashboardPage(
  header,
  sidebar,#dashboardSidebar(disable = F),
  body
)



