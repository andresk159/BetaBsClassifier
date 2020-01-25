#
# SHINY APP CREATED BY *****ANDRES CAMILO MENDEZ****
# THIS APP COMES WITH ABSOLUTLEY NO WARRANTY
# FEEL FREE TO DISTRIBUTE
#
# APP TO CLASSIFY IN A FRIENDLY WAY SICK DENGUE INDIVIDUALS  
#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

suppressMessages(if(!require(yaml)){install.packages('yaml'); library(yaml)} else {library(yaml)})
suppressMessages(if(!require(shiny)){install.packages("shiny");library(shiny)}else{library(shiny)})
suppressMessages(if(!require(shinydashboard)){install.packages("shinydashboard");library(shinydashboard)}else{library(shinydashboard)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(gtools)){install.packages('gtools'); library(gtools)} else {library(gtools)})
suppressMessages(if(!require(caTools)){install.packages('caTools'); library(caTools)} else {library(caTools)})
suppressMessages(if(!require(shinyjs)){install.packages('shinyjs'); library(shinyjs)} else {library(shinyjs)})



###### SPACE TO DEFINE IMPORTANT FUNCTIONS ############
#calculate hemogram variable using the tresshods forleucos and plaq 


source("www/helpers.R", local = TRUE)
source("r_functions.R", local = TRUE)

############### END FUNCTION SPACE ######################


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
   
  #upload data bases and check some possible errors
  bd <- reactive ({
  
   
  req(input$data_in)
    
    tryCatch(
      {
        df.raw <- read.csv(input$data_in$datapath, header = TRUE, sep = ",") 
        
        df <- read.csv(input$data_in$datapath, header = TRUE, sep = ",") %>% 
          filter(complete.cases( dplyr::select(., 1:dxdengue) ) )
          
        vals <- df %>% 
          dplyr::select(., 1:dxdengue) %>% 
          apply(., 2, function(x){ length(unique(x))  })
         
        vals <- if_else(all( vals == 2 ) , "There is Not problems", "Error: There was found different values from  0 and 1." )
        
        nsimt <- df %>% 
          dplyr::select(., 1:dxdengue) %>% 
          dplyr::select(., 1:(ncol(.)-1)) %>%
          ncol(.)
        
        nsimt <- ifelse(nsimt <= 12 & nsimt >= 5,"There is Not problems with the number of Symptoms", "Error: The number of symptoms must be between 5 and 12")
     
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
  
   
  return(list(bd.filtered = df, bd.raw = df.raw  ,logs= data_frame(valores = vals, sintomas = nsimt) ))  
  } )
  
  #render database
  output$data_in_prev <- renderDataTable({ 
    bd()[[1]]
  } ,options = list(
    pageLength = 5,
    scrollX = TRUE
  ))
  
  #write and show the status of database
  output$logfile <- renderText({
    
    if( length(grep("Error:",  bd()$logs$valores)) != 0 ) {
      
     x1 <- paste('<li> Values:', '<font color = "red">', bd()$logs$valores , '</font></li> ')
      
    }else{ x1 <-  paste('<li> Values:', '<font color = "green">', bd()$logs$valores , '</font></li> ') }

    if( length(grep("Error:",  bd()$logs$sintomas)) != 0 ) {
      
     x2 <-  paste('<li> Number of Symptoms:', '<font color = "red">', bd()$logs$sintomas , '</font></li>  ')
      
    }else{ x2 <-paste('<li> Number of Symptoms:', '<font color = "green">', bd()$logs$sintomas , '</font></li>') }
    
    return(paste(x1,"\n", x2))

  })
 
  #Calculate Bayesian probabilities
  
 arbol <- eventReactive(input$cal_probs, {

       datosim <- bd()[[1]]
   
   withBusyIndicatorServer("cal_probs", {
    arbol <- calculate_probs(datosim = datosim  , 
                             method = input$prType,
                             prev = input$prev,
                             hemogram = input$quest , 
                             pt.cut = c(input$Leucos_pc, input$plaq_pc), 
                             hmethod = input$hmog)
   })
   updateButton(session, "cal_probs",label = "Calculated", block = F, style = "success") 
   
    return(arbol)
     
  })
 
 output$arbol_prev <- renderDataTable({ 
   
   arbol()
   
 } ,options = list(
   pageLength = 5,
   scrollX = TRUE
 ))

 #calculate the summary of performance measures
 summa<- eventReactive(input$cal_probs,{
   c.roc <- ROC_curves(arbol = arbol(), 
                       n.sint = bd()$logs$sintomas , 
                       hemogram = input$quest)
   
   msg <- summary_function(data =  c.roc, 
                           hemogram = input$quest)
   return( msg )
   
 })
 
croc <- eventReactive(input$cal_probs,{
  
  c.roc <- ROC_curves(arbol = arbol(), 
                      n.sint = bd()$logs$sintomas , 
                      hemogram = input$quest)
  
  return(c.roc)

})
 #show the summary table
 output$summ <-  renderTable({
  summa()
 }, digits = 4)
 
 #download probabilities
 output$save_probs <- downloadHandler(
   filename = function(){
     paste0("dengue_probabilities_",Sys.Date(), ".csv")
     
   }
   ,content = function(file) {
     write.csv(arbol(), file, row.names = FALSE)
   })
 
#download summary
 
 output$save_summ <- downloadHandler(
   filename = function(){
     paste0("summary_probabilities_",Sys.Date(), ".csv")
     
   }
   ,content = function(file) {
     write.csv(summa(), file, row.names = FALSE)
   })
 #Plot roc curves
 output$roc_curve <- renderPlot({
  
    df <- croc()
    ## ROC curve without hemgoram
    if(input$quest == 1){
     
     plot(1-df$especi, df$sensi, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve")
     lines(seq(0,1,0.01),seq(0,1,0.01))  
     AUC <- trapz(1-df$especi, df$sensi)
     legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
     
   }
   
   if(input$quest == 2 ){
     
     par(mfrow = c(1,2))
     plot(1-df$especi, df$sensi, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve")
     lines(seq(0,1,0.01),seq(0,1,0.01))  
     AUC <- trapz(1-df$especi, df$sensi)
     legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
     
     ## ROC curve with hemgoram
     
  
     plot(1-df$especi.hem, df$sensi.hem, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve  with Hemogram")
     lines(seq(0,1,0.01),seq(0,1,0.01))  
     AUC <- trapz(1-df$especi.hem, df$sensi.hem)
     legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
     
     
   }
   
 })
 
 #OPTIONAL download the data frame to cosntruct ROC curve
 output$save_roc <- downloadHandler(
   filename = function(){
     paste0("ROC_data",Sys.Date(), ".csv")
     
   }
   ,content = function(file) {
     write.csv(croc(), file, row.names = FALSE)
   })
 
 
 #################################################
 ################ ADD NEW RECORDS ################
 #################################################
 
 #upload probailities data base (arbol)
 probs_arbol <- reactive ({
   
   req(input$arbol_in)
   tryCatch(
     {
       df<- read.csv(input$arbol_in$datapath, header = TRUE, sep = ",") 
       
       nsimt <- df %>% 
         dplyr::select(., 1:bayes) %>% 
         dplyr::select(., 1:(ncol(.)-1)) %>%
         ncol(.)
     },
     error = function(e) {
       # return a safeError if a parsing error occurs
       stop(safeError(e))
     }
   )

   return(list(probs_db = df, sint = nsimt  ))
 } )
 
 #previsualize the database
 output$arbol_in_prev <- renderDataTable({ 
   probs_arbol()[[1]]
 } ,options = list(
   pageLength = 5,
   scrollX = TRUE
 ))
 
 
new_records <- reactive ({
   
   req(input$new_records)
   tryCatch(
     {
       df<- read.csv(input$new_records$datapath, header = TRUE, sep = ",") 
       
       nsimt <- df %>% 
         dplyr::select(., 1:dxdengue) %>% 
         dplyr::select(., 1:(ncol(.)-1)) %>%
         ncol(.)
       nsimt <- ifelse(nsimt  == probs_arbol()[[2]]  ,"There is Not problems with the number of Symptoms", "Error: The number of symptoms differ from the probabilities files")
     
       vals <- df %>% 
         dplyr::select(., 1:dxdengue) %>% 
         apply(., 2, function(x){ length(unique(x))  })
       
       vals <- if_else(all( vals == 2 ) , "There is Not problems", "Error: There was found different values from  0 and 1." )
       
       hyper <-  if_else(length(grep("alpha", names(probs_arbol()[[1]]))) != 0, "Hyperparameters are Ok", "Error: Hyperparameters not found"  )
       
     },
     error = function(e) {
       # return a safeError if a parsing error occurs
       stop(safeError(e))
     }
   )
   
   return(list(newr_db = df, logs= data_frame(valores = vals, sintomas = nsimt, hyper = hyper) ))
 } )
 
 
#write and show the status of database
output$logfile2 <- renderText({
  
  if( length(grep("Error:",  new_records()$logs$valores)) != 0 ) {
    
    x1 <- paste('<li> Values:', '<font color = "red">', new_records()$logs$valores , '</font></li> ')
    
  }else{ x1 <-  paste('<li> Values:', '<font color = "green">', new_records()$logs$valores , '</font></li> ') }
  
  if( length(grep("Error:",  new_records()$logs$sintomas)) != 0 ) {
    
    x2 <-  paste('<li> Number of Symptoms:', '<font color = "red">', new_records()$logs$sintomas , '</font></li>  ')
    
  }else{ x2 <-paste('<li> Number of Symptoms:', '<font color = "green">', new_records()$logs$sintomas , '</font></li>') }
  
  if_else(length(grep("Error:",  new_records()$logs$hyper)) != 0, 
          x3 <- paste('<li> Hyperparameters:', '<font color = "red">', new_records()$logs$hyper , '</font></li>  '),
          x3 <- paste('<li> Hyperparameters:', '<font color = "green">', new_records()$logs$hyper , '</font></li>  '))
  
  return(paste(x1,"\n", x2, "\n", x3))
  
})

output$newr_in_prev <- renderDataTable({ 
  new_records()[[1]]
} ,options = list(
  pageLength = 5,
  scrollX = TRUE
))
 
 
#assign probabilitites to each new individual
probs_assigned <- eventReactive(input$assign_probs, {
 arbolx <- probs_arbol()[[1]]
 #number of symptoms
 n.sint <- probs_arbol()$sint
  
  withBusyIndicatorServer("assign_probs", {
    
    assignation <- assigns_ind(new_records = new_records()[[1]], 
                               arbol = arbolx, 
                               n.sint = n.sint , 
                               hemogram = input$quest2, 
                               pt.cut = c(input$Leucos_pc2, input$plaq_pc2),  
                               hmethod = input$hmog2)

    })
  updateButton(session, "assign_probs",label = "Probs assigned", block = F, style = "success") 
  
  return(assignation)
  
})


output$assigned_in_prev <- renderDataTable({ 
  probs_assigned()
} ,options = list(
  pageLength = 5,
  scrollX = TRUE
))

output$save_assignation <- downloadHandler(
  filename = function(){
    paste0("probs_assigned",Sys.Date(), ".csv")
    
  }
  ,content = function(file) {
    write.csv(probs_assigned(), file, row.names = FALSE)
  })
  
#calculate preformance measures when button is clicked
summ2 <- eventReactive(input$calc_summ2, {
  

 dfx <- performance_measure(assigned_records =  probs_assigned(), 
                      hemogram = input$quest2, 
                      cut.points = c(input$bayes_pc, input$bayes.hemo.pos_pc, input$bayes.hemo.neg_pc))
  
 updateButton(session, "calc_summ2",label = "Calculated", block = F, style = "success") 
  return(dfx)
})
#print the performance measures 
output$summ2 <-  renderTable({
  summ2()
}, digits = 4)

probs_updated <- eventReactive(input$update_prs, {

  #number of symptoms
  n.sint <- probs_arbol()$sint
  
  withBusyIndicatorServer("update_prs", {
    
    updated_prs <- update_probs(arbol = probs_arbol()[[1]], 
                                  assigned_records = probs_assigned(), 
                                  n.sint =  n.sint,
                                  hemogram = input$quest2)

    
  })
  updateButton(session, "update_prs",label = "Probs updated", block = F, style = "success") 
  
  return(updated_prs)
  
})
#show probabilities updated
output$updated_probs_prev <- renderDataTable({ 
  probs_updated()
} ,options = list(
  pageLength = 5,
  scrollX = TRUE
))
#download updated probabilities
output$save_updated_probs <- downloadHandler(
  
  filename = function(){
   
    paste0("updated_probs_", Sys.Date(), ".csv")
    
  }
  ,content = function(file) {
    write.csv( probs_updated(), file, row.names = FALSE)
  })


#caculate the performance measires for updated probabilities
summ_uptdated_probs<- eventReactive(input$update_prs,{
  c.roc <- ROC_curves(arbol = probs_updated(), 
                      n.sint = probs_arbol()$sint , 
                      hemogram = input$quest2)
  
  msg <- summary_function(data =  c.roc, 
                          hemogram = input$quest2)
  return( msg )
  
})

#calcultate ROC curves for updated probabilities
croc_updated_probs <- eventReactive(input$calc_croc2,{
  withBusyIndicatorServer("calc_croc2", {
  c.roc <- ROC_curves(arbol = probs_updated(), 
                      n.sint = probs_arbol()$sint , 
                      hemogram = input$quest2)
  })
  updateButton(session, "calc_croc2",label = "CROC calculated", block = F, style = "success") 
  
  return(c.roc)
  
})
#show the summary table
output$summ3 <-  renderTable({
  summ_uptdated_probs()
}, digits = 4)

output$save_summ3 <- downloadHandler(
  
  filename = function(){
    
    paste0("Performance_measures_", Sys.Date(), ".csv")
    
  }
  ,content = function(file) {
    write.csv( summ_uptdated_probs(), file, row.names = FALSE)
  })

#plot ROC curves for updated probabilities
output$roc_curve2 <- renderPlot({
  
  df <- croc_updated_probs()
  ## ROC curve without hemgoram
  if(input$quest2 == 1){
    
    plot(1-df$especi, df$sensi, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve")
    lines(seq(0,1,0.01),seq(0,1,0.01))  
    AUC <- trapz(1-df$especi, df$sensi)
    legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
    
  }
  
  if(input$quest2 == 2 ){
    
    par(mfrow = c(1,2))
    plot(1-df$especi, df$sensi, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve")
    lines(seq(0,1,0.01),seq(0,1,0.01))  
    AUC <- trapz(1-df$especi, df$sensi)
    legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
    
    ## ROC curve with hemgoram
    
    
    plot(1-df$especi.hem, df$sensi.hem, xlab="1 - Specificity",ylab="Sensitivity",type="l", main = "ROC curve  with Hemogram")
    lines(seq(0,1,0.01),seq(0,1,0.01))  
    AUC <- trapz(1-df$especi.hem, df$sensi.hem)
    legend("bottomright",paste("AUC=",round(AUC*100,1),"%"),cex=1.5,bty="n",inset=0.000000001)
    
    
  }
  
})


#download the ROc curve data frame for updated probabilities
output$save_roc2 <- downloadHandler(
  filename = function(){
    paste0("ROC_data_updated_probs",Sys.Date(), ".csv")
    
  }
  ,content = function(file) {
    write.csv(croc_updated_probs(), file, row.names = FALSE)
  })


})#end all
