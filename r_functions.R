#
# SCRIPT  CREATED BY: ANDRES CAMILO MENDEZ, DIANA CAICEDO, RAFAEL TOVAR.
# NECESSARY FUNCTIONS TO CALCULATE BAYESIAN PROBABILITIES AND 
# CLASSIFY SICK INDIVIDUALS OF DENGUE 
# 
# THIS SCRIPT COMES WITH ABSOLUTLEY NO WARRANTY
# FEEL FREE TO DISTRIBUTE IT
#





#importing packages
suppressMessages(if(!require(yaml)){install.packages('yaml'); library(yaml)} else {library(yaml)})
suppressMessages(if(!require(shiny)){install.packages("shiny");library(shiny)}else{library(shiny)})
suppressMessages(if(!require(shinydashboard)){install.packages("shinydashboard");library(shinydashboard)}else{library(shinydashboard)})
suppressMessages(if(!require(dplyr)){install.packages('dplyr'); library(dplyr)} else {library(dplyr)})
suppressMessages(if(!require(gtools)){install.packages('gtools'); library(gtools)} else {library(gtools)})
suppressMessages(if(!require(caTools)){install.packages('caTools'); library(caTools)} else {library(caTools)})
suppressMessages(if(!require(shinyjs)){install.packages('shinyjs'); library(shinyjs)} else {library(shinyjs)})


#calculate hemogram variable using the tressholds for leucos and plaq 
hemo_class <- function( data , method = c(1,2,3,4)  ) {
  # method is the way of build the variable hemogram
  # 1 => H+ = (leuco+ and plaq+) everything else will be Hemogram negative
  # 2 => H+ = (leuco+ and plaq+) also (leuco+ and plaq-)  everything else will be Hemogram negative
  # 3 => H+ = (leuco+ and plaq+) also (leuco- and plaq+)  everything else will be Hemogram negative
  # 4 => H- = (luecos- and plaqu-)   everything else will be H+
  
  if(method == 1){
   hemogram <- ifelse(data$leuco.recode == 1 & data$plaq.recode == 1, 1, 0 )
  }
  if(method == 2){
    hemogram <- ifelse(data$leuco.recode == 1 & data$plaq.recode == 1 | data$leuco.recode == 1 & data$plaq.recode == 0  , 1, 0 )
    
  }
  if(method == 3){
    hemogram <-  ifelse(data$leuco.recode == 1 & data$plaq.recode == 1 | data$leuco.recode == 0 & data$plaq.recode == 1  , 1, 0 )
    
  }
  if(method == 4){
    hemogram <-  ifelse(data$leuco.recode == 1 & data$plaq.recode == 1 | data$leuco.recode == 1 & data$plaq.recode == 0 | data$leuco.recode == 0 & data$plaq.recode == 1 , 1, 0 )
    
  }
  
  return(hemogram)
}




#function to generate the vectors with all permutations depending of the number of symptoms
vectors <- function(n.sint){
  require(gtools)
  vectors <- gtools::permutations(n=2, r=n.sint, v=c(0,1), repeats.allowed = T)
  vectors <- vectors[nrow(vectors):1, ]
 return(data.frame(vectors))
}


#function to calculate the number of sick  and healthy individual
contIndi<-function(datosim, n.sint, x){
  nm<-names(datosim)
  name <- "datosim"
  enfer<-c()
  san<-c()
  
  for(i in 1:nrow(x)){
    
    pos<- paste0(name,"$",nm[1:n.sint],"==",x[i,]," ","&"," ",collapse="")
    pos1<-substr(pos,1,nchar(pos)-3)
    pos2<-parse(text=pos1)
    enfer[i]<- nrow(subset(datosim,eval(pos2) & datosim$dxdengue == 1))
    san[i]<- nrow(subset(datosim,eval(pos2) & datosim$dxdengue == 0))
    
  }
  enfer.y.san<-cbind(enfer,san)
  return(data.frame(enfer.y.san))
}

# function to calculate the number of sick and healty individual with hemogram pos or neg
contHemo<-function(datosim, x, n.sint){
  nm<-names(datosim)
  
  hemo.pos<-c()
  hemo.neg<-c()
  hemo.pos.ND<-c()
  hemo.neg.ND<-c()
  for(i in 1:nrow(x)){
    pos<-paste0("datosim","$",nm[1:n.sint],"==",x[i,]," ","&"," ",collapse="");pos
    pos1<-substr(pos,1,nchar(pos)-3)
    pos2<-parse(text=pos1)
    
    hemo.pos[i]<-nrow(subset(datosim,eval(pos2) & datosim$dxdengue == 1 & datosim$hemogram == 1, selec = 1:n.sint) )
    hemo.neg[i]<-nrow(subset(datosim,eval(pos2) & datosim$dxdengue == 1 & datosim$hemogram == 0, selec = 1:n.sint) )
    hemo.pos.ND[i]<-nrow(subset(datosim,eval(pos2)& datosim$dxdengue == 0 & datosim$hemogram == 1, selec = 1:n.sint) )
    hemo.neg.ND[i]<-nrow(subset(datosim,eval(pos2)& datosim$dxdengue == 0 & datosim$hemogram == 0, selec = 1:n.sint) )
    
    
  }
  hemogramas<-data.frame(hemo.pos.D = hemo.pos,hemo.neg.D = hemo.neg,hemo.pos.ND = hemo.pos.ND,hemo.neg.ND = hemo.neg.ND)
  return(hemogramas)
}



# function to calculate discrete probabilities
prob.indep<-function(datosim, x, n.sint){
  
  Xi_D <- apply(datosim[, 1:n.sint], 2, function(x){ 
    
    sum( x == 1 & datosim$dxdengue == 1, na.rm = TRUE)
    
  }) #vector para almacenar los conteos de individuos con el sintoma dado que tienen dengue
  Xi_NoD <- apply(datosim[, 1:n.sint], 2, function(x){ 
    
    sum( x == 1 & datosim$dxdengue == 0, na.rm = TRUE)
    
  })#vector  para almacenar los conteos de individuos con el sintoma dado que NO tienen dengue
  
  
  PXi_D <- Xi_D/sum(datosim$dxdengue == 1, na.rm = TRUE) #se divide por el total de individuos con dengue y se almacena en un nuevo vector
  PXi_NoD <- Xi_NoD/sum(datosim$dxdengue == 0, na.rm = TRUE)#s#se divide por el total de individuos SIN dengue y se almacena en un nuevo vector
  
  z<-c() 
  z2<-c() # vectores para multiplicar las probabilidades
  for(i in 1:length(x)){
    if(x[i]==1){
      z[i]<-PXi_D[i] #se adiciona la probabilidad al vector dependiendo de las coordenadas de x
      z2[i]<-PXi_NoD[i]
    }else{
      z[i]<-1-PXi_D[i] 
      z2[i]<-1-PXi_NoD[i]
    }
  }
  
  m<-z[1] 
  m2<-z2[1]
  for(i in 2:length(z)){#multiplicar los elementos del vector
    m<-m*z[i]
    m2<-m2*z2[i]
  }
  m
  m2#multiplicacion de la prob independientes
  return(c(m,m2))
}


# function to calculate bayesian probs both case discrete and continuos
calculate_probs <- function(datosim, prev = 0.1 ,method = c(1,2), hemogram = c(1,2), pt.cut = c(4200, 165000), hmethod = c(1,2,3,4)){
  #datosim parameter is the database with individuals containing the symptoms, dxdengue and hemogram cuantitative vars
  
  #prev parameter is the prevalence of dengue disease
  
  #method parameter control the type of probability to be calculate
  # 1 = discrete probabilities
  # 2 = Continuos probabilities
  
  #hemogram parameter control if you wish to add hemogram var to the probabilities calculation process
  # 1 = No hemogram information avaliable
  # 2 = Hemogram information avaliable
  
  # pt.cut IF hemogram is == 2 then you must set two values to  the thresholds for leucos and plaqu
  # the first element of the vector correspond to leucos
  # the second element of the vectr correspond to plaque
  
  # hmethod is the way of build the variable hemogram
  # 1 => H+ = (leuco+ and plaq+) everything else will be Hemogram negative
  # 2 => H+ = (leuco+ and plaq+) also (leuco+ and plaq-)  everything else will be Hemogram negative
  # 3 => H+ = (leuco+ and plaq+) also (leuco- and plaq+)  everything else will be Hemogram negative
  # 4 => H- = (luecos- and plaqu-)   everything else will be H+
  
  
  require(dplyr)
  
  #Define number of symptoms
  n.sint <- datosim %>% 
    dplyr::select(., 1:dxdengue) %>% 
    dplyr::select(., 1:(ncol(.)-1)) %>%
    ncol(.)
  
  #generate permutations
  x <- vectors(n.sint)
  names(x) <- names(datosim)[1:n.sint]
 
    
  #Calculte  discrete probs
  if(method == 1){
    
    
    
    d<-c() #vecotres vacios para almacenar el producto de las multiplicaciones de los vectores realizados 
           #por la funcion "Prob.indep"
    d2<-c()
    bayes<-c()
    prevD<- prev #Prevalencias de Dengue
    prevND<-1- prevD # 1-Prevlanecia de dengue
    
    
    bayes <-  apply(x, 1, function(i){
      pr <- prob.indep(datosim, x = i, n.sint)
      d  <- pr[1] #se asigna al vector vacio en la posicion i el valor de la funcion
      d2 <- pr[2] #se asigna al vector vacio en la posicion i el valor de la funcion
      
      #Bayes theorem for discrete case
      return((prevD*d)/((prevD*d)+(prevND*d2)))
    })
    
    #count the number of sick and healty individual in the database
    enfer.y.san <- contIndi(datosim = datosim, n.sint = n.sint, x = x)
    
    #A la matrix de permutaciones x se le combina la columna con todas las probabilidades
    arbol <- data.frame(x, bayes, enfer.y.san)
    
    
    # if hemogram == 2 then we can use te secuential Bayes Theorem
    if(hemogram == 2){
      
      datosim2 <- datosim %>% dplyr::filter(complete.cases(.))
      #####calculate hemogram varaible ####
      #find the position for each variable in the hemogram
      leucom <- as.numeric(grep("leuco",names(datosim2)))
      plaq <-  as.numeric(grep("plaq",names(datosim2)))
      #recode leucos and plaq taking into account the thresholds in pt.cut and hmethod
      datosim2 <- datosim2 %>% dplyr::mutate(., leuco.recode = if_else(.[ , leucom] <= pt.cut[1], 1, 0 ),
                                                 plaq.recode = if_else(.[ , plaq ] <= pt.cut[2], 1, 0)) %>%
        dplyr::mutate(.,hemogram = hemo_class(data = ., method = hmethod)) %>% 
        dplyr::select(., -leuco.recode, -plaq.recode, -leuco, -plaq)
      
      hemogram.counts <- contHemo(datosim = datosim2, x ,n.sint)
      
      arbol <- data.frame(arbol, hemogram.counts)
      
      #secuencial nature of Bayes Theorem
      bayes.secu<-c()
      bayes.secu.neg<-c()
     
       for(i in 1:nrow(arbol)){
        set.seed(1000)
         #bayes secuential formula for both when the hemogram is positive and negative respectively
        bayes.secu[i]<-(arbol$bayes[i]*mean(rbeta(1000,1+arbol$hemo.pos.D[i],1+arbol$hemo.neg.D[i])))/((arbol$bayes[i]*mean(rbeta(1000,1+arbol$hemo.pos.D[i],1+arbol$hemo.neg.D[i])))+((1-arbol$bayes[i])*(mean(rbeta(1000,1+arbol$hemo.pos.ND[i],1+arbol$hemo.neg.ND[i])))))
        bayes.secu.neg[i]<-(arbol$bayes[i]*mean(rbeta(1000,1+arbol$hemo.neg.D[i],1+arbol$hemo.pos.D[i])))/((arbol$bayes[i]*mean(rbeta(1000,1+arbol$hemo.neg.D[i],1+arbol$hemo.pos.D[i])))+((1-arbol$bayes[i])*(mean(rbeta(1000,1+arbol$hemo.neg.ND[i],1+arbol$hemo.pos.ND[i])))))
        
        
      }
     arbol <- data.frame(arbol, "bayes.secu.pos" = bayes.secu , "bayes.secu.neg" = bayes.secu.neg)

     
    }
    
    
    
  }
  
  #Calculte  continuos probs
  if(method == 2){
    
    #count the sick indiviudal and healty in each vector and calculate continuos probabilities
    en.san.elici<- data.frame(x,contIndi(datosim = datosim, n.sint = n.sint, x = x)) %>% 
      dplyr::mutate(.,alpha=1+enfer, beta=1+san ) %>% 
      #bayes formula for the continuos case
      mutate(., bayes = apply(dplyr::select(., alpha, beta), 1, function(x){
        mean(rbeta(1000, x[1], x[2]))
      }) )
    # final table
    arbol <- data.frame(x,bayes=en.san.elici$bayes,enfer=en.san.elici$enfer,san=en.san.elici$san, alpha = en.san.elici$alpha, beta = en.san.elici$beta) 
    
    # if hemogram == 2 then we can use te secuential Bayes Theorem
    if(hemogram == 2){
      datosim2 <- datosim %>% dplyr::filter(complete.cases(.))
      
      leucom <- as.numeric(grep("leuco",names(datosim2)))
      plaq <-  as.numeric(grep("plaq",names(datosim2)))
      #recode leucos and plaq taking into account the thresholds in pt.cut and hmethod
      datosim2 <- datosim2 %>% dplyr::mutate(., leuco.recode = if_else(.[ , leucom] <= pt.cut[1], 1, 0 ),
                                             plaq.recode = if_else(.[ , plaq ] <= pt.cut[2], 1, 0)) %>%
        dplyr::mutate(.,hemogram = hemo_class(data = ., method = hmethod)) %>% 
        dplyr::select(., -leuco.recode, -plaq.recode, -leuco, -plaq)
      
      
      hemogram.counts <- contHemo(datosim = datosim2, x ,n.sint) 
      #create variables to calculate scuential bayes 
      arbol <- arbol %>% dplyr::mutate(.,hemo.pos.D = hemogram.counts[,1],
                                         hemo.neg.D = hemogram.counts[,2],
                                         hemo.pos.ND = hemogram.counts[,3],
                                         hemo.neg.ND = hemogram.counts[,4],
                                         alpha.secu.1 = 1 + hemogram.counts[,1],
                                         beta.secu.1 = 1 + hemogram.counts[,2],
                                         alpha.secu.2 = 1 + hemogram.counts[,3],
                                         beta.secu.2 = 1 + hemogram.counts[,4]) %>%
        # bayes secuential formula
        mutate(., bayes.secu.pos = apply(dplyr::select(., bayes, alpha.secu.1, beta.secu.1, alpha.secu.2, beta.secu.2), 1, function(x){
          ( x[1] * mean(rbeta(1000, x[2], x[3])) )/(  ( x[1] * mean(rbeta(1000, x[2], x[3])) ) + ((1-x[1])* (mean(rbeta(1000, x[4], x[5] ))) )  )
        }) ) %>%
        mutate(., bayes.secu.neg = apply(dplyr::select(., bayes, alpha.secu.1, beta.secu.1, alpha.secu.2, beta.secu.2), 1, function(x){
          ( x[1] * mean(rbeta(1000, x[3], x[2])) )/(  ( x[1] * mean(rbeta(1000, x[3], x[2])) ) + ((1-x[1])* (mean(rbeta(1000, x[5], x[4] ))) )  )
        }) ) 
        
   
    }
    
    
  }#end method 2
  
  return(arbol)
}



#function to calculate ROC curve, sensitivity, specificity and youden index
ROC_curves <- function(arbol, n.sint, hemogram = c(1,2)){
  #arbol parameter refers to the dataset with the probabilities
  
  #n.sint is the number of symptoms in de database
  
  #hemogram parameter control if you wish to add hemogram var to the probabilities calculation process
  # 1 = No hemogram information avaliable
  # 2 = Hemogram information avaliable
  
  require(dplyr)
  

   c.roc <- arbol%>% dplyr::select(., bayes, enfer, san) 
    tresh <- sort(unique(arbol$bayes), decreasing = TRUE)
   #calculate sensi and speci for the whole set of thresholds
    df <- sapply(tresh, function(i){
      
      cont.a<- sum(as.matrix(c.roc[which(c.roc$bayes >= i), "enfer"]))
      cont.b<- sum(as.matrix(c.roc[which(c.roc$bayes >= i), "san"]))
      cont.c<-sum(as.matrix(c.roc[which(c.roc$bayes < i), "enfer"]))
      cont.d<-sum(as.matrix(c.roc[which(c.roc$bayes < i), "san"]))
      
      sensi<-cont.a/(cont.a+cont.c)
      especi<-cont.d/(cont.b+cont.d)
      prev <- (cont.a+cont.b)/(cont.a+cont.b+cont.c+cont.d)
      vpp <- (prev*sensi)/(prev*sensi + (1-prev)*(1 - especi))  
      vpn <-  ((1- prev)*especi)/((1-prev)*especi + (1 - sensi)*prev)
      rp_pos <- sensi/(1-especi)
      rp_neg <- (1-sensi)/especi
      return(c(sensi,especi,prev, vpp, vpn, rp_pos, rp_neg))
    }, simplify =  FALSE ) %>% unlist() %>% matrix(., ncol = 7, byrow = TRUE) %>% cbind(tresh, .) %>% as.data.frame() 
    
    names(df) <- c("prob","sensi", "especi", "prev", "vpp", "vpn", "rp_pos", "rp_neg")
    df <- df %>% dplyr::mutate(., youden = sensi + especi - 1)
    
  
  # if the hemogran information is available then the ROC curve is calculate like follow
  if(hemogram == 2){
    
    s <- c()
    e <- c()
    prev <- c()
    vpp <- c()
    vpn <- c()
    rp_pos <- c()
    rp_neg <- c() 
    ptc1 <- sort((arbol$bayes.secu.pos),decreasing = T)
    ptc2 <- sort((arbol$bayes.secu.neg),decreasing = T)
    for(i in 1:nrow(arbol)){
      
      a<-sum(as.matrix(arbol[which(arbol$bayes.secu.pos >= ptc1[i]),"hemo.pos.D"]))+sum(as.matrix(arbol[which(arbol$bayes.secu.neg[i] >= ptc2[i]),"hemo.neg.D"]))
      b<-sum(as.matrix(arbol[which(arbol$bayes.secu.pos >= ptc1[i]),"hemo.pos.ND"]))+sum(as.matrix(arbol[which(arbol$bayes.secu.neg >= ptc2[i]),"hemo.neg.ND"]))
      c<-sum(as.matrix(arbol[which(arbol$bayes.secu.pos < ptc1[i]),"hemo.pos.D"]))+sum(as.matrix(arbol[which(arbol$bayes.secu.neg < ptc2[i]),"hemo.neg.D"]))
      d<-sum(as.matrix(arbol[which(arbol$bayes.secu.pos < ptc1[i]),"hemo.pos.ND"]))+sum(as.matrix(arbol[which(arbol$bayes.secu.neg < ptc2[i]),"hemo.neg.ND"]))
      s[i] <- a/(a+c)
      e[i] <- d/(b+d)
      prev[i] <- (a+c)/(a+b+c+d)
      vpp[i] <- (prev[i]*s[i])/(prev[i]*s[i] + (1-prev[i])*(1 - e[i]))  
      vpn[i] <-  ((1- prev[i])*e[i])/((1-prev[i])*e[i] + (1 - s[i])*prev[i])
      rp_pos[i] <- s[i]/(1-e[i])
      rp_neg[i] <- (1-s[i])/e[i]
    }
    #join all ROC curves in one single file
    df <- data.frame(df,"ptc.hem.pos" = ptc1, 
                     "ptc.hem.neg" = ptc2,
                     "sensi.hem" = s, 
                     "especi.hem" = e, 
                     "youden.hem" = s+e-1, 
                     "prev.hem" = prev,
                     "vpp.hem" = vpp,
                     "vpn.hem" = vpn,
                     "rp_pos.hem" = rp_pos,
                     "rp_neg.hem" = rp_neg)
    
  }
    
 return(df)
    
}# end ROC curves function


## function to show you the summary or performance measures for the algorithm
summary_function <- function(data, hemogram = c(1,2)){
  #data parameter refers to the dataset with the probabilities, sentivities, specificities,... called c.roc
  
  #hemogram parameter control if you wish to add hemogram var to the probabilities calculation process
  # 1 = No hemogram information avaliable
  # 2 = Hemogram information avaliable
  
  
 require(caTools)
  
 
    auc <- data %>% dplyr::select(., sensi, especi) %>% dplyr::mutate(., "1-especi" = 1-especi)%>% dplyr::select(., sensi, "1-especi") 
    auc <- caTools::trapz( auc$`1-especi`, auc$sensi  ) %>% round(., 4)
    
    res <-  data %>% dplyr::filter(., youden == max(youden)) %>% dplyr::filter(., prob == min(prob))  %>% round(., 4)
   
    # join all results in one table
   summ <- data.frame( vars = c("Results without hemogram ", "AUC", "Sensitivity", "Specificity", "Youden index", "Threshold", "Prevalence",
                                "VPP", "VPN", "Likelihood Ratio +", "Likelihood Ratio -"),
   values = c(NA ,auc, res$sensi, res$especi, res$youden, res$prob, res$prev, res$vpp, res$vpn, res$rp_pos, res$rp_neg))
   
  
  
  if( hemogram == 2 ){
    
    auc2 <- data %>% dplyr::select(., sensi.hem, especi.hem) %>% dplyr::mutate(., "1-especi.hem" = 1-especi.hem)%>% dplyr::select(., sensi.hem, "1-especi.hem") 
    auc2 <- caTools::trapz( auc2$`1-especi.hem`, auc2$sensi.hem )  %>% round(., 4)
                           
    res2 <-  data %>% dplyr::filter(., youden.hem == max(youden.hem)) %>% dplyr::filter(., ptc.hem.pos == min(ptc.hem.pos))  %>% round(., 4)
    
    # join all results in one table
    summ2 <- data.frame(vars = c("Results with hemogram ", "AUC", "Sensitivity", "Specificity", "Youden index", 
                                 "Threshold for positive hemogram", "Threshold for negative hemogram", 
                                 "Prevalence", "VPP", "VPN", "likelihood Ratio +", "Likelihood Ratio -"),
    values = c(NA ,auc2, res2$sensi.hem, res2$especi.hem, res2$youden.hem, res2$ptc.hem.pos, res2$ptc.hem.neg, res2$prev.hem, res2$vpp.hem, res2$vpn.hem, res2$rp_pos.hem, res2$rp_neg.hem))
    summ <- rbind(summ, summ2) 

          
  } 

  return(summ)
}

#############################################################
########### HOW TO CALCULATE A NEW ALGORITHM ################
############################################################

#Dengue_calculator process (Run from here)

#cargar la base de datos
# datosim <- read.csv("C:/Users/acmendez/Google Drive (andres.mendez@correounivalle.edu.co)/TESIS/TESIS/shiny_app/dengue_claculator/BD_prototipo_1.csv", header = TRUE)
#remove NA values
# datosim <- datosim %>% filter(complete.cases( dplyr::select(., 1:dxdengue) ) )

#First step
#Generate combiantions od symptoms and calculate probs
#x <- vectors(n.sint = 12)
#arbol <- calculate_probs(datosim = datosim , prev = 0.1 , method = 2, hemogram = 2 , pt.cut = c(4200, 165000), hmethod = 1)

#Second step 
#calcualte ROC curve
#c.roc <- ROC_curves(arbol, n.sint = 12, hemogram = 2)

#fourth step
#calculate summary 
# summ <- summary_function(data = c.roc, hemogram = 2) 

#fifth step
#to save the results use:
#write.csv(arbol, "c:/User/.../arbol_25_jun_12sint.csv" , header = TRUE, row.names = FALSE)




#------------------------------------------------------------------------------------------------------------

#=*=*=*=*=*=**=*=*=*=*=*=*=*=*=*=*==**=*=*=*=*==*=*=*=*=*=*=*==*=*=*=*=*=*==*=*=*=*
#=*=*=*=*=*=*=*=*=**==*ADD NEW RECORDS AND RECALCULATE PROBS =*=*=*=*=*=*=*=*=*=*=*
#=*=*=*==*=*=*=*=*=*=*==*==*=*=*=*=*=*=*=*=*=*=*=*=*=**=*=*=*=*==*=*=*=*=*=*=*=*=*=

################################ !!!!!IMPORTANT !!!!!  ############################
#----------- THIS SECTION ONLY WORKS WITH CONTINUOS BAYESIAN ROBABILITIES -------#



#Function to assign probabilities to each individual in the database
assigns_ind <- function(new_records, arbol, n.sint, hemogram = c(1,2), pt.cut = c(4200, 165000),  hmethod = c(1,2,3,4)){ 
  #new_records parameter is a dataframe with the new indidivuald to be clasiffied
  #arbol parameter is a dataframe with the probabilities for each combination of symptoms
  # n.sint parameter is the total number of Symptons 
  #hemogram parameter control if you wish to add hemogram var to the probabilities calculation process
  # 1 = No hemogram information avaliable
  # 2 = Hemogram information avaliable
  
  # pt.cut IF hemogram is == 2 then you must set two values to  the thresholds for leucos and plaqu
  # the first element of the vector correspond to leucos
  # the second element of the vectr correspond to plaque
  
  # hmethod is the way of build the variable hemogram
  # 1 => H+ = (leuco+ and plaq+) everything else will be Hemogram negative
  # 2 => H+ = (leuco+ and plaq+) also (leuco+ and plaq-)  everything else will be Hemogram negative
  # 3 => H+ = (leuco+ and plaq+) also (leuco- and plaq+)  everything else will be Hemogram negative
  # 4 => H- = (luecos- and plaqu-)   everything else will be H+
  
  #create data frame with the whole combinations between the symptoms
  x <- vectors(n.sint)
  
  
  arbolx <- arbol %>% na.omit() %>% dplyr::mutate(., code = apply(dplyr::select(., 1:n.sint),1, function(x){  
    paste0(x,collapse = "")
  } )) %>% as_tibble()
  
  #assign probabilities to each new individual
  new_records <- new_records %>% na.omit() %>% dplyr::mutate(., code = apply(dplyr::select(., 1:n.sint),1, function(x){  
    paste0(x,collapse = "")
  } )) %>% as_tibble() %>% 
    dplyr::mutate(., bayes = unlist(apply(dplyr::select(., code), 1, function(x){
      
      dplyr::filter(arbolx, code == x) %>% dplyr::select(., bayes) 
      
    })  )) %>% as_tibble()
  #If the new individuals have hemogram information then
  if(hemogram == 2){
    
    
    leucom <- as.numeric(grep("leuco",names(new_records)))
    plaq <-  as.numeric(grep("plaq",names(new_records)))
    
    #recode leucos and plaq taking into account the thresholds in pt.cut and hmethod
    new_records <- new_records %>% dplyr::mutate(., leuco.recode = if_else(.[ , leucom] <= pt.cut[1], 1, 0 ),
                                                 plaq.recode = if_else(.[ , plaq ] <= pt.cut[2], 1, 0)) %>%
      dplyr::mutate(.,hemogram = hemo_class(data = ., method = hmethod)) %>% 
      dplyr::select(., -leuco.recode, -plaq.recode, -leuco, -plaq) 
    
    #assign hemogram probabilities to each new individual
    new_records <- new_records %>%
      dplyr::mutate(., bayes.secu.pos = unlist(apply(dplyr::select(., code), 1, function(x){
        
        dplyr::filter(arbolx, code == x) %>% dplyr::select(., bayes.secu.pos) 
        
      })  )) %>%
      dplyr::mutate(., bayes.secu.neg = unlist(apply(dplyr::select(., code), 1, function(x){
        
        dplyr::filter(arbolx, code == x) %>% dplyr::select(., bayes.secu.neg) 
        
      })  )) 
    
    
    
  }
  
  new_records <- new_records %>% dplyr::select(., -code)  
  
  return(new_records)
}#end function to assign new individuals


#assigned_records <- assigns_ind(new_records = new_records, arbol = arbol, n.sint = 12 , hemogram = 2, pt.cut = c(4200, 165000),  hmethod = 1)

#function to calculate the performance measures for the assigned probabilities to the new records
performance_measure <- function(assigned_records, hemogram = c(1, 2), cut.points = c(0.5, 0.5, 0.5)){

  a.b <- assigned_records %>% dplyr::filter(., bayes >= cut.points[1] ) %>% dplyr::select(., dxdengue) %>% table() %>% .[2:1]
  c.d <- assigned_records %>% dplyr::filter(., bayes < cut.points[1] ) %>% dplyr::select(., dxdengue) %>% table()%>% .[2:1]
  CM <- matrix(c(a.b, c.d), 2, 2, byrow = TRUE)  
  
 sensi <- rbeta(10000, 1 + CM[1,1], 1 + CM[2,1] ) %>% mean()     
 speci <- rbeta(10000, 1 + CM[2,2], 1 + CM[1,2] ) %>% mean()  
 prev <- (CM[1,1]+CM[1,2])/sum(CM, na.rm =TRUE)
 vpp <- (prev*sensi)/(prev*sensi + (1-prev)*(1 - speci))  
 vpn <-  ((1- prev)*speci)/((1-prev)*speci + (1 - sensi)*prev)
 rp_pos <- sensi/(1-speci)
 rp_neg <- (1-sensi)/speci
 
 df <- tibble( variables = c("Confusion Matrix", "a", "b", "c", "d", "Sensitivity", "Specificity","Prevalence", "VPP", "VPN", "Likelihood Ratio +", "Likelihood Ratio -"),
               values = c(NA, CM[1,1], CM[1,2], CM[2,1], CM[2,2], sensi, speci, prev, vpp,  vpn, rp_pos, rp_neg)  )
 if( hemogram == 2){
   #conunts to calculate sensitivity and specificity when the hemogram is avaliable
   assigned_records <- assigned_records %>%  dplyr::mutate(., hemo.pos.D = if_else(dxdengue == 1 & hemogram == 1 , 1, 0)
                                                           ,hemo.pos.ND = if_else(dxdengue == 0 & hemogram == 1, 1, 0 )
                                                           ,hemo.neg.D = if_else(dxdengue == 1 & hemogram == 0, 1, 0)
                                                           ,hemo.neg.ND = if_else(dxdengue == 0 & hemogram == 0, 1, 0))
   
   a<-sum(assigned_records[which(assigned_records$bayes.secu.pos >= cut.points[[2]]),"hemo.pos.D"])+sum(assigned_records[which(assigned_records$bayes.secu.neg >= cut.points[[3]]),"hemo.neg.D"])
   b<-sum(assigned_records[which(assigned_records$bayes.secu.pos >= cut.points[[2]]),"hemo.pos.ND"])+sum(assigned_records[which(assigned_records$bayes.secu.neg >= cut.points[[3]]),"hemo.neg.ND"])
   c<-sum(assigned_records[which(assigned_records$bayes.secu.pos < cut.points[[2]]),"hemo.pos.D"])+sum(assigned_records[which(assigned_records$bayes.secu.neg < cut.points[[3]]),"hemo.neg.D"])
   d<-sum(assigned_records[which(assigned_records$bayes.secu.pos < cut.points[[2]]),"hemo.pos.ND"])+sum(assigned_records[which(assigned_records$bayes.secu.neg < cut.points[[3]]),"hemo.neg.ND"])
   
   sensi.hem <- rbeta(10000, 1 + a, 1 + c ) %>% mean()     
   speci.hem <- rbeta(10000, 1 + d, 1 + b ) %>% mean()  
   prev.hem <- (a+c)/(a+b+c+d)
   vpp.hem <- (prev.hem*sensi.hem)/(prev.hem*sensi.hem + (1-prev.hem)*(1 - speci.hem))  
   vpn.hem <-  ((1- prev.hem)*speci.hem)/((1-prev.hem)*speci.hem + (1 - sensi.hem)*prev.hem)
   rp_pos.hem <- sensi.hem/(1-speci.hem)
   rp_neg.hem <- (1-sensi.hem)/speci.hem
 
   df <- bind_rows(df,  tibble( variables = c("Confusion Matrix hemogram","a", "b", "c", "d", "Sensitivity", "Specificity","Prevalence", "VPP", "VPN", "Likelihood Ratio +", "Likelihood Ratio -" ),
                                values = c(NA,a, b, c, d, sensi.hem, speci.hem, prev.hem, vpp.hem,  vpn.hem, rp_pos.hem, rp_neg.hem)) )
   
 }
return(df)
 }


update_probs <- function(arbol, assigned_records, n.sint ,hemogram = c(1, 2)){
  # arbol parameter is a dataframe with the continuos bayesian probabilities and counts of sick and healty individuals
  #with or without hemogram
  #also arbol dataframe must have the hyperparameters of the beta distribution (alpha, beta)
  
  #assigned_records parameter is the data frame obtained by the function assigns_ind  
  
  #hemogram parameter control if you wish to add hemogram var to the probabilities calculation process
  # 1 = No hemogram information avaliable
  # 2 = Hemogram information avaliable
  
  x <- vectors(n.sint)
  
  
  arbolx <- arbol %>% dplyr::mutate(., code = apply(dplyr::select(., 1:n.sint),1, function(x){  
    paste0(x,collapse = "")
  } )) %>% as_tibble()
  
 
res_new_indi <- assigned_records %>% na.omit() %>% dplyr::mutate(., code = apply(dplyr::select(., 1:n.sint),1, function(x){  
    paste0(x,collapse = "")
  } ))  %>% as_tibble()  %>% 
  dplyr::group_by(., code)  %>% 
  summarise(enfer.y = sum(dxdengue == 1), 
            san.y = sum(dxdengue == 0))
  
#uptade probabilitites using the new records
arbolx <- left_join(arbolx, res_new_indi, by = "code")%>% 
  dplyr::mutate_at(., vars(contains(".y")), funs(replace(., is.na(.), 0)) )  %>% 
  dplyr::mutate(., enfer = enfer + enfer.y,  
                san = san + san.y, 
                alpha = alpha + enfer.y, 
                beta = beta + san.y, 
                #Bayes secuential formula
                bayes = apply(dplyr::select(., alpha, beta),1, function(x){
                  mean(rbeta(10000,x[1], x[2] ), na.rm= TRUE)
                }))  %>% dplyr::select(., -enfer.y , - san.y)


if(hemogram == 2){
  
  #conunts to calculate the beta poterior distribution when the hemogram is avaliable
  assigned_records <- assigned_records %>%  dplyr::mutate(., hemo.pos.D = if_else(dxdengue == 1 & hemogram == 1 , 1, 0)
                                      ,hemo.pos.ND = if_else(dxdengue == 0 & hemogram == 1, 1, 0 )
                                      ,hemo.neg.D = if_else(dxdengue == 1 & hemogram == 0, 1, 0)
                                      ,hemo.neg.ND = if_else(dxdengue == 0 & hemogram == 0, 1, 0))
  
  
  res_new_indi <- assigned_records %>% na.omit() %>% dplyr::mutate(., code = apply(dplyr::select(., 1:n.sint),1, function(x){  
    paste0(x,collapse = "")
  } ))  %>% as_tibble()  %>% 
    dplyr::group_by(., code)  %>% 
    summarise(hemo.pos.D.y = sum(hemo.pos.D == 1), 
              hemo.neg.D.y = sum(hemo.neg.D == 1),
              hemo.pos.ND.y = sum(hemo.pos.ND == 1),
              hemo.neg.ND.y = sum(hemo.neg.ND == 1))
  
  #uptade probabilitites using the new records when hemogram is available
  arbolx <- left_join(arbolx, res_new_indi, by = "code")%>% 
    dplyr::mutate_at(., vars(contains(".y")), funs(replace(., is.na(.), 0)) ) %>% 
    dplyr::mutate(., hemo.pos.D = hemo.pos.D + hemo.pos.D.y,  
                  hemo.neg.D = hemo.neg.D + hemo.neg.D.y, 
                  hemo.pos.ND = hemo.pos.ND + hemo.pos.ND.y, 
                  hemo.neg.ND = hemo.neg.ND + hemo.neg.ND.y,
                  alpha.secu.1 = alpha.secu.1 + hemo.pos.D.y,
                  beta.secu.1 = beta.secu.1 + hemo.neg.D.y,
                  alpha.secu.2 = alpha.secu.2 + hemo.pos.ND.y,
                  beta.secu.2 =  beta.secu.2 + hemo.neg.ND.y,
                  #Bayes secuential formula
                  bayes.secu.pos = apply(dplyr::select(., bayes, alpha.secu.1, beta.secu.1, alpha.secu.2, beta.secu.2), 1, function(x){
                    ( x[1] * mean(rbeta(1000, x[2], x[3])) )/(  ( x[1] * mean(rbeta(1000, x[2], x[3])) ) + ((1-x[1])* (mean(rbeta(1000, x[4], x[5] ))) )  )
                  }),
                  bayes.secu.neg = apply(dplyr::select(., bayes, alpha.secu.1, beta.secu.1, alpha.secu.2, beta.secu.2), 1, function(x){
                    ( x[1] * mean(rbeta(1000, x[3], x[2])) )/(  ( x[1] * mean(rbeta(1000, x[3], x[2])) ) + ((1-x[1])* (mean(rbeta(1000, x[5], x[4] ))) )  )
                  }) 
                  ) %>% dplyr::select(.,  -hemo.pos.D.y, -hemo.neg.D.y, -hemo.pos.ND.y, -hemo.neg.ND.y )
  
  
  
  
}

arbolx <- arbolx %>% dplyr::select(., -code)
return(arbolx)
  
}#end function update probs


#################################################################
######## HOW TO ADD NEW RECORDS TO THE ALGORITHM ################
################################################################

#data base with the new individuals to assign probabilities
#new_records <- read.csv("C:/Users/.../new_records.csv")
#Upload the probabities file called normally arbol
#arbol <- read.csv("C:/Users/.../arbol.csv")

#first step
#assign probabilities to the new individuals
#assigned_records <- assigns_ind(new_records = new_records, arbol = arbol, n.sint = 12 , hemogram = 2, pt.cut = c(4200, 165000),  hmethod = 1)

#Second step
#calculate performance measures (sensitivity, specificity)
# summary <- performance_measure(assigned_records =  assigned_records, hemogram = 2, cut.points = c(0.5, 0.05, 0.05))

#Third step 
#Update the algorithm with the new records
# arbol_new <- update_probs(arbol = arbol, assigned_records = assigned_records, n.sint = 12 ,hemogram = 2)

#Fourth step 
#Calculate ROC curves and performance measures for the updated algorithm
#c.roc <- ROC_curves(arbol = arbol2, n.sint = 12, hemogram = 2)

















