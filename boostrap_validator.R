suppressMessages(if(!require(tidyverse)){install.packages("tidyverse");library(tidyverse)}else{library(tidyverse)})

#pts.corte <- c(pt.bayes, pt.bayes.secu.pos, pt.bayes.secu.neg)
bosstrap_sampler <- function(arbol , pts.corte = c(NA, NA , NA), use.hemogram = TRUE, n.sample = 1000){

  n.sint <- which(names(arbol) == "bayes") -1 
  cat("Creating data base from Arbol \n")
  if(use.hemogram){
    
   bd <- apply(arbol, 1, function(x){
     
     if(any(is.na(pts.corte))   ){stop("One of the cut points is missing")}
     
      x <- as.numeric(x)
      enfer <- x[n.sint +2]
      san <-   x[n.sint +3]
      #print(x[n.sint+1])
      hemo.pos.D  <- x[n.sint + 6]
      hemo.neg.D  <- x[n.sint + 7]
      hemo.pos.ND <- x[n.sint + 8]
      hemo.neg.ND <- x[n.sint + 9]
      
      
      tot <- enfer + san
      
      pob <- base::matrix(data = x[1:n.sint] , nrow = tot, ncol = n.sint, byrow = F)
      pob <- as.data.frame(pob)
      names(pob) <- names(arbol)[1:n.sint]
      
      pob$dxdengue <- c(rep(1, enfer) , rep(0, san))
      
      
      hemo.D <- c(rep(1, hemo.pos.D), rep(0, hemo.neg.D))
      hemo.ND <- c(rep(1, hemo.pos.ND), rep(0, hemo.neg.ND))
      
      if(length(hemo.D) < enfer){
        hemo.D <- c(hemo.D, rep(NA, enfer - length(hemo.D)))
      }
      if(length(hemo.ND) < san){
        hemo.ND <- c(hemo.ND, rep(NA, san - length(hemo.ND)))
      }
      
      pob$hemogra  <-c(hemo.D, hemo.ND)
      
      
      if(length(pob$dxdengue) != 0){
        pob$bayes <- x[n.sint +1]  
        pob$bayes.secu.pos <- x[length(x)-1]
        pob$bayes.secu.neg <- x[length(x)]
      }else{
        pob$bayes <- numeric(0)
        pob$bayes.secu.pos <- numeric(0)
        pob$bayes.secu.neg <- numeric(0)
        
        }
      
      
      return(pob)
    })
    
   
    
  }else{
    
    bd <- apply(arbol, 1, function(x){
      x <- as.numeric(x)
      enfer <- x[n.sint +2]
      san <-   x[n.sint +3]

      tot <- enfer + san
      
      pob <- base::matrix(data = x[1:n.sint] , nrow = tot, ncol = n.sint, byrow = F)
      pob <- as.data.frame(pob)
      names(pob) <- names(arbol)[1:n.sint]
      
      pob$dxdengue <- c(rep(1, enfer) , rep(0, san))

      if(length(pob$dxdengue) != 0){
        pob$bayes <- x[n.sint +1]  
      }else{
        pob$bayes <- numeric(0)
        
      }
      
      
      return(pob)
    })
    

  }
 
  
  bd <- do.call(rbind, bd)
  
  cat("Creating Resamples \n")
  boostrap_sampler <- tibble(resample = 1:n.sample, sampled_db =  lapply(1:n.sample, function(i){
    sample_n(as_tibble(bd), size = nrow(bd), replace = TRUE)
  }) )
  cat("Calculating performance measures for each resample \n")
 
   boostrap_sampler <- boostrap_sampler %>% dplyr::mutate(., metrics_no_hemo = purrr::map(.x = sampled_db, .f = function(.x){
    
    
    a <- .x %>% dplyr::filter(., dxdengue == 1 & bayes >= pts.corte[1]) %>% dplyr::select(., dxdengue) %>% nrow(.)
    b <- .x %>% dplyr::filter(., dxdengue == 0 & bayes >= pts.corte[1]) %>% dplyr::select(., dxdengue) %>% nrow(.)
    c <- .x %>% dplyr::filter(., dxdengue == 1 & bayes < pts.corte[1]) %>% dplyr::select(., dxdengue) %>% nrow(.)
    d <- .x %>% dplyr::filter(., dxdengue == 0 & bayes < pts.corte[1]) %>% dplyr::select(., dxdengue) %>% nrow(.)
    
    sensi <- a/(a+c)
    especi <- d/(b+d)
    prev <- (a+b)/(a+b+c+d)
    vpp <- (prev*sensi)/(prev*sensi + (1-prev)*(1 - especi))  
    vpn <-  ((1- prev)*especi)/((1-prev)*especi + (1 - sensi)*prev)
    rp_pos <- sensi/(1-especi)
    rp_neg <- (1-sensi)/especi
    
    res <- data.frame("sensi"=sensi, "especi"= especi, "prev"=prev,"vpp"=  vpp, "vpn" = vpn, "lr_pos" = rp_pos, "lr_neg" = rp_neg)
    return(res)
    
  }))
  boostrap_metrics <- do.call(rbind, boostrap_sampler$metrics_no_hemo)
  boostrap_metrics <- data.frame('Metric' = names(boostrap_metrics), 
                                           'Mean'   = colMeans(boostrap_metrics),
                                           'Ic_low' = apply(boostrap_metrics, 2, function(x){quantile(x, probs = c(.025))}),
                                           'Ic_upp' = apply(boostrap_metrics, 2, function(x){quantile(x, probs = c(.975))}))

   
   if(use.hemogram){
     boostrap_sampler <- boostrap_sampler %>% dplyr::mutate(., metrics_si_hemo =  purrr::map(.x = sampled_db, .f = function(.x){
       
       a1 <- .x %>% dplyr::filter(.x$bayes.secu.pos >= pts.corte[2] & .x$dxdengue == 1 & .x$hemogra == 1) %>% nrow()
       a2 <- .x %>% dplyr::filter(.x$bayes.secu.neg >= pts.corte[3] & .x$dxdengue == 1 & .x$hemogra == 0) %>% nrow()
       a <- a1+ a2
       
       b1 <- .x %>% dplyr::filter(.x$bayes.secu.pos >= pts.corte[2] & .x$dxdengue == 0 & .x$hemogra == 1) %>% nrow()
       b2 <- .x %>% dplyr::filter(.x$bayes.secu.neg >= pts.corte[3] & .x$dxdengue == 0 & .x$hemogra == 0) %>% nrow()
       b <- b1+b2
       
       c1 <- .x %>% dplyr::filter(.x$bayes.secu.pos < pts.corte[2] & .x$dxdengue == 1 & .x$hemogra == 1) %>% nrow()
       c2 <- .x %>% dplyr::filter(.x$bayes.secu.neg < pts.corte[3] & .x$dxdengue == 1 & .x$hemogra == 0) %>% nrow()
       c <- c1 + c2
       
       d1 <- .x %>% dplyr::filter(.x$bayes.secu.pos < pts.corte[2] & .x$dxdengue == 0 & .x$hemogra == 1) %>% nrow()
       d2 <- .x %>% dplyr::filter(.x$bayes.secu.neg < pts.corte[3] & .x$dxdengue == 0 & .x$hemogra == 0) %>% nrow()
       d <- d1 + d2

       sensi <- a/(a+c)
       especi <- d/(b+d)
       prev <- (a+b)/(a+b+c+d)
       vpp <- (prev*sensi)/(prev*sensi + (1-prev)*(1 - especi))  
       vpn <-  ((1- prev)*especi)/((1-prev)*especi + (1 - sensi)*prev)
       rp_pos <- sensi/(1-especi)
       rp_neg <- (1-sensi)/especi
       
       res <- data.frame("sensi"=sensi, "especi"= especi, "prev"=prev,"vpp"=  vpp, "vpn" = vpn, "lr_pos" = rp_pos, "lr_neg" = rp_neg)
       return(res)
     }) 
   
     )
     
     boostrap_metrics_hemograma <- do.call(rbind, boostrap_sampler$metrics_si_hemo)
     boostrap_metrics_hemograma <- data.frame('Metric' = names(boostrap_metrics_hemograma), 
                                              'Mean'   = colMeans(boostrap_metrics_hemograma),
                                              'Ic_low' = apply(boostrap_metrics_hemograma, 2, function(x){quantile(x, probs = c(.025))}),
                                              'Ic_upp' = apply(boostrap_metrics_hemograma, 2, function(x){quantile(x, probs = c(.975))}))
   }else{
     boostrap_metrics_hemograma <- NA
   }
   
  return(list("bostrap_metrics_no_hemograma" = boostrap_metrics, "boostrap_metrics_si_hemograma" = boostrap_metrics_hemograma) )
   
}



######## HOW TO USE

#load Data
arbol <- read.csv("C:/Users/acmendez/Downloads/dengue_probabilities_2019-01-12_12sint_retro_artre_ratio_arinorrea46000186000.csv")
#cut points for the algorithm
pts.corte = c(0.4, 0.45 , 0.3)

results <- bosstrap_sampler(arbol = arbol , 
                            pts.corte = pts.corte,#cut points c(cut point bayes, cut point bayes.hemo.pos, cut point bayes.hemo.neg) 
                            use.hemogram = TRUE, #control if the hemogram has to be proccesed 
                            n.sample = 1000)










pacman::p_load(raster, rasterVis, maptools, tidyverse, xlsx, latticeExtra, sp)


mThm <- rasterTheme(region = brewer.pal(9, "YlGn"),
                    panel.background = list(col = "#708090"))
data(wrld_simpl)

# Creating reference lines
grat <<- sp::gridlines(wrld_simpl, easts = seq(-180, 180, by = 10), norths = seq(-90, 90, by = 15))

ext_ctm <- c(-123.017 , -34.442,  -41.280,   34.291)#extent para que salgan bien los graficos


sdm_meso <- raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/andean/americas/prj_models/andean_prj_median.tif")


cost_dist <- raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/andean/americas/gap_models/gap_class_final.tif")

cost_dist <- raster::mask(cost_dist, sdm_meso)
r_temp <- raster::crop(cost_dist, ext_ctm)
#r_temp <- r_temp/max(na.omit(r_temp[]))

scale <- c(quantile(r_temp[], probs = seq(0,1, 0.045), na.rm = T),1)
p <- rasterVis::levelplot(r_temp,
                          at = scale,
                          margin = F,
                          par.settings = mThm,
                          maxpixels = ncell(r_temp)) + 
  latticeExtra::layer(sp.lines(grat, lwd = 0.5, lty = 2, col = "white")) +
  latticeExtra::layer(sp.polygons(wrld_simpl, lwd = 0.8, col = "black", fill = "#4B4B4B"), under = T)



png("D:/OneDrive - CGIAR/Desktop/mapas/andean_gap_score_delaunay.png", height = 7, width = 10, units = "in", res = 300)
print(p)
dev.off()



cost_dist <- raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/andean/americas/gap_models/gap_class_final.tif")

rsin <- raster::crop(cost_dist, ext_ctm)
rsin <- ratify(rsin)
rat <- levels(rsin)[[1]]
rat$level <- c('Well conserved', 'Gap by one approach', 'Gap by both approaches')
levels(rsin) <- rat

#figure details
ht <- 12
fct <- (rsin@extent@xmin-rsin@extent@xmax)/(rsin@extent@ymin-rsin@extent@ymax)
wt <- ht*(fct+.1)

# My theme
mThm <- rasterTheme(region = brewer.pal(9, "YlGn"),
                    panel.background = list(col = "#708090"))


cols <- c('grey70', 'goldenrod3', 'red2')

p <- rasterVis:::levelplot(rsin,
                           att = 'level',
                           margin = F,
                           par.settings = mThm,
                           col.regions = cols,
                           maxpixels = ncell(rsin)) + 
  latticeExtra::layer(sp.lines(grat, lwd = 0.5, lty = 2, col = "white")) +
  latticeExtra::layer(sp.polygons(wrld_simpl, lwd = 0.8, col = "black", fill = "#4B4B4B"), under = T)


png("D:/OneDrive - CGIAR/Desktop/mapas/andean_gap_class_final.png", height = 7, width = 10, units = "in", res = 300)
print(p)
dev.off()



gc_andean <- raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/andean/americas/gap_models/gap_class_final.tif")
gc_meso   <- raster("Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/mesoamerican/americas/gap_models/gap_class_final.tif")

gc_andean[which(gc_andean[] == 1)]<- 0
gc_andean[which(gc_andean[] == 2)]<- 1
gc_andean[which(gc_andean[] == 0)]<- 999


gc_meso[which(gc_meso[] == 1)] <- 0
gc_meso[which(gc_meso[] == 0)] <- 999

stk <- stack(gc_meso, gc_andean)
f_rast <- raster::calc(stk, sum, na.rm = TRUE)
f_rast[which(f_rast[] == 0)]<- NA
f_rast[which(f_rast[] >= 999)]<- 0


rsin <- raster::crop(f_rast, ext_ctm)
rsin <- ratify(rsin)
rat <- levels(rsin)[[1]]
rat$level <- c('Well conserved', 'Andean Gaps', 'Mesoamerican Gaps', "Overlap")
levels(rsin) <- rat


writeRaster(rsin, "Z:/gap_analysis_landraces/runs/results/common_bean/lvl_1/common_bean_final_gap_map.tif", overwrite = TRUE)


k2 <- raster("Z:/gap_analysis_landraces/runs/results/rice_african/lvl_1/K2/africa/gap_models/gap_class_final.tif") 
k4 <- raster("Z:/gap_analysis_landraces/runs/results/rice_african/lvl_1/K4/africa/gap_models/gap_class_final.tif")
k5 <- raster("Z:/gap_analysis_landraces/runs/results/rice_african/lvl_1/K5/africa/gap_models/gap_class_final.tif")



lista <- list(k2, k4,k5)
names(lista ) <- c("k2", "k4", "k5")

l <- lapply(1:length(lista), function(i){
  x <- list[[i]]
  names(x)<- names(list[[i]])
  df<- as.data.frame(x, xy = T)
  colnames(df) <- c("lon", "lat", names(x))
  return(df)
  
})

df <- reduce(l, left_join, by = c("lon", "lat"))

##### gap richness ######

crop <- "wheat_durum"
level <- paste0(c(2,3,4))
region <- "world"

msk <- raster(paste0("//dapadfs/Workspace_cluster_9/gap_analysis_landraces/runs/input_data/mask/mask_", region, ".tif"))

gp_finals <- lapply(1:length(level), function(i){
  
  gap_dir <- paste0("Z:/gap_analysis_landraces/runs/results/", crop, "/lvl_1/", level[i], "/", region, "/gap_models/")  
  
  gp_final <- raster(paste0(gap_dir,"gap_class_final.tif"))
  gp_final[which(gp_final[] != 2)] <- 0
  gp_final[which(gp_final[] == 2)] <- 1
  
  
  df<- as.data.frame(gp_final, xy = T)
  colnames(df) <- c("lon", "lat", level[i])
  
  return(df)
})

df <- reduce(gp_finals, left_join, by = c("lon", "lat"))

df$gap_richness <- apply(df[, 3:ncol(df)], 1, function(x){
  
  ifelse( all(is.na(x)), sum <- NA, sum <- sum(x, na.rm = TRUE) )
  return(sum)
})

to_rast <- df %>% dplyr::select(., lon ,lat, gap_richness)
g_richness_map <- rasterFromXYZ(to_rast, res = res(msk), crs = crs(msk))

rsin <- (g_richness_map)
rat <- levels(rsin)[[1]]
rat$level <- c('Well conserved', 'Gaps in one race', "Gaps in two races", "Gaps in all races") # "Gaps in two races", "Gaps in three races",
levels(rsin) <- rat

writeRaster(rsin, paste0("Z:/gap_analysis_landraces/runs/results/", crop, "/lvl_1/", crop,"_gap_richness.tif"))


wmask <- raster::raster(paste0(mask_dir, "/mask_world", ".tif" ))


crop<-c("banana", "sorghum", "potato", "barley", "maize", "african_maize", "wheat_durum", "wheat_bread", "rice_african", "common_bean")

level <- "lvl_1"
crops<-lapply(1:length(crop), function(j){
  cat(j, "\n")
  
  x <- paste0(results_dir, "/", crop[j], "/", level, "/", crop[j], "_gap_richness.tif")
  rst <- raster(x)
  rst <- raster::crop(rst, extent(wmask))
  
  df <- as.data.frame(rst, xy = T) 
  colnames(df) <- c("lon", "lat", names(rst))
  df[,3] <- as.numeric(df[,3])
  df <- df[!is.na(df[,3]),]
  
  #rst[which(rst[] == 0)] <- 999
  return(df)
  
})

df <- purrr::reduce(crops, full_join, by = c("lon", "lat"))


df$gap_richness <- apply(df[, 3:ncol(df)], 1, function(x){
  
  ifelse( all(is.na(x)), sum <- NA, sum <- sum(x, na.rm = TRUE) )
  return(sum)
})

to_rast <- df %>% dplyr::select(., lon ,lat, gap_richness)
g_richness_map <- rasterFromXYZ(to_rast, res = res(wmask), crs = crs(wmask))

writeRaster(g_richness_map, "Z:/gap_analysis_landraces/runs/gap_richness_map.tif", overwrite = TRUE)

#### para el landrace richness

stk1 <- raster::stack("Z:/gap_analysis_landraces/runs/models.tif")

sdms_rich <- lapply(1:nlayers(stk1), function(i){
  cat(i, "\n")
  x <- stk1[[i]]
  x[!is.na(x[])] <- 1
  df <- as.data.frame(x, xy = T) 
  colnames(df) <- c("lon", "lat", names(x))
  #df[,3] <- as.numeric(df[,3])
  df <- df[!is.na(df[,3]),]
  return(df)
})

df <- purrr::reduce(sdms_rich, full_join, by = c("lon", "lat"))

df2 <- apply(df[, -c(1,2)], 2, function(x){
  x[which(!is.na(x))] <- 1
  return(x)
})
df_final <- data.frame(df[,1:2], df2)

df_final$gap_richness <- apply(df_final[, 3:ncol(df_final)], 1, function(x){
  
  ifelse( all(is.na(x)), sum <- NA, sum <- sum(x, na.rm = TRUE) )
  return(sum)
})

to_rast <- df_final %>% dplyr::select(., lon ,lat, gap_richness)
g_richness_map <- rasterFromXYZ(to_rast, res = res(wmask), crs = crs(wmask))

writeRaster(g_richness_map, "Z:/gap_analysis_landraces/runs/landrace_richness_map.tif", overwrite = TRUE)

