rm(list = ls())
require(shiny)

require(ade4)
require(maptools)
require(spdep)
require(sp) 
require(rgdal) 
require(tripack) 
require(RColorBrewer) 
require(DCluster) 
require(MASS) 
require(dispmod) 
require(SpatialPack) 
require(nlme) 
require(mgcv) 
require(sf) 
require(ggplot2) 
require(viridis) 
require(tmap) 
require(cartogram)


boxmap=function(vnam,df,legtitle=NA,mtitle="Box Map",mult=1.5){
  # box map
  # arguments:
  #   vnam: variable name (as character, in quotes)
  #   df: simple features polygon layer
  #   legtitle: legend title
  #   mtitle: map title
  #   mult: multiplier for IQR
  # returns:
  #   a tmap-element (plots a map)
  var <- get.var(vnam,df)
  bb <- boxbreaks(var)
  tm_shape(df) +
    tm_fill(vnam,title=legtitle,breaks=bb,palette=rev(brewer.pal(6,"RdBu")),
            labels = c("lower outlier", "< 25%", "25% - 50%", "50% - 75%","75%", "upper outlier"))  +
    tm_borders() +
    tm_layout(title = mtitle, title.position = c("right","bottom"))
}

boxbreaks=function(v,mult=1.5) {
  # break points for box map
  # arguments:
  #   v: vector with observations
  #   mult: multiplier for IQR (default 1.5)
  # returns:
  #   bb: vector with 7 break points
  # compute quartile and fences
  qv <- unname(quantile(v))
  iqr <- qv[4] - qv[2]
  upfence <- qv[4] + mult * iqr
  lofence <- qv[2] - mult * iqr
  # initialize break points vector
  bb <- vector(mode="numeric",length=7)
  # logic for lower and upper fences
  if (lofence < qv[1]) {  # no lower outliers
    bb[1] <- lofence
    bb[2] <- floor(qv[1])
  } else {
    bb[2] <- lofence
    bb[1] <- qv[1]
  }
  if (upfence > qv[5]) { # no upper outliers
    bb[7] <- upfence
    bb[6] <- ceiling(qv[5])
  } else {
    bb[6] <- upfence
    bb[7] <- qv[5]
  }
  bb[3:5] <- qv[2:4]
  return(bb)
}

load("muncol.RData")

ui <- fluidPage(h4("Spatial modelling of infant leukaemia in Colombia"),tabsetPanel(
  tabPanel(h4("Exploratory spatial data analysis"),   
  fluidRow(
  column(width = 4),
  column(width = 4,
         selectInput('xcol', 'Variable', choices=list("Standardised morbidity ratio"="SMR_NCancer",
                                                      "Armed conflict incidence index"="ACII",
                                                      "Unsatisfied basic needs index"="UBN",
                                                      "Percentage of people living in rural areas"="Rurality",
                                                      "Percentage of people with health insurance coverage"="Insurance"),
                     selected = "SMR_NCancer")),
  column(width = 4)
  ),
  fluidRow(
    column(width = 6,         plotOutput('plot1')),
    column(width = 6,         plotOutput('plot2')       
    )
 ),
 fluidRow(
   column(width = 6,          plotOutput('plot3')
   ),
   column(width = 6,          plotOutput('plot4'))
 )
  )))

server = function(input, output, session) { 
  output$plot1 <- renderPlot({ 
    muncol_st=st_as_sf(muncol)
    ggplot() + geom_sf(data=muncol_st, aes(fill=input$xcol)) +
      theme_minimal()
  })
  
  
  output$plot2 <- renderPlot({ 
    muncol_st=st_as_sf(muncol)
    boxmap(input$xcol,muncol_st,mtitle = " ")
  })
  

  
  output$plot3 <- renderPlot({
   moran.plot(muncol[[input$xcol]], matrices[[input$xcol]],xlab=input$xcol,ylab = paste("spatially lagged",input$xcol))
  })
 
  
  output$plot4 <- renderPlot({ 
    
  locm <- localmoran(muncol[[input$xcol]], matrices[[input$xcol]])  #calculate the local moran's I
  
  # manually make a moran plot standarize variables
  muncol[[paste0("s",input$xcol)]] <- scale(muncol[[input$xcol]])  #save to a new column
  
  # create a lagged variable
  muncol[[paste0("lag_s",input$xcol)]] <- lag.listw(matrices[[input$xcol]], muncol[[paste0("s",input$xcol)]])
  
  muncol$quad_sig <- NA
  muncol@data[(muncol[[paste0("s",input$xcol)]] >= 0 & muncol[[paste0("lag_s",input$xcol)]] >= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 1
  muncol@data[(muncol[[paste0("s",input$xcol)]] <= 0 & muncol[[paste0("lag_s",input$xcol)]] <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 2
  muncol@data[(muncol[[paste0("s",input$xcol)]] >= 0 & muncol[[paste0("lag_s",input$xcol)]] <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 3
  muncol@data[(muncol[[paste0("s",input$xcol)]] >= 0 & muncol[[paste0("lag_s",input$xcol)]] <= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 4
  muncol@data[(muncol[[paste0("s",input$xcol)]] <= 0 & muncol[[paste0("lag_s",input$xcol)]] >= 0) & (locm[, 5] <= 0.05), "quad_sig"] <- 5 
  
  # Set the breaks for the thematic map classes
  breaks <- seq(1, 5, 1)
  
  # Set the corresponding labels for the thematic map classes
  labels <- c("high-High", "low-Low", "High-Low", "Low-High", "Not Signif.")
  
  # see ?findInterval - This is necessary for making a map
  np <- findInterval(muncol$quad_sig, breaks)
  
  # Assign colors to each map class
  colors <- c("red", "blue", "lightpink", "skyblue2", "white")
  plot(muncol, col = colors[np])  #colors[np] manually sets the color for each county
  mtext("Local Moran's I", cex = 1.5, side = 3, line = 1)
  legend("topleft", legend = labels, fill = colors, bty = "n")
  
  
  })
  
  
  
}

shinyApp(ui = ui, server = server)
