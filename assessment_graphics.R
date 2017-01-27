###############################################################################
# LIBRARIES
###############################################################################
library(extrafont)
library(ggplot2)
library(plyr)
library(latticeExtra)
library(gridExtra)
library(Hmisc)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(coda)


#------------------------------------------------------------------------------
# Graphics defaults
#------------------------------------------------------------------------------
#font_import()
#loadfonts(device="win")
windowsFonts(Times=windowsFont("TT Calibri"))
theme_set(theme_bw(base_size=18,base_family="Calibri")+
            #theme(legend.position = c(0.25,0.9))+
            theme(legend.direction = "vertical")+
            theme(legend.margin = unit(1,"cm"))+
            theme(axis.title.y=element_text(vjust=1))+
            theme(axis.title.x=element_text(vjust=-0.5))+
            theme(legend.title = element_blank())+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))


#------------------------------------------------------------------------------
# Data calls
#------------------------------------------------------------------------------
.REP <- "model" 
source("globals.R")


#------------------------------------------------------------------------------
# Model comparisons
#
# 'model'   = global model with Sex-specific structure
# 'model_M' = global model with Sex-Unified structure
#------------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# ADMB report files - output defined in FINAL_CALCS section of .tpl file
#
# NOTE: you must have all four ADMB files written here: .par, .std, .cor, .rep
#
#------------------------------------------------------------------------------
model    <- read.admb("model_output/model")
modelM   <- read.admb("model_output/modelM")


# ADMB .std files (contain parameter values and standard deviation)
#       - all defined parameters with whatever 'sdreport_vector / matrix'
#         calls are defined in the .tpl file
par   <- read.delim("model_output/model.std",header=TRUE,sep="")$value
var   <- read.delim("model_output/model.std",header=TRUE,sep="")$std
parM  <- read.delim("model_output/modelM.std",header=TRUE,sep="")$value
varM  <- read.delim("model_output/modelM.std",header=TRUE,sep="")$std

# MCMC - primary parameters to check distributions and bounds
# BOOT - primary parameters to check distributions and bounds
# These can be either from the bootstrap or MCMC
theta    <- read.delim("model_output/boot.txt",header=FALSE,sep="")
thetaM   <- read.delim("model_output/bootM.txt",header=FALSE,sep="")

# 
# # MCMC - derived quantities
if (file.exists("MR.ps"))
MR            <- read.delim("model_output/MR.ps",header=FALSE,sep="")
if (file.exists("MRM.ps"))
MRM           <- read.delim("model_output/MRM.ps",header=FALSE,sep="")
if (file.exists("MR_srv.ps"))
MR_srv        <- read.delim("model_output/MR_srv.ps",header=FALSE,sep="")
if (file.exists("MR_srvM.ps"))
MR_srvM       <- read.delim("model_output/MR_srvM.ps",header=FALSE,sep="")
if (file.exists("recruit.ps"))
recruit       <- read.delim("model_output/recruit.ps",header=FALSE,sep="")
if (file.exists("recruitM.ps"))
recruitM      <- read.delim("model_output/recruitM.ps",header=FALSE,sep="")
if (file.exists("spbm.ps"))
spbiomass     <- read.delim("model_output/spbm.ps",header=FALSE,sep="")
if (file.exists("spbmM.ps"))
spbiomassM    <- read.delim("model_output/spbmM.ps",header=FALSE,sep="")



#------------------------------------------------------------------------------
#
# RETROSPECTIVE ANALYSES 
#   You'll have to run these individually and manually
#
#-------------------------------------------------------------------------------
# par_R1m   <- read.delim("model_output/model_r1m.std",header=TRUE,sep="")$value
# var_R1m   <- read.delim("model_output/model_r1m.std",header=TRUE,sep="")$std
# par_R2m   <- read.delim("model_output/model_r2m.std",header=TRUE,sep="")$value
# var_R2m   <- read.delim("model_output/model_r2m.std",header=TRUE,sep="")$std
# par_R4m   <- read.delim("model_output/model_r4m.std",header=TRUE,sep="")$value
# var_R4m   <- read.delim("model_output/model_r4m.std",header=TRUE,sep="")$std
# par_R5m   <- read.delim("model_output/model_r5m.std",header=TRUE,sep="")$value
# var_R5m   <- read.delim("model_output/model_r5m.std",header=TRUE,sep="")$std
# par_R9m   <- read.delim("model_output/model_r9m.std",header=TRUE,sep="")$value
# var_R9m   <- read.delim("model_output/model_r9m.std",header=TRUE,sep="")$std
# par_R10m  <- read.delim("model_output/model_r10m.std",header=TRUE,sep="")$value
# var_R10m  <- read.delim("model_output/model_r10m.std",header=TRUE,sep="")$std
# 
# par_R1   <- read.delim("model_output/model_r1.std",header=TRUE,sep="")$value
# var_R1   <- read.delim("model_output/model_r1.std",header=TRUE,sep="")$std
# par_R2   <- read.delim("model_output/model_r2.std",header=TRUE,sep="")$value
# var_R2   <- read.delim("model_output/model_r2.std",header=TRUE,sep="")$std
# par_R6   <- read.delim("model_output/model_r6.std",header=TRUE,sep="")$value
# var_R6   <- read.delim("model_output/model_r6.std",header=TRUE,sep="")$std
# par_R9   <- read.delim("model_output/model_r9.std",header=TRUE,sep="")$value
# var_R9   <- read.delim("model_output/model_r9.std",header=TRUE,sep="")$std
# par_R10  <- read.delim("model_output/model_r10.std",header=TRUE,sep="")$value
# var_R10  <- read.delim("model_output/model_r10.std",header=TRUE,sep="")$std


#------------------------------------------------------------------------------
# Parse and sort MCMC for variances
#   These already have every 100th saved, and if I recall correctly, as ADMB
#   begins the draws from the approximately multivariate normal of the Hessian
#   matrix, a deleted burn-in isn't stricly necessary, but as I am a
#   statistical dunce, I still do it.
#
# MCMC outputs of 1,000,000 draws w/100th saved =  10,000
#  Trim first 20% = 8,000 remaining
#------------------------------------------------------------------------------
# theta      <- theta[c(2001:10000),]
# thetaM     <- thetaM[c(2001:10000),]

# MR         <- MR[c(2001:10000),]
# MRM        <-  MRM[c(2001:10000),]
# MR_srv     <- MR_srv[c(2001:10000),]
# MR_srvM    <-  MR_srvM[c(2001:10000),]
# recruit    <- recruit[c(2001:10000),]
# recruitM   <- recruitM[c(2001:10000),]
# spbiomass  <- spbiomass[c(2001:10000),]
# spbiomassM <- spbiomassM[c(2001:10000),]
# 
# # Sort each MCMC matrix by column and select 5% and 95%
# # ASSUMES ALL MCMC MATRICES ARE OF THE SAME DIMENSIONS
# min<-0.05*length(theta[,1])
# max<-0.95*length(theta[,1])
# 
# for(i in 1:10)
# {
#   theta[,i]<-sort(theta[,i]) 
#   thetaM[,i]<-sort(thetaM[,i]) 
# }
# 
# for(i in 1:36)
# {
#   MR[,i]         <- sort(MR[,1])
#   MRM[,i]        <- sort(MRM[,1])
#   MR_srv[,i]     <- sort(MR_srv[,1])
#   MR_srvM[,i]    <- sort(MR_srvM[,1])
#   recruit[,i]    <- sort(recruit[,1])
#   recruitM[,i]   <- sort(recruitM[,1])
#   spbiomass[,i]  <- sort(spbiomass[,1])
#   spbiomassM[,i] <- sort(spbiomassM[,1])
# }
# 
# 
# 
# #------------------------------------------------------------------------------
# # THETA PARAMETERS - these are currently from the bootstrap - kvk 1/26/2017
# #------------------------------------------------------------------------------
parMcmc <-data.frame(Par  = c(unlist(theta[,1:6]),
                              unlist(thetaM[,1:6])),
                      Name = c(rep(c(rep("Mean Recruitment",10000),
                               rep("Mean Year 1 N",10000),rep("Mean F",10000),
                               rep("Catchability (commercial)",10000),
                               rep("Catchability (survey 1)",10000),
                               rep("Catchability (survey 2)",10000)),2)),
                      Model =  c(rep("Sex-specific",60000),
                                 rep("Unified",60000)))

plot.parmc <- function(df) {
  ggplot(subset(df, Name %in% c("Mean Recruitment", "Mean Year 1 N",
                                      "Mean F",
                                      "Catchability (commercial)",
                                      "Catchability (survey 1)",
                                      "Catchability (survey 2)",
                                      "Age 50% sel.")),
         aes(Par,group=interaction(Name,Model))) + geom_density(aes(Par, y=..scaled.., fill = Model),alpha=.2)+
         facet_wrap(~Name,scales="free",ncol=2)+
         labs(y="")+
         labs(x="")+
         theme(legend.margin=unit(-0.25, "cm"))+
         guides(fill = guide_legend(ncol=3))+
         theme(legend.position = 'bottom')+
         theme(legend.key.size = unit(0.5, "cm"))+
         scale_fill_manual(name='', values=c('Sex-specific' = 'black',
                                             'Unified' = 'blue'))
}


#------------------------------------------------------------------------------
# Mark-recapture fishery abundance
#------------------------------------------------------------------------------
dfmr <- data.frame(Year  = rep(c(2003,2004,2005,2006,2007,2008,2009,2010,2012,2013,2015),2), 
                  Density = c(par[382:392],parM[273:283]),
                  Stdev   = c(var[382:392],varM[273:283]),
                  Model = c(rep("Sex-specific",11), rep("Unified",11)),
                  MR    = c(model$obs_MR_abund,modelM$obs_MR_abund),
                  Var   = c(model$obs_MR_var,modelM$obs_MR_var))

plot.mr <- function(df) {
  df <- df %>%
    mutate(p.var = 1.96 * sqrt(Var)) %>%
    mutate(p.lower = MR - p.var, p.upper = MR + p.var)
  df <- df %>% 
    mutate(err= 1.96 * Stdev) %>%
    mutate(lower = Density - err, upper = Density + err)
  
  
  ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    geom_point(aes(x=Year,y= MR, colour=Model),size=3, alpha = 0.3)+
    geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper, colour=Model),lwd=0.8, alpha = 0.3,width=0.4)+
    facet_grid(~ Model)+
    scale_x_continuous(breaks = df$Year)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.6))+
     scale_colour_manual(name='', values=c('Sex-specific' = 'black', 
                                           'Unified' = "blue"))+
     scale_fill_manual(name='', values=c('Sex-specific' = 'black', 
                                         'Unified' = "blue"))+
    labs(x="Year")+
    labs(y="Mark-Recapture abundance - fishery")+
    theme(legend.position = "null")+
    #theme(legend.direction = "vertical")+
    theme(legend.key = element_rect(fill = "white"))
  
}


#------------------------------------------------------------------------------
# Mark-recapture survey abundance
#------------------------------------------------------------------------------
dfmrs <- data.frame(Year  = rep(c(2003,2004,2005,2006,2007,2008,2009,2010,2012,2013,2015),2), 
                  Density = c(par[393:403],parM[284:294]),
                  Stdev   = c(var[393:403],varM[284:294]),
                  Model = c(rep("Sex-specific",11), rep("Unified",11)),
                  MR    = c(model$obs_MR_abund_SRV,modelM$obs_MR_abund_SRV),
                  Var   = c(model$obs_MR_var_SRV,modelM$obs_MR_var_SRV))

plot.mrs <- function(df) {
  df <- df %>%
    mutate(p.var = 1.96 * sqrt(Var)) %>%
    mutate(p.lower = MR - p.var, p.upper = MR + p.var)
  df <- df %>% 
    mutate(err= 1.96 * Stdev) %>%
    mutate(lower = Density - err, upper = Density + err)
  
  
  ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    geom_point(aes(x=Year,y= MR, colour=Model),size=3, alpha = 0.3)+
    geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper, colour=Model),lwd=0.8, alpha = 0.3,width=0.4)+
    facet_grid(~ Model)+
    scale_x_continuous(breaks = df$Year)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.6))+
    scale_colour_manual(name='', values=c('Sex-specific' = 'black', 
                                          'Unified' = "blue"))+
    scale_fill_manual(name='', values=c('Sex-specific' = 'black', 
                                        'Unified' = "blue"))+
    labs(x="Year")+
    labs(y="Mark-Recapture abundance - survey")+
    theme(legend.position = "null")+
    theme(legend.key = element_rect(fill = "white"))
  
}




#------------------------------------------------------------------------------
# Spawning Biomass (tons)
#------------------------------------------------------------------------------
dfA <- data.frame(Year  = rep(seq(1980,2015,by=1),2), 
                  Density = c(par[310:345],parM[201:236]),
                  Stdev   = c(var[310:345],varM[201:236]),
                  Model = c(rep("Sex-specific",36), rep("Unified",36)))

plot.A <- function(df) {
  df <- df %>% 
    mutate(err= 1.96 * Stdev/1000) %>%
    mutate(lower = Density/1000 - err, upper = Density/1000 + err)
  
  
  ggplot(df,aes(Year,Density/1000,group=Model)) + geom_line(aes(Year,Density/1000,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
    scale_colour_manual(name='', values=c('Sex-specific' = 'black', 
                                          'Unified' = "blue"))+
    scale_fill_manual(name='', values=c('Sex-specific' = 'black', 
                                        'Unified' = "blue"))+
    labs(x="Year")+
    labs(y="Biomass (metric tons)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))
  
}


#------------------------------------------------------------------------------
# RECRUITMENT
#------------------------------------------------------------------------------
dfR <- data.frame(Year  = rep(seq(1980,2015,by=1),2),
                  Density = c(par[346:381],parM[237:272]),
                  Stdev   = c(var[346:381],varM[237:272]),
                  Model = c(rep("Sex-specific",36),rep("Unified",36)))

plot.rec <- function(df) {
  df <- df %>%
    mutate(err= 1.96 * Stdev/1000) %>%
    mutate(lower = Density/1000 - err, upper = Density/1000 + err)


  ggplot(df,aes(Year,Density/1000,group=Model)) + geom_line(aes(Year,Density/1000,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) +
    scale_colour_manual(name='', values=c('Sex-specific' = 'black',
                                          'Unified' = "blue"))+
    scale_fill_manual(name='', values=c('Sex-specific' = 'black',
                                        'Unified' = "blue"))+
    labs(x="Year")+
    labs(y="Recruitment (thousands)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    theme(legend.key.size = unit(1., "cm"))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))

}


#------------------------------------------------------------------------------
# SPAWNING BIOMASS
#------------------------------------------------------------------------------
dfS <- data.frame(Year  = rep(seq(1980,2015,by=1),2),
                  Density = c(par[310:345],parM[201:236]),
                  Stdev   = c(var[310:345],varM[201:236]),
                  Model = c(rep("Sex-specific",36),rep("Unified",36)))

plot.spawners <- function(df) {
  df <- df %>%
    mutate(err= 1.96 * Stdev/1000) %>%
    mutate(lower = Density/1000 - err, upper = Density/1000 + err)


  ggplot(df,aes(Year,Density/1000,group=Model)) + geom_line(aes(Year,Density/1000,colour=Model),size=1)+
    geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) +
    scale_colour_manual(name='', values=c('Sex-specific' = 'black',
                                          'Unified' = "blue"))+
    scale_fill_manual(name='', values=c('Sex-specific' = 'black',
                                        'Unified' = "blue"))+
    labs(x="Year")+
    labs(y="Spawning biomass (metric tons)")+
    theme(legend.position = c(0.8,0.88))+
    theme(legend.direction = "vertical")+
    #theme(legend.position = c(0.8,0.8))+
    theme(legend.key.size = unit(1., "cm"))+
    #ggtitle("Total density (individuals per square kilometer)")+
    theme(legend.key = element_rect(fill = "white"))

}

# 


#------------------------------------------------------------------------------
# Abundance at age
#------------------------------------------------------------------------------
# yearT <- seq(1980,2015,by=1)
# yearX <- c('1988','1989',seq(1991,2005,by=1),seq(2008,2015,by=1))
age <- seq(2,42,by=1)


SSf <- data.frame(Year = rep(as.numeric(model$yy),length(model$aa)),
                 Age = rep(as.numeric(model$aa), each=length(model$yy)),
                 p = as.vector(model$natagef))

SSm <- data.frame(Year = rep(as.numeric(model$yy),length(model$aa)),
                 Age = rep(as.numeric(model$aa), each=length(model$yy)),
                 p = as.vector(model$natagem))

SU <- data.frame(Year = rep(as.numeric(model$yy),length(model$aa)),
                 Age = rep(as.numeric(model$aa), each=length(model$yy)),
                 p = as.vector(modelM$natage))


png(file='figures/age.png', res=300, width=7, height=9, units ="in", bg="transparent")
par(mfrow = c(3,1))
symbols(SSf$Year, SSf$Age, SSf$p, inches=0.2, xlim=c(1980,2015), ylim=c(2,42),
        xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Abundance at age females")
symbols(SSm$Year, SSm$Age, SSm$p, inches=0.2, xlim=c(1980,2015), ylim=c(2,42),
        xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Abundance at age males")
symbols(SU$Year, SU$Age, SU$p, inches=0.2, xlim=c(1980,2015), ylim=c(2,42),
        xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Abundance at age unified")
  dev.off()




  #------------------------------------------------------------------------------
  # Catch at age (looking for recruitment indicators)
  #------------------------------------------------------------------------------

  f <- data.frame(Year = rep(as.numeric(model$yrs_fish_age),length(model$aa)),
                   Age = rep(as.numeric(model$aa), each=length(model$yrs_fish_age)),
                   p = as.vector(model$oac_fishf))
  
  m <- data.frame(Year = rep(as.numeric(model$yrs_fish_age),length(model$aa)),
                  Age = rep(as.numeric(model$aa), each=length(model$yrs_fish_age)),
                  p = as.vector(model$oac_fishm))
  
  u <- data.frame(Year = rep(as.numeric(modelM$yrs_fish_age),length(modelM$aa)),
                  Age = rep(as.numeric(modelM$aa), each=length(modelM$yrs_fish_age)),
                  p = as.vector(modelM$oac_fish))



  png(file='figures/catch_age_dataf.png', res=300, width=9, height=7, units ="in", bg="transparent")
  symbols(f$Year, f$Age, f$p, inches=0.2, xlim=c(2002,2015), ylim=c(2,42),
          xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition females")
  dev.off()

  png(file='figures/catch_age_datam.png', res=300, width=9, height=7, units ="in", bg="transparent")
  symbols(m$Year, m$Age, m$p, inches=0.2, xlim=c(2002,2015), ylim=c(2,42),
          xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition males")
  dev.off()

  png(file='figures/catch_age_datau.png', res=300, width=9, height=7, units ="in", bg="transparent")
  symbols(u$Year, u$Age, u$p, inches=0.2, xlim=c(2002,2015), ylim=c(2,42),
          xlab="Year", ylab = "Observed age", cex.axis=1.5, cex.lab=1.5, cex.main=1.,main="Observed catch composition unified")

  dev.off()
  #




  
  #------------------------------------------------------------------------------
  # Selectivity
  #------------------------------------------------------------------------------
  dfSel <- data.frame(Year =  rep(seq(2,25,by=1),5),
                      Sel  =  c(model$fish_sel_f[1:24],model$fish_sel_m[1:24],
                                modelM$fish_sel[1:24],model$vfish_sel_f[1:24],
                                model$vfish_sel_m[1:24]),
                      Sex   =  c(rep("Female",24),rep("Male",24),rep("Unified",24),
                                 rep("NOAA - female",24),rep("NOAA - male",24)))

  plot.sel <- function(df) {

    ggplot(df,aes(Year,Sel, group=Sex, colour=Sex)) + geom_line(aes(Year,Sel,group=Sex, colour=Sex),size=1)+
      scale_colour_manual(name='', values=c('Female' = 'red',
                                            'Male' = 'blue',
                                            'Unified' = 'black',
                                            "NOAA - female" = 'firebrick',
                                            "NOAA - male" = "green"))+
      labs(x="Age")+
      labs(y="Selectivity")+
      theme(legend.position = c(0.75,0.8))+
      theme(legend.direction = "vertical")+
      theme(legend.key = element_rect(fill = "white"))

  }
  
  
  #------------------------------------------------------------------------------
  # Fishing mortality
  #------------------------------------------------------------------------------
  dfF <- data.frame(Year =   rep(seq(1980,2015,by=1),3),
                    Mort  =  c(model$Fmortf,model$Fmortm,modelM$Fmort),
                    Model =  c(rep("Female",36),rep("Male",36),rep("Unified",36)))

  plot.Fmt <- function(df) {

    ggplot(df,aes(Year,Mort, group=Model, colour=Model)) +
      geom_line(aes(Year,Mort,group=Model, colour=Model))+
      geom_hline(yintercept=0.07)+
      scale_colour_manual(name='', values=c('Female' = 'red',
                                            'Male' = 'blue',
                                            'Unified' = 'black'))+
      labs(x="Year")+
      labs(y="Mortality")+
      theme(legend.position = c(0.25,0.85))+
      theme(legend.direction = "vertical")+
      theme(legend.key = element_rect(fill = "white"))

  }


  #------------------------------------------------------------------------------
  # Catch age residuals
  #------------------------------------------------------------------------------
#   yearT <- c('1988','1989',seq(1991,2005,by=1),seq(2008,2015,by=1))
#   age <- seq(8,75,by=1)
#   
# 
#   dfs <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
#                     Age = rep(as.numeric(age), each=length(yearT)),
#                     p = as.vector(model$res_age))
#   
#   dfm <- data.frame(Year = rep(as.numeric(yearT),length(age)), 
#                     Age = rep(as.numeric(age), each=length(yearT)),
#                     p = as.vector(modelM$res_age))
#   
#   
#   plot.age <- function(dfs,dfm) {
#     
#  A <-   ggplot(df,aes(Year,Age,size=p))+
#       geom_point(aes(colour=p))+
#       scale_size_area(max_size = 5, guide=FALSE)+
#       scale_colour_gradient2(limits=c(-0.04,0.08),low="red", high="blue", mid = "grey", midpoint = 0.0)+
#       xlab("Age")+
#       ylab("Residuals (o - p)")+
#       theme(legend.margin=unit(0.2, "cm"))+
#       theme(legend.position = "right")
#     
# B <-    ggplot(dfm,aes(Year,Age,size=p))+
#       geom_point(aes(colour=p))+
#       scale_size_area(max_size = 5, guide=FALSE)+
#       scale_colour_gradient2(limits=c(-0.04,0.08),low="red", high="blue", mid = "grey", midpoint = 0.0)+
#       xlab("Age")+
#       ylab("Residuals (o - p)")+
#       theme(legend.margin=unit(0.2, "cm"))+
#       theme(legend.position = "right")
# 
# grid.arrange(A,B)
#     
#   }
#   
  
  #------------------------------------------------------------------------------
  # COMMERCIAL CPUE 
  #------------------------------------------------------------------------------
  dfcp <- data.frame(Year  =  as.factor(c(rep(seq(1980,2015,by=1),3))),
                     CPUE   =  c(model$obs_fshy_cpue,model$pred_fshy_cpue,modelM$pred_fshy_cpue),
                     Var    =  c(model$obs_fshy_cpue_var,var[432:467],varM[323:358]),
                     Model  =  c(rep("Observed",36),rep("Sex-specific",36),rep("Unified",36)))


  plot.cpue <- function(df) {
    df <- df %>%
      mutate(err= 1.96 * Var) %>%
      mutate(lower = CPUE - err, upper = CPUE + err)

    ggplot(df,aes(Year,CPUE, group=Model, colour=Model)) +
      geom_line(data=df,aes(Year,CPUE,group=Model, colour=Model),lwd=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model, colour=NA),alpha=0.2, show.legend=FALSE) +
      geom_point(data=df, aes(Year,CPUE,group=Model),size=1)+
      labs(y="Fishery CPUE (pounds per hook)")+
      scale_x_discrete(breaks=seq(1980,2015,by=4))+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
      theme(legend.position = c(0.4,0.85))+
      theme(legend.direction="vertical")+
      scale_colour_manual(name='', values=c('Sex-specific' = 'black', "Observed" = 'firebrick', 'Unified' = 'blue'))+
      scale_fill_manual(name='', values=c('Sex-specific' = 'black', "Observed" = 'firebrick', 'Unified' = 'blue'))

  }
  # 
  # 
  # 
  # #------------------------------------------------------------------------------
  # # IPHC CPUE
  # #------------------------------------------------------------------------------
  dfcp2<- data.frame(Year  =  as.factor(c(rep(seq(1997,2015,by=1),3))),
                     CPUE   =  c(model$obs_cpue2_biom,model$pred_cpue2,modelM$pred_cpue2),
                     Var    =  c(model$obs_cpue2_var,var[413:431],varM[304:322]),
                     Model  =  c(rep("Observed",19),rep("Sex-specific",19),
                                 rep("Unified",19)))

  plot.cpue2 <- function(df) {
    df <- df %>%
      mutate(err= 1.96 * Var) %>%
      mutate(lower = CPUE - err, upper = CPUE + err)

    ggplot(df,aes(Year,CPUE, group=Model, colour=Model)) +
      geom_line(data=df,aes(Year,CPUE,group=Model, colour=Model), lwd=1)+
      geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model, colour=NA),alpha=0.2, show.legend=FALSE) +
      geom_point(data=df,aes(Year,CPUE,group=Model),size=1)+
      labs(y=("Survey CPUE (individuals per hook)"))+
      scale_x_discrete(breaks=seq(1997,2015,by=4))+
      theme(legend.position = c(0.7,0.85))+
      theme(legend.background = element_blank())+
      theme(legend.margin = unit(0.5,"cm"))+
      theme(legend.direction="vertical")+
      theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
      scale_colour_manual(name='', values=c('Sex-specific' = 'black', "Observed" = 'firebrick', "Unified" = "blue"))+
      scale_fill_manual(name='', values=c('Sex-specific' = 'black', "Observed" = 'firebrick', "Unified" = 'blue'))

  }


  
# #------------------------------------------------------------------------------
# # Spawning biomass projections 
# #------------------------------------------------------------------------------
#   dfproj <- data.frame(Year  = rep(seq(2016,2030,by=1),2), 
#                        Density = c(par[367:381],parM[366:380]),
#                        Error  =   c(var[367:381],varM[366:380]),
#                        Model = c(rep("Recommended F = 0.022 (F_65) Sex-specific",15),
#                                  rep("Recommended F = 0.022 (F_65) Unified",15)))
#   
#   
#   plot.proj<- function(df) {
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       scale_colour_manual(name='', values=c('Recommended F = 0.022 (F_65) Sex-specific' = 'black', 
#                                             'Recommended F = 0.022 (F_65) Unified' = 'blue'))+
#       scale_fill_manual(name='', values=c('Recommended F = 0.022 (F_65) Estiated M' = 'black',  
#                                           'Recommended F = 0.022 (F_65) Unified' = 'blue'))+
#       labs(x="Year")+
#       ylim(0,6000)+
#       labs(y="Projected spawning biomass")+
#       theme(legend.position = c(0.5,0.2))+
#       theme(legend.direction = "vertical")+
#       theme(legend.key.size = unit(1., "cm"))+
#       theme(legend.background = element_blank())+
#       #ggtitle("Projected spawning biomass (mt)")+
#       theme(legend.key = element_rect(fill = "white"))
#     
#   }
#   
#   
#   #------------------------------------------------------------------------------
#   # DENSITY GLOBAL RETROSPECTIVE
#   #------------------------------------------------------------------------------
#   dfDR <- data.frame(Year  = rep(seq(1980,2015,by=1),6), 
#                      Density = c(par[283:363],
#                                  par_R1[217:246],NA,
#                                  par_R2[213:241],NA,NA,
#                                  par_R6[197:221],NA,NA,NA,NA,NA,NA,
#                                  par_R9[185:206],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                  par_R10[181:201],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                      Error = c(var[283:363],
#                                var_R1[217:246],NA,
#                                var_R2[213:241],NA,NA,
#                                var_R6[197:221],NA,NA,NA,NA,NA,NA,
#                                var_R9[185:206],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                var_R10[181:201],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                      Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),rep("retro 6",36),
#                                rep("retro 9",36), rep("retro 10",36)),
#                      Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),6),
#                      Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),6))
#   
#   dfDR$Model <- factor(dfDR$Model, as.character(dfDR$Model))
#   
#   
#   plot.retrodensity <- function(df) {
#     df <- df %>%
#       mutate(p.var = 1.96 * Var) %>%
#       mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       geom_point(aes(x=Year,y= Rov), colour="black",size=3)+
#       geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                             'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                           'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       labs(x="Year")+
#       labs(y="Density (ind. per square kilometer)")+
#       theme(legend.position = c(0.87,0.75))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#   }  
#   
#   
#   #------------------------------------------------------------------------------
#   # DENSITY GLOBAL RETROSPECTIVE Unified
#   #------------------------------------------------------------------------------
#   dfDRm <- data.frame(Year  = rep(seq(1980,2015,by=1),7), 
#                      Density = c(parM[282:362],
#                                  par_R1m[216:245],NA,
#                                  par_R2m[212:240],NA,NA,
#                                  par_R4m[204:230],NA,NA,NA,NA,
#                                  par_R5m[200:225],NA,NA,NA,NA,NA,
#                                  par_R9m[184:205],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                  par_R10m[180:200],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                      Error = c(varM[282:362],
#                                var_R1m[216:245],NA,
#                                var_R2m[212:240],NA,NA,
#                                var_R4m[204:230],NA,NA,NA,NA,
#                                var_R5m[200:225],NA,NA,NA,NA,NA,
#                                var_R9m[184:205],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                var_R10m[180:200],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                      Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),rep("retro 4",36),rep("retro 5",36),
#                                rep("retro 9",36), rep("retro 10",36)),
#                      Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),7),
#                      Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),7))
#   
#   dfDRm$Model <- factor(dfDRm$Model, as.character(dfDRm$Model))
#   
#   
#   plot.retrodensitym <- function(df) {
#     df <- df %>%
#       mutate(p.var = 1.96 * Var) %>%
#       mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       geom_point(aes(x=Year,y= Rov), colour="black",size=3)+
#       geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
#                                             'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
#                                           'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       labs(x="Year")+
#       labs(y="Density (ind. per square kilometer)")+
#       theme(legend.position = c(0.87,0.75))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#   }
#   
#   
#   
#   #------------------------------------------------------------------------------
#   # SPAWNING BIOMASS GLOBAL RETROSPECTIVE
#   #------------------------------------------------------------------------------
#   dfSR<-data.frame(Year  = rep(seq(1980,2015,by=1),6), 
#                    Density = c(par[190:220],
#                                par_R1[187:216],NA,
#                                par_R2[184:212],NA,NA,
#                                par_R6[172:196],NA,NA,NA,NA,NA,NA,
#                                par_R9[163:184],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                par_R10[160:180],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                    Error = c(var[190:220],
#                              var_R1[187:216],NA,
#                              var_R2[184:212],NA,NA,
#                              var_R6[172:196],NA,NA,NA,NA,NA,NA,
#                              var_R9[163:184],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                              var_R10[160:180],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                    Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),
#                              rep("retro 6",36),rep("retro 9",36),
#                              rep("retro 10",36)))
#   
#   dfSR$Model <- factor(dfSR$Model, as.character(dfSR$Model))
#   
#   plot.retrospawn <- function(df) {
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                             'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                           'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       #coord_cartesian(xlim=c(1994,2015))+
#       #scale_y_continuous(limits=c(500,6000))+
#       labs(x="Year")+
#       labs(y="Spawning biomass (tons)")+
#       theme(legend.position = c(0.85,0.70))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#     
#   }
#   
# 
#   
#   #------------------------------------------------------------------------------
#   # SPAWNING BIOMASS GLOBAL RETROSPECTIVE Unified
#   #------------------------------------------------------------------------------
#   dfSRm <- data.frame(Year  = rep(seq(1980,2015,by=1),7), 
#                       Density = c(parM[189:219],
#                                   par_R1m[186:215],NA,
#                                   par_R2m[183:211],NA,NA,
#                                   par_R4m[177:203],NA,NA,NA,NA,
#                                   par_R5m[174:199],NA,NA,NA,NA,NA,
#                                   par_R9m[162:183],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                   par_R10m[159:179],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                       Error = c(varM[189:219],
#                                   var_R1m[186:215],NA,
#                                   var_R2m[183:211],NA,NA,
#                                   var_R4m[177:203],NA,NA,NA,NA,
#                                   var_R5m[174:199],NA,NA,NA,NA,NA,
#                                   var_R9m[162:183],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                   var_R10m[159:179],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                       Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),rep("retro 4",36),rep("retro 5",36),
#                                 rep("retro 9",36), rep("retro 10",36)))
#   
#   dfSRm$Model <- factor(dfSRm$Model, as.character(dfSRm$Model))
#   
#   
#   plot.retrospawnm <- function(df) {
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
#                                             'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
#                                           'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       labs(x="Year")+
#       labs(y="Spawning biomass (tons)")+
#       theme(legend.position = c(0.87,0.75))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#   }
#   
#   
#   
#   
#   
#   #------------------------------------------------------------------------------
#   # Recruitment GLOBAL RETROSPECTIVE
#   #------------------------------------------------------------------------------
#   dfRR<-data.frame(Year  = rep(seq(1980,2015,by=1),6), 
#                    Density = c(par[159:189],
#                                par_R1[157:186],NA,
#                                par_R2[155:183],NA,NA,
#                                par_R6[147:171],NA,NA,NA,NA,NA,NA,
#                                par_R9[141:162],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                par_R10[139:159],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                    Error= c(var[159:189],
#                             var_R1[157:186],NA,
#                             var_R2[155:183],NA,NA,
#                             var_R6[147:171],NA,NA,NA,NA,NA,NA,
#                             var_R9[141:162],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                             var_R10[139:159],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                    Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),
#                              rep("retro 6",36),rep("retro 9",36),
#                              rep("retro 10",36)))
#   
#   
#   dfRR$Model <- factor(dfRR$Model, as.character(dfRR$Model))
#   
#   
#   plot.retrorec <- function(df) {
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                             'retro 6' = 'green',  'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 
#                                           'retro 6' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       labs(x="Year")+
#       labs(y="Age 8 recruits (thousands)")+
#       theme(legend.position = c(0.8,0.7))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#     
#   }
#   
#   
#   #------------------------------------------------------------------------------
#   # RECRUITMENT GLOBAL RETROSPECTIVE Unified
#   #------------------------------------------------------------------------------
#   dfRRm <- data.frame(Year  = rep(seq(1980,2015,by=1),7), 
#                       Density = c(parM[158:188],
#                                   par_R1m[156:185],NA,
#                                   par_R2m[154:182],NA,NA,
#                                   par_R4m[150:176],NA,NA,NA,NA,
#                                   par_R5m[148:173],NA,NA,NA,NA,NA,
#                                   par_R9m[140:161],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                   par_R10m[138:158],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                       Error= c(varM[158:188],
#                                   var_R1m[156:185],NA,
#                                   var_R2m[154:182],NA,NA,
#                                   var_R4m[150:176],NA,NA,NA,NA,
#                                   var_R5m[148:173],NA,NA,NA,NA,NA,
#                                   var_R9m[140:161],NA,NA,NA,NA,NA,NA,NA,NA,NA,
#                                   var_R10m[138:158],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
#                       Model = c(rep("Full",36),rep("retro 1",36),rep("retro 2",36),rep("retro 4",36),rep("retro 5",36),
#                                 rep("retro 9",36), rep("retro 10",36)))
#   
#   dfRRm$Model <- factor(dfRRm$Model, as.character(dfRRm$Model))
#   
#   
#   plot.retrorecm <- function(df) {
#     df <- df %>% 
#       mutate(err= 1.96 * Error) %>%
#       mutate(lower = Density - err, upper = Density + err)
#     
#     
#     ggplot(df,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill=Model),alpha=0.2, show.legend=FALSE) + 
#       scale_colour_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange', 'retro 4' = 'gold',
#                                             'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       scale_fill_manual(name='', values=c('Full' = '#000000', 'retro 1' = 'red', 'retro 2' = 'orange','retro 4' = 'gold', 
#                                           'retro 5' = 'green', 'retro 9' = 'darkblue', 'retro 10' = 'violet'))+
#       labs(x="Year")+
#       labs(y="Age 8 recruits (thousands)")+
#       theme(legend.position = c(0.87,0.75))+
#       theme(legend.direction = "vertical")+
#       guides(fill=guide_legend(ncol=2))
#   }
#   
#   
  
#------------------------------------------------------------------------------
# SDNR - Global
#------------------------------------------------------------------------------
  # dfsdnr <- data.frame(Year  = rep(yearT,4), 
  #                      Result= sdnr$Result,
  #                      Iteration = sdnr$Iteration,
  #                      Metric = sdnr$Metric)
  # 
  # plot.sample <- function(df) {
  #   ggplot(df, aes(Year,Result,group=Iteration)) + geom_line(aes(Year,Result,colour=Iteration),size=1)+
  #     facet_wrap(~Metric, scales = "free", ncol = 1)+
  #     scale_colour_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                           'Third pass' = 'green'))+
  #     scale_fill_manual(name='', values=c('Data' = '#000000', 'First pass' = 'purple', 'Second pass' = 'blue', 
  #                                         'Third pass' = 'green'))+
  #     labs(x="Year")+
  #     labs(y="")+
  #     theme(legend.margin=unit(0.2, "cm"))+
  #     theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 12))+
  #     theme(legend.direction = "vertical")+
  #     guides(fill=guide_legend(ncol=2))
  # }
  # 
  

  
#------------------------------------------------------------------------------
# FUNCTION CALLS
#------------------------------------------------------------------------------
d1 <- plot.parmc(df = parMcmc) 
d2 <- plot.mr(df = dfmr) 
d3 <- plot.mrs(df = dfmrs) 
d4 <- plot.A(df = dfA) 
d5 <- plot.rec(df = dfR) 
d6 <- plot.spawners(df = dfS) 
d7 <- plot.sel(df = dfSel) 
d8 <- plot.Fmt(df = dfF) 
# d7 <- plot.age(dfs = dfs, dfm = dfm)
d9 <- plot.cpue(df = dfcp)
d10 <- plot.cpue2(df = dfcp2)
# d10 <- plot.proj(df = dfproj)
# d11 <- plot.retrodensity(df = dfDR)
# d12 <- plot.retrospawn(df = dfSR)
#d13 <- plot.retrorec(df = dfRR)
#d14 <- plot.sample(df = dfsdnr)
#d15 <- plot.sample(df = dfsdnrM)
# d16 <- plot.retrodensitym(df = dfDRm)
# d17 <- plot.retrospawnm(df = dfSRm)
# d18 <- plot.retrorecm(df = dfRRm)
#d19 <- plot.age(df = dfx)

#------------------------------------------------------------------------------
# EXPORT
#------------------------------------------------------------------------------
ggsave(d1,file='figures/theta.png', dpi=300, width=7, height=8, units ="in", bg="transparent") 
ggsave(d2,file='figures/mark_recapture.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d3,file='figures/mark_recapture_S.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d4,file='figures/abundance.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d5,file='figures/recruitment.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d6,file='figures/spawning.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d7,file='figures/selectivity.png', dpi=300, width=9, height=7, units ="in", bg="transparent") 
ggsave(d8,file='figures/Fmort.png', dpi=300, width=7, height=4., units ="in", bg="transparent") 
# ggsave(d7,file='figures/age_res.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d9,file='figures/cpue.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
ggsave(d10,file='figures/cpue2.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d10,file='figures/spbm_proj.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d11,file='figures/retro_den.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d12,file='figures/retro_spawn.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d13,file='figures/retro_rec.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
#ggsave(d14,file='../../figures/sample.png', dpi=300, width=7, height=4, units ="in", bg="transparent") 
#ggsave(d15,file='../../figures/sample2.png', dpi=300, width=7, height=4, units ="in", bg="transparent") 
# ggsave(d16,file='figures/retro_denm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d17,file='figures/retro_spawnm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
# ggsave(d18,file='figures/retro_recm.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 
#ggsave(d19,file='../../figures/age_res_uncorrected.png', dpi=300, width=7, height=5, units ="in", bg="transparent") 










# Ds <- read.delim("../model_output/dens_sdnr.ps",header=FALSE,sep="")
# Dg <- read.delim("../model_output/dens_g.ps",header=FALSE,sep="")
# Do <- read.delim("../model_output/dens_4.ps",header=FALSE,sep="")
# 
# 
# Ds<-Ds[c(10001:25000),]
# Dg<-Dg[c(430358:445357),]
# Do<-Do[c(5001:20000),]
# 
# for(i in 1:36)
# {
#   Ds[,i]<-sort(Ds[,i]) 
#   Dg[,i]<-sort(Dg[,i]) 
#   Do[,i]<-sort(Do[,i]) 
# }
# 
# MenD<-as.matrix(c(colMeans(Ds),colMeans(Dg),colMeans(Do)))
# MinD<-as.matrix(c(Ds[min,],Dg[min,],Do[min,]))
# MaxD<-as.matrix(c(Ds[max,],Dg[max,],Do[max,]))
# 
# dfD <- data.frame(Year  = rep(seq(1980,2015,by=1),3), 
#                   Density = MenD,
#                   MinD  = unlist(MinD),
#                   MaxD  = unlist(MaxD),
#                   Model = c(rep("S",36),rep("G",36),rep("O",36)),
#                   Rov   = rep(c(rep(NA,10),2820,NA,2103.5,NA,1980,rep(NA,3),2839,NA,2357,NA,1050,NA,1930,NA,NA,752,986,NA,1641),3),
#                   Var   = rep(c(rep(NA,10),549.5,NA,474.5,NA,380,rep(NA,3),417.5,NA,424,NA,126.3,NA,320,NA,NA,97,217,NA,288),3))
# 
# plot.density <- function(dfD) {
#   dfD <- dfD %>%
#     mutate(p.var = 1.96 * Var) %>%
#     mutate(p.lower = Rov - p.var, p.upper = Rov + p.var)
#   
#   
#   ggplot(dfD,aes(Year,Density,group=Model)) + geom_line(aes(Year,Density,colour=Model),size=1)+
#     geom_ribbon(aes(ymin = MinD, ymax = MaxD, fill=Model),alpha=0.2, show.legend=FALSE) + 
#     geom_point(aes(x=Year,y= Rov), colour="black",size=4)+
#     geom_errorbar(aes(x = Year,ymin = p.lower, ymax = p.upper),colour="black",lwd=0.8, width=0.5)+
#     scale_colour_manual(name='', values=c('S' = '#000000', 'G' = '#009E73', 'O' = '#e79f00'))+
#     scale_fill_manual(name='', values=c('S' = '#000000', 'G' = '#009E73', 'O' = '#e79f00'))+
#     labs(x="Year")+
#     labs(y="Density (ind. per square kilometer)")+
#     theme(legend.position = c(0.8,0.8))+
#     theme(legend.direction = "vertical")+
#     #theme(legend.position = c(0.8,0.8))+
#     theme(legend.key.size = unit(1., "cm"))+
#     #ggtitle("Total density (individuals per square kilometer)")+
#     theme(legend.key = element_rect(fill = "white"))
#   
#}

# dev_mcmc <-data.frame(Par  = c(unlist(mcmc_Gx[,11:41]),unlist(mcmc_2x[,11:41])),
#                       Time =   c(rep(c(rep(seq(1,36,by=1),each=15000)),2)),
#                       Model =  c(rep("Global",465000),rep("Uncorrected global",465000)))
# 
# plot.devmc <- function(df) {
#   ggplot(dev_mcmc,
#          aes(Par,group = interaction(Time,Model))) + geom_density(aes(Par, y=..scaled.., fill = Model),alpha=.2)+
#     facet_wrap(~Time,ncol=5)+
#     labs(y="")+
#     scale_fill_manual(name='', values=c('Global' = '#009E73',
#                                         'Uncorrected global' = '#e79f00'))
# }
