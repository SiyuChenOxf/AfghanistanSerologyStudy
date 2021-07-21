##Afghanistan exposure by region##
rm(list = ls())
set.seed(100)

library(rstan)
library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(dplyr)
library(tidyr)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggpubr)
library(egg)
library(RColorBrewer)
library(gridExtra)
library(bayesplot)
library(rstanarm)
setwd("D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes")

# #Extracting the whole posterior distribution of seroreversion rate
# SeroModelConstantIFR <- readRDS("SeroModelConstantIFR.rds")
# alpha <- rstan::extract(SeroModelConstantIFR)$beta
# alpha<-sample(alpha,size=10000)

##Common values for figure generation
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 14
font_size_title = 18
lwd = 1
pt_size = 0.3
right_margin=1

##Serolory test performance
kse<-1         #Sensitivity of serology test((IgM+ or IgG+; Total) https://www.fda.gov/media/138438/download)
ksp<-0.975     #Specificity of serology test((IgM- / IgG-; Total) https://www.fda.gov/media/138438/download)

##MCMC parameterss
niter <- 10000
chains<-4

##Serology positivity ratios by region derived from serology survey directly (unadjusted by the serology test performance)
seropositivity_Kabul<-0.53
seropositivity_East<-0.429
seropositivity_Central<-0.363
seropositivity_West<-0.341
seropositivity_NorthEast<-0.324
seropositivity_SouthEast<-0.322
seropositivity_North<-0.307
seropositivity_South<-0.258
seropositivity_CentralHighland<-0.211

#Serology survey sample sizes by region
seroN_Kabul<-1104                                                        
seroN_East<-1233                                                        
seroN_Central<-1056                                                      
seroN_West<-1176                                                         
seroN_NorthEast<-1265                                                    
seroN_SouthEast<-969                                                     
seroN_North<-1071                                                        
seroN_South<-738                                                         
seroN_CentralHighland<-902                                               

#This section is to use Baysian inference to estimate the seroprevalence by region adjusted by serology test performance
posi_nega_Kabul <- c(rep(1,round(seroN_Kabul*seropositivity_Kabul)),rep(0,seroN_Kabul-round(seroN_Kabul*seropositivity_Kabul)))                     
sero_data_Kabul <- list(N = length(posi_nega_Kabul), y = posi_nega_Kabul,kse=kse,ksp=ksp) 

posi_nega_East <- c(rep(1,round(seroN_East*seropositivity_East)),rep(0,seroN_Kabul-round(seroN_East*seropositivity_East)))                     
sero_data_East <- list(N = length(posi_nega_East), y = posi_nega_East,kse=kse,ksp=ksp) 

posi_nega_Central <- c(rep(1,round(seroN_Central*seropositivity_Central)),rep(0,seroN_Kabul-round(seroN_Central*seropositivity_Central)))                     
sero_data_Central <- list(N = length(posi_nega_Central), y = posi_nega_Central,kse=kse,ksp=ksp) 

posi_nega_West <- c(rep(1,round(seroN_West*seropositivity_West)),rep(0,seroN_Kabul-round(seroN_West*seropositivity_West)))                     
sero_data_West <- list(N = length(posi_nega_West), y = posi_nega_West,kse=kse,ksp=ksp) 

posi_nega_NorthEast <- c(rep(1,round(seroN_NorthEast*seropositivity_NorthEast)),rep(0,seroN_Kabul-round(seroN_NorthEast*seropositivity_NorthEast)))                     
sero_data_NorthEast <- list(N = length(posi_nega_NorthEast), y = posi_nega_NorthEast,kse=kse,ksp=ksp) 

posi_nega_SouthEast <- c(rep(1,round(seroN_SouthEast*seropositivity_SouthEast)),rep(0,seroN_Kabul-round(seroN_SouthEast*seropositivity_SouthEast)))                     
sero_data_SouthEast <- list(N = length(posi_nega_SouthEast), y = posi_nega_SouthEast,kse=kse,ksp=ksp) 

posi_nega_North <- c(rep(1,round(seroN_North*seropositivity_North)),rep(0,seroN_Kabul-round(seroN_North*seropositivity_North)))                     
sero_data_North <- list(N = length(posi_nega_North), y = posi_nega_North,kse=kse,ksp=ksp) 

posi_nega_South <- c(rep(1,round(seroN_South*seropositivity_South)),rep(0,seroN_Kabul-round(seroN_South*seropositivity_South)))                     
sero_data_South <- list(N = length(posi_nega_South), y = posi_nega_South,kse=kse,ksp=ksp) 

posi_nega_CentralHighland <- c(rep(1,round(seroN_CentralHighland*seropositivity_CentralHighland)),rep(0,seroN_Kabul-round(seroN_CentralHighland*seropositivity_CentralHighland)))                     
sero_data_CentralHighland <- list(N = length(posi_nega_CentralHighland), y = posi_nega_CentralHighland,kse=kse,ksp=ksp) 

# AdjSeroprevalence_Kabul <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#            data = sero_data_Kabul,
#            iter = niter,
#            chains = chains,
#            control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_East <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                          data = sero_data_East,
#                          iter = niter,
#                          chains = chains,
#                          control = list(adapt_delta = 0.99)) 
# 
# 
# AdjSeroprevalence_South <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                               data = sero_data_South,
#                               iter = niter,
#                               chains = chains,
#                               control = list(adapt_delta = 0.99))
# 
# AdjSeroprevalence_West <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                               data = sero_data_West,
#                               iter = niter,
#                               chains = chains,
#                               control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_North <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                                data = sero_data_North,
#                                iter = niter,
#                                chains = chains,
#                                control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_Central <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                          data = sero_data_Central,
#                          iter = niter,
#                          chains = chains,
#                          control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_NorthEast <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                          data = sero_data_NorthEast,
#                          iter = niter,
#                          chains = chains,
#                          control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_SouthEast <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                          data = sero_data_SouthEast,
#                          iter = niter,
#                          chains = chains,
#                          control = list(adapt_delta = 0.99)) 
# 
# AdjSeroprevalence_CentralHighland <-stan(file="D:/OxfDPhil/Research/Afganistan/SeroEpipaper/Codes/AdjSero.stan",
#                          data = sero_data_CentralHighland,
#                          iter = niter,
#                          chains = chains,
#                          control = list(adapt_delta = 0.99)) 
# saveRDS(AdjSeroprevalence_CentralHighland, file="AdjSeroprevalence_CentralHighland.rds")
# saveRDS(AdjSeroprevalence_SouthEast, file="AdjSeroprevalence_SouthEast.rds")
# saveRDS(AdjSeroprevalence_NorthEast, file="AdjSeroprevalence_NorthEast.rds")
# saveRDS(AdjSeroprevalence_Central, file="AdjSeroprevalence_Central.rds")
# saveRDS(AdjSeroprevalence_West, file="AdjSeroprevalence_West.rds")
# saveRDS(AdjSeroprevalence_East, file="AdjSeroprevalence_East.rds")
# saveRDS(AdjSeroprevalence_Kabul, file="AdjSeroprevalence_Kabul.rds")
# saveRDS(AdjSeroprevalence_North, file="AdjSeroprevalence_North.rds")
# saveRDS(AdjSeroprevalence_South, file="AdjSeroprevalence_South.rds")


AdjSeroprevalence_CentralHighland<-readRDS("AdjSeroprevalence_CentralHighland.rds")
AdjSeroprevalence_SouthEast<-readRDS("AdjSeroprevalence_SouthEast.rds")
AdjSeroprevalence_NorthEast<-readRDS("AdjSeroprevalence_NorthEast.rds")
AdjSeroprevalence_Central<-readRDS("AdjSeroprevalence_Central.rds")
AdjSeroprevalence_West<-readRDS("AdjSeroprevalence_West.rds")
AdjSeroprevalence_East<-readRDS("AdjSeroprevalence_East.rds")
AdjSeroprevalence_North<-readRDS("AdjSeroprevalence_North.rds")
AdjSeroprevalence_Kabul<-readRDS("AdjSeroprevalence_Kabul.rds")
AdjSeroprevalence_South<-readRDS("AdjSeroprevalence_South.rds")

posx_CentralHighland <- rstan::extract(AdjSeroprevalence_CentralHighland)$p
posx_North <- rstan::extract(AdjSeroprevalence_North)$p
posx_South <- rstan::extract(AdjSeroprevalence_South)$p
posx_SouthEast <- rstan::extract(AdjSeroprevalence_SouthEast)$p
posx_NorthEast <- rstan::extract(AdjSeroprevalence_NorthEast)$p
posx_Kabul <- rstan::extract(AdjSeroprevalence_Kabul)$p
posx_Central <- rstan::extract(AdjSeroprevalence_Central)$p
posx_West <- rstan::extract(AdjSeroprevalence_West)$p
posx_East <- rstan::extract(AdjSeroprevalence_East)$p

# posx_Kabul<-sample(posx_Kabul,size = 10000)
# posx_North<-sample(posx_North,size = 10000)
# posx_South<-sample(posx_South,size = 10000)
# posx_SouthEast<-sample(posx_SouthEast,size = 10000)
# posx_NorthEast<-sample(posx_NorthEast,size = 10000)
# posx_Central<-sample(posx_Central,size = 10000)
# posx_West<-sample(posx_West,size = 10000)
# posx_East<-sample(posx_East,size = 10000)
# posx_CentralHighland<-sample(posx_East,size = 10000)

posx=as.data.frame(cbind(posx_Kabul,posx_East,posx_NorthEast,posx_West,posx_Central,posx_North,posx_SouthEast,posx_South,posx_CentralHighland))
data<-data.frame(label=factor(c("Kabul","East","Northeast","West","Central","North","Southeast","South","Central highland"), level = c("Kabul","East","Northeast","West","Central","North","Southeast","South","Central highland")),val=apply(posx, 2, median),upp=apply(posx, 2, function(x) quantile(x, probs = 0.95)),low=apply(posx, 2, function(x) quantile(x, probs = 0.025)))

##Saving seroprevalence by region 
p<-ggplot(data, aes(x=label, y=100*val, ymin = 100*low, ymax = 100*upp)) +
  geom_pointrange()+ggtitle("Seroprevalence on 21/7/2020 adjusted by serology test")+
  xlab(" ")+
  ylab("Percentage (%)")+
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  scale_y_continuous(breaks = 100*c(0,0.10,0.20,0.30,0.40,0.50,0.6), limit = 100*c(0, 0.6))+
  theme_minimal() +
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = 16,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    # legend.position = "None",
    legend.title = element_blank(),
    axis.text.y = element_text(size=14),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

tiff(file="AdjSeroprevalence.tiff",
     width=27, height=15, units="cm", res=300)
ggarrange(p)
dev.off()

##This section is to calculate the total exposure in the population by region
AFG_data<-read.csv("RegionalDeath.csv",header = TRUE)        #Load regional daily death data
serosurvey_date<-"2020-07-21"                                #Define the serology survey date

delta_epsilon<-21                                            #fixed time lag between exposure and death(seroconversion)
kd<-1/193                                                    #Assumed seroreversion rate

daily_death_Kabul<-AFG_data$Kabul
daily_death_Central<-AFG_data$Central
daily_death_East<-AFG_data$East
daily_death_NorthEast<-AFG_data$NorthEast 
daily_death_CentralHighland<-AFG_data$CentralHighland 
daily_death_Southeast<-AFG_data$Southeast 
daily_death_North<-AFG_data$North 
daily_death_South<-AFG_data$South
daily_death_West<-AFG_data$West

n<-length(daily_death_Kabul)                                  
t<-seq(1,n,by=1)
n2 <- as.integer(as.Date(serosurvey_date)-as.Date("2020-01-01")+1)                                                  
t2<-seq(1,n2,by=1)
exposure_Kabul<-seroprevalence_Kabul<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_Kabul[i,]<-posx_Kabul[i]*exp(kd*n2)*cumsum(daily_death_Kabul)/sum(exp(kd*t2)*daily_death_Kabul[1:n2])
  seroprevalence_Kabul[i,]<-posx_Kabul[i]/exp(kd*(t-n2))*cumsum(daily_death_Kabul*exp(kd*t))/sum(exp(kd*t2)*daily_death_Kabul[1:n2]) 
}


exposure_Kabul<-exposure_Kabul[,(delta_epsilon+1):n]


data1_Kabul = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                      t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                      median = c(100*apply(exposure_Kabul, 2, function(x) quantile(x, probs = 0.5)), 
                                 100*apply(seroprevalence_Kabul, 2, function(x) quantile(x, probs = 0.5))), 
                      lower1 = c(100*apply(exposure_Kabul, 2, function(x) quantile(x, probs = 0.025)), 
                                 100*apply(seroprevalence_Kabul, 2, function(x) quantile(x, probs = 0.025))), 
                      upper1 = c(100*apply(exposure_Kabul, 2, function(x) quantile(x, probs = 0.975)),
                                 100*apply(seroprevalence_Kabul, 2, function(x) quantile(x, probs = 0.975))))
data2_Kabul = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_Kabul), upper= 100*quantile(posx_Kabul, probs = 0.975), lower = 100*quantile(posx_Kabul, probs = 0.025))
ymax_exposure<-80#round(max(data1_Kabul$upper1*10))/10

p1Kabul<-ggplot(data1_Kabul, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("Kabul")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_Kabul, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1Kabul <- p1Kabul +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%)") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#South
exposure_South<-seroprevalence_South<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_South[i,]<-posx_South[i]*exp(kd*n2)*cumsum(daily_death_South)/sum(exp(kd*t2)*daily_death_South[1:n2])
  seroprevalence_South[i,]<-posx_South[i]/exp(kd*(t-n2))*cumsum(daily_death_South*exp(kd*t))/sum(exp(kd*t2)*daily_death_South[1:n2]) 
}
exposure_South<-exposure_South[,(delta_epsilon+1):n]


data1_South = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                         t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                         median = c(100*apply(exposure_South, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(seroprevalence_South, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(exposure_South, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(seroprevalence_South, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(exposure_South, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(seroprevalence_South, 2, function(x) quantile(x, probs = 0.975))))
data2_South = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_South), upper= 100*quantile(posx_South, probs = 0.975), lower = 100*quantile(posx_South, probs = 0.025))

p1South<-ggplot(data1_South, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("South")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_South, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1South <- p1South +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#North
exposure_North<-seroprevalence_North<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_North[i,]<-posx_North[i]*exp(kd*n2)*cumsum(daily_death_North)/sum(exp(kd*t2)*daily_death_North[1:n2])
  seroprevalence_North[i,]<-posx_North[i]/exp(kd*(t-n2))*cumsum(daily_death_North*exp(kd*t))/sum(exp(kd*t2)*daily_death_North[1:n2]) 
}
exposure_North<-exposure_North[,(delta_epsilon+1):n]


data1_North = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                         t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                         median = c(100*apply(exposure_North, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(seroprevalence_North, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(exposure_North, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(seroprevalence_North, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(exposure_North, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(seroprevalence_North, 2, function(x) quantile(x, probs = 0.975))))
data2_North = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_North), upper= 100*quantile(posx_North, probs = 0.975), lower = 100*quantile(posx_North, probs = 0.025))

p1North<-ggplot(data1_North, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("North")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_North, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1North <- p1North +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.y= element_blank(),
    axis.text.x= element_text(size=font_size),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#West
exposure_West<-seroprevalence_West<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_West[i,]<-posx_West[i]*exp(kd*n2)*cumsum(daily_death_West)/sum(exp(kd*t2)*daily_death_West[1:n2])
  seroprevalence_West[i,]<-posx_West[i]/exp(kd*(t-n2))*cumsum(daily_death_West*exp(kd*t))/sum(exp(kd*t2)*daily_death_West[1:n2]) 
}
exposure_West<-exposure_West[,(delta_epsilon+1):n]

data1_West = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                        t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                        median = c(100*apply(exposure_West, 2, function(x) quantile(x, probs = 0.5)), 
                                   100*apply(seroprevalence_West, 2, function(x) quantile(x, probs = 0.5))), 
                        lower1 = c(100*apply(exposure_West, 2, function(x) quantile(x, probs = 0.025)), 
                                   100*apply(seroprevalence_West, 2, function(x) quantile(x, probs = 0.025))), 
                        upper1 = c(100*apply(exposure_West, 2, function(x) quantile(x, probs = 0.975)),
                                   100*apply(seroprevalence_West, 2, function(x) quantile(x, probs = 0.975))))
data2_West = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_West), upper= 100*quantile(posx_West, probs = 0.975), lower = 100*quantile(posx_West, probs = 0.025))

p1West<-ggplot(data1_West, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_West, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1West <- p1West +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%)") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.title.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#NorthEast
exposure_NorthEast<-seroprevalence_NorthEast<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_NorthEast[i,]<-posx_NorthEast[i]*exp(kd*n2)*cumsum(daily_death_NorthEast)/sum(exp(kd*t2)*daily_death_NorthEast[1:n2])
  seroprevalence_NorthEast[i,]<-posx_NorthEast[i]/exp(kd*(t-n2))*cumsum(daily_death_NorthEast*exp(kd*t))/sum(exp(kd*t2)*daily_death_NorthEast[1:n2]) 
}
exposure_NorthEast<-exposure_NorthEast[,(delta_epsilon+1):n]


data1_NorthEast = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                             t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                             median = c(100*apply(exposure_NorthEast, 2, function(x) quantile(x, probs = 0.5)), 
                                        100*apply(seroprevalence_NorthEast, 2, function(x) quantile(x, probs = 0.5))), 
                             lower1 = c(100*apply(exposure_NorthEast, 2, function(x) quantile(x, probs = 0.025)), 
                                        100*apply(seroprevalence_NorthEast, 2, function(x) quantile(x, probs = 0.025))), 
                             upper1 = c(100*apply(exposure_NorthEast, 2, function(x) quantile(x, probs = 0.975)),
                                        100*apply(seroprevalence_NorthEast, 2, function(x) quantile(x, probs = 0.975))))
data2_NorthEast = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_NorthEast), upper= 100*quantile(posx_NorthEast, probs = 0.975), lower = 100*quantile(posx_NorthEast, probs = 0.025))

p1NorthEast<-ggplot(data1_NorthEast, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("Northeast")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_NorthEast, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1NorthEast <- p1NorthEast +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("Percentage (%)") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.title.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#Southeast
exposure_Southeast<-seroprevalence_Southeast<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_Southeast[i,]<-posx_SouthEast[i]*exp(kd*n2)*cumsum(daily_death_Southeast)/sum(exp(kd*t2)*daily_death_Southeast[1:n2])
  seroprevalence_Southeast[i,]<-posx_SouthEast[i]/exp(kd*(t-n2))*cumsum(daily_death_Southeast*exp(kd*t))/sum(exp(kd*t2)*daily_death_Southeast[1:n2]) 
}
exposure_Southeast<-exposure_Southeast[,(delta_epsilon+1):n]


data1_Southeast = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                             t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                             median = c(100*apply(exposure_Southeast, 2, function(x) quantile(x, probs = 0.5)), 
                                        100*apply(seroprevalence_Southeast, 2, function(x) quantile(x, probs = 0.5))), 
                             lower1 = c(100*apply(exposure_Southeast, 2, function(x) quantile(x, probs = 0.025)), 
                                        100*apply(seroprevalence_Southeast, 2, function(x) quantile(x, probs = 0.025))), 
                             upper1 = c(100*apply(exposure_Southeast, 2, function(x) quantile(x, probs = 0.975)),
                                        100*apply(seroprevalence_Southeast, 2, function(x) quantile(x, probs = 0.975))))
data2_Southeast = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_SouthEast), upper= 100*quantile(posx_SouthEast, probs = 0.975), lower = 100*quantile(posx_SouthEast, probs = 0.025))

p1Southeast<-ggplot(data1_Southeast, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("Southeast")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_Southeast, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1Southeast <- p1Southeast +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#Central
exposure_Central<-seroprevalence_Central<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_Central[i,]<-posx_Central[i]*exp(kd*n2)*cumsum(daily_death_Central)/sum(exp(kd*t2)*daily_death_Central[1:n2])
  seroprevalence_Central[i,]<-posx_Central[i]/exp(kd*(t-n2))*cumsum(daily_death_Central*exp(kd*t))/sum(exp(kd*t2)*daily_death_Central[1:n2]) 
}
exposure_Central<-exposure_Central[,(delta_epsilon+1):n]


data1_Central = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                           t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                           median = c(100*apply(exposure_Central, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(seroprevalence_Central, 2, function(x) quantile(x, probs = 0.5))), 
                           lower1 = c(100*apply(exposure_Central, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(seroprevalence_Central, 2, function(x) quantile(x, probs = 0.025))), 
                           upper1 = c(100*apply(exposure_Central, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(seroprevalence_Central, 2, function(x) quantile(x, probs = 0.975))))
data2_Central = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_Central), upper= 100*quantile(posx_Central, probs = 0.975), lower = 100*quantile(posx_Central, probs = 0.025))

p1Central<-ggplot(data1_Central, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("Central")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_Central, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1Central <- p1Central +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%)") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.x = element_text(size = font_size),
    axis.text.y = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#Central Highland
exposure_CentralHighland<-seroprevalence_CentralHighland<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_CentralHighland[i,]<-posx_CentralHighland[i]*exp(kd*n2)*cumsum(daily_death_CentralHighland)/sum(exp(kd*t2)*daily_death_CentralHighland[1:n2])
  seroprevalence_CentralHighland[i,]<-posx_CentralHighland[i]/exp(kd*(t-n2))*cumsum(daily_death_CentralHighland*exp(kd*t))/sum(exp(kd*t2)*daily_death_CentralHighland[1:n2]) 
}
exposure_CentralHighland<-exposure_CentralHighland[,(delta_epsilon+1):n]


data1_CentralHighland = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                                   t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                                   median = c(100*apply(exposure_CentralHighland, 2, function(x) quantile(x, probs = 0.5)), 
                                              100*apply(seroprevalence_CentralHighland, 2, function(x) quantile(x, probs = 0.5))), 
                                   lower1 = c(100*apply(exposure_CentralHighland, 2, function(x) quantile(x, probs = 0.025)), 
                                              100*apply(seroprevalence_CentralHighland, 2, function(x) quantile(x, probs = 0.025))), 
                                   upper1 = c(100*apply(exposure_CentralHighland, 2, function(x) quantile(x, probs = 0.975)),
                                              100*apply(seroprevalence_CentralHighland, 2, function(x) quantile(x, probs = 0.975))))
data2_CentralHighland = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_CentralHighland), upper= 100*quantile(posx_CentralHighland, probs = 0.975), lower = 100*quantile(posx_CentralHighland, probs = 0.025))

p1CentralHighland<-ggplot(data1_CentralHighland, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("Central Highland")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_CentralHighland, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1CentralHighland <- p1CentralHighland +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

exposure_East<-seroprevalence_East<-matrix(0,nrow = dim(posx)[1],ncol = n)

for (i in 1:dim(posx)[1]) {
  exposure_East[i,]<-posx_East[i]*exp(kd*n2)*cumsum(daily_death_East)/sum(exp(kd*t2)*daily_death_East[1:n2])
  seroprevalence_East[i,]<-posx_East[i]/exp(kd*(t-n2))*cumsum(daily_death_East*exp(kd*t))/sum(exp(kd*t2)*daily_death_East[1:n2]) 
}
exposure_East<-exposure_East[,(delta_epsilon+1):n]


data1_East = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)), 
                        t=c(as.Date(AFG_data$Date)[1:(n-delta_epsilon)],as.Date(AFG_data$Date)[1:n]), 
                        median = c(100*apply(exposure_East, 2, function(x) quantile(x, probs = 0.5)), 
                                   100*apply(seroprevalence_East, 2, function(x) quantile(x, probs = 0.5))), 
                        lower1 = c(100*apply(exposure_East, 2, function(x) quantile(x, probs = 0.025)), 
                                   100*apply(seroprevalence_East, 2, function(x) quantile(x, probs = 0.025))), 
                        upper1 = c(100*apply(exposure_East, 2, function(x) quantile(x, probs = 0.975)),
                                   100*apply(seroprevalence_East, 2, function(x) quantile(x, probs = 0.975))))
data2_East = data.frame( t=as.Date(serosurvey_date), value=100*median(posx_East), upper= 100*quantile(posx_East, probs = 0.975), lower = 100*quantile(posx_East, probs = 0.025))

p1East<-ggplot(data1_East, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = 100*c(0,0.2,0.4,0.6,0.8), limit = c(0,ymax_exposure ))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01")), labels=c("Jan","Mar", "May", "Jul","Sep"), limit = as.Date(c("2020-01-01","2020-08-04")))+
  geom_pointrange(data=data2_East, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1East <- p1East +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size_title),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 2, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#Exposure map
course.url <- "https://kinglab.eeb.lsa.umich.edu/480"
destfile <- file.path(tempdir(),'afghan_province.zip')
download.file(paste0(course.url,'/data/afghan_province.zip'),
              destfile=destfile,mode='wb')
unzip(destfile,exdir=tempdir())
options(stringsAsFactors=FALSE)
library(rgdal)
afg <- readOGR(dsn=tempdir(),layer='afghan_province')
AFG<-sf::st_as_sf(afg)

#Northeast
AFG$expos[1]<-median(exposure_NorthEast[,n-delta_epsilon])
AFG$expos[3]<-median(exposure_NorthEast[,n-delta_epsilon])
AFG$expos[19]<-median(exposure_NorthEast[,n-delta_epsilon])
AFG$expos[32]<-median(exposure_NorthEast[,n-delta_epsilon])

#West
AFG$expos[12]<-median(exposure_West[,n-delta_epsilon])
AFG$expos[2]<-median(exposure_West[,n-delta_epsilon])
AFG$expos[10]<-median(exposure_West[,n-delta_epsilon])
AFG$expos[7]<-median(exposure_West[,n-delta_epsilon])

#South
AFG$expos[11]<-median(exposure_South[,n-delta_epsilon])
AFG$expos[15]<-median(exposure_South[,n-delta_epsilon])
AFG$expos[24]<-median(exposure_South[,n-delta_epsilon])
AFG$expos[33]<-median(exposure_South[,n-delta_epsilon])
AFG$expos[34]<-median(exposure_South[,n-delta_epsilon])

#East
AFG$expos[18]<-median(exposure_East[,n-delta_epsilon])
AFG$expos[20]<-median(exposure_East[,n-delta_epsilon])
AFG$expos[23]<-median(exposure_East[,n-delta_epsilon])
AFG$expos[25]<-median(exposure_East[,n-delta_epsilon])

#North
AFG$expos[4]<-median(exposure_North[,n-delta_epsilon])
AFG$expos[8]<-median(exposure_North[,n-delta_epsilon])
AFG$expos[13]<-median(exposure_North[,n-delta_epsilon])
AFG$expos[30]<-median(exposure_North[,n-delta_epsilon])
AFG$expos[31]<-median(exposure_North[,n-delta_epsilon])

#Cenral
AFG$expos[16]<-median(exposure_Central[,n-delta_epsilon])
AFG$expos[21]<-median(exposure_Central[,n-delta_epsilon])
AFG$expos[22]<-median(exposure_Central[,n-delta_epsilon])
AFG$expos[28]<-median(exposure_Central[,n-delta_epsilon])
AFG$expos[29]<-median(exposure_Central[,n-delta_epsilon])

#Southeast
AFG$expos[9]<-median(exposure_Southeast[,n-delta_epsilon])
AFG$expos[26]<-median(exposure_Southeast[,n-delta_epsilon])
AFG$expos[17]<-median(exposure_Southeast[,n-delta_epsilon])
AFG$expos[27]<-median(exposure_Southeast[,n-delta_epsilon])

#CH
AFG$expos[5]<-median(exposure_CentralHighland[,n-delta_epsilon])
AFG$expos[6]<-median(exposure_CentralHighland[,n-delta_epsilon])

AFG$expos[14]<-median(exposure_Kabul[,n-delta_epsilon])

a<-data.frame(province=afg$prov,exps=AFG$expos)

a %>%
  left_join(as.data.frame(afg),by=c("province"="prov")) -> a


AFG_map<-ggplot(data=a,aes(fill=exps))+
  geom_map(aes(map_id=id),map=fortify(afg))+
  coord_map()+
  expand_limits(x=c(60,75),y=c(29,38.5))+
  scale_fill_distiller(palette = "Spectral",labels = scales::percent_format(accuracy = 5))+
  theme_void()+ 
  theme(
    legend.position = c(0.85,0.2),
    legend.text=element_text(size=font_size), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = font_size_title),
    legend.title = element_blank(),
    legend.key.width = unit(1, "cm"), legend.key.height = unit(1, "cm"))+
  ggtitle("Predicted exposure on 21/7/2020")
  

##Saving figures##
lay <- rbind(c(1,1,2,3,4),
             c(1,1,5,6,7),
             c(1,1,8,9,10))

tiff(file="ExposureSeroprevalence.tiff",
     width=34, height=20, units="cm", res=300)
grid.arrange(AFG_map,styled1Kabul,styled1North,styled1West,styled1NorthEast,styled1South,styled1Southeast,styled1Central,styled1CentralHighland,styled1East, layout_matrix = lay, widths=c(0.75,1.5,1.5,1.5,1.5))
dev.off()
