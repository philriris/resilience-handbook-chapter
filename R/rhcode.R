### The biogeography of population resilience in lowland South America
### Philip Riris
### Institute of Modelling Socio-Environmental Transitions
### Bournemouth University

## 0. Setup & data loading

require(rcarbon)
require(sp)
require(ggplot2)
require(rgdal)
require(sf)
require(ggsn)
require(rworldmap)
require(patchwork)
require(cowplot)
require(dplyr)
require(lme4)
require(effects)
require(cAIC4)
require(scales)

dataset <- read.csv("rhdata.csv")
dataset <- subset(dataset, dataset$calcurve == "shcal20") # Only terrestrial dates
dataset$SiteCode <- paste("S", as.numeric(as.factor(dataset$SiteName)), sep="")

calDates = calibrate(dataset$Age,
                     dataset$Error,
                     normalised = FALSE,
                     dateDetails = dataset$SiteCode,
                     calCurves = dataset$calcurve)

bins = binPrep(sites = dataset$SiteCode, ages = dataset$Age, 
               h = 50)

## 1. Data processing and display

# 1.1 SPDs

e.spd <- spd(x = calDates, timeRange = c(8000,500), spdnormalised = FALSE)
s.spd <- spd(x = calDates, timeRange = c(8000,500), runm = 100, spdnormalised = FALSE)

# 1.2 Base maps

base=getMap(resolution="low")
locations = unique(data.frame(SiteID = dataset$Age, 
                              Longitude = dataset$Lon, 
                              Latitude = dataset$Lat,
                              Where = dataset$L2Desc))

locations= locations[,-1]
coordinates(locations) <-c("Longitude", "Latitude")
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84")

base2 <- st_as_sf(base)
loc.sf <- st_as_sf(locations)
loc.bb <- st_as_sfc(st_bbox(loc.sf))


ecoregions <- readOGR(dsn=".", layer = "ecoregions") 
ecoregions <- st_as_sf(ecoregions)
ecoregions$LEVEL2[ecoregions$LEVEL2 == 17.1] <- 16.3
ecoregions$LEVEL2[ecoregions$LEVEL2 == 21.1] <- 20.3

cols <- data.frame(cols = c("#B5D12F", "#B4D030", "#7E8D54", 
          "#B4DDA1",
          "#F3CF3B", "#8DB465", "#69BD5D", "#9FBC4E", "#9DD357",
          "#AAA676", "#D9C5A0", "#879E1C"),
          LEVEL2 = c(16.1, 16.2, 16.3,
                     17.1,
                     20.1, 20.2, 20.3, 20.4, 20.5,
                     21.2, 21.3, 21.4))

ecoregions$col <- cols$cols[match(ecoregions$LEVEL2, cols$LEVEL2)]

ylims <- c(st_bbox(loc.sf)[2]-2, st_bbox(loc.sf)[4]+2)

box <- st_bbox(loc.sf)
box[4] <- box[4]+2
box[2] <- box[2]-2
box[3] <- box[3]+1.5

## 2. Radiocarbon analysis

# 2.1 Mark-permutation test

pt <- permTest(x = calDates, marks = as.factor(dataset$L2Desc), 
               timeRange = c(8000,500), nsim = 1000, bins = bins, runm = 50, 
               datenormalised = FALSE, spdnormalised = FALSE)

# 2.2 Resilience metrics

# Function for extracting resilience metrics from a permTest object
resmet <- function (x, focalm=1){
  name <- names(x$observed)[focalm]
  obs <- x$observed[[focalm]]
  envelope <- x$envelope[[focalm]]
  colnames(obs)<-c("calBP","PrDens")
  colnames(envelope)<-c("lo","hi")
  obs$Years <- obs$calBP
  busts <- which(obs$PrDens<envelope[,1])
  years <- obs$calBP[busts]
  ranges <- split(years, cumsum(c(1, diff(years) != -1)))
  #ranges <- ranges[-which(lapply(ranges, length)<cutoff)]
  
  mm <- data.frame(Tstart=numeric(), Tend=numeric(),
                   Duration = numeric(), Resistance = numeric(),
                   Resilience = numeric(), Lag = numeric(), LD = numeric(),
                   Recoveries = numeric(), Region = character())
  
  steps <- length(ranges)
  
  for (i in 1:steps){
    mm[i,1] = max(ranges[[i]]) #starts
    mm[i,2] = min(ranges[[i]]) #ends
    mm[i,3] = mm[i,1] - mm[i,2] #total time
    index1 = match(mm[i,1], obs$calBP)
    index2 = match(mm[i,2], obs$calBP)
    
    p1 = obs$PrDens[index1] # PrDens beginning
    p2 = obs$PrDens[index2] # PrDens end
    p3 = obs[index1:index2, 1][which.min(obs[index1:index2, 2])] # min point between p1:p2
    
    index3 = match(p3, obs$calBP) # index of middle
    p3 = obs$PrDens[index3] # PrDens middle
    
    d <- abs(p1 - p3) #resistance
    dd <- 2 * d
    p <- abs(p1) + d # abs() because negative growth rates
    dp <- dd/p
    mm[i,4] <- 1-dp
    
    e <- abs(p1 - p2) #resilience
    e <- d + e
    e <- dd/e
    mm[i,5] <- e - 1
    
    mm[i,6] = mm[i,1] - obs$calBP[index3]
    mm[i,7] = mm[i,6] / mm[i,3]
    mm[i,8] = steps
    mm[i,9] = name
  }
  
  return(mm)
}

metrics <- data.frame()

# This subsets out the Orinoco Llanos (no downturns)
#for (i in 1:9){
#  out = resmet(pt, focalm=i)
#  metrics = rbind(metrics, out)
#}
#
#for (i in 11){
#  out = resmet(pt, focalm=i)
#  metrics = rbind(metrics, out)
#}
#
#metrics <- subset(metrics, metrics$Resistance != 1) # Keeps only true downturns
#write.csv(x = metrics, file = "metrics_regular.csv")

metrics <- read.csv("metrics_regular.csv")
metrics$stdevents2 <- log10((metrics$Cumulative/metrics$Duration)*1000) # Cumulative number of events/1000 yrs accounting for duration

# 2.3 Plotting preparation
metsum <- metrics %>%  
  group_by(Region) %>% 
  tally() %>%  
  arrange(desc(n)) 

metrics$Rebound <- metrics$Duration - metrics$Lag #Rebound is the time to completion recovery

binbreaks <- c(-1,50,100,500,1000,2000)

metrics$rbins <- cut(metrics$Rebound, breaks=binbreaks)

metrics2 <- metrics %>% 
  group_by(rbins) %>% 
  summarise(rec.count = sum(Cumulative))

metsub <- metrics %>% #Subsetting out regions where obs = 1
  filter(Region != "Bahamas") %>% 
  filter(Region != "Atlantic Forests")

## 3. Statistical Modelling

# 3.1 Resistance

fits <- lmer(Resistance ~ (1|Region), data=metrics, REML=T)

cAIC4::stepcAIC(fits, fixEfCandidates=c("LD", "stdevents2", "Lag", "Duration"), 
                data=metrics, trace=TRUE, 
                direction="forward", returnResult=TRUE)

#Best model is as follows:
fit.resistance <- lmer(Resistance ~ (1|Region) + stdevents2, data=metrics, REML=T)

# 3.2 Resilience

fits <- lmer(Resilience ~ (1|Region), data=metrics, REML=T)

cAIC4::stepcAIC(fits, fixEfCandidates=c("LD", "stdevents2", "Lag", "Duration"), 
                data=metrics, trace=TRUE, 
                direction="forward", returnResult=TRUE)

#Best model is as follows:
fit.resilience <- lmer(Resilience ~ (1|Region) + LD, data=metrics, REML=T)

# 3.3 Effect sizes

effects_rt <- effect(term="stdevents2", mod=fit.resistance)

effects_rt <- as.data.frame(effects_rt)

effects_rs <- effect(term="LD", mod=fit.resilience)

effects_rs <- as.data.frame(effects_rs)

## 4. output

# Figure 1
m <- ggplot() + 
  geom_sf(data=base2) +
  xlim(st_bbox(loc.sf)[c(1, 3)]) +
  ylim(ylims) +
  geom_sf(data=ecoregions, aes(col = as.factor(ecoregions$LEVEL2), fill=as.factor(LEVEL2)))  +
  geom_sf(data=loc.sf, col="grey30", size=0.9) +
  geom_sf(data=base2, col="black", fill=NA) +
  theme_bw() +
  scale_fill_manual(values = c("#B5D12F", "#B4D030", "#7E8D54", #16
                               "#B4DDA1", #17
                               "#F3CF3B", "#8DB465", "#69BD5D", "#9FBC4E", "#9DD357", #20
                               "#AAA676", "#D9C5A0", "#879E1C")) + #21
  scale_color_manual(values = c("#B5D12F", "#B4D030", "#7E8D54", #16
                                "#B4DDA1", #17
                                "#F3CF3B", "#8DB465", "#69BD5D", "#9FBC4E", "#9DD357", #20
                                "#AAA676", "#D9C5A0", "#879E1C")) + #21
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank(),
        legend.position = "none") +
  north(symbol = 10,
        x.min = box[1], x.max = box[3], y.min = box[2], y.max = box[4])

inset <- ggplot(data=e.spd$grid, aes(x = calBP, y = PrDens)) +
  geom_area(fill="grey70") +
  geom_line(data = s.spd$grid, (aes(x = calBP, y = PrDens))) +
  xlab("cal years BP") +
  ylab("Probability density") +
  theme_bw() 


  cowplot::ggdraw() +
  cowplot::draw_plot(m) +
  cowplot::draw_plot(inset, x = 0.5, y = 0.7, width = 0.35, height = 0.25)


# Figure 2

covars1 <- ggplot(metsub) + 
  geom_boxplot(aes(y=reorder(as.factor(Region), -stdevents2), x=stdevents2), fill="grey50") +
  geom_jitter(aes(y=reorder(as.factor(Region), -stdevents2), x=stdevents2), shape="o") +
  theme_bw() +
  geom_vline(xintercept = median(metsub$stdevents2), linetype=2) +
  scale_y_discrete(labels=abbreviate) + ylab("Region") + 
  xlab("Standardised rate of disturbance/1000 yrs")

covars2 <- ggplot(metsub) + 
  geom_boxplot(aes(y=reorder(as.factor(Region), -LD), x=LD), fill="grey50") +
  geom_jitter(aes(y=reorder(as.factor(Region), -LD), x=LD), shape="o") +
  theme_bw() +
  geom_vline(xintercept = median(metsub$LD), linetype=2) +
  ylab("Region") + xlab("Normalised speed of downturn (Lag:Duration)") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

covars1+covars2

# Figure 3

#Effect size resilience
eff_rs <- ggplot() + 
  geom_point(data=metrics, aes(x=LD, y=Resilience)) + 
  geom_line(data=effects_rs, aes(x=LD, y=fit), color="purple", size=1) +
  geom_ribbon(data= effects_rs, aes(x=LD, ymin=lower, ymax=upper), alpha= 0.3, fill="purple") +
  theme_bw() + xlab("Lag:Duration") 

#Effect size resistance
eff_rt <- ggplot() + 
  geom_point(data=metrics, aes(x=stdevents2, y=Resistance)) + 
  geom_line(data=effects_rt, aes(x=stdevents2, y=fit), color="darkorange", size=1) +
  geom_ribbon(data= effects_rt, aes(x=stdevents2, ymin=lower, ymax=upper), alpha= 0.3, fill="darkorange") +
  theme_bw() + xlab("Standardised rate of downturn/1000 yrs") +
  coord_cartesian(ylim=c(0,1))

#Rebound time
rebound.cum <- ggplot(metrics2, aes(x=rbins, y=rec.count, group=rbins)) + geom_bar(stat="identity") +
  scale_x_discrete(labels = c("0-50", "50-100", "100-500", "500-1000", ">1000")) + 
  xlab("Rebound time (yrs)") + ylab("Cumulative # of recoveries") + theme_bw() +
  theme(axis.title.y = element_blank())

#Domesticates

edi <- read.csv("domesticates.csv") 

edi.cum <- ggplot(edi, aes(Median, weight=Sum)) + stat_ecdf(geom="step", col="#29AF7FFF", lwd=1.5) +
  scale_y_continuous(name = "Cumulative density", trans="reverse", labels = seq(1,0, -0.25)) + 
  theme_bw() + xlab("") +
  geom_vline(xintercept = 8200, linetype = "dashed", col="grey60") +
  geom_vline(xintercept = 4200, linetype = "dashed", col="grey60") +
  geom_rug()

edi.dens <- ggplot(edi, aes(Median, weight=Sum)) + 
  geom_density(kernel="epanechnikov", 
               bw="SJ-dpi", fill="#29AF7FFF",
               col="#1f968bff", lwd=1.5) + 
  theme_bw() + xlab("cal years BP") + ylab("Density") +
  theme(axis.text.y = element_text(angle = 90, size=7, hjust = .4)) + 
  scale_y_continuous(labels = scientific) +
  geom_vline(xintercept = 8200, linetype = "dashed", col="grey60") +
  geom_vline(xintercept = 4200, linetype = "dashed", col="grey60") +
  geom_rug()

(eff_rt + eff_rs) / (rebound.cum + (edi.cum/edi.dens)) + plot_annotation(tag_levels = 'A')

## 5. Supplementary output

# 5.1 Plots

par(mfrow=c(3,4))
for (i in 1:11){
  plot(pt, focalm=i, main = names(pt$observed.roc)[i])
}

# 5.2 Model diagnostics

#Resistance
sjPlot::tab_model(fit.resistance, file="resistance_fit.html")

plot(fit.resistance, resid(., scaled=TRUE) ~fitted(.), abline=0, pch=16,col=as.factor(metrics$Region))
plot(fit.resistance, as.factor(Region) ~ resid(., scaled=TRUE),abline=0,pch=16,xlab="Standardised residuals",ylab="Regions")
plot(fit.resistance, resid(., scaled=TRUE) ~ fitted(.)| Region, abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
car::qqPlot(resid(fit.resistance), col=as.factor(metrics$Region))

#Resilience
sjPlot::tab_model(fit.resilience,  file="resilience_fit.html")

plot(fit.resilience, resid(., scaled=TRUE) ~fitted(.), abline=0, pch=16,col=as.factor(metrics$Region))
plot(fit.resilience, as.factor(Region) ~ resid(., scaled=TRUE),abline=0,pch=16,xlab="Standardised residuals",ylab="Regions")
plot(fit.resilience, resid(., scaled=TRUE) ~ fitted(.)| Region, abline = 0,pch=16,xlab="Fitted values",ylab="Standardised residuals")
car::qqPlot(resid(fit.resilience), col=as.factor(metrics$Region))
