# Ben Sullender
# sullendb@uw.edu
# Oct 8, 2024
# Apex Predators x Snow
# R Script Part 3: Cougar Kills

library(tidyverse)
library(mgcv)

cougKills <- read.csv("/Users/bensullender/Documents/UW/Research/ChapterTwo/Final_Code_for_Github/cougarKills.csv")

# ID_wyr is concatenated cougar ID and water year (e.g. Oct 2023 thru Apr 2024 = water year2024)
# case 1 == kill/cluster site, 0 == randomly generated
# snod = snow depth (m)
# dens = snow density (kg/m3)
# Shr5 = percent shrubs (0-1), within 5 cell radius of focal location
# deer = deer index, 0 is lowest and 1 is highest (see Supplement)
# tri = terrain ruggedness index
# cc = canopy cover (0-10000, with 10000 = 100% canopy cover)


# Step 1: Use AIC to determine most parsimonious combination of snow variables
killGAMfull <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + te(snod,dens,k=5),
                   method = 'REML',family=binomial(),data=cougKills)
killGAMsnod <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(snod,k=5),
                   method = 'REML',family=binomial(),data=cougKills)
killGAMdens <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(dens,k=5),
                   method = 'REML',family=binomial(),data=cougKills)
killGAMti <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + ti(snod,dens,k=5),
                 method = 'REML',family=binomial(),data=cougKills)
killGAMnosnow <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5),
                     method = 'REML',family=binomial(),data=cougKills)
killGAMdenssnod <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(dens,k=5) + s(snod,k=5),
                       method = 'REML',family=binomial(),data=cougKills)

# Compare the AICs of each model
AIC(killGAMfull)
AIC(killGAMsnod)
AIC(killGAMdens)
AIC(killGAMti)
AIC(killGAMnosnow)

# Use the full tensor product in future steps
killGAM <- gam(case ~ s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + te(snod,dens,k=5),
               method = 'REML',family=binomial(),data=cougKills)


# Step 2: Generate relative selection strength plots (RSS)

# First, create prediction data frame
killPred <- cougKills
# create new version of kill predictions (="kp"), with all the same other covariate values, except fixed density at low (200 kg/m3) and high densities (350 kg/m3)
kpLite <- killPred %>% mutate(dens = 200)
kpDense <- killPred %>% mutate(dens=350)  
# same, but with fixed depth at shallow (5cm) and deep (30cm)
kpShal <- killPred %>% mutate(snod = 0.05)
kpDeep <- killPred %>% mutate(snod = 0.3)



# custom function to create relative selection strength.
#     Takes input GAM ("gammod"), input predictions data frame ("df"), focal variable ("variable"),
#     and value to compare predictions to ("valu")
calcRSS <- function(gammod,df,variable,valu = min(variable)){
  avgDF <- df
  # remove first variable (=y variable)
  allVars <- names(gammod$model)[-1]
  # remove focal variable, since we'll set that to a specific value
  allVars <- allVars[allVars!=variable]
  # for all other variables, set them at average value
  #     because we'll feed in new dataframes with fixed depths (kpShal and kpDeep) and densities (kpLite and kpDense),
  #     the average of those fixed values is still the fixed value.
  for (i in 1:length(allVars)){
    # !! and := allow tidyverse to interject a variable as df column name
    avgDF <- avgDF %>% mutate(!!allVars[i]:=mean(get(allVars[i]),na.rm=T))
  }
  avgDF$mod <- predict.gam(gammod,newdata=avgDF,type="response")
  avgDF$se <- predict.gam(gammod,newdata=avgDF,type="response",se.fit=T)$se.fit
  minDF <- avgDF %>% mutate(!!variable:=valu)
  avgDF$baseline <- predict.gam(gammod,newdata=minDF,type="response")
  avgDF$rss <- avgDF$mod/avgDF$baseline
  avgDF$rssMin <- c(avgDF$mod-avgDF$se)/avgDF$baseline
  avgDF$rssMax <- c(avgDF$mod+avgDF$se)/avgDF$baseline
  return(avgDF)
}

# Run RSS for each of the four fixed depths/densities, relative to the minimum observed depth (1cm) or density (56kg/m3)
snodLite <- calcRSS(killGAM,kpLite,"snod",0.01)
snodDense <- calcRSS(killGAM,kpDense,"snod",0.01)
densShal <- calcRSS(killGAM,kpShal,"dens",56)
densDeep <- calcRSS(killGAM,kpDeep,"dens",56)

# Now plot these for snow depth
ggplot()+
  geom_line(aes(x=snod,y=rss,color="Dense (350kg/m^3)"),size=1.5,alpha=0.7,data=snodDense)+
  geom_ribbon(aes(x=snod,ymin=rssMin,ymax=rssMax,fill="Dense (350kg/m^3)"),alpha=0.5,data=snodDense)+
  geom_line(aes(x=snod,y=rss,color="Light (200kg/m^3)"),size=1.5,alpha=0.7,data=snodLite[snodLite$snod<1,])+
  geom_ribbon(aes(x=snod,ymin=rssMin,ymax=rssMax,fill="Light (200kg/m^3)"),alpha=0.5,data=snodLite[snodLite$snod<1,])+
  geom_hline(yintercept=1,linetype=2,size=1.5)+
  coord_cartesian(ylim=c(0,3),expand=F)+
  scale_color_manual("Snow Density",values=c("Dense (350kg/m^3)" ="#FC7200","Light (200kg/m^3)"="black"))+
  scale_fill_manual("Snow Density",values=c("Dense (350kg/m^3)" ="#FC7200","Light (200kg/m^3)"="black"))+
  theme_bw(base_size = 24)+
  labs(x = 'Snow depth (m)', y = 'Kill site relative probability')

# And for snow density
ggplot()+
  geom_line(aes(x=dens,y=rss,color="Deep (30cm)"),size=1.5,alpha=0.7,data=densDeep)+
  geom_ribbon(aes(x=dens,ymin=rssMin,ymax=rssMax,fill="Deep (30cm)"),alpha=0.5,data=densDeep)+
  geom_line(aes(x=dens,y=rss,color="Shallow (5cm)"),size=1.5,alpha=0.7,data=densShal)+
  geom_ribbon(aes(x=dens,ymin=rssMin,ymax=rssMax,fill="Shallow (5cm)"),alpha=0.5,data=densShal)+
  coord_cartesian(ylim=c(0,3),expand=F)+
  geom_hline(yintercept=1,linetype=2,size=1.5)+
  scale_color_manual("Snow Depth",values=c("Deep (30cm)" ="#FC7200","Shallow (5cm)"="black"))+
  scale_fill_manual("Snow Depth",values=c("Deep (30cm)" ="#FC7200","Shallow (5cm)"="black"))+
  theme_bw(base_size = 24)+
  labs(x = bquote('Snow density'~(kg/m^3)), y = 'Kill site relative probability')


