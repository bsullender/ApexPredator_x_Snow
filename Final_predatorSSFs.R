# Ben Sullender
# sullendb@uw.edu
# Oct 8, 2024
# Apex Predators x Snow
# R Script Part 2: Predation SSF

#
# Part 0: Set-up
#

# Set up libraries
library(tidyverse)
library(mgcv)
library(gratia)

# Read in data
coug <- readRDS("/Users/bensullender/Documents/UW/Research/ChapterTwo/Final_Code_for_Github/cougarSteps.rds")
wolf <- readRDS("/Users/bensullender/Documents/UW/Research/ChapterTwo/Final_Code_for_Github/wolfSteps.rds")


#
# Part 1: Select most parsimonious of colinear variables for both predators
#


# If you run this code block, each model will take 20-45 minutes to run, for a total of 6-15 hours processing time. 
# Instead, I recommend skipping ahead to lines 206 and 208, which read in the already-run R model objects.

#     Rename so that columns match across species
wComb <- wolf %>% rename(forest = For19,shr = Shr19,open = Open19) %>% mutate(spp = "Wolf")
cComb <- coug %>% rename(forest = For5,shr = Shr5,open = Open5) %>% mutate(spp = "Coug")
wc <- rbind(wComb,cComb)
wc$spp <- as.factor(wc$spp)

# Check for combinations of variables that are too highly correlated to use in same models
wcC <- wc[,c(9:17,22)]
cor(wcC)

# First, check whether slope or tri is a better fit; each model takes ~45 mins to run
#     For this and all other models in Part 1, cbind(dummy, individ_step_f) conditions each used (case = 1) location against other steps that were not taken.
#     s() indicates a smooth term, which will be split into k-1 basis functions. We have chosen k=5 based on a sensitivity analysis (not shown here) indicating that
#     k=5 is a reasonable compromise of precision and computational cost.
#     The by=spp is a interaction of the predictor variable (slope in this first model) and species, which is either "Wolf" or "Coug".
system.time(sloMod <- gam(cbind(dummy, indiv_step_f) ~  s(slo,by=spp,k=5),
                          weights = case, # 1 0 reponse
                          method = 'REML', # recommended
                          select = TRUE, # penalty to favor parsimonious model
                          family=cox.ph(),
                          data=wc))

triMod <- gam(cbind(dummy, indiv_step_f) ~  s(tri,by=spp,k=5),
                          weights = case, # 1 0 reponse
                          method = 'REML', # recommended
                          select = TRUE, # penalty to favor parsimonious model
                          family=cox.ph(),
                          data=wc)
AIC(sloMod)
AIC(triMod)

# Second, compare Open, Shrub, Forest, Canopy Cover, Shrub and Canopy Cover, Forest and Shrub, and Shrub and Open
ccMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,by=spp,k=5),
                         weights = case, # 1 0 reponse
                         method = 'REML', # recommended
                         select = TRUE, # penalty to favor parsimonious model
                         family=cox.ph(),
                         data=wc)
AIC(ccMod)
forMod <- gam(cbind(dummy, indiv_step_f) ~  s(forest,by=spp,k=5),
              weights = case, # 1 0 reponse
              method = 'REML', # recommended
              select = TRUE, # penalty to favor parsimonious model
              family=cox.ph(),
              data=wc)
AIC(forMod)
shrMod <- gam(cbind(dummy, indiv_step_f) ~  s(shr,by=spp,k=5),
                          weights = case, # 1 0 reponse
                          method = 'REML', # recommended
                          select = TRUE, # penalty to favor parsimonious model
                          family=cox.ph(),
                          data=wc)
AIC(shrMod)
openMod <- gam(cbind(dummy, indiv_step_f) ~  s(open,by=spp,k=5),
                           weights = case, # 1 0 reponse
                           method = 'REML', # recommended
                           select = TRUE, # penalty to favor parsimonious model
                           family=cox.ph(),
                           data=wc)
AIC(openMod)
openShrMod <- gam(cbind(dummy, indiv_step_f) ~  s(open,by=spp,k=5) + s(shr,by=spp,k=5),
                              weights = case, # 1 0 reponse
                              method = 'REML', # recommended
                              select = TRUE, # penalty to favor parsimonious model
                              family=cox.ph(),
                              data=wc)
AIC(openShrMod)
forShrMod <- gam(cbind(dummy, indiv_step_f) ~  s(forest,by=spp,k=5) + s(shr,by=spp,k=5),
                             weights = case, # 1 0 reponse
                             method = 'REML', # recommended
                             select = TRUE, # penalty to favor parsimonious model
                             family=cox.ph(),
                             data=wc)
AIC(forShrMod)
shrCCMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,by=spp,k=5) + s(shr,by=spp,k=5),
                            weights = case, # 1 0 reponse
                            method = 'REML', # recommended
                            select = TRUE, # penalty to favor parsimonious model
                            family=cox.ph(),
                            data=wc)
AIC(shrCCMod)
# The model with both Shrub and Canopy Cover performs best and is therefore used throughout the rest of the code.

# Now that the best covariates for both predators are known, we'll run the same AIC-based comparisons to see if we can pare down the snow terms.
# For wolves:
wDepthDensMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + s(snod,k=5) + s(dens,k=5),
                                 weights = case, # 1 0 reponse
                                 method = 'REML', # recommended
                                 select = TRUE, # penalty to favor parsimonious model
                                 family=cox.ph(),
                                 data=wolf)
wDensMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + s(dens,k=5),
                            weights = case, # 1 0 reponse
                            method = 'REML', # recommended
                            select = TRUE, # penalty to favor parsimonious model
                            family=cox.ph(),
                            data=wolf)
wDepthMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + s(snod,k=5),
                             weights = case, # 1 0 reponse
                             method = 'REML', # recommended
                             select = TRUE, # penalty to favor parsimonious model
                             family=cox.ph(),
                             data=wolf)
wTIMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + ti(snod,dens,k=5),
                          weights = case, # 1 0 reponse
                          method = 'REML', # recommended
                          select = TRUE, # penalty to favor parsimonious model
                          family=cox.ph(),
                          data=wolf)
wFinalMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + te(snod,dens,k=5),
                 weights = case, # 1 0 reponse
                 method = 'REML', # recommended
                 select = TRUE, # penalty to favor parsimonious model
                 family=cox.ph(),
                 data=wolf)
AIC(wDepthMod)
AIC(wDensMod)
AIC(wDepthDensMod)
AIC(wTIMod)
AIC(wFinalMod)
# Unsurprisingly, the one labeled "FinalMod" works best (with a full tensor product of snow depth and density).

# Same for cougars:
cDepthDensMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(snod,k=5) + s(dens,k=5),
                                 weights = case, # 1 0 reponse
                                 method = 'REML', # recommended
                                 select = TRUE, # penalty to favor parsimonious model
                                 family=cox.ph(),
                                 data=coug)
cDensMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(dens,k=5),
                            weights = case, # 1 0 reponse
                            method = 'REML', # recommended
                            select = TRUE, # penalty to favor parsimonious model
                            family=cox.ph(),
                            data=coug)
cDepthMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + s(snod,k=5),
                             weights = case, # 1 0 reponse
                             method = 'REML', # recommended
                             select = TRUE, # penalty to favor parsimonious model
                             family=cox.ph(),
                             data=coug)
cTIMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + ti(snod,dens,k=5),
                          weights = case, # 1 0 reponse
                          method = 'REML', # recommended
                          select = TRUE, # penalty to favor parsimonious model
                          family=cox.ph(),
                          data=coug)
AIC(cDepthMod)
AIC(cDensMod)
AIC(cDepthDensMod)
AIC(cTIMod)
AIC(cFinalMod)
# Unsurprisingly, the one labeled "FinalMod" works best (with a full tensor product of snow depth and density).


#
# Part 2: Model Predator Step-Selection Functions
#


# Run actual final models; similarly computationally expensive to run as above models, each taking ~20-45 minutes to run.
wFinalMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr19,k=5) + s(tri,k=5) + s(deer,k=5) + te(snod,dens,k=5),
                             weights = case, # 1 0 reponse
                             method = 'REML', # recommended
                             select = TRUE, # penalty to favor parsimonious model
                             family=cox.ph(),
                             data=wolf)

cFinalMod <- gam(cbind(dummy, indiv_step_f) ~  s(cc,k=5) + s(Shr5,k=5) + s(tri,k=5) + s(deer,k=5) + te(snod,dens,k=5),
                             weights = case, # 1 0 reponse
                             method = 'REML', # recommended
                             select = TRUE, # penalty to favor parsimonious model
                             family=cox.ph(),
                             data=coug)

# Rather than take all the time to run these models, you can simply load the GAMs here:
#     wolf final SSF model:
wFinalMod <- readRDS("/Users/bensullender/Documents/UW/Research/ChapterTwo/Final_Code_for_Github/wFinalMod.rds")
#     cougar final SSF model:
cFinalMod <- readRDS("/Users/bensullender/Documents/UW/Research/ChapterTwo/Final_Code_for_Github/cFinalMod.rds")

# Visually investigate what these models look like:
gratia::draw(wFinalMod,rug=F)
gratia::draw(cFinalMod,rug=F)


#
# Part 3: Partial Effects Plots
#

wPartFX <- wolf
cPartFX <- coug

# Not very efficient, but relatively concise. Takes input animal steps and generate new predictions based on only one term at a time.
#       This shows partial effects over all observed values of that covariate.
wPartFX <- wPartFX %>% mutate(triP = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(tri)")$fit,
                              triSE = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(tri)")$se.fit,
                              ccP = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(cc)")$fit,
                              ccSE = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(cc)")$se.fit,
                              shrP = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(Shr19)")$fit,
                              shrSE = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(Shr19)")$se.fit,
                              deerP = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(deer)")$fit,
                              deerSE = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "s(deer)")$se.fit,
                              snowP = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "te(snod,dens)")$fit,
                              snowSE = predict.gam(wFinalMod,newdata=wPartFX,se.fit=T,terms = "te(snod,dens)")$se.fit)

cPartFX <- cPartFX %>% mutate(triP = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(tri)")$fit,
                              triSE = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(tri)")$se.fit,
                              ccP = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(cc)")$fit,
                              ccSE = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(cc)")$se.fit,
                              shrP = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(Shr5)")$fit,
                              shrSE = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(Shr5)")$se.fit,
                              deerP = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(deer)")$fit,
                              deerSE = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "s(deer)")$se.fit,
                              snowP = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "te(snod,dens)")$fit,
                              snowSE = predict.gam(cFinalMod,newdata=cPartFX,se.fit=T,terms = "te(snod,dens)")$se.fit)

# Create Partial Effects Plot for TRI
ggplot()+
  geom_ribbon(aes(x=tri,ymin=triP-1.96*triSE,ymax=triP+1.96*triSE,fill="wolf"),alpha=0.3,data=wPartFX)+
  geom_ribbon(aes(x=tri,ymin=triP-1.96*triSE,ymax=triP+1.96*triSE,fill="cougar"),alpha=0.3,data=cPartFX)+
  geom_line(aes(x=tri,y=triP,color="wolf"),data=wPartFX)+
  geom_line(aes(x=tri,y=triP,color="cougar"),data=cPartFX)+
  labs(x = "Terrain Ruggedness Index", y = 'Partial Effect', fill = "", colour = '')+
  scale_color_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  coord_cartesian(ylim=c(-5,5),expand=F)+
  theme_minimal(base_size = 24)

# Create Partial Effects Plot for Canopy Cover
ggplot()+
  geom_ribbon(aes(x=cc*100,ymin=ccP-1.96*ccSE,ymax=ccP+1.96*ccSE,fill="wolf"),alpha=0.3,data=wPartFX)+
  geom_ribbon(aes(x=cc*100,ymin=ccP-1.96*ccSE,ymax=ccP+1.96*ccSE,fill="cougar"),alpha=0.3,data=cPartFX)+
  geom_line(aes(x=cc*100,y=ccP,color="wolf"),data=wPartFX)+
  geom_line(aes(x=cc*100,y=ccP,color="cougar"),data=cPartFX)+
  labs(x = "Canopy Cover (%)", y = 'Partial Effect', fill = "", colour = '')+
  scale_color_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  coord_cartesian(ylim=c(-1,1),expand=F)+
  theme_minimal(base_size = 24)

# Create Partial Effects Plot for Percent Shrub
#       Rename shrub variables so we can plot them together.
wPartFX$shr <- wPartFX$Shr19
cPartFX$shr <- cPartFX$Shr5
ggplot()+
  geom_ribbon(aes(x=shr*100,ymin=shrP-1.96*shrSE,ymax=shrP+1.96*shrSE,fill="wolf"),alpha=0.3,data=wPartFX)+
  geom_ribbon(aes(x=shr*100,ymin=shrP-1.96*shrSE,ymax=shrP+1.96*shrSE,fill="cougar"),alpha=0.3,data=cPartFX)+
  geom_line(aes(x=shr*100,y=shrP,color="wolf"),data=wPartFX)+
  geom_line(aes(x=shr*100,y=shrP,color="cougar"),data=cPartFX)+
  labs(x = "Shrub Cover (%)", y = 'Partial Effect', fill = "", colour = '')+
  scale_color_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  coord_cartesian(ylim=c(-1,1),expand=F)+
  theme_minimal(base_size = 24)

# Create Partial Effects Plot for Deer Index
ggplot()+
  geom_ribbon(aes(x=deer,ymin=deerP-1.96*deerSE,ymax=deerP+1.96*deerSE,fill="wolf"),alpha=0.3,data=wPartFX)+
  geom_ribbon(aes(x=deer,ymin=deerP-1.96*deerSE,ymax=deerP+1.96*deerSE,fill="cougar"),alpha=0.3,data=cPartFX)+
  geom_line(aes(x=deer,y=deerP,color="wolf"),data=wPartFX)+
  geom_line(aes(x=deer,y=deerP,color="cougar"),data=cPartFX)+
  labs(x = "Deer Index", y = 'Partial Effect', fill = "", colour = '')+
  scale_color_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('wolf'="#9E0142",'cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  coord_cartesian(ylim=c(-1,1),expand=F)+
  theme_minimal(base_size = 24)

# Create tensor product plots showing partial effects of interactive continuous variables (snow depth and density)

# First, create predictions that we can interpolate from
wSnowTE <- gratia::smooth_estimates(wFinalMod, smooth = 'te(snod,dens)',n=500)
cSnowTE <- gratia::smooth_estimates(cFinalMod, smooth = 'te(snod,dens)',n=500)
min(wSnowTE$.estimate)
# Tensor product plot for wolves
ggplot(wSnowTE, aes(x=snod,y=dens,z=.estimate)) +
  geom_raster(aes(fill = .estimate), interpolate = T)+
  geom_contour(aes(z=.estimate),alpha=0.1,color="black",breaks=c(-0.5,-.4,-.3,-.2,-.1,0,.1,.2))+
  scale_fill_gradient2(low='black',mid="white",high='#9E0142',midpoint=0,name="Partial Effect")+
  #scale_fill_gradient(colours = c('black', "white",'#9E0142'),values=c(-0.5,1),name="Partial Effect")+
  #  geom_text_contour(aes(z=.estimate),alpha=0.5)+
  coord_cartesian(expand = FALSE,xlim=c(0,1.25))+
  theme_minimal(base_size = 24)+
  theme(axis.ticks.x = element_line(),axis.ticks.y=element_line())+
  labs(y = bquote('Snow Density'~(kg/m^3)), x = 'Snow Depth (m)')

# Tensor product plot for cougars
ggplot(cSnowTE, aes(x=snod,y=dens,z=.estimate)) +
  geom_raster(aes(fill = .estimate), interpolate = T)+
  geom_contour(aes(z=.estimate),alpha=0.1,color="black",breaks=c(-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1))+
  scale_fill_gradient2(low='black',mid="white",high='#FC7200',midpoint=0,name="Partial Effect")+
  #  scale_fill_gradientn(colours = c('black', "white",'#FC7200'),values=c(,1.5),name="Partial Effect")+
  #  geom_text_contour(aes(z=.estimate),alpha=0.5,breaks=c(-2,-1.75,-1.5,-1.25,-1,-.75,-.5,-.25,0,.25,.5,.75,1))+
  coord_cartesian(expand = FALSE,xlim=c(0,1.25))+
  theme_minimal(base_size = 24)+
  theme(axis.ticks.x = element_line(),axis.ticks.y=element_line())+
  labs(y = bquote('Snow Density'~(kg/m^3)), x = 'Snow Depth (m)')


#
# Part 4: Average Effects Plots
#


# For average effects plots, use only values of covariates at the available locations (Avgar et al. 2017)
# Create new data frame for each species, with model predictions based on all observed combinations of covariate values
wAvgDF <- wolf[wolf$case==0,]
wAvgDF$mod <- predict.gam(wFinalMod,newdata=wAvgDF)

cAvgDF <- coug[coug$case==0,]
cAvgDF$mod <- predict.gam(cFinalMod,newdata=cAvgDF)

# Next, we'll need to compare the distributions of snow depth and density for cougars and wolves. 
#     To do this, we'll make sure all covariates have the same names.

cRenamed <- cAvgDF %>% mutate(For = For5,Shr = Shr5,spp="Cougar") %>% dplyr::select(-For5,-Shr5,-Open5)
wRenamed <- wAvgDF %>% mutate(For = For19,Shr = Shr19,spp="Wolf") %>% dplyr::select(-For19,-Shr19,-Open19)
allRenamed <- rbind(cRenamed,wRenamed)
allRenamed$spp <- as.factor(allRenamed$spp)

ggplot()+
  geom_histogram(aes(x=dens,fill="Cougar"),bins=50,data=allRenamed[allRenamed$spp=="Cougar",])+
  geom_histogram(aes(x=dens,fill="Wolf"),bins=50,data=allRenamed[allRenamed$spp=="Wolf",])+
  labs(x = "Snow Density (kg/m^3)", y = 'Count')+
  #scale_color_manual(values = c('Wolf'="#9E0142",'Cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('Wolf'="#9E0142",'Cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  theme_minimal(base_size = 16)

ggplot()+
  geom_histogram(aes(x=snod,fill="Cougar"),bins=50,data=allRenamed[allRenamed$spp=="Cougar",])+
  geom_histogram(aes(x=snod,fill="Wolf"),bins=50,data=allRenamed[allRenamed$spp=="Wolf",])+
  labs(x = "Snow Depth (m)", y = 'Count')+
  #scale_color_manual(values = c('Wolf'="#9E0142",'Cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  scale_fill_manual(values = c('Wolf'="#9E0142",'Cougar'='#FC7200'),labels=c("Cougar","Wolf"))+
  labs(fill = "Species",color="Species")+
  coord_cartesian(xlim=c(0,1.5))+
  theme_minimal(base_size = 16)  

# It looks like wolves tend to travel in slightly more dense and deep snow.
#       However, it would get unnecessarily confusing to have four different breakpoints, so we will pool cougars and wolves.


# Suppress scientific notation
options(scipen=999)
quantile(allRenamed$snod,probs=seq(0,1,0.25))
# We will round these to 5cm (25%) and 30cm (75%)
quantile(allRenamed$dens,probs=seq(0,1,0.25))
# We will round these to 200kg/m3 (25%) and 350kg/m3 (75%)


# Now we'll split the data, for each species, using these quartiles
# dens = high density (>350kg/m3), lite = low density (<200kg/m3), shal = shallow depth (<5cm), deep = deep depth (>30cm)
# Because this is still an overwhelming number of points, we will add a "draw" variable that will let us randomly omit half of the points
cDens <- cAvgDF[cAvgDF$dens>=350,] %>% mutate(draw = round(runif(nrow(.))))
cLite <- cAvgDF[cAvgDF$dens<=200,] %>% mutate(draw = round(runif(nrow(.))))
cDeep <- cAvgDF[cAvgDF$snod>=.30,] %>% mutate(draw = round(runif(nrow(.))))
cShal <- cAvgDF[cAvgDF$snod<=.05,] %>% mutate(draw = round(runif(nrow(.))))

wDens <- wAvgDF[wAvgDF$dens>=350,] %>% mutate(draw = round(runif(nrow(.))))
wLite <- wAvgDF[wAvgDF$dens<=200,] %>% mutate(draw = round(runif(nrow(.))))
wDeep <- wAvgDF[wAvgDF$snod>=.30,] %>% mutate(draw = round(runif(nrow(.))))
wShal <- wAvgDF[wAvgDF$snod<=.05,] %>% mutate(draw = round(runif(nrow(.))))

# Wolves, snow depth x relative use
ggplot()+
  geom_point(aes(x=snod,y=mod),color="#CE025F",size=1,alpha=0.1,data=wDens[wDens$draw==1,])+
  geom_point(aes(x=snod,y=mod),color="gray30",size=1,alpha=0.1,data=wLite[wLite$draw==1,])+
  geom_smooth(aes(x=snod,y=mod,color="Light (<200kg/m^3)"),data=wLite)+
  geom_smooth(aes(x=snod,y=mod,color="Dense (>350kg/m^3)"),data=wDens)+
  scale_color_manual("Snow Density",values=c("Dense (>350kg/m^3)" ="#82023C","Light (<200kg/m^3)"="black"))+
  labs(x = "Snow depth (m)", y = 'Relative use', fill = "", linetype = "Snow Depth")+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(-2.5,2),xlim=c(0,1.5))+
  guides(linetype=guide_legend(keywidth=2.5,keyheight=1),color=guide_legend(keywidth=2.5,keyheight=1))

# Cougars, snow depth x relative use
ggplot()+
  geom_point(aes(x=snod,y=mod),color="#F98732",size=1,alpha=0.1,data=cDens[cDens$draw==1,])+
  geom_point(aes(x=snod,y=mod),color="gray30",size=1,alpha=0.1,data=cLite[cLite$draw==1,])+
  geom_smooth(aes(x=snod,y=mod,color="Light (<200kg/m^3)"),data=cLite)+
  geom_smooth(aes(x=snod,y=mod,color="Dense (>350kg/m^3)"),data=cDens)+
  scale_color_manual("Snow Density",values=c("Dense (>350kg/m^3)" ="#E56200","Light (<200kg/m^3)"="black"))+
  labs(x = "Snow depth (m)", y = 'Relative use', fill = "", linetype = "Snow Depth")+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(-2.5,2),xlim=c(0,1.5))+
  guides(linetype=guide_legend(keywidth=2.5,keyheight=1),color=guide_legend(keywidth=2.5,keyheight=1))

# Cougars, snow density x relative use
ggplot()+
  geom_point(aes(x=dens,y=mod),color="#F98732",size=1,alpha=0.1,data=cDeep)+
  geom_point(aes(x=dens,y=mod),color="gray30",size=1,alpha=0.1,data=cShal)+
  geom_smooth(aes(x=dens,y=mod,color="Shallow (<5cm)"),data=cShal)+
  geom_smooth(aes(x=dens,y=mod,color="Deep (>30cm)"),data=cDeep)+
  scale_color_manual("Snow Depth",values=c("Deep (>30cm)" ="#E56200","Shallow (<5cm)"="black"))+
  labs(x = bquote('Snow density'~(kg/m^3)), y = 'Relative use', fill = "", linetype = "Snow Depth")+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(-2,2),xlim=c(50,550))+
  guides(linetype=guide_legend(keywidth=2.5,keyheight=1),color=guide_legend(keywidth=2.5,keyheight=1))

# Wolves, snow density x relative use
ggplot()+
  geom_point(aes(x=dens,y=mod),color="#CE025F",size=1,alpha=0.1,data=wDeep)+
  geom_point(aes(x=dens,y=mod),color="gray30",size=1,alpha=0.1,data=wShal)+
  geom_smooth(aes(x=dens,y=mod,color="Shallow (<5cm)"),data=wShal)+
  geom_smooth(aes(x=dens,y=mod,color="Deep (>30cm)"),data=wDeep)+
  scale_color_manual("Snow Depth",values=c("Deep (>30cm)" ="#82023C","Shallow (<5cm)"="black"))+
  labs(x = bquote('Snow density'~(kg/m^3)), y = 'Relative use', fill = "", linetype = "Snow Depth")+
  theme_bw(base_size = 24)+
  coord_cartesian(ylim = c(-2,2),xlim=c(50,550))+
  guides(linetype=guide_legend(keywidth=2.5,keyheight=1),color=guide_legend(keywidth=2.5,keyheight=1))
