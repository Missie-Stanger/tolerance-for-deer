# In this code I first used Exploratory Factor Analysis to determine the 
# latent variable tolerance. Then I used Confirmatory Factor Analysis to 
# confirm this analysis as suggested by Costello and Osborne (2005). 
# Then last of all I created a linear model explaining tolerance. 

landowner <- read.csv("C:/Users/Andre/Desktop/Missie/Thesis/Oh Deer/old stuff/Private Land and Deer in Ohio All Data.csv", header=T)
names(landowner)

####### get rid of 9 in acceptance capacity; 9=no opinion
landowner$D2.acceptance_capacity[landowner$D2.acceptance_capacity == 9] <- "NA"
## switch back to numeric
landowner$D2.acceptance_capacity <- as.numeric(landowner$D2.acceptance_capacity)
summary(landowner$D2.acceptance_capacity)
# reverse code some variables
landowner$rA1b.worry_problems<- 6-landowner$A1b.worry_problems
landowner$rA1c.nuisance<- 6-landowner$A1c.nuisance
landowner$rD2.acceptance_capacity<- 4-landowner$D2.acceptance_capacity
#write.csv(landowner, "landowner.csv")
down<-na.omit(landowner[,c(128,129)])
plot(down$rD2.acceptance_capacity, down$affect)
library(psych)
library(MVN)
library(ggplot2)
library(car)
library(stargazer)
library(lmtest)
library(gvlma)
#library(sjPlot)
#library(GPArotation)
citation("psych")

## CREATE LATENT VARIABLE TOLERANCE
landowner$Std_A1a.enjoy_seeing <- scale(landowner$A1a.enjoy_seeing)
landowner$Std_A1d.comfort <- scale(landowner$A1d.comfort)
landowner$Std_rA1b.worry_problems <- scale(landowner$rA1b.worry_problems)
landowner$Std_rA1c.nuisance <- scale(landowner$rA1c.nuisance)
landowner$Std_rD2.acceptance_capacity <- scale(landowner$rD2.acceptance_capacity)
landowner$tolerance <- (landowner$Std_A1a.enjoy_seeing + landowner$Std_A1d.comfort + landowner$Std_rA1b.worry_problems
                        + landowner$Std_rA1c.nuisance + landowner$Std_rD2.acceptance_capacity) / 5
###create variable for conservation stewardship
landowner$steward <- (landowner$E2a.volunteer+landowner$E2b.rec_activities
                      +landowner$E2c.plant+landowner$E2d.money
                      +landowner$E2e.practices)/5

### Standardize all variabless for Factor Analysis
landowner$Std_A3a.did_you_hunt <- scale(landowner$A3a.did_you_hunt)
landowner$Std_C3.lease <- scale(landowner$C3.lease)
landowner$Std_D5.garden_damage <- scale(landowner$D5.garden_damage)
landowner$Std_D4.commercial_damage <- scale(landowner$D4.commercial_damage)
landowner$Std_E3.vehicle_collision <- scale(landowner$E3.vehicle_collision)
landowner$Std_B1.operate_ag <- scale(landowner$B1.operate_ag)
landowner$Std_E5a.farmer.rancher <- scale(landowner$E5a.farmer.rancher)
landowner$Std_E5b.hunter <- scale(landowner$E5b.hunter)
landowner$Std_E5d.enviromentalist <- scale(landowner$E5d.enviromentalist)
landowner$Std_E5e.animal_rights_advocate <- scale(landowner$E5e.animal_rights_advocate)
landowner$Std_E5d.enviromentalist <- scale(landowner$E5d.enviromentalist)
landowner$Std_E5f.property_rights_advocate <- scale(landowner$E5f.property_rights_advocate)
head(landowner)

###create variable for poaching social norm
landowner$rD3d.account_community <- 6 - landowner$D3d.account_community
landowner$rD3e.account_officer <- 6 - landowner$D3e.account_officer
landowner$poach <- (landowner$D3a.poaching_common+landowner$D3b.people_accept+landowner$D3c.area_accept+landowner$rD3d.account_community+landowner$rD3e.account_officer)/5

### create a smaller data frame to work with
names(landowner)
owner<-na.omit(landowner[,c(9,68,96,97,104,18,134,135,2,106:111,129:133,136:146,149,122)])## All Landowners
describe(owner)
names(owner)
head(owner)

#####################################
### Testing Assumptions for Factor Analysis
########################################

###Henze-Zirkler's MVN test
result <- mvn(data = owner[,c(16:26)], mvnTest = "hz")
result$multivariateNormality

#Doornik-Hansen's MVN test
result <- mvn(data = owner[,c(16:26)], mvnTest = "dh")
result$multivariateNormality

##Kernel density plot (from Anderson PCA in R)
mu<-colMeans(owner[ ,16:26]);sig<-cov(owner[ ,16:26])
stopifnot(mahalanobis(owner[ ,16:26],0,diag(ncol(owner[ ,16:26])))==rowSums(owner[ ,16:26]*owner[ ,16:26]))
D2<-mahalanobis(owner[ ,16:26],mu,sig)
n=dim(owner[ ,16:26])[1]
p=dim(owner[ ,16:26])[2]
par(mfrow=c(1,2))
plot(density(D2,bw=.5),main=paste("Mahalanobis Distances"))
qqplot(qchisq(ppoints(n),df=p),D2,main="Q-Q Plot of Mahalanobis D Squared vs. Quantiles of Chi Squared")
abline(0,1,col='gray')


###################
### EFA ############
################
names(owner)
describe(owner)
head(owner)
#fa.owner_IVs<-factanal(owner[ ,c(16:26)], factors=2, rotation="promax") # use an orthogonal rotation (varimax for example) since the factors correlate under .32 (Tabachnick et al., book entitled "Using Multivatiate Statistics" 5th edition)
fa.owner_IVs<-factanal(owner[ ,c(16:26)], factors=2, rotation="varimax")
fa.owner_IVs
fa.owner_IVs$loadings
print(loadings(fa.owner_IVs), digits = 2, cutoff = .2, sort = TRUE)

#plot factor 1 by factor 2
load <- fa.owner_IVs$loadings[,1:2]
par(mfrow=c(1,1))
plot(load,type="n") # set up plot
text(load,labels=names(owner),cex=1) # add variable names

#Function to make a scree plot (from PCA and EFA, an example):
screeplot.factanal <- function(fa.fit,xlab="factor",ylab="eigenvalue",...) {
  sosq <- function(v) {sum(v^2)}
  my.loadings <- as.matrix(fa.fit$loadings)
  evalues <- apply(my.loadings,2,sosq)
  plot(evalues,xlab=xlab,ylab=ylab,...)
}

screeplot.factanal(factanal(owner[ ,c(16:26)],factors=2,rotation="varimax"))
# I have identifyed the latent variable tolerance on the first axis of the EFA.
# The variable 
#Interpret factors (Steiger EFA with R):
#1. Examine variables with high loadings on each factor.
#2. Decide what construct is common to the variables that load high on each factor.
#3. Name the each factor after the constructs.
#Factor1 = "Tolerance"


#############
### CFA #########
###############
names(owner)
describe(owner)
fa.owner_IVs<-factanal(owner[ ,c(16:20)], factors=1, rotation="varimax")
fa.owner_IVs
print(loadings(fa.owner_IVs), digits = 2, cutoff = .2, sort = TRUE)

#Function to make a scree plot (from PCA and EFA, an example):
screeplot.factanal <- function(fa.fit,xlab="factor",ylab="eigenvalue",...) {
  sosq <- function(v) {sum(v^2)}
  my.loadings <- as.matrix(fa.fit$loadings)
  evalues <- apply(my.loadings,2,sosq)
  plot(evalues,xlab=xlab,ylab=ylab,...)
}

screeplot.factanal(factanal(owner[ ,c(16:20)],factors=1,rotation="varimax"))

#Interpret factors (Steiger EFA with R):
#1. Examine variables with high loadings on each factor.
#2. Decide what construct is common to the variables that load high on each factor.
#3. Name the each factor after the constructs.
#Factor1 = "Tolerance"


###############################
############ GLM ############
#########################


mod <- glm(owner$tolerance ~ owner$D5.garden_damage + owner$D4.commercial_damage 
             + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist 
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$DMURes , family=gaussian(link=identity))

summary(mod)
par(mfrow=c(2,2))
plot(mod)

#Testing fit  
#Diagnostics--Goodness of Fit for mod:
1-pchisq(deviance(mod),df.residual(mod))
#P=0.00 so strong evidence for lack of fit, dispesion parameter(want it near 1)
sum(residuals(mod,type="pearson")^2)/mod$df.res
#dispersion parameter indicates severe overdispersion; r^2 (deviance explained)
(mod$null.deviance-mod$deviance)/mod$null.deviance

############ Test for multicollinearity
### first print out HTML code or text for a correlation matrix
correlation.matrix <- cor(owner[,c(1:8,10:15)])
stargazer(correlation.matrix, type="text")
#for Varience Inflation Factor, anything over 2.5 may be problematic (some literature suggests 5 as a cutoff)
vif(mod)

# Test for Heteroscedasticity (if p-value is more the 0.05 then you are good)
ncvTest(mod)
bptest(mod)

# Test for Outliers
outlierTest(mod)
crPlots(mod)
describe(owner)
names(owner)

#### Extra cool trick:
# Testing Linear Model Assumptions
gvlma(mod15)


#############################################
#### Code for candidate model set ###########
#############################################
mod1 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
           + owner$E3.vehicle_collision)
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)


mod2 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
           + owner$E3.vehicle_collision + owner$B1.operate_ag)
summary(mod2)
par(mfrow=c(2,2))
plot(mod2)


mod4 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
           + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
           + owner$E5b.hunter + owner$E5c.conservationist + owner$E5d.enviromentalist
           + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate)
summary(mod4)
par(mfrow=c(2,2))
plot(mod4)
describe(owner)
head(owner)

mod5 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
           + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
           + owner$E5b.hunter + owner$E5c.conservationist + owner$E5d.enviromentalist
           + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
           + owner$steward)
summary(mod5)
par(mfrow=c(2,2))
plot(mod5)


mod8 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
           + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
           + owner$E5b.hunter + owner$E5d.enviromentalist
           + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
           + owner$steward)
summary(mod8)
par(mfrow=c(2,2))
plot(mod8)


mod10 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage)
summary(mod10)
par(mfrow=c(2,2))
plot(mod10)

mod11 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage + owner$B1.operate_ag)
summary(mod11)
par(mfrow=c(2,2))
plot(mod11)

mod12 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage 
            + owner$DMURes)
summary(mod12)
par(mfrow=c(2,2))
plot(mod12)

mod13 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$B1.operate_ag
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage 
            + owner$DMURes)
summary(mod13)
par(mfrow=c(2,2))
plot(mod13)

mod14 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$B1.operate_ag*owner$E5a.farmer.rancher
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage 
            + owner$DMURes)
summary(mod14)
par(mfrow=c(2,2))
plot(mod14)

mod15 <- lm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
            + owner$E3.vehicle_collision + owner$B1.operate_ag+owner$E5a.farmer.rancher
            + owner$E5b.hunter + owner$E5d.enviromentalist
            + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
            + owner$steward + owner$D4.commercial_damage 
            + owner$DMURes)
summary(mod15)
par(mfrow=c(2,2))
plot(mod15)

mod16 <- glm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$B1.operate_ag+owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage 
             + owner$DMURes, family=gaussian(link=identity))
summary(mod16)
par(mfrow=c(2,2))
plot(mod16)

#Testing fit
#Diagnostics--Goodness of Fit for modelb:
1-pchisq(deviance(mod16),df.residual(mod16))
#P=0.00 so strong evidence for lack of fit, dispesion parameter(want it near 1)
sum(residuals(mod16,type="pearson")^2)/mod16$df.res
#dispersion parameter indicates severe overdispersion; r^2 (deviance explained)
(mod16$null.deviance-mod16$deviance)/mod16$null.deviance
#pseudo R-squared not bad??

mod17 <- glm(owner$tolerance ~ owner$A3a.did_you_hunt + owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$B1.operate_ag+owner$E5a.farmer.rancher
             + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage
             + owner$DMURes, family=gaussian(link=identity))
summary(mod17)
par(mfrow=c(2,2))
plot(mod17)
#Testing fit
#Diagnostics--Goodness of Fit for modelb:
1-pchisq(deviance(mod17),df.residual(mod17))
#P=0.00 so strong evidence for lack of fit, dispesion parameter(want it near 1)
sum(residuals(mod17,type="pearson")^2)/mod17$df.res
#dispersion parameter indicates severe overdispersion; r^2 (deviance explained)
(mod17$null.deviance-mod17$deviance)/mod17$null.deviance
#pseudo R-squared not bad??

mod18 <- glm(owner$tolerance ~ owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$B1.operate_ag+owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage
             + owner$DMURes, family=gaussian(link=identity))
summary(mod18)
par(mfrow=c(2,2))
plot(mod18)
#Testing fit
#Diagnostics--Goodness of Fit for modelb:
1-pchisq(deviance(mod18),df.residual(mod18))
#P=0.00 so strong evidence for lack of fit, dispesion parameter(want it near 1)
sum(residuals(mod18,type="pearson")^2)/mod18$df.res
#dispersion parameter indicates severe overdispersion; r^2 (deviance explained)
(mod18$null.deviance-mod18$deviance)/mod18$null.deviance
#pseudo R-squared not bad??

mod19 <- glm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage 
             + owner$DMURes, family=gaussian(link=identity))
summary(mod19)
par(mfrow=c(2,2))
plot(mod19)
names(owner)
#Testing fit
#Diagnostics--Goodness of Fit for modelb:
1-pchisq(deviance(mod19),df.residual(mod19))#P=0.00 so strong evidence for lack of fit
#dispesion parameter(want it near 1)
sum(residuals(mod19,type="pearson")^2)/mod19$df.res#dispersion parameter indicates underdispersion
#r^2 (deviance explained)
(mod19$null.deviance-mod19$deviance)/mod19$null.deviance
#pseudo R-squared not bad??

mod99 <- glm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage + owner$poach
             + owner$DMURes, family=gaussian(link=identity))
summary(mod99)
par(mfrow=c(2,2))
plot(mod99)

mod999 <- glm(owner$tolerance ~ owner$A3a.did_you_hunt  + owner$D5.garden_damage 
             + owner$E3.vehicle_collision + owner$E5a.farmer.rancher
             + owner$E5b.hunter + owner$E5d.enviromentalist
             + owner$E5e.animal_rights_advocate + owner$E5f.property_rights_advocate
             + owner$steward + owner$D4.commercial_damage + owner$poach + owner$E7.gender
             + owner$DMURes, family=gaussian(link=identity))

summary(mod999)
par(mfrow=c(2,2))
plot(mod999)


AIC(mod,mod1,mod2,mod4,mod5,mod8,mod10,mod11,mod12,mod13,mod14,mod15,mod16,mod17
    ,mod18,mod19,mod99,mod999)

