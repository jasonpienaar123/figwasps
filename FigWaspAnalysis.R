#figwasps

library(geiger)
library(ggplot2)
library(ape)
library(mvSLOUCH) # preliminary analysis
library(slouch)
library(olsrr)
library(car)

# set working directory to ("../manuscript/analysis") once downloaded
# e.g. in Rstudio, once downloaded, set working directory to r source file location 

# Read in tree 
WaspChrono<-read.tree("chronogram_sexratio_1myr.nwk")
plot(WaspChrono, cex=0.4)
class(WaspChrono)

# Read in Data set
Dat<-read.csv("CleanData.csv")
head(Dat)

# rerun analysis from here after statistical outliers are removed (see below)

# Prune tree to include only wasps which have data for analysis 
tree_wasps_names<-geiger::name.check(WaspChrono,Dat$voucher, data.names=Dat$voucher)
tree_wasps_names
tips_todrop<-tree_wasps_names$tree_not_data
PrunedWaspChrono<-ape::drop.tip(WaspChrono,tips_todrop)
PrunedWaspChrono<-ape::di2multi(PrunedWaspChrono)
rownames(Dat)<-Dat$voucher
Dat<-Dat[PrunedWaspChrono$tip.label,]
mvSLOUCH::phyltree_paths(PrunedWaspChrono)$tree_height
tree_height<-mvSLOUCH::phyltree_paths(PrunedWaspChrono)$tree_height
ScaledPrunedWaspChrono<-PrunedWaspChrono
ScaledPrunedWaspChrono$edge.length<-ScaledPrunedWaspChrono$edge.length/tree_height
mvSLOUCH::phyltree_paths(ScaledPrunedWaspChrono)$tree_height

geiger::name.check(ScaledPrunedWaspChrono, Dat)
rownames(Dat)<-Dat$voucher
Dat<-Dat[ScaledPrunedWaspChrono$tip.label,]

# Speeds up multivariate estimation (not run here - for future analysis)
# mvStree<-mvSLOUCH::phyltree_paths(ScaledPrunedWaspChrono)

#_______________________________________________________________________________
# SONS 

# Brownian motion model
brown.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
          response = Dat$mean_males, mv.response = Dat$mv_males)
# global regime model
OU_global<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                       response = Dat$mean_males,mv.response = Dat$mv_males)
OU_global

# Regressions on randomly evolving niche variables 
# Multi foundress frequency 
OU_RandMFP<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                        response = Dat$mean_males, random.cov = (1-Dat$single_foundress_frequency), mv.response = Dat$mv_males) 
OU_RandMFP

# Flower numbers
OU_RandFN<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                       response = Dat$mean_males, random.cov = Dat$mean_figflowernumber, mv.response = Dat$mv_males) 
OU_RandFN              
                
# Daughters and multifoundress frequency                
randcov2<-cbind(Dat$mean_females, (1-Dat$single_foundress_frequency))
                      colnames(randcov2)<-c("Females", "MFP")
                      mv_randcov2<-cbind(Dat$mv_females, rep(0, times=nrow(Dat)))
OU_RandFemMFP<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,response = Dat$mean_males, random.cov = randcov2, mv.response = Dat$mv_males, mv.random.cov=mv_randcov2)                 
OU_RandFemMFP                
                
# Daughters and flower number                 
randcov3<-cbind(Dat$mean_females, Dat$mean_figflowernumber)               
                      colnames(randcov3)<-c("Females", "Flowers")
                      # use mv_randcov2 for measurement variance 
OU_RandFemFlowers<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,response = Dat$mean_males, random.cov = randcov3, mv.response = Dat$mv_males, mv.random.cov = mv_randcov2)              
OU_RandFemFlowers               
                
# Flower number and multi foundress probablity 
randcov4<-cbind(Dat$mean_figflowernumber, (1-Dat$single_foundress_frequency))
                colnames(randcov4)<-c("Flowers", "MFP")
OU_RandMFPFlowers<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,response = Dat$mean_males, random.cov = randcov4, mv.response = Dat$mv_males)                
OU_RandMFPFlowers               

# Daughters, multi foundress probability and flower number
randcov5<-cbind(Dat$mean_females, (1-Dat$single_foundress_frequency), Dat$mean_figflowernumber)
                colnames(randcov5)<-c("Females", "MFP", "Flowers")
                mv_randcov5<-cbind(Dat$mv_females, rep(0, times=nrow(Dat)), rep(0, times=nrow(Dat)))
OU_RandFemMFPFlowers<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,response = Dat$mean_males, random.cov = randcov5, mv.response = Dat$mv_males, mv.random.cov = mv_randcov5)               
OU_RandFemMFPFlowers               
summary(OU_RandFemMFPFlowers)               

# Grid search for best Sons model                 
h<-1:100*0.01
vy<-1:100               
OU_RandFemMFP_gs<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher, 
                              hl_values = h, vy_values = vy,response = Dat$mean_males, random.cov = randcov2, mv.response = Dat$mv_males, mv.random.cov = mv_randcov2, hillclimb = FALSE)
plot(OU_RandFemMFP_gs)              

#_______________________________________________________________________________
# DAUGHTERS

# Brownian motion model                
brown.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                          response = Dat$mean_females, mv.response = Dat$mv_females)

# Global regime model
OU_global_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                           response = Dat$mean_females, mv.response = Dat$mv_females)
OU_global_fem

# Regressions on Randomly evolving niche variables 
# Multi foundress probablity 
OU_RandMFP_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                            response = Dat$mean_females, random.cov = (1-Dat$single_foundress_frequency), mv.response = Dat$mv_females)                 
OU_RandMFP_fem            

# Flower number                
OU_RandFN_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                           response = Dat$mean_females, random.cov = Dat$mean_figflowernumber, mv.response = Dat$mv_females)                
OU_RandFN_fem               

# Sons                
OU_RandSons_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                            response = Dat$mean_females, random.cov = Dat$mean_males, mv.response = Dat$mv_females, mv.random.cov=Dat$mv_males) 
OU_RandSons_fem
                
# Sons and multi foundress probablity                
randcov2<-cbind(Dat$mean_males, (1-Dat$single_foundress_frequency))
                 colnames(randcov2)<-c("Males", "MFP")  
                 mv_randcov2<-cbind(Dat$mv_males, rep(0, times=nrow(Dat)))
OU_RandSonsMFP_fem<-slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                                  response = Dat$mean_females, random.cov = randcov2, mv.response = Dat$mv_females, mv.random.cov = mv_randcov2)                 
OU_RandSonsMFP_fem                 

# Sons and flower number               
randcov3<-cbind(Dat$mean_males, Dat$mean_figflowernumber)
                 colnames(randcov3)<-c("Males", "Flowers")
                # use mv_randcov2 for measurement variance
OU_RandSonsFlowers_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                                    response = Dat$mean_females, random.cov = randcov3, mv.response = Dat$mv_females, mv.random.cov = mv_randcov2)                  
OU_RandSonsFlowers_fem                
         
# multi foundress probablity and flower number                 
randcov4<-cbind((1-Dat$single_foundress_frequency), Dat$mean_figflowernumber)
                colnames(randcov4)<-c("MFP", "Flowers")                
OU_RandMFPFlowers_fem<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                                   response = Dat$mean_females, random.cov = randcov4, mv.response = Dat$mv_females)
OU_RandMFPFlowers_fem

# Sons,  multi foundress probablity and flower number 
randcov5<-cbind(Dat$mean_males, (1-Dat$single_foundress_frequency), Dat$mean_figflowernumber)
                colnames(randcov5)<-c("Males", "MFP", "Flowers")                
                mv_randcov5<-cbind(Dat$mv_males, rep(0, times=nrow(Dat)), rep(0, times=nrow(Dat)))
OU_RandSonsMFPFlowers<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher,
                                   response = Dat$mean_females, random.cov = randcov5, mv.response = Dat$mv_females, mv.random.cov=mv_randcov5) 
OU_RandSonsMFPFlowers

# Grid search for best model
h<-1:125*0.001
vy<-1:50*100
OU_RandSons_fem_gs<- slouch.fit(phy = ScaledPrunedWaspChrono, species = Dat$voucher, hl_values = h, vy_values = vy, response = Dat$mean_females, random.cov = Dat$mean_males, mv.response = Dat$mv_females, mv.random.cov = Dat$mv_males, hillclimb = FALSE) 
plot(OU_RandSons_fem_gs)

#_______________________________________________________________________________
# Oulier analysis
# Sons
MFP <-(1-Dat$single_foundress_frequency)
model <- lm(Dat$mean_males ~ Dat$mean_females + MFP + Dat$mean_figflowernumber)
plot(model)
ols_test_normality(model) # normality 
bp_test <- ncvTest(model)
print(bp_test) # heteroscedasticity

Dat1<-Dat[-26,]
MFP <-(1-Dat1$single_foundress_frequency)
model1 <- lm(Dat1$mean_males ~ Dat1$mean_females + MFP + Dat1$mean_figflowernumber)
plot(model1)
ols_test_normality(model1)
bp_test <- ncvTest(model1)
print(bp_test) # heteroscedasticity

Dat2<-Dat1[-13,]
MFP <-(1-Dat2$single_foundress_frequency)
model2 <- lm(Dat2$mean_males ~ Dat2$mean_females + MFP + Dat2$mean_figflowernumber)
plot(model2)
ols_test_normality(model2)
bp_test <- ncvTest(model2)
print(bp_test)

# rerun models with outliers removed (go back to top where indicated)
Dat<-Dat2

#_______________________________________________________________________________
      
                         