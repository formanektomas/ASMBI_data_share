#
rm(list=ls())
DF <- read.csv("ASMBIdf.csv")
DF <- DF[order(DF$time, DF$geo),]
DF$ERR <- DF$ERR*100 # percentages: 0-100 scale
DF$ERR_1 <- DF$ERR_1*100
#
# Connectivity matrix 
require(lme4)
require(splm)
require(spmoran)
coords <- DF[DF$time==2010,c("geo","long", "lat")]
IDs <- coords$geo
coords <- coordinates(coords[,2:3])
#
##############################################################################################
#
# Main estimation loop for plot
#
##############################################################################################
#
#
#
# Stability evaluation
s2.df <- NULL
s2.df <- data.frame(max.dist=0, LL=0, 
                    Unem=0, Unem.SE = 0, 
                    lRD_1=0, lRD_1.SE = 0, 
                    GDP_b=0, GDP_b.SE = 0, 
                    Emp_B_E=0, Emp_B_E.SE = 0,
                    logTotEn=0, logTotEn.SE=0, 
                    REng=0, REng.SE = 0)
####################################################################
set.seed(1300)
for(j in 16:70) {
  nb.km <- dnearneigh(coords, d1=0, d2=j*10, longlat=T, row.names = IDs)
  listw.km <- nb2listw(nb.km)
  C1 <- nb2mat(nb.km,style = "B")
  E1 <- meigen(cmat=C1)
  MEM <- kronecker(diag(7),E1$sf[,1:4])
  lme.2 <- lmer(lGDPpc ~ Unem+lRD_1+Y09GDPHab+Rel_Emp_B_E+logTotEn+ERR+NUTS0+MEM + (1|geo), data=DF)
  O1 <- summary(lme.2)
  loopresids <-  O1$residuals
  #
  #
  s2.df <- rbind(s2.df, c(j*10, O1$logLik, 
                          O1$coefficients[2,1], O1$coefficients[2,2], # Unem
                          O1$coefficients[3,1], O1$coefficients[3,2], # lRD_1
                          O1$coefficients[4,1], O1$coefficients[4,2], # Y09GDPHab
                          O1$coefficients[5,1], O1$coefficients[5,2], # Rel_Emp_B_E
                          O1$coefficients[6,1], O1$coefficients[6,2], # logTotEn
                          O1$coefficients[7,1], O1$coefficients[7,2])) # RenEnegyRate
} 
s2.df <- s2.df[-1,]
which(s2.df$LL==max(s2.df$LL))
s2.df[35,]
#
colnames(s2.df)
library(tidyr)
s2.df.long <- s2.df %>% gather(key="key", value="value", 
                               colnames(s2.df)[c(2:3,5,7,9,11,13)]) 
s2.df.long$SE <- with(s2.df.long, 
                 ifelse(key == "Unem", Unem.SE,
                 ifelse(key == "lRD_1", lRD_1.SE,
                 ifelse(key == "GDP_b", GDP_b.SE,
                 ifelse(key == "Emp_B_E", Emp_B_E.SE,    
                 ifelse(key == "logTotEn", logTotEn.SE,
                 ifelse(key == "REng", REng.SE,
                 NA)))))))
#
s2.df.long$key <- as.factor(s2.df.long$key)
levels(s2.df.long$key)
levels(s2.df.long$key) <- c("NACE2 B-E Employment (%)", 
                            "Base GDP pc (2009)", "LogLik", "log Energy Consumption","log R&D Exp (t-1)", 
                             "Renewable Energy Share", "Unemployment (%)")
s2.df.long$key <- factor(s2.df.long$key, levels=c("log Energy Consumption", "Renewable Energy Share", 
                                                  "NACE2 B-E Employment (%)", 
                                                  "Base GDP pc (2009)", "log R&D Exp (t-1)", "Unemployment (%)",   
                                                  "LogLik"))
# #
require(ggplot2)
ggplot(s2.df.long)+
  geom_line(aes(max.dist,value), size=0.8)+
  geom_ribbon(aes(ymin=value-(2*SE), ymax=value+(2*SE), x=max.dist), alpha=0.3, colour="blue", fill="grey")+
  geom_vline(aes(xintercept=500), linetype = 2)+
  facet_wrap(~key, scales="free", nrow=3)+
  xlab("Maximum neighbor distance (km)")+
  theme_minimal()
#
##############################################################################################