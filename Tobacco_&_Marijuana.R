# library
library(tidyverse)
library(survival)
library(survminer)
library(gregmisc)
#==================================================
# Balloonplot
#==================================================
LTx$Marijuana <- with(LTx, ifelse(m.history==0,"Never",
                                  ifelse(m.history==1,"Former",
                                         ifelse(m.history==2,"Current",NA))))
LTx$Marijuana <- factor(LTx$Marijuana, levels=c("Never","Former","Current"))

LTx$Tobacco<- with(LTx, ifelse(Cig_N_F_C==0,"Never",
                               ifelse(Cig_N_F_C==1,"Former",
                                      ifelse(Cig_N_F_C==2,"Current",NA))))
LTx$Tobacco <- factor(LTx$Tobacco, levels=c("Never","Former","Current"))

TM <- count(LTx,Tobacco,Marijuana)

balloonplot(TM$Marijuana,TM$Tobacco,TM$n,label.size=1.5,dotsize=2.5/max(strwidth(19),strheight(19)),dotcolor=8,
            colmar=0.5,rowmar=0.5, main="Classification of smoking history",xlab="Marijuana",ylab="Tobacco",text.size=1.2,cum.margins=F,show.margins=F,sorted=T,label.lines=T)
#==================================================
# category
#==================================================
LTx$Category <- with(LTx, ifelse(Cig_N_F_C==0&m.history==0,"NS donor",
                                 ifelse(Cig_N_F_C==0 & m.history==1|Cig_N_F_C==1 & m.history==0|Cig_N_F_C==1 & m.history==1,"FS donor",
                                        ifelse(Cig_N_F_C==2|m.history==2,"CS donor",NA))))
LTx$Category <- factor(LTx$Category,levels=c("NS donor","FS donor","CS donor"))
#==================================================
# Survival curve
#==================================================
LTx <- LTx %>% mutate(surv3y_year=surv3y/365)
sf_graft_3y<-survfit(Surv(surv3y_year,graft_3y)~Category, data=LTx)
ggsurv <- ggsurvplot(sf_graft_3y, break.time.by=1,title = "",legend.title = "",legend=c(.8,.4),legend.labs=c("NS donor","FS donor","CS donor"),font.legend=list(color = "black"),
                     risk.table=T,palette=c("lightcoral" ,"green3" ,"cornflowerblue"),linetype = c(1,4,3),risk.table.y.text = T,
                     xlab="Time (years)",xlim=c(0,3),ylim=c(0,1))
ggpar(ggsurv,
      font.tickslab = 17,
      font.x        = 17,
      font.y        = 17)
#==================================================
# stepwise for multivariate covariates
#==================================================
library(MASS)
COX <- coxph(Surv(surv3y,graft_3y)~Age+Sex+donor.Last.PO2+Aspiration_checked+GIT+
                R_Age+R_sex+Dx_PH_CLAD+
                SLT+CPB, data=LTx)
stepAIC(COX)
LTx <- LTx %>% mutate(GIT.60=GIT/60)

# main analysis; + GIT + Dx_PH_CLAD + SLT
COX <- coxph(Surv(surv3y,graft_3y)~current_cig+current_m+GIT.60 + Dx_PH_CLAD + SLT, data=LTx)
summary(COX)

# sensitivity analysis
## donor
COX <- coxph(Surv(surv3y,graft_3y)~current_cig+current_m+Age_more_55 + GIT.60 + PF_300_CoLog, data=LTx)
summary(COX)

## Recipient
COX <- coxph(Surv(surv3y,graft_3y)~current_cig+current_m+Dx_PH_CLAD + SLT+CPB, data=LTx)
summary(COX)

#==================================================
# Cumulative incidence
#==================================================
Dataset <- read.csv("Review.csv")
Dataset$Category <- with(Dataset,ifelse(classification_0_NS_1_FS_2_CS==0, "NS donor",
                                        ifelse(classification_0_NS_1_FS_2_CS==1,"FS donor","CS donor")))
Dataset$Category <- factor(Dataset$Category,levels=c("NS donor","FS donor","CS donor"))

# result_setting (including re-LTx)
Dataset$all_result <- with(Dataset,ifelse(reLTx==1,4,C_death_3y))

# Setup for analysis
DummyEventForCI <- Dataset$all_result
DummyEventForCI <- as.factor(DummyEventForCI) # Required for Surv() with mstate option
res <- NULL
ci <- NULL
ci.summary.table <- NULL
ci <- survfit(Surv((surv3y_year/1), DummyEventForCI, type="mstate")~Category, data=Dataset)
if(is.null(ci$surv) & is.null(ci$prev)) ci$surv <- 1-ci$pstate
library(cmprsk)
res <- with(Dataset, cuminc((surv3y_year/1), all_result, Category, cencode=0, na.action = na.omit))
compevents <- levels(factor(Dataset$all_result))
nevents <- length(compevents)
if (compevents[1]=="0") {compevents <- compevents[2:nevents]; nevents <- nevents - 1}
len <- nchar("Category")
groups <- substring(names(ci$strata), len+2)
ngroups <- length(groups)
legend <- groups

# CI plot
par(lwd=2)
plot(ci[,1],main="A. CLAD-Cause Mortality", bty="l", col=c("lightcoral" ,"green3" ,"cornflowerblue"), lty=c(1,4,3), lwd=5,ylim=c(0, 0.3), xaxp=c(0,3,3),yaxp=c(0,0.3,3),
     cex.lab=2.6,cex.axis=2.6,cex.main=3.0,xlab="Time (years)", ylab="", mark.time=TRUE,las=1)

plot(ci[,2],main="B. Infection-Cause Mortality", bty="l", col=c("lightcoral" ,"green3" ,"cornflowerblue"), lty=c(1,4,3), lwd=5, ylim=c(0, 0.3), xaxp=c(0,3,3),yaxp=c(0,0.3,3),
     cex.lab=2.6,cex.axis=2.6,cex.main=3.0,xlab="Time (years)", ylab="", mark.time=TRUE,las=1)

plot(ci[,3],main="C. Others-Cause Mortality", bty="l", col=c("lightcoral" ,"green3" ,"cornflowerblue"), lty=c(1,4,3), lwd=5, ylim=c(0, 0.3), xaxp=c(0,3,3),yaxp=c(0,0.3,3),
     cex.lab=2.6,cex.axis=2.6,cex.main=3.0,xlab="Time (years)", ylab="", mark.time=TRUE,las=1)