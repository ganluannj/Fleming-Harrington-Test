library(FHtest)
library(KMsurv)
library(survival)

#VA lung cancer data 

#Patients with advanced, inoperable lung cancer were treated with 
#chemotherapy. 
#N = 137 
#  
# Variables 
# Treatment 1=standard, 2=test 
# Cell type 1=squamous, 2=small cell, 3=adeno, 4=large 
# Survival in days 
# Status 1=dead, 0=censored 
# Karnofsky score (measure of general performance, 100=best) 
# Months from Diagnosis 
# Age in years 
# Prior therapy 0=no, 10=yes 



VA<-read.table("D:\\Googledrive\\courses\\Math 659 Survival Analysis\\project\\VA lung cancer.txt")
names(VA)<-c("treatment", "Cell", "Time", "Status",
                "Score", "Diagnosis", "Age", "Prior")
VA<-as.data.frame(VA)

hist(VA$Age, xlab="patients' age", cex.lab=1.3,
     main="Histogram of patients' age")
abline(v=60, col="red", lwd=3)
text (40, 30, "Group I", cex=1.3, col="red")
text (80, 30, "Group II", cex=1.3, col="red")

AgeI<-as.numeric(factor(VA$Age>59))
VA<-cbind(VA, AgeI)

fit1<-survfit(Surv(Time, Status)~1, data=VA, subset=(AgeI==1))
fit2<-survfit(Surv(Time, Status)~1, data=VA, subset=(AgeI==2))
#xmax<-max(VA$Time)
plot(fit1, col = 'red', conf.int=FALSE,xlim=c(0, 500), 
     main="Kaplan-Meier survival curves", xlab="Days",
     ylab="Estimated Survival", cex.lab=1.3)
par(new=TRUE)
plot(fit2, col ='blue', conf.int=FALSE, xlim=c(0, 500))
legend(300, 0.9, legend=c("Group I", "Group II"),
       col=c("red", "blue"), lty=c(1,1), cex=1.1)

# Log-rank test
lr<-survdiff(Surv(Time, Status)~AgeI, data=VA)
lrp<-1-pchisq(lr$chisq,1) # 0.09
lrp

# FH p=0, q=1
fhl<-FHtestrcc(Surv(Time, Status)~AgeI, data=VA, rho=0, lambda=1)
fhlp<-fhl$pvalue # 0.014
fhlp

# FH p=0, q=3
fhl<-FHtestrcc(Surv(Time, Status)~AgeI, data=VA, rho=0, lambda=3)
fhlp<-fhl$pvalue # 0.02
fhlp

# FH p=0, q=2
fhl<-FHtestrcc(Surv(Time, Status)~AgeI, data=VA, rho=0, lambda=2)
fhlp<-fhl$pvalue # 0.015
fhlp


