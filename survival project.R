library(FHtest)
library(KMsurv)
#### first we trying to use simulatino to show that log-rank test or Fleming-Harrington 
#### test with predifined (p=3,q=0; or p=0, q=3) to compare two group, then the 
#### type I error will be around 0.05  
#### we use exponential distribution (lambda = 0.5) to generate survival data
#### we use binomial distribution (p=0.5) we generate censering indicator
#### we generate 40 data data for each group 
#### we will run 5000 times for total simulation study
#### Please all these parameters can be chaged to better illustrate our points


# create a dataframe to store the p values from the test
# LR is for log-rank test
# FHE (early departure) is for Fleming-Harrington test with p=3, q=0
# FHL (late departure) is for Fleming-Harrington test with p=0, q=3
pvalues1<-data.frame(LR=numeric(), FHE=numeric(), FHL=numeric())

# write for loop for simulation
for (i in 1:5000){
  print (i) 
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df<-cbind(time=rweibull(160, shape=1.5, scale=1), status=rbinom(160, 1, 0.7), group=rep(c(0,1),each=80))
  # log-rank test
  lr<-survdiff(Surv(time, status)~group, data=as.data.frame(df))
  lrp<-1-pchisq(lr$chisq,1)
  # Fleming-Harrington test with p=3, q=0 (early departure)
  fhe<-FHtestrcc(Surv(time, status)~group, data=as.data.frame(df), rho=1, lambda=0)
  fhep<-fhe$pvalue
  # Fleming-Harrington test with p=0, q=3 (late departure)
  fhl<-FHtestrcc(Surv(time, status)~group, data=as.data.frame(df), rho=0, lambda=1)
  fhlp<-fhl$pvalue
  #add these newly generated p values to the dataframe pvalues
  temp<-data.frame(LR=lrp, FHE=fhep, FHL=fhlp)
  pvalues1<-rbind(pvalues1, temp)
}

# count the percentage of p vlues less than 0.05 in each of these three tests
sum(pvalues1$LR<0.05)/5000  
# 0.0532 (binom p=0.5, rho=lambda=3); 
# 0.0538 (binom p=0.7, rho=lambda=3)
# 0.0566 (binom p=0.7, rho=lambda=1)

sum(pvalues1$FHE<0.05)/5000 
# 0.051 (binom p=0.5, rho=lambda=3);  
# 0.0524 (binom p=0.7, rho=lambda=3)
# 0.0518 (binom p=0.7, rho=lambda=1)

sum(pvalues1$FHL<0.05)/5000 
# 0.0756 (binom p=0.5, rho=lambda=3);
# 0.05812 (binom p=0.7, rho=lambda=3)
# 0.0624 (binom p=0.7, rho=lambda=1)

# we can see that for all these three tests have type I error near 0.05
# note Fleming-Harrington test with p=0, q=3 have a slightly higher type I error

# next we try to use simulation to see the type I error in the following situation
# we will use fleming-harrington test, however, instead of determing p and q values 
# beforing looking at the data, we determine the p and q values after obsering data
# if we see early departure we use p=3, q=0; if we see late departure we use p=0, q=3
# Still same set up as above, survival data generated from exponential distribution
# with lambda be 0.5, generate censoring indicator with binomial distribution with
# p be 0.5. Each sample has 40 data points. 
# here to determine early or late departure, we directly perform two tests, one with p=3, q=0; 
# the other with p=0, q=3, then use the one with smaller p values as the final p value

# create a dataframe to store the p values from the test
pvalues3<-c()
# create a dataframe to store the choice between these two FH tests
# with 0 indicate early departure and 1 indicate late departure
#FHtest3<-c()

# write a loop for simulation
for (i in 1:5000){
  print (i)
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df<-cbind(time=rexp(80, 0.5), status=rbinom(80, 1, 0.7), group=rep(c(0,1),each=40))
  df<-as.data.frame(df)
  
  # fit the early departure model
  fitearly<-FHtestrcc(Surv(time, status)~group, data=df, rho=4, lambda=0)
  p1<-fitearly$pvalue
  
  # fit the late departure model
  fitlate<-FHtestrcc(Surv(time, status)~group, data=df, rho=0, lambda=4)
  p2<-fitlate$pvalue
  
  if (p1<=p2){
    #FHtest3<-c(FHtest3,0)
    pvalues3<-c(pvalues3,p1)
  }
  if (p2<p1){
    #FHtest3<-c(FHtest3,1)
    pvalues3<-c(pvalues3,p2)
  }
}

sum(pvalues3<0.05)/5000 
# 0.1104 (binom p=0.5, rho=lambda=3);
# 0.12 (binom p=0.7, rho=lambda=3)
# 0.0952 (binom p=0.7, rho=lambda=1)
# we can see that the type I error is bigger than 0.05

# generate table for the choice of two tests
table(FHtest3) 
# 2362 for early departure, 2638 for late departure (binom p=0.5, rho=lambda=3)
# 2326 for early departure, 2674 for late departure (binom p=0.7, rho=lambda=3)
# 2472 for early departure, 2528 for late departure (binom p=0.7, rho=lambda=1)
# percentage of p vlaues less than 0.05


#### next we try to see the power of these three test: log-rank test, FH test for early (p=3, q=0) 
# and late departure (p=0,q=3). Here we generate two survival dataset, both from exponential distribution,
# one with lambda as 0.5, and the other with lambda as 1. The censoring indicator was generating by 
# binomial distribution with p=0.7
# the percentage that p values less than 0.05 is treated as power

# create a vector to store p values
# LR is for log-rank test
# FHE (early departure) is for Fleming-Harrington test with p=3, q=0
# FHL (late departure) is for Fleming-Harrington test with p=0, q=3
pvalues4<-data.frame(LR=numeric(), FHE=numeric(), FHL=numeric())

# write for loop for simulation
for (i in 1:5000){
  print (i) 
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df1<-cbind(time=rweibull(80, shape=1.5, scale=1), status=rbinom(80, 1, 0.7), group=rep(0,80))
  df2<-cbind(time=rweibull(80, shape=1.5, scale=1.5), status=rbinom(80, 1, 0.7), group=rep(1,80))
  df<-as.data.frame(rbind(df1, df2))
  # log-rank test
  lr<-survdiff(Surv(time, status)~group, data=as.data.frame(df))
  lrp<-1-pchisq(lr$chisq,1)
  # Fleming-Harrington test with p=3, q=0 (early departure)
  fhe<-FHtestrcc(Surv(time, status)~group, data=as.data.frame(df), rho=1, lambda=0)
  fhep<-fhe$pvalue
  # Fleming-Harrington test with p=0, q=3 (late departure)
  fhl<-FHtestrcc(Surv(time, status)~group, data=as.data.frame(df), rho=0, lambda=1)
  fhlp<-fhl$pvalue
  #add these newly generated p values to the dataframe pvalues
  temp<-data.frame(LR=lrp, FHE=fhep, FHL=fhlp)
  pvalues4<-rbind(pvalues4, temp)
}

# count the percentage of p vlues less than 0.05 in each of these three tests
sum(pvalues4$LR<0.05)/5000  
# 0.6944 (binom p=0.7, rho=lambda=1)
# 0.7112 (binom p=0.7, rho=lambda=3)

sum(pvalues4$FHE<0.05)/5000 
# 0.6332 (binom p=0.7, rho=lambda=1)
# 0.4692 (binom p=0.7, rho=lambda=3)

sum(pvalues4$FHL<0.05)/5000 
# 0.5922 (binom p=0.7, rho=lambda=1)
# 0.4448 (binom p=0.7, rho=lambda=3)


### next we decide if we determine if we determine the p,q values after seeing the data what is the 
# power, we use the same set up as above 

# create a dataframe to store the p values from the test
pvalues5<-c()
# create a dataframe to store the choice between these two FH tests
# with 0 indicate early departure and 1 indicate late departure
FHtest5<-c()

for (i in 1:5000){
  print (i)
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df1<-cbind(time=rweibull(40, shape=1.5, scale=1), status=rbinom(40, 1, 0.7), group=rep(0,40))
  df2<-cbind(time=rweibull(40, shape=2, scale=1), status=rbinom(40, 1, 0.7), group=rep(1,40))
  df<-as.data.frame(rbind(df1, df2))
  
  # fit the early departure model
  fitearly<-FHtestrcc(Surv(time, status)~group, data=df, rho=3, lambda=0)
  p1<-fitearly$pvalue
  
  # fit the late departure model
  fitlate<-FHtestrcc(Surv(time, status)~group, data=df, rho=0, lambda=3)
  p2<-fitlate$pvalue
  
  if (p1<=p2){
    FHtest5<-c(FHtest5,0)
    pvalues5<-c(pvalues5,p1)
  }
  if (p2<p1){
    FHtest5<-c(FHtest5,1)
    pvalues5<-c(pvalues5,p2)
  }
}

# generate table for the choice of two tests
table(FHtest5) 
# 2571 for early departure, 2429 for late departure (binom p=0.7, rho=lambda=1)
# 2585 for early departure, 2415 for late departure (binom p=0.7, rho=lambda=3)


# percentage of p vlaues less than 0.05
sum(pvalues5<0.05)/5000 
# 0.7664 (binom p=0.7, rho=lambda=1)
# 0.6902 (binom p=0.7, rho=lambda=3)

#### we try the method to split data and use 20% of the data to determine p and q and then use the rest 80% 
# data to perform the test
# create a dataframe to store the p values from the test
# pvalues6_part to store p values if we only use the rest 80% data for real analysis
# pvalues6_all to store p values if we use all data for real analysis
pvalues6_part<-c()
pvalues6_all<-c()


# write a loop for simulation
for (i in 1:5000){
  print (i)
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df1<-cbind(time=rweibull(80, shape=1.5, scale=1), status=rbinom(80, 1, 0.7), group=rep(0,80))
  df2<-cbind(time=rweibull(80, shape=1.5, scale=1.5), status=rbinom(80, 1, 0.7), group=rep(1,80))
  df<-as.data.frame(rbind(df1, df2))
  # generate one vector to randomly select 20% of the data to determine p and q
  V_20<-sample(1:80, 16, replace=FALSE)
  # V_80 is the vector for the rest 80% of the data
  V_80<-(1:80)[-V_20]
  
  # select the 20% data in each group
  df1_20<-df1[V_20,]
  df2_20<-df2[V_20,]
  df_20<-as.data.frame(rbind(df1_20, df2_20))
  
  # generate the rest 80% data in each group
  df1_80<-df1[V_80,]
  df2_80<-df2[V_80,]
  df_80<-as.data.frame(rbind(df1_80, df2_80))
  
  # fit the early departure model
  fitearly<-FHtestrcc(Surv(time, status)~group, data=df_20, rho=1, lambda=0)
  p1<-fitearly$pvalue
  
  # fit the late departure model
  fitlate<-FHtestrcc(Surv(time, status)~group, data=df_20, rho=0, lambda=1)
  p2<-fitlate$pvalue
  
  # firse use rest 80% data
  if (p1<=p2){
    fit_80<-FHtestrcc(Surv(time, status)~group, data=df_80, rho=1, lambda=0)
  }
  if (p2<p1){
    fit_80<-FHtestrcc(Surv(time, status)~group, data=df_80, rho=0, lambda=1)
  }
  pvalues6_part<-c(pvalues6_part, fit_80$pvalue)
  
  # if we use all data
  if (p1<=p2){
    fit_all<-FHtestrcc(Surv(time, status)~group, data=df, rho=1, lambda=0)
  }
  if (p2<p1){
    fit_all<-FHtestrcc(Surv(time, status)~group, data=df, rho=0, lambda=1)
  }
  pvalues6_all<-c(pvalues6_all, fit_all$pvalue)
}

# percentage of p vlaue less than 0.05 (type I error)
sum(pvalues6_part<0.05)/5000
# 0.0626 (binom p=0.7, rho=lambda=1)
# 0.0674 (binom p=0.7, rho=lambda=3)

sum(pvalues6_all<0.05)/5000
# 0.0694 (binom p=0.7, rho=lambda=1)
# 0.0786 (binom p=0.7, rho=lambda=3)

#### we try the method to split data and use 20% of the data to determine p and q and then use the rest 80% 
# data or all data to perform the final test and determine the power
# we generate data from exponential distribution with lambda = 0.5 or 1 

# create a dataframe to store the p values from the test
# pvalues7_part to store p values if we only use the rest 80% data for real analysis
# pvalues7_all to store p values if we use all data for real analysis
pvalues7_part<-c()
pvalues7_all<-c()

# write a loop for simulation
for (i in 1:5000){
  print (i)
  # generate survival data, with 40 data in group 0 and the other 40 data in group 1
  df1<-cbind(time=rexp(80, 0.5), status=rbinom(80, 1, 0.7), group=rep(0,80))
  df2<-cbind(time=rexp(80, 1), status=rbinom(80, 1, 0.7), group=rep(1,80))
  df<-as.data.frame(rbind(df1, df2))
  # generate one vector to randomly select 20% of the data to determine p and q
  V_20<-sample(1:80, 16, replace=FALSE)
  # V_80 is the vector for the rest 80% of the data
  V_80<-(1:80)[-V_20]
  
  # select the 20% data in each group
  df1_20<-df1[V_20,]
  df2_20<-df2[V_20,]
  df_20<-as.data.frame(rbind(df1_20, df2_20))
  
  # generate the rest 80% data in each group
  df1_80<-df1[V_80,]
  df2_80<-df2[V_80,]
  df_80<-as.data.frame(rbind(df1_80, df2_80))
  
  # fit the early departure model
  fitearly<-FHtestrcc(Surv(time, status)~group, data=df_20, rho=1, lambda=0)
  p1<-fitearly$pvalue
  
  # fit the late departure model
  fitlate<-FHtestrcc(Surv(time, status)~group, data=df_20, rho=0, lambda=1)
  p2<-fitlate$pvalue
  
  # firse use rest 80% data
  if (p1<=p2){
    fit_80<-FHtestrcc(Surv(time, status)~group, data=df_80, rho=1, lambda=0)
  }
  if (p2<p1){
    fit_80<-FHtestrcc(Surv(time, status)~group, data=df_80, rho=0, lambda=1)
  }
  pvalues7_part<-c(pvalues7_part, fit_80$pvalue)
  
  # if we use all data
  if (p1<=p2){
    fit_all<-FHtestrcc(Surv(time, status)~group, data=df, rho=1, lambda=0)
  }
  if (p2<p1){
    fit_all<-FHtestrcc(Surv(time, status)~group, data=df, rho=0, lambda=1)
  }
  pvalues7_all<-c(pvalues7_all, fit_all$pvalue)
}

# percentage of p vlaue less than 0.05 (power)
sum(pvalues7_part<0.05)/5000
# 0.5152 (binom p=0.7, rho=lambda=1)
# 0.3864 (binom p=0.7, rho=lambda=3)

sum(pvalues7_all<0.05)/5000
# 0.628 (binom p=0.7, rho=lambda=1)
# 0.4884 (binom p=0.7, rho=lambda=3)