rm(list=ls())
PATH = getwd()
load(paste0(PATH, "/Data/sourse/aids.RData"))
aids$death.time = aids$time
aids$time = aids$obstime
aids$Y = aids$CD4

n = length(table(aids$time))
N = length(unique(aids$id))
aids.mat = aids.time = aids.D = NULL
kk = 0
for(i in 1:N)
{
  Y = aids[aids$id==unique(aids$id)[i],]$Y
  time = aids[aids$id==unique(aids$id)[i],]$time
  drug = unique(aids[aids$id==unique(aids$id)[i],]$drug)
  aids.mat = rbind(aids.mat, c(Y, rep(NA, n - length(Y)), drug))
  aids.time = rbind(aids.time, c(time, rep(NA, n- length(time))))
  aids.D = rbind(aids.D, c(max(time), drug))
  if(max(aids[aids$id==unique(aids$id)[i],]$time)==6) kk = kk+1
}
aids.mat[is.na(aids.mat[,5]),]
aids.time[is.na(aids.time[,5]),]

table(aids.D[,2], aids.D[,1])
table(aids.D[,1])/N
N - cumsum(colSums(table(aids.D[,2], aids.D[,1])))
table(aids$time)/N
1-table(aids$time)/N
### ### ### ### ### ### 
### ### ### ### ### ###  trajectory plot
### ### ### ### ### ### 
id = unique(aids$id)
n = length(id)
nj = numeric(n)
for(i in 1: n) nj[i] = length(aids$time[aids$id == id[i]])
nj == table(aids$id)
aids$Subject = rep(1:n, nj)
aids$Dropout = rep(aids.D[,1], nj)

library(ggplot2)
ya1 = ggplot(data = aids, aes(x=time, y=Y, col=drug, group=as.factor(id))) + 
  geom_line(aes(group=as.factor(id))) + geom_point(aes(shape = drug)) + xlab("Time") + ylab("Score") + ggtitle("aids") +
  theme(legend.position="top") +
  theme(text = element_text(size=20))
ya1

ya3 = ggplot(data = aids, aes(x=time, y=Y, col=drug, group=as.factor(id))) + 
  geom_line() + geom_point(aes(shape = drug)) + xlab("Time") + ylab("Score") + ggtitle("aids") +
  scale_linetype_manual(values = c(2,1), breaks = c(0,1), labels = c("complete case", "dropout"), name = NULL) +
  facet_wrap(~drug, ncol = 4) +
  geom_line(data = aids[which(aids$Dropout != 18),], color = "black")  +
  geom_point(data = aids[which(aids$Dropout != 18),], color = "black", shape = 1)  +
  theme(legend.position="top") +
  theme(text = element_text(size=20))
ya3


######## ######## ######## 
######## ######## 
######## ######## 
library(nlme)
levels(aids$gender) = c(0,1)
levels(aids$prevOI) = c(0,1)
levels(aids$AZT) = c(0,1)
aids = data.frame(aids$Subject, aids$time, aids$Y, aids$drug, aids$gender, aids$prevOI, aids$AZT)
colnames(aids) = c("Subject", "week", "Var1", "drug", "gender", "prevOI", "AZT")
aids = groupedData(Var1 ~ week|Subject, data=aids)
aids$drug = as.character(aids$drug)
aids$drug[which(aids$drug=="ddI")] = 0
aids$drug[which(aids$drug=="ddC")] = 1
aids$drug = as.numeric(aids$drug)
aids$D = rep(aids.D[,1], nj)
aids$gender = as.character(aids$gender)
aids$gender = as.numeric(aids$gender)
aids$prevOI = as.character(aids$prevOI)
aids$prevOI = as.numeric(aids$prevOI)
aids$AZT = as.character(aids$AZT)
aids$AZT = as.numeric(aids$AZT)
aids$week = aids$week

table(aids$drug, aids$week)
table(aids$drug, aids$D)

cumsum.nj = cumsum(nj)

###########################################################################
###########################################################################
##############################
##############################  Selection model
##############################
##############################  load("D:/code_joineR_AIDS_113_0504/code_joineR_AIDS_113_0504/Untitled3.RData")
#############################
############################# contruct missing data
D.max = max(aids$week)
aids.miss = NULL
for(i in 1:n)
{
  if(unique(aids[aids$Subject==i,]$D) != D.max){
    aids.i = aids[aids$Subject==i,]
    if(aids.i$D[1]==0)
    {
      week = aids.i$week
      aids.i = rbind(aids.i, aids.i[dim(aids.i)[1],])
      aids.i$week = c(week, 2)
      aids.i$Var1[which(aids.i$week==2)] = NA
    }
    if(aids.i$D[1]==2)
    {
      week = aids.i$week
      aids.i = rbind(aids.i, aids.i[dim(aids.i)[1],])
      aids.i$week = c(week, 6)
      aids.i$Var1[which(aids.i$week==6)] = NA
    }
    if(aids.i$D[1]==6)
    {
      week = aids.i$week
      aids.i = rbind(aids.i, aids.i[dim(aids.i)[1],])
      aids.i$week = c(week, 12)
      aids.i$Var1[which(aids.i$week==12)] = NA
    }
    if(aids.i$D[1]==12)
    {
      week = aids.i$week
      aids.i = rbind(aids.i, aids.i[dim(aids.i)[1],])
      aids.i$week = c(week, 18)
      aids.i$Var1[which(aids.i$week==18)] = NA
    }
    aids.i = cbind(aids.i, miss = 1)
    aids.miss = rbind(aids.miss, aids.i)
  }else{
    aids.i = aids[aids$Subject==i,]
    aids.i = cbind(aids.i, miss = 0)
    aids.miss = rbind(aids.miss, aids.i)
  }
}

na.ind = is.na(as.vector(t(aids.miss$Var1)))
aids.miss[aids.miss$Subject==25,]
aids.miss[aids.miss$Subject==17,]

aids.miss = cbind(aids.miss, R = 0)
aids.miss$R[is.na(aids.miss$Var1)] = 1


aids.miss$Subject = as.numeric(as.character(aids.miss$Subject))
aids.miss$week = aids.miss$week
aids.miss$Var1 = aids.miss$Var1


cs1 <- corCAR1(0.2, form = ~ Time | Mare)
cs1AR1 <- corAR1(0.8, form = ~ 1 | Subject)
cs1AR1. <- Initialize(cs1AR1, data = Orthodont)
corMatrix(cs1AR1.)

cs1CAR1 <- corCAR1(0.8, form = ~ 1 | Subject)
cs1CAR1. <- Initialize(cs1CAR1, data = Orthodont)
corMatrix(cs1CAR1.)


Data = aids.miss
Data = groupedData(Var1 ~ week|Subject, data=Data)
fm1 = lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R!=1),], random  = ~ week|Subject, method = "ML")
summary(fm1)
intervals(fm1)
fm2 = lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R!=1),], random  = ~ week|Subject, correlation = corCAR1(), method = "ML")
summary(fm2)
intervals(fm2)

Data$Subject = as.numeric(as.character(Data$Subject))
class(Data)
Data = as.data.frame(Data)


#############################
#############################  Initialization
X = cbind(1, Data$week, Data$prevOI, Data$drug, Data$gender, Data$AZT)
# Z = as.matrix(X[,1])     # RI
Z = X[,c(1,2)]            # RIS
p = ncol(X)
q = ncol(Z)
g =1

Beta = matrix(c(fm2$coefficients$fixed), ncol=1)
DD = matrix(0, q,q)
DD[1,1] = as.numeric(VarCorr(fm2)[1,1])
DD[2,2] = as.numeric(VarCorr(fm2)[2,1])
DD[1,2] = DD[2,1] = as.numeric(VarCorr(fm2)[2,3])*sqrt(as.numeric(VarCorr(fm2)[1,1])*as.numeric(VarCorr(fm2)[2,1]))
sigma = c(fm2$sigma)
Phi = runif(g)
ga = 1
cor.type = "UNC"
tol = 1e-6
max.iter=100
per=100
M = 1000 ## Metropolis algorithm
M.LL = 500 ## importance sampling

## prediction of missing response 
source('D:/code_joineR_AIDS_113_0504/code_joineR_AIDS_113_0504/funciton.r')
alpha = as.numeric(prediction_ym(Data, fm2, X, Z, "ARp")$alpha.hat)
Data$yc = prediction_ym(Data, fm2, X, Z, "ARp")$yc


##### run MNAR
source('D:/code_joineR_AIDS_113_0504/code_joineR_AIDS_113_0504/LMMmissingSEM.r')
source('D:/code_joineR_AIDS_113_0504/code_joineR_AIDS_113_0504/tLMMmissingSEM.r')
# source('D:/code_joineR_AIDS_113_0504/code_joineR_AIDS_113_0504/LMMmissingSEM1.r')
############################################
###########  ARp
set.seed(1)
alpha[4] = 1
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha, nu=20)
mechanism='MNAR'
cor.type = "ARp"
est.MNAR.ARp = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 1000, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=per,
                             mechanism='MNAR')

