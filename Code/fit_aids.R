rm(list=ls())
load("/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/aids.RData")
# ddI = didanosine and ddC = zalcitabine.

# aids$drug[which(aids$drug=="ddI")] = "didanosine"
# aids$drug[which(aids$drug=='ddC')] = "zalcitabine"

aids$death.time = aids$time
aids$time = aids$obstime
aids$Y = aids$CD4

library(ggplot2)
kk = ggplot(data = aids, aes(x=time, y=Y, col=as.factor(drug), grop=as.factor(id)))+
  geom_line() + xlab("Time") + ylab("CD4") + ggtitle("AIDS") +
  theme(legend.position="top") +
  theme(text = element_text(size=20))
kk


## change to matrix
rm(list=ls())
load("/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/aids.RData")
aids$death.time = aids$time
aids$time = aids$obstime
aids$Y = aids$CD4

# n = max(aids$time)
n = 5
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

colSums(table(aids.D[,2], aids.D[,1]))/N

N - cumsum(colSums(table(aids.D[,2], aids.D[,1])))

### 
id = unique(aids$id)
n = length(id)
nj = numeric(n)
for(i in 1: n) nj[i] = length(aids$time[aids$id == id[i]])
nj == table(aids$id)
aids$Subject = rep(1:n, nj)

aids$Dropout = rep(aids.D[,1], nj)
# aids$drug = as.factor(aids$drug)
# aids$id = as.factor(aids$id)
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


########  fit model
library(nlme)
levels(aids$gender) = c(0,1)
levels(aids$prevOI) = c(0,1)
levels(aids$AZT) = c(0,1)
aids = data.frame(aids$Subject, aids$time, aids$Y, aids$drug, aids$gender, aids$prevOI, aids$AZT)
colnames(aids) = c("Subject", "week", "Var1", "drug", "gender", "prevOI", "AZT")
aids = groupedData(Var1 ~ week|Subject, data=aids)
aids$week = sqrt(aids$week)
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

table(aids$drug, aids$week)
table(aids$drug, aids$D)

# aids$drug1 = aids$drug2 = rep(NA, dim(aids)[1])
# for(i in 1:(dim(aids)[1]))
# {
#   if(aids$drug[i] == 1)
#   {
#     aids$drug1[i] = 1
#     aids$drug2[i] = 0
#   }
#   if(aids$drug[i] == 2)
#   {
#     aids$drug1[i] = 0
#     aids$drug2[i] = 0
#   }
#   if(aids$drug[i] == 3)
#   {
#     aids$drug1[i] = 0
#     aids$drug2[i] = 1
#   }
# }

# Zi = matrix(rep(1, Nj), ncol=q)
cumsum.nj = cumsum(nj)
g=1

### Initial values
fm1 = lme(Var1 ~ week + prevOI, data = aids, random  = ~ 1 + week, method = "ML")
summary(fm1)
acf(aids$Var1)
fm2 = lme(Var1 ~ week + prevOI, data = aids, random  = ~ 1, correlation = corAR1(), method = "ML")
intervals(fm2)

# source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/LMM.r')
# source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/function.micture.mix.model.r')
### UNC
# fit100 = LMM.AECM(Data, X, Z, g=1, init.para, cor.type = "UNC", tol = 1e-4, max.iter=1000, per=1)
# fit100$para.est$Beta # aids$week, aids$drug, aids$week * aids$drug
# fit100$model.inf$loglik


# source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/LMM.missing.r')
### UNC
# init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, alpha = c(-1.3460, -0.4689))
# fit100 = LMM.AECM(Data, X, Z, g=1, init.para, cor.type = "UNC", tol = 1e-3, max.iter=1000, per=1)
# fit100$para.est$Beta # aids$week, aids$drug, aids$week * aids$drug

###########################################################################
###########################################################################
##############################
##############################  Selection model
##############################
##############################  load("/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/Untitled3.RData")
#############################
############################# contruct missing data
aids$week = round(aids$week^2)
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

# which(aids.D==5)
# i=329
# aids.mat[which(aids.D==5),]
# aids.time[which(aids.D==5),]

na.ind = is.na(as.vector(t(aids.miss$Var1)))
# aids.miss[is.na(aids.miss$Var1),]

aids.miss[aids.miss$Subject==25,]
aids.miss[aids.miss$Subject==17,]
Data = aids.miss

aids.miss = cbind(aids.miss, R = 0)
aids.miss$R[is.na(aids.miss$Var1)] = 1


aids.miss$Subject = as.numeric(as.character(aids.miss$Subject))
aids.miss$week = aids.miss$week
aids.miss$Var1 = aids.miss$Var1

# nlme(R ~ week + week:drug + gender + prevOI + AZT, data=aids.miss)

Data = aids.miss
Data = groupedData(Var1 ~ week|Subject, data=Data)
fm2 = lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R!=1),], random  = ~ week|Subject, correlation = corAR1(), method = "ML")
summary(fm2)
intervals(fm2)

fm1 = lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R!=1),], random  = ~ week|Subject, method = "ML")
summary(fm1)
intervals(fm1, which = "fixed")
Data$Subject = as.numeric(as.character(Data$Subject))
class(Data)
Data = as.data.frame(Data)


#############################
############################# run missing model
X = cbind(1, Data$week, Data$prevOI, Data$drug, Data$gender, Data$AZT)
# Z = as.matrix(X[,1])     # RI
Z = X[,c(1,2)]            # RIS
p = ncol(X)
q = ncol(Z)

Beta = matrix(c(fm1$coefficients$fixed), ncol=1)
DD = matrix(0, q,q)
DD[1,1] = as.numeric(VarCorr(fm1)[1,1])
DD[2,2] = as.numeric(VarCorr(fm1)[2,1])
DD[1,2] = DD[2,1] = as.numeric(VarCorr(fm1)[2,3])*sqrt(as.numeric(VarCorr(fm1)[1,1])*as.numeric(VarCorr(fm1)[2,1]))
sigma = c(fm1$sigma)
alpha = c(-1.22,-0.05,-0.05,-0.05,-0.05, 0.038, -0.02)  ## Intercept, yij-1, yij
Phi = runif(g)
ga = 1
cor.type = "UNC"
tol = 1e-6
max.iter=1000
per=20
M = 1000
M.LL = 1000

V.fun = function(y, Data.miss)
{
  V=NULL
  for(i in 1:N)
  {
    if(i == 1) idx1 = 1: cumsum.ni[1]
    else idx1 = (cumsum.ni[i-1]+1): cumsum.ni[i]
    Data.i = Data.miss[which(Data.miss$Subject==i),]
    week = Data.i$week
    k = length(week)
    fData.i = y[1:(k-1)]
    fData.i = y[idx1][-k]
    fData.i = c(0, fData.i)
    bData.i = y[idx1]
    # V.i = cbind(1, Data.i$treat[1:k], fData.i, bData.i)
    V.i = cbind(1, Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
    V = rbind(V, V.i)
  }
  return(V)
}

fm2 = glm(R ~ prevOI, data = Data, family = binomial())
summary(fm2)
alpha = c(fm2$coefficients, 0.38, -0.2)

##### run MNAR
source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/funciton.r')
source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/LMMmissingSEM.r')
source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/tLMMmissingSEM.r')
# source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/LMMmissingSEM1.r')
############################################
###########  UNC
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha, nu=20)
mechanism='MNAR'
# init.para = fit.MNAR$para.est
# init.para = list(Beta=fit.MNAR$para.est$Beta, DD=fit.MNAR$para.est$D, sigma=fit.MNAR$para.est$sigma, Phi=Phi, alpha=fit.MNAR$para.est$alpha)
fit.MNAR = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                        mechanism='MNAR')

#### run MAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha = alpha[-length(alpha)])
fit.MAR = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                       mechanism='MAR')

#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=c(-5.3))
fit.MCAR = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                        mechanism='MCAR')


fit.MCAR$model.inf$loglik; fit.MAR$model.inf$loglik; fit.MNAR$model.inf$loglik
# -2463.935 -2463.769 -2466.402
fit.MCAR$model.inf$aic; fit.MAR$model.inf$aic; fit.MNAR$model.inf$aic
fit.MCAR$model.inf$bic; fit.MAR$model.inf$bic; fit.MNAR$model.inf$bic

############################################
###########  ARp
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha)
mechanism='MNAR'
cor.type = "ARp"
# init.para = fit.MNAR$para.est
# init.para = list(Beta=fit.MNAR$para.est$Beta, DD=fit.MNAR$para.est$D, sigma=fit.MNAR$para.est$sigma, Phi=Phi, alpha=fit.MNAR$para.est$alpha)
fit.MNAR.ARp = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                            mechanism='MNAR')

#### run MAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha = alpha[-length(alpha)])
fit.MAR.ARp = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                           mechanism='MAR')

#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=c(-5.3))
fit.MCAR.ARp = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                            mechanism='MCAR')

############################################
###########  "BAND1"
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0.2, ga = ga, alpha=alpha, nu=20)
mechanism='MNAR'
cor.type = c("BAND1")
fit.MNAR.BAND1 = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                              mechanism='MNAR')

#### run MAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0, ga = ga, alpha = alpha[-length(alpha)])
fit.MAR.BAND1 = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                             mechanism='MAR')
# fit.MAR.treat = LMM.AECM.miss(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, P=1, tol = 1e-2, max.iter=max.iter, per=per)

#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0, ga = ga, alpha=c(-2.15), nu=20)
fit.MCAR.BAND1 = LMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                              mechanism='MCAR')

########################################################################################
##### run Student t
########################################################################################
##### run MNAR
source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/funciton.r')
source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/tLMMmissingSEM.r')
# source('/Users/yangyucheng/Desktop/Boss/NLME-missing/Run_Data/code_joineR_AIDS_112_1221_2/code_joineR_AIDS_112_1221_2/LMMmissingSEM1.r')
############################################
###########  UNC
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha, nu=20)
mechanism='MNAR'
est.MNAR = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=1,
                         mechanism='MNAR')

set.seed(20240128)
#### run MAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha[-length(alpha)], nu=20)
est.MAR = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=per,
                        mechanism='MAR')
# fit.MAR.drug = LMM.AECM.miss(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, P=1, tol = 1e-2, max.iter=max.iter, per=1)

set.seed(20240119)
#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=c(-5.3), nu=20)
est.MCAR = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=per,
                         mechanism='MCAR')


############################################
###########  ARp
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha, nu=20)
mechanism='MNAR'
cor.type = "ARp"
est.MNAR.ARp = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=1,
                             mechanism='MNAR')

set.seed(20240103)
#### run MAR
mechanism='MNAR'
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=alpha[-length(alpha)], nu=20)
est.MAR.ARp = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=per,
                            mechanism='MAR')
# fit.MAR.drug = LMM.AECM.miss(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, P=1, tol = 1e-2, max.iter=max.iter, per=1)

set.seed(20240157)
#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, ga=ga, alpha=c(-5.3), nu=20)
est.MCAR.ARp = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = tol, max.iter=max.iter, per=per,
                             mechanism='MCAR')

############################################
###########  "BAND1"
set.seed(20240101)
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0.2, ga = ga, alpha=alpha, nu=20)
mechanism='MNAR'
cor.type = c("BAND1")
est.MNAR.BAND1 = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                               mechanism='MNAR')

#### run MAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0.2, ga = ga, alpha=alpha[-length(alpha)], nu=20)
est.MAR.BAND1 = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                              mechanism='MAR')
# est.MAR.treat = LMM.AECM.miss(Data, X, Z, V, g, init.para, cor.type = c("UNC"), M = 100, P=1, tol = 1e-2, max.iter=max.iter, per=per)

#### run MCAR
init.para = list(Beta=Beta, DD=DD, sigma=sigma, Phi=0.2, ga = ga, alpha=c(-2.15), nu=20)
est.MCAR.BAND1 = tLMM.miss.SEM(Data, X, Z, V, g, init.para, cor.type = c("BAND1"), M = 100, M.LL = M.LL, P=1, tol = tol, max.iter=max.iter, per=per,
                               mechanism='MCAR')


c(fit.MNAR$model.inf$bic, fit.MAR$model.inf$bic, fit.MCAR$model.inf$bic)
c(fit.MNAR.ARp$model.inf$bic, fit.MAR.ARp$model.inf$bic, fit.MCAR.ARp$model.inf$bic)
c(fit.MNAR.DEC$model.inf$bic, fit.MAR.DEC$model.inf$bic, fit.MCAR.DEC$model.inf$bic)
c(fit.MNAR.CS$model.inf$bic, fit.MAR.CS$model.inf$bic, fit.MCAR.CS$model.inf$bic)
c(fit.MNAR.BAND1$model.inf$bic, fit.MAR.BAND1$model.inf$bic, fit.MCAR.BAND1$model.inf$bic)


c(est.MNAR$model.inf$bic, est.MAR$model.inf$bic, est.MCAR$model.inf$bic)
c(est.MNAR.ARp$model.inf$bic, est.MAR.ARp$model.inf$bic, est.MCAR.ARp$model.inf$bic)
c(est.MNAR.DEC$model.inf$bic, est.MAR.DEC$model.inf$bic, est.MCAR.DEC$model.inf$bic)
c(est.MNAR.CS$model.inf$bic, est.MAR.CS$model.inf$bic, est.MCAR.CS$model.inf$bic)
c(est.MNAR.BAND1$model.inf$bic, est.MAR.BAND1$model.inf$bic, est.MCAR.BAND1$model.inf$bic)



