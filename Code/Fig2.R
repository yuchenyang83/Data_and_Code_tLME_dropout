set.seed(14)
source(paste0(PATH, '/function/simulation/funciton.r'))

n=100

# Parameter setting:
p = 4; si = 6; q = 2; g=1
Beta = matrix(c(5.35, -0.3, -0.1, -0.65), ncol=g) ## intercept; time; treat; time x treat
DD = matrix(0, q, q)
DD[1,1] = 0.4
DD[2,2] = 0.2
DD[2,1] = DD[1,2] = 0.1
sigma = 0.5
nu = 5
Phi = 1e-6

### 75 dropout
# 
alpha = c(-0.72, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
para = list(alpha=alpha, Beta=Beta, DD=DD, sigma=sigma, Phi=Phi, nu=nu)

cor.type = "UNC"
gen.Data = gen.tlmm(n, para, cor.type='UNC', si, q)
Data = gen.Data$Data
Ymat = gen.Data$Ymat

gen.Data.miss = add.MNAR.tLMM(gen.Data$Data, gen.Data$Ymat, si, para$alpha)
Data.miss = gen.Data.miss$Data.sub
Data.miss$R = rep(0, nrow(Data.miss))
Data.miss$R[is.na(Data.miss$Var1)] = 1
Data.miss$week = sqrt(Data.miss$Time)
X = gen.Data.miss$X
Z = gen.Data.miss$Z
Data = Data.miss


library(ggplot2)
k1 = ggplot(data = Data, aes(x=Time, y=Var1, col=as.factor(treat),shape = as.factor(treat), grop=as.factor(Subject)))+
  geom_line() + geom_point(size = 3) + xlab("Time") + ylab(NULL) + 
  ylab(expression(paste("Measurements of ", y[ij]^"o"))) +
  theme(legend.position="top")+ 
  scale_color_discrete(labels = c("Placebo", "Therapy")) +
  scale_shape_discrete(labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size=20)) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5, "cm"))+
  theme(legend.title=element_blank())


k2 = ggplot(data = Data, aes(x=as.factor(Time), y=Var1, col=as.factor(treat),shape = as.factor(treat), grop=as.factor(Time)))+
  geom_boxplot(na.rm = T) + xlab(NULL) + ylab(NULL) + 
  ylab(expression(paste("Measurements of ", y[ij]^"o"))) +
  theme(legend.position="top")+
  scale_color_discrete(labels = c("Placebo", "Therapy")) +
  scale_shape_discrete(labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size=20)) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5, "cm"))+
  theme(legend.title=element_blank())

k3 = ggplot(data = Data[-which(is.na(Data$Var1)),], aes(x=as.factor(Time), fill = as.factor(treat))) +
  geom_bar(position = position_dodge()) + xlab("Time") + ylab("Number of observed responses") +
  geom_text(stat='count', aes(label = paste0(round(..count../n, 3)*100, "%")), vjust=1.6, color="black", size=5, 
            position = position_dodge(0.9))+
  theme(legend.position="top")+
  scale_fill_discrete(name = NULL, labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size=20))

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(1,2,1,3), nrow = 2, byrow = TRUE)
multiplot(plotlist = list(k1, k2, k3), layout = layout)