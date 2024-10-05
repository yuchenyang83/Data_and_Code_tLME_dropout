### ### ### ### ### ### ### ### ### 
load(paste0(PATH, "/Data/fit_result.RData"))
Data.t = as.data.frame(cbind(rbind(est.MNAR.ARp$para.est$Beta, fit.MNAR.ARp$para.est$Beta), rep(c("t", "n"), each = 6)))


### ### ### ### ### ### ### ### ### 
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha1.RData"))
alpha[4]
fit.beta = est.beta = NULL
fix_alpha = alpha[4]
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))

load(paste0(PATH, "/Data/fixed_alpha/fix_alpha2.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha3.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha4.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha5.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha6.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha7.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha8.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha9.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM2$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha10.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM2$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha13.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha14.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha15.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha16.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha17.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha18.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha19.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha20.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM2$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha21.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))
load(paste0(PATH, "/Data/fixed_alpha/fix_alpha22.RData"))
fix_alpha = c(fix_alpha, alpha[4])
est.beta = rbind(est.beta, t(est.MNAR$IM$out[,1:6]))



beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), 
                expression(paste(beta[3])), expression(paste(beta[4])), expression(paste(beta[5])))

kk = data.frame(est.beta = est.beta[,1],
                est.upper = est.beta[,1] + 1.96*est.beta[,2],
                est.lower = est.beta[,1] - 1.96*est.beta[,2],
                alpha = rep(rep(fix_alpha, each = 6)),
                beta = c(as.character(gl(6, 1, labels = beta_labels))))
kk$beta = as.factor(kk$beta)


Data.t$beta = c(as.character(gl(6, 1, labels = beta_labels)))
Data.t$V1 = as.numeric(Data.t$V1)

ya1 = ggplot(kk, aes(x=alpha, y=est.beta)) + 
  geom_line(data = kk, aes(x=alpha, y=est.upper, color = "red"), linetype = 2, size = 1) +
  geom_line(data = kk, aes(x=alpha, y=est.lower, color = "red"), linetype = 2, size = 1) +
  # geom_hline(data = Data.t[which(Data.t$V2=="t"),], aes(yintercept = V1), color= "red", size = 1.5) +
  geom_point(size = 2) + geom_line(size = 1) +
  facet_wrap(.~beta, labeller = "label_parsed", scales = "free") + 
  xlab(expression(paste(alpha[2]))) + 
  scale_linetype_manual(values = c(2,1), labels = c("LME", "tLME")) +
  scale_x_continuous(breaks = c(-16, -8, -4, 0, 4, 8, 16)) +
  theme(legend.key.width = unit(2, "cm")) +
  guides(shape = guide_legend(reverse=F,title=NULL),
         linetype = guide_legend(reverse=F,title=NULL),
         color = guide_legend(reverse=F,title=NULL),
         fill = guide_legend(reverse=F,title=NULL)) +
  theme(legend.position="none") +
  scale_y_continuous(name = "Estimated fixed effects")+
  theme(text = element_text(size=25), 
        strip.text.x = element_text(size = 25), 
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 25)) + 
  theme(axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_text(color = "red", vjust = 2))
print(ya1)


