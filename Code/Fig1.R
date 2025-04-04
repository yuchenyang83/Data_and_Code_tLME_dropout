load(paste0(PATH, "/Data/sourse/aids.RData"))

aids$death.time = aids$time
aids$time = aids$obstime
aids$Y = aids$CD4
aids$obstime = factor(aids$obstime)
N = length(unique(aids$id))

ggplot()+
  geom_boxplot(data = aids, aes(y = Y, x=time, group = time)) +
  facet_grid(.~drug) + 
  xlab("Month") + ylab(expression(sqrt("CD4"))) +
  scale_x_continuous(breaks = c(0,2,6,12,18)) + 
  theme(legend.position="none") +
  theme(text = element_text(size=20))

library(ggplot2)
kk = ggplot()+
  geom_line(data = aids, aes(x=time, y=Y, col=as.factor(drug), group=as.factor(id))) + 
  geom_point(data = aids, aes(x=time, y=Y, col=as.factor(drug), group=as.factor(id))) + 
  geom_boxplot(data = aids, aes(x=time, y=Y, group = time, fill = NULL), alpha = 0, 
               outlier.colour = "black", outlier.shape = 2, outlier.size = 3, outlier.alpha = 1, linewidth = 0.9) +
  facet_grid(.~drug) + 
  xlab("Month") + ylab(expression(sqrt("CD4"))) +
  scale_x_continuous(breaks = c(0,2,6,12,18)) + 
  theme(legend.position="none") +
  theme(text = element_text(size=20))
kk

ya6 = ggplot(data = aids, aes(x=as.factor(time), fill = drug)) +
  geom_bar(position = position_dodge()) + xlab("Month") + #ylab(str2expression(paste("Number","of"," sqrt(CD4)", sep = "~" ))) +
  ylab("Number of\n observed responses") +
  geom_text(stat='count', aes(label = paste0(round(..count../N, 3)*100, "%")), vjust=1.6, color="black", size=5, 
            position = position_dodge(0.9)) +
  scale_fill_discrete(name = "") + 
  theme(legend.position="top") +
  theme(text = element_text(size=20), axis.text.y.left = element_text(size=1))
ya6

source(paste0(PATH, "/function/multiplot.R"))
layout <- matrix(c(1,2), ncol = 1, byrow = TRUE)
multiplot(plotlist = list(kk, ya6), layout = layout)
