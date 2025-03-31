### ### ### ### ### ### ### ### ### 
kk = as.data.frame(read.table(paste0(PATH, "/Data/fixed_alpha.txt"),na.strings="NA",sep=""))


beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), 
                expression(paste(beta[3])), expression(paste(beta[4])), expression(paste(beta[5])))

kk$beta = c(as.character(gl(6, 1, labels = beta_labels))) 
kk$beta = as.factor(kk$beta)


ya1 = ggplot(kk, aes(x=alpha, y=est.beta, color = "black")) + 
  geom_point(size = 2) + geom_line(size = 1) +
  facet_wrap(.~beta, labeller = "label_parsed", scales = "free") + 
  xlab(expression(paste(alpha[2]))) + 
  geom_line(data = kk, aes(x=alpha, y=est.upper, color = "red"), linetype = 2, size = 1) +
  geom_line(data = kk, aes(x=alpha, y=est.lower, color = "red"), linetype = 2, size = 1) +
  scale_x_continuous(breaks = c(-16, -8, -4, 0, 4, 8, 16)) +
  theme(legend.key.width = unit(2, "cm")) +
  guides(shape = guide_legend(reverse=F,title=NULL),
         linetype = guide_legend(reverse=F,title=NULL),
         color = guide_legend(reverse=F,title=NULL),
         fill = guide_legend(reverse=F,title=NULL)) +
  theme(legend.position="top") +
  scale_y_continuous(name = "Estimated fixed effects") +
  scale_color_manual(values = c("black", "red"),labels = c("Estimated value", "95% Confidence bound")) + 
  theme(text = element_text(size=25), 
        strip.text.x = element_text(size = 25), 
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 25)) + 
  theme(axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(), 
        axis.title.y.right = element_text(color = "red", vjust = 2))
print(ya1)

