element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
  }
  NextMethod()
}
############################################################################
############################################################################  sim 1
############################################################################
PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t25/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
MAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

names(realpara) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")
colnames(MNAR.est1) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")

realt = realpara
colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
value = c(apply((t(MAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est1)[1],
          apply((t(MAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est2)[1],
          apply((t(MAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est3)[1],
          apply((t(MAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est4)[1],
          apply((t(MAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est5)[1],
          apply((t(MCAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est1)[1],
          apply((t(MCAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est2)[1],
          apply((t(MCAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est3)[1],
          apply((t(MCAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est4)[1],
          apply((t(MCAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est5)[1],
          apply((t(MNAR.est1) - realt)^2, 1, sum)/ dim(MNAR.est1)[1],
          apply((t(MNAR.est2) - realt)^2, 1, sum)/ dim(MNAR.est2)[1],
          apply((t(MNAR.est3) - realt)^2, 1, sum)/ dim(MNAR.est3)[1],
          apply((t(MNAR.est4) - realt)^2, 1, sum)/ dim(MNAR.est4)[1],
          apply((t(MNAR.est5) - realt)^2, 1, sum)/ dim(MNAR.est5)[1]
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))

rep(gl(9, 1, labels = beta_labels), 10)
rep(gl(13, 1, labels = beta_labels2), 5)
ppc = data.frame(value = value,
                 sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 9), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 13), 1)),
                 beta = factor(c(as.character(rep(gl(9, 1, labels = beta_labels), 10)),
                                 as.character(rep(gl(13, 1, labels = beta_labels2), 5)))),
                 G = c(rep(c("MAR", "MCAR"), each = 9*5), rep(c("MNAR"), each = 13*5))
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))


ppc1.1 = rbind(cbind(ppc, model = "tLME"))



ppc1.1$sample = factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc1.1$beta = factor(ppc1.1$beta, levels = c("paste(beta[0])", "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
                                             "paste(sigma^2)",
                                             "paste(d[11])",     "paste(d[12])", "paste(d[21])", "paste(d[22])", "paste(nu)",
                                             "paste(alpha[0][0])",  "paste(alpha[0][1])", "paste(alpha[1])" , "paste(alpha[2])"))
ppc1.1$G = factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc1.1$model = factor(ppc1.1$model, levels = c("tLME"))
ppc1.1$GG = c(c(rep(c("MAR t 25%", "MCAR t 25%"), each = 9*5), rep(c("MNAR t 25%"), each = 13*5)))
ppc1.1 = cbind(ppc1.1, missing = "5% missing")

############################################################################
############################################################################  sim 1
############################################################################
PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t50/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
MAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

names(realpara) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")
colnames(MNAR.est1) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")

realt = realpara

colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
value = c(apply((t(MAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est1)[1],
          apply((t(MAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est2)[1],
          apply((t(MAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est3)[1],
          apply((t(MAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est4)[1],
          apply((t(MAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est5)[1],
          apply((t(MCAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est1)[1],
          apply((t(MCAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est2)[1],
          apply((t(MCAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est3)[1],
          apply((t(MCAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est4)[1],
          apply((t(MCAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est5)[1],
          apply((t(MNAR.est1) - realt)^2, 1, sum)/ dim(MNAR.est1)[1],
          apply((t(MNAR.est2) - realt)^2, 1, sum)/ dim(MNAR.est2)[1],
          apply((t(MNAR.est3) - realt)^2, 1, sum)/ dim(MNAR.est3)[1],
          apply((t(MNAR.est4) - realt)^2, 1, sum)/ dim(MNAR.est4)[1],
          apply((t(MNAR.est5) - realt)^2, 1, sum)/ dim(MNAR.est5)[1]
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))

rep(gl(9, 1, labels = beta_labels), 10)
rep(gl(13, 1, labels = beta_labels2), 5)
ppc = data.frame(value = value,
                 sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 9), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 13), 1)),
                 beta = factor(c(as.character(rep(gl(9, 1, labels = beta_labels), 10)),
                                 as.character(rep(gl(13, 1, labels = beta_labels2), 5)))),
                 G = c(rep(c("MAR", "MCAR"), each = 9*5), rep(c("MNAR"), each = 13*5))
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))


ppc2.1 = rbind(cbind(ppc, model = "tLME"))



ppc2.1$sample = factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc2.1$beta = factor(ppc1.1$beta, levels = c("paste(beta[0])", "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
                                             "paste(sigma^2)",
                                             "paste(d[11])",     "paste(d[12])", "paste(d[21])", "paste(d[22])", "paste(nu)",
                                             "paste(alpha[0][0])",  "paste(alpha[0][1])", "paste(alpha[1])"  , "paste(alpha[2])"))
ppc2.1$G = factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc2.1$model = factor(ppc1.1$model, levels = c("tLME"))
ppc2.1$GG = c(c(rep(c("MAR t 50%", "MCAR t 50%"), each = 9*5), rep(c("MNAR t 50%"), each = 13*5)))
ppc2.1 = cbind(ppc2.1, missing = "10% missing")


############################################################################
############################################################################  sim 3
############################################################################
PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t75/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
MAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.est5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.para.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

names(realpara) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")
colnames(MNAR.est1) = c("Beta0","Beta1","Beta2","Beta3", "sigma", "d11", "d12","d22", "nu", "alpha00", "alpha01","alpha1","alpha2")

realt = realpara

colMeans(MCAR.est2)
colMeans(MAR.est2)
colMeans(MNAR.est2)
value = c(apply((t(MAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est1)[1],
          apply((t(MAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est2)[1],
          apply((t(MAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est3)[1],
          apply((t(MAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est4)[1],
          apply((t(MAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MAR.est5)[1],
          apply((t(MCAR.est1[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est1)[1],
          apply((t(MCAR.est2[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est2)[1],
          apply((t(MCAR.est3[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est3)[1],
          apply((t(MCAR.est4[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est4)[1],
          apply((t(MCAR.est5[,1:9]) - realt[1:9])^2, 1, sum)/ dim(MCAR.est5)[1],
          apply((t(MNAR.est1) - realt)^2, 1, sum)/ dim(MNAR.est1)[1],
          apply((t(MNAR.est2) - realt)^2, 1, sum)/ dim(MNAR.est2)[1],
          apply((t(MNAR.est3) - realt)^2, 1, sum)/ dim(MNAR.est3)[1],
          apply((t(MNAR.est4) - realt)^2, 1, sum)/ dim(MNAR.est4)[1],
          apply((t(MNAR.est5) - realt)^2, 1, sum)/ dim(MNAR.est5)[1]
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))

rep(gl(9, 1, labels = beta_labels), 10)
rep(gl(13, 1, labels = beta_labels2), 5)
ppc = data.frame(value = value,
                 sample = c(rep(rep(c("25", "50", "100", "200", "400"), each = 9), 2), rep(rep(c("25", "50", "100", "200", "400"), each = 13), 1)),
                 beta = factor(c(as.character(rep(gl(9, 1, labels = beta_labels), 10)),
                                 as.character(rep(gl(13, 1, labels = beta_labels2), 5)))),
                 G = c(rep(c("MAR", "MCAR"), each = 9*5), rep(c("MNAR"), each = 13*5))
)

beta_labels = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])))

beta_labels2 = c(expression(paste(beta[0])), expression(paste(beta[1])), expression(paste(beta[2])), expression(paste(beta[3])),
                 expression(paste(sigma^2)), expression(paste(d[11])), expression(paste(d[12])), expression(paste(d[22])), expression(paste(nu)),
                 expression(paste(alpha[0][0])), expression(paste(alpha[0][1])), expression(paste(alpha[1])), expression(paste(alpha[2])))


ppc3.1 = rbind(cbind(ppc, model = "tLME"))



ppc3.1$sample = factor(ppc1.1$sample, c("25", "50", "100", "200", "400"))
ppc3.1$beta = factor(ppc1.1$beta, levels = c("paste(beta[0])", "paste(beta[1])", "paste(beta[2])", "paste(beta[3])",
                                             "paste(sigma^2)",
                                             "paste(d[11])",     "paste(d[12])", "paste(d[21])", "paste(d[22])", "paste(nu)",
                                             "paste(alpha[0][0])",  "paste(alpha[0][1])" , "paste(alpha[1])" , "paste(alpha[2])"))
ppc3.1$G = factor(ppc1.1$G, c("MCAR", "MAR", "MNAR"))
ppc3.1$model = factor(ppc1.1$model, levels = c("tLME"))
ppc3.1$GG = c(c(rep(c("MAR t 75%", "MCAR t 75%"), each = 9*5), rep(c("MNAR t 75%"), each = 13*5)))
ppc3.1 = cbind(ppc3.1, missing = "20% missing")

ppcTotal = rbind(ppc1.1, ppc2.1, ppc3.1)
ppcTotal = ppcTotal[ppcTotal$model=="tLME",]
ppcTotal = ppcTotal[ppcTotal$G=="MNAR",]
ppcTotal$missing = factor(ppcTotal$missing, levels = c("5% missing", "10% missing", "20% missing"))
levels(ppcTotal$missing) =  c("25% dropout", "50% dropout", "75% dropout")
ppcTotal$missing = factor(ppcTotal$missing)

ppcTotal = ppcTotal[-which(ppcTotal$beta=="paste(alpha[0][0])"),]

mse.sim2 = ggplot(data=ppcTotal, aes(x=sample, y=value, colour = missing, group = missing)) +
  geom_line(aes(linetype = missing), size=1) + geom_point(aes(shape=missing), size=4) +  xlab("sample size N") + ylab("MSE") +
  facet_wrap(~beta, labeller = "label_parsed",scales="free_y", ncol = 3, dir = "h") + 
  guides(shape = guide_legend(reverse=F,title=NULL),
         linetype = guide_legend(reverse=F,title=NULL),
         color = guide_legend(reverse=F,title=NULL)) +
  scale_shape_manual(values=c(2,0,1)) + 
  scale_linetype_manual(values=c(1, 2, 3)) + 
  theme(legend.position="top")+ 
  theme(text = element_text(size=20)) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit(1.5, "cm"))+
  ylab("MSE")  + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

print(mse.sim2)
