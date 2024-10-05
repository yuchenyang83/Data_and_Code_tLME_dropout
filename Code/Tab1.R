PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t25/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
names(realpara) = c(rep("Beta", 4), "sigma", rep("DD", 3), "nu", rep("alpha", 4))

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

MAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

MAR.se1	=	MAR.se1[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se1	=	MCAR.se1[,c(1:4, 8, 5:7)]
MNAR.se1	=	MNAR.se1[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se2	=	MAR.se2[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se2	=	MCAR.se2[,c(1:4, 8, 5:7)]
MNAR.se2	=	MNAR.se2[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se3	=	MAR.se3[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se3	=	MCAR.se3[,c(1:4, 8, 5:7)]
MNAR.se3	=	MNAR.se3[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se4	=	MAR.se4[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se4	=	MCAR.se4[,c(1:4, 8, 5:7)]
MNAR.se4	=	MNAR.se4[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se5	=	MAR.se5[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se5	=	MCAR.se5[,c(1:4, 8, 5:7)]
MNAR.se5	=	MNAR.se5[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]


sum.table = rbind(round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow=1), 3),
                  round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow=1), 3),
                  round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow=1), 3),
                  round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow=1), 3),
                  round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow=1), 3))

##########################################################################################
########################################################################
########################################################################  50%
PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t50/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
names(realpara) = c(rep("Beta", 4), "sigma", rep("DD", 3), "nu", rep("alpha", 4))

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

MAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

MAR.se1	=	MAR.se1[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se1	=	MCAR.se1[,c(1:4, 8, 5:7)]
MNAR.se1	=	MNAR.se1[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se2	=	MAR.se2[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se2	=	MCAR.se2[,c(1:4, 8, 5:7)]
MNAR.se2	=	MNAR.se2[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se3	=	MAR.se3[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se3	=	MCAR.se3[,c(1:4, 8, 5:7)]
MNAR.se3	=	MNAR.se3[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se4	=	MAR.se4[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se4	=	MCAR.se4[,c(1:4, 8, 5:7)]
MNAR.se4	=	MNAR.se4[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se5	=	MAR.se5[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se5	=	MCAR.se5[,c(1:4, 8, 5:7)]
MNAR.se5	=	MNAR.se5[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]


sum.table = rbind(sum.table,
                  rbind(round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow=1), 3)))


##########################################################################################
########################################################################
########################################################################  75%
PATH1 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t75/", sep = "")

realpara = colMeans(as.matrix(read.table(paste(PATH1,'SIM1/realpara.txt',sep=""),na.strings="NA",sep="")))
names(realpara) = c(rep("Beta", 4), "sigma", rep("DD", 3), "nu", rep("alpha", 4))

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

MAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se1 = as.matrix(read.table(paste(PATH1,'SIM1/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se2 = as.matrix(read.table(paste(PATH1,'SIM2/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se3 = as.matrix(read.table(paste(PATH1,'SIM3/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se4 = as.matrix(read.table(paste(PATH1,'SIM4/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MCAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MCAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]
MNAR.se5 = as.matrix(read.table(paste(PATH1,'SIM5/MNAR.se.est.txt',sep=""),na.strings="NA",sep=""))[,-1]

MAR.se1	=	MAR.se1[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se1	=	MCAR.se1[,c(1:4, 8, 5:7)]
MNAR.se1	=	MNAR.se1[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se2	=	MAR.se2[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se2	=	MCAR.se2[,c(1:4, 8, 5:7)]
MNAR.se2	=	MNAR.se2[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se3	=	MAR.se3[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se3	=	MCAR.se3[,c(1:4, 8, 5:7)]
MNAR.se3	=	MNAR.se3[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se4	=	MAR.se4[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se4	=	MCAR.se4[,c(1:4, 8, 5:7)]
MNAR.se4	=	MNAR.se4[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]
MAR.se5	=	MAR.se5[,c(1:4, 8, 5:7, 9:dim(MAR.se1)[2])]
MCAR.se5	=	MCAR.se5[,c(1:4, 8, 5:7)]
MNAR.se5	=	MNAR.se5[,c(1:4, 8, 5:7, 9:dim(MNAR.se1)[2])]


sum.table = rbind(sum.table,
                  rbind(round(colMeans(MNAR.se1, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est1, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se2, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est2, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se3, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est3, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se4, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est4, 2, sd))), nrow=1), 3),
                        round(colMeans(MNAR.se5, na.rm = T), 3), round(matrix(c(rbind(apply(MNAR.est5, 2, sd))), nrow=1), 3)))


print(sum.table)

