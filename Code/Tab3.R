#################
#################  25 %; nu=5
#################
PATH0 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t25/", sep = "")
PATH1 = paste(PATH0, "SIM1/", sep = "")
PATH2 = paste(PATH0, "SIM2/", sep = "")
PATH3 = paste(PATH0, "SIM3/", sep = "")
PATH4 = paste(PATH0, "SIM4/", sep = "")
PATH5 = paste(PATH0, "SIM5/", sep = "")

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))


MSE1 = t(MSE1)
sum.tabel = rbind(colMeans(t(MSE1))[-1][c(1,4)],
                  apply(MSE1, 1, sd)[-1][c(1,4)],
                  c(sum(apply(MSE1[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE1[-1,][c(1,4),], 2, order)[2,]==1)))

MSE2 = t(MSE2)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE2))[-1][c(1,4)],
                        apply(MSE2, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE2[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE2[-1,][c(1,4),], 2, order)[2,]==1))))

MSE3 = t(MSE3)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE3))[-1][c(1,4)],
                        apply(MSE3, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE3[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE3[-1,][c(1,4),], 2, order)[2,]==1))))

MSE4 = t(MSE4)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE4))[-1][c(1,4)],
                        apply(MSE4, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE4[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE4[-1,][c(1,4),], 2, order)[2,]==1))))

MSE5 = t(MSE5)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE5))[-1][c(1,4)],
                        apply(MSE5, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE5[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE5[-1,][c(1,4),], 2, order)[2,]==1))))

round(sum.tabel, 3)



for(i in c(3,6,9,12,15))
{
  if(sum(sum.tabel[i,][1:2])>100 ) sum.tabel[i,] = round((sum.tabel[i,]/sum(sum.tabel[i,][1:2])) * 100)
}

sum.tabel.TT = NULL
sum.tabel.TT = cbind(sum.tabel.TT, round(sum.tabel, 3))

#################
#################  50 %; nu=5
#################
PATH0 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t50/", sep = "")
PATH1 = paste(PATH0, "SIM1/", sep = "")
PATH2 = paste(PATH0, "SIM2/", sep = "")
PATH3 = paste(PATH0, "SIM3/", sep = "")
PATH4 = paste(PATH0, "SIM4/", sep = "")
PATH5 = paste(PATH0, "SIM5/", sep = "")

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))


MSE1 = t(MSE1)
sum.tabel = rbind(colMeans(t(MSE1))[-1][c(1,4)],
                  apply(MSE1, 1, sd)[-1][c(1,4)],
                  c(sum(apply(MSE1[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE1[-1,][c(1,4),], 2, order)[2,]==1)))

MSE2 = t(MSE2)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE2))[-1][c(1,4)],
                        apply(MSE2, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE2[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE2[-1,][c(1,4),], 2, order)[2,]==1))))

MSE3 = t(MSE3)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE3))[-1][c(1,4)],
                        apply(MSE3, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE3[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE3[-1,][c(1,4),], 2, order)[2,]==1))))

MSE4 = t(MSE4)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE4))[-1][c(1,4)],
                        apply(MSE4, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE4[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE4[-1,][c(1,4),], 2, order)[2,]==1))))

MSE5 = t(MSE5)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE5))[-1][c(1,4)],
                        apply(MSE5, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE5[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE5[-1,][c(1,4),], 2, order)[2,]==1))))

round(sum.tabel, 3)



for(i in c(3,6,9,12,15))
{
  if(sum(sum.tabel[i,][1:2])>100 ) sum.tabel[i,] = round((sum.tabel[i,]/sum(sum.tabel[i,][1:2])) * 100)
}


sum.tabel.TT = cbind(sum.tabel.TT, round(sum.tabel, 3))


#################
#################  75%; nu=5
#################
PATH0 = paste(PATH, "/Data/Simulation/SS-simulationSEM-t75/", sep = "")
PATH1 = paste(PATH0, "SIM1/", sep = "")
PATH2 = paste(PATH0, "SIM2/", sep = "")
PATH3 = paste(PATH0, "SIM3/", sep = "")
PATH4 = paste(PATH0, "SIM4/", sep = "")
PATH5 = paste(PATH0, "SIM5/", sep = "")

loglik1 =    as.matrix(read.table(paste(PATH1,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE1 =    as.matrix(read.table(paste(PATH1,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara1 = colMeans(as.matrix(read.table(paste(PATH1,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik2 =    as.matrix(read.table(paste(PATH2,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE2 =    as.matrix(read.table(paste(PATH2,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara2 = colMeans(as.matrix(read.table(paste(PATH2,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik3 =    as.matrix(read.table(paste(PATH3,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE3 =    as.matrix(read.table(paste(PATH3,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara3 = colMeans(as.matrix(read.table(paste(PATH3,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik4 =    as.matrix(read.table(paste(PATH4,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE4 =    as.matrix(read.table(paste(PATH4,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara4 = colMeans(as.matrix(read.table(paste(PATH4,'realpara.txt',sep=""),na.strings="NA",sep="")))

loglik5 =    as.matrix(read.table(paste(PATH5,'loglik.txt',sep=""),na.strings="NA",sep=""))
MSE5 =    as.matrix(read.table(paste(PATH5,'MSE.txt',sep=""),na.strings="NA",sep=""))
realpara5 = colMeans(as.matrix(read.table(paste(PATH5,'realpara.txt',sep=""),na.strings="NA",sep="")))


MSE1 = t(MSE1)
sum.tabel = rbind(colMeans(t(MSE1))[-1][c(1,4)],
                  apply(MSE1, 1, sd)[-1][c(1,4)],
                  c(sum(apply(MSE1[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE1[-1,][c(1,4),], 2, order)[2,]==1)))

MSE2 = t(MSE2)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE2))[-1][c(1,4)],
                        apply(MSE2, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE2[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE2[-1,][c(1,4),], 2, order)[2,]==1))))

MSE3 = t(MSE3)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE3))[-1][c(1,4)],
                        apply(MSE3, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE3[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE3[-1,][c(1,4),], 2, order)[2,]==1))))

MSE4 = t(MSE4)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE4))[-1][c(1,4)],
                        apply(MSE4, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE4[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE4[-1,][c(1,4),], 2, order)[2,]==1))))

MSE5 = t(MSE5)
sum.tabel = rbind(sum.tabel,
                  rbind(colMeans(t(MSE5))[-1][c(1,4)],
                        apply(MSE5, 1, sd)[-1][c(1,4)],
                        c(sum(apply(MSE5[-1,][c(1,4),], 2, order)[1,]==1), sum(apply(MSE5[-1,][c(1,4),], 2, order)[2,]==1))))

round(sum.tabel, 3)

sum.tabel.TT = cbind(sum.tabel.TT, round(sum.tabel, 3))

print(sum.tabel.TT[-c(2,5,8,11,14),c(2,1,4,3,6,5)])
