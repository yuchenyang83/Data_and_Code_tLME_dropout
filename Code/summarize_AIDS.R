################################################################################
#  Filename   : summarize_AIDS.R
#   Project   : BiomJ article "Extending t linear mixed models for longitudinal 
#               data with non-ignorable dropout applied to AIDS studies"                                                           
#  Authors    : Yu-Chen Yang, Wan-Lun Wang, Luis M. Castro, Tsung-I Lin
#  Date       : 19.10.2025
#  Purpose    : Reproduce the Section 2 descriptive summaries for the AIDS dataset:
#               (1) total number of patients;
#               (2) male/female counts;
#               (3) treatment proportions (didanosine, zalcitabine);
#               (4) baseline AIDS status (proportion diagnosed at study entry);
#               (5) number of dropouts and dropout percentage.
#
#  PATH       : Typically points to the Data_and_Code folder.
#               Example (Windows): "D:/Data_and_Code"
#
#  Input file : Data_and_Code/Data/source/aids.RData
#               (expected object with variables: id, obstime, gender, drug, prevOI, …)
#
#  Output file: results/summarize_AIDS.csv
#
#
#  R Version  : R-4.3.1  
################################################################################
PATH <- getwd()
load(paste0(PATH, "/Data/source/aids.RData"))

N <- length(unique(aids$id))
N

aids
table(aids[aids$obstime == 0, ]$gender)

round(round(table(aids[aids$obstime == 0, ]$drug)[1] / N * 100, 2), 1)
round(round(table(aids[aids$obstime == 0, ]$drug)[2] / N * 100, 2), 1)

round(table(aids[aids$obstime == 0, ]$prevOI) / N * 100, 1)

sum(aids$obstime == 18)
N - sum(aids$obstime == 18)
round((N - sum(aids$obstime == 18)) / N * 100, 1)

sum.tabel.1 <- matrix(c("total", NA, N, NA), 2, byrow = T)
sum.tabel.2 <- matrix(c("male", "female", rev(table(aids[aids$obstime == 0, ]$gender))), 2, byrow = T)
sum.tabel.3 <- matrix(c("didanosine", "zalcitabine", round(round(table(aids[aids$obstime == 0, ]$drug)[2] / N * 100, 2), 1), round(round(table(aids[aids$obstime == 0, ]$drug)[1] / N * 100, 2), 1)), 2, byrow = T)
sum.tabel.4 <- matrix(c("AIDS", "noAIDS", rev(round(table(aids[aids$obstime == 0, ]$prevOI) / N * 100, 1))), 2, byrow = T)
sum.tabel.5 <- matrix(c("Number of dropouts", NA, N - sum(aids$obstime == 18), NA), 2, byrow = T)
sum.tabel.6 <- matrix(c("Dropout rate", NA, round((N - sum(aids$obstime == 18)) / N * 100, 1), NA), 2, byrow = T)

sum.tabel.T <- rbind(sum.tabel.1, sum.tabel.2, sum.tabel.3, sum.tabel.4, sum.tabel.5, sum.tabel.6)

write.csv(sum.tabel.T, paste0(PATH, "/Result/summarize_AIDS.csv"))
