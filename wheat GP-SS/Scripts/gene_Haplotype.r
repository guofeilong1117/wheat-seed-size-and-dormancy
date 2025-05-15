###以TaMFT-3A的Haopotype分析
setwd("F:/wheat/GWAS-小麦/Haplotype")
library(geneHapR)
genotype <- import_vcf("TaMFT-3A.recode.vcf")
pheno <- read.csv("Phenotype.csv", row.names = 1)
hapResult <- vcf2hap(genotype,
                     hapPrefix = "H",
                     pad = 3,
                     hetero_remove =T,
                     na_drop = T)
hapResult
hapSummary <- hap_summary(hapResult)
plotHapTable(hapSummary)
write.csv(hapResult, "TaMFT-3A.csv")
results <-hapVsPheno(hapResult, hapPrefix = "H", mergeFigs = TRUE, pheno = pheno, phenoName = "GP3D", minAcc = 3)
plot(results$figs)
results <-hapVsPheno(hapResult, hapPrefix = "H", mergeFigs = TRUE, pheno = pheno, phenoName = "GP7D", minAcc = 3)
plot(results$figs)
results <-hapVsPheno(hapResult, hapPrefix = "H", mergeFigs = TRUE, pheno = pheno, phenoName = "TKW", minAcc = 3)
plot(results$figs)
results <-hapVsPheno(hapResult, hapPrefix = "H", mergeFigs = TRUE, pheno = pheno, phenoName = "KL", minAcc = 3)
plot(results$figs)
results <-hapVsPheno(hapResult, hapPrefix = "H", mergeFigs = TRUE, pheno = pheno, phenoName = "KW", minAcc = 3)
plot(results$figs)
q()
