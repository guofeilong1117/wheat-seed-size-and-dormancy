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
plotHapTable(hapSummary)#绘制单倍型表格图
write.csv(hapResult, "TaMFT-3A.csv")
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "GP3D",
                     minAcc = 3)
plot(results$figs)#绘制SNP与表型关联图
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "GP7D",
                     minAcc = 3)
plot(results$figs)#绘制SNP与表型关联图
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "TKW",
                     minAcc = 3)
plot(results$figs)#绘制SNP与表型关联图
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "KL",
                     minAcc = 3)
plot(results$figs)#绘制SNP与表型关联图
results <-hapVsPheno(hapResult,
                     hapPrefix = "H",
                     mergeFigs = TRUE,
                     pheno = pheno,
                     phenoName = "KW",
                     minAcc = 3)
plot(results$figs)#绘制SNP与表型关联图
