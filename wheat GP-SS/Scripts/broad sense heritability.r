###Calculating broad-sense heritability
library(lme4)
library(lsmeans)
dat <- read.table("F:/wheat/phenotype/phe.txt",header = T)
trait <- phe$Phe
year <- phe$Year
line <- phe$Line
loc <- phe$Location
rep <- phe$Rep
TRAIT <- as.numeric(trait)
LINE <- as.factor(line)
LOC <- as.factor(loc)
YEAR <- as.factor(year)
REP <- as.factor(rep)
## Modeling
blup <-lmer(TRAIT~(1|LINE)+(1|LOC)+(1|YEAR)+(1|REP%in%LOC:YEAR)+(1|LINE:LOC)+(1|LINE:YEAR))
summary(blup)

###Calculating h2
h2 = (σG2) / (σG2 + (σGL2 ) / L + (σGY 2 ) / Y +(σE 2) / YRL)
h2
