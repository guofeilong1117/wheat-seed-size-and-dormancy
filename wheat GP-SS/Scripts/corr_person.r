###correlation plot
library(corrplot)
td <-read.table("F:/wheat/phenotype/phe.txt.txt", header=T)
cor (td, method="pearson")
tdc <- cor (td, method="pearson")
corrplot(tdc)

tdc <- cor(td, method = "pearson")
png("F:/wheat/phenotype/correlation_plot.png", width = 800, height = 800, res = 100)
corrplot(tdc, 
         method = "color", 
         type = "upper",  
         order = "hclust", 
         tl.col = "black",
         tl.srt = 45,
         addCoef.col = "black", 
         diag = FALSE) 
dev.off()
