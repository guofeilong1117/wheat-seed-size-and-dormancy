###计算Fst
vcftools --gzvcf wheat-c.vcf.gz --weir-fst-pop landrace.txt --weir-fst-pop cultivar.txt --out Lan_Cul --fst-window-size 200000 --fst-window-step 100000
sed "1d" Lan_Cul.windowed.weir.fst|awk '{if($5<0)print $1"\t"$2"\t0";else print $1"\t"$2"\t"$5}' > Lan_Cul_fst.txt

####Fst可视化
rm(list=ls())
library(qqman)
library(Cairo)
pifile<-read.table("F:/wheat/Lan_Cul_fst.txt", header=F, stringsAsFactors=F)
SNP<-paste(pifile[,1],pifile[,2],sep = ":")
pifile=cbind(SNP,pifile)
colnames(pifile)<-c("SNP","CHR","POS","fst")
outfile<-"fst_wheat"
filePNG<-paste(outfile,"manhattan.png",sep=".")
CairoPNG(file=filePNG,width=1000,height=500)
colorset<-c("#e64a37","#2E8B57","#6495ED")
manhattan(pifile,chr="CHR",bp="POS",p="fst",snp="SNP",col=colorset,logp=F,suggestiveline = F,genomewideline = F,ylab="pi",ylim=c(0,1),font.lab=4,cex.lab=1.2,main="LvsC",cex=0.8)
dev.off()

####Fst_top5%_gene
sort -k 5 -g Lan_Cul.windowed.weir.fst > Lan_Cul.sorted.fst
#窗口计数
wc -l Lan_Cul.sorted.fst #假设为20000
#取前5%
tail -n 1000 Lan_Cul.sorted.fst >Lan_Cul.sorted.5%.fst
cut -f 1-3 Lan_Cul.sorted.5%.fst | sort -k1,1n -k2,2n > Fst-5%-1.pos
bedtools intersect -a Fst-5%-1.pos -b wheat.gff3 -wb > Fst-5%-gene.txt


###计算Pi
vcftools --gzvcf wheat-c.vcf.gz --keep landrace.txt --window-pi 200000 --window-pi-step 100000 --out landrace_Pi
vcftools --gzvcf wheat-c.vcf.gz --keep cultivar.txt --window-pi 200000 --window-pi-step 100000 --out cultivar_Pi
###计算Pi_radio
sed -n '5p' landrace_pi.windowed.pi > Lan_Pi.txt
sed -n '5p' Cultivar_pi.windowed.pi > Cul_Pi.txt
awk 'FNR==NR {a=$1; next} {print $1/a}' Cul_Pi.txt Lan_Pi.txt > Lan_Cul_Pi.txt
awk 'BEGIN{OFS="\t"} {print $1, $2}' landrace_pi.windowed.pi > pos_Pi.txt
paste pos_Pi.txt Lan_Cul_Pi.txt > Pi_radio.txt

####pi可视化
rm(list=ls())
library(qqman)
library(Cairo)
pifile<-read.table("F:/wheat/Pi_radio.txt", header=F, stringsAsFactors=F)
SNP<-paste(pifile[,1],pifile[,2],sep = ":")
pifile=cbind(SNP,pifile)
colnames(pifile)<-c("SNP","CHR","POS","pi")
outfile<-"Pi_radio_wheat"
filePNG<-paste(outfile,"manhattan.png",sep=".")
CairoPNG(file=filePNG,width=1500,height=500)
colorset<-c("#e64a37","#2E8B57","#6495ED")
manhattan(pifile,chr="CHR",bp="POS",p="pi",snp="SNP",col=colorset,logp=F,suggestiveline = F,genomewideline = F,ylab="pi ratio",ylim=c(0,100),font.lab=4,cex.lab=1.2,main="LvsC",cex=0.8)
dev.off()

####Pi_top5%_gene
sort -k 5 -g Pi_radio.txt > Lan_Cul.sorted.pi
#窗口计数
wc -l Lan_Cul.sorted.Pi #假设为20000
#取前5%
tail -n 1000 Lan_Cul.sorted.Pi >Lan_Cul.sorted.5%.pi
cut -f 1-3 Lan_Cul.sorted.5%.Pi | sort -k1,1n -k2,2n > Pi-5%-1.pos
bedtools intersect -a Pi-5%-1.pos -b wheat.gff3 -wb > Pi-5%-gene.txt

###vcftools
Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, et al. The variant call format and VCFtools. Bioinformatics, 2011, 27(15), 2156–2158.
###bedtools
Quinlan AR & Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010, 26, 6, pp. 841–842.
