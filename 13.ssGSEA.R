library(GSVA)
library(limma)
library(GSEABase)
expFile=".txt"     
gmtFile="h.all.v7.5.1.symbols.gmt"   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat=avereps(mat)
mat=normalizeBetweenArrays(mat)
mat=mat[rowMeans(mat)>0,]

geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())

ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
	return((x-min(x))/(max(x)-min(x)))}
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut, file="ssgseaOut.txt", sep="\t", quote=F, col.names=F)



library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)


expFile="ssgseaOut.txt"      
typeFile="sample.txt"          
C="C"                    
P="P"                    
Ccol="blue"           
Pcol="orange"             
afmethod="t.test"     

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

type=read.table(typeFile, sep="\t", header=F, check.names=F, row.names=1)
colnames(type)="data"
sameSample=intersect(row.names(data),row.names(type))
rt1=cbind(data[sameSample,],type[sameSample,])
colnames(rt1)[ncol(rt1)]="data"
rt1=as.data.frame(rt1)
rt1[,1:(ncol(rt1)-1)]=lapply(rt1[,1:(ncol(rt1)-1)],as.numeric)
rt1=melt(rt1,id.vars=c("data"))
colnames(rt1)=c("data","Gene","Expression")
group=levels(factor(rt1$data))
rt1$data=factor(rt1$data, levels=c(C,P))
comp=combn(group,2)
my_comparisons=list()
for(j in 1:ncol(comp)){my_comparisons[[j]]<-comp[,j]}

boxplot=ggboxplot(rt1, x="Gene", y="Expression", fill ="data",
				  xlab="",
				  ylab="Score",
				  legend.title="Type",
				  width=0.8,
				  palette = c(Ccol,Pcol) )+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=data),
	method=afmethod,
	symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 0.2,1), symbols=c("***", "**", "*", "#","ns")), label="p.signif")+
  theme(axis.text=element_text(size=7,face = "bold"))+rotate()

pdf(file="h.all.diff.pdf",  width=7, height=12)
print(boxplot)
dev.off()



library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)


h.all="ssgseaOut.txt"     
exp=".txt"   
hub=".txt"    

rt=read.table(exp,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
immune=read.table(h.all, header=T, sep="\t", check.names=F, row.names=1)
immune=t(immune)
data=immune
risk=rt
af=read.table(file = hub,sep = "\t",header = F,check.names = F)
risk=risk[af[,1],]
sameSample=intersect(colnames(risk),rownames(data))
data=data[sameSample,]
risk=risk[,sameSample]
risk=t(risk)
data1=data
outTab=data.frame()
for(immune in colnames(data)){
	for(gene in colnames(risk)){
		x=as.numeric(data[,immune])
		y=as.numeric(risk[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}
outTab$cor=as.numeric(outTab$cor)
pdf(file="cor.pdf", width=8, height=15)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(high="blue", mid = "white",low= "orange") + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +   
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  
        axis.text.y = element_text(size = 8, face = "bold")) +      
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) + 
  scale_x_discrete(position = "bottom")     
dev.off()

