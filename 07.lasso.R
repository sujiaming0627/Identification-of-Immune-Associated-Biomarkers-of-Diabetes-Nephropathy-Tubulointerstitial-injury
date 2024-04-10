library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)
library(patchwork)
library(limma)

inputFile=".txt"     
C="C"                        


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


afcon=sum(sample[,1]==C)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])

set.seed(123)
cvfit = cv.glmnet(x, y,family = "binomial", nlambda=100, alpha=1,nfolds = 10) 


fit <- glmnet(x,y,family = "binomial")
cvfit$lambda.min

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(geneCoef, file="geneCoef.xls", sep="\t", quote=F, row.names=F)
write.table(file="lassoset.txt",lassoGene,sep="\t",quote=F,col.names=F,row.names=F) 

#######################################################
pdf("lasso.pdf",height = 5,width = 7)
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))   
#pdf("lambda.pdf")
plot(fit,xvar = 'lambda')
#dev.off()
#pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
#dev.off()
dev.off()
############################################################
