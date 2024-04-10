library(limma)
library(ConsensusClusterPlus)

rt=read.table("matrix.txt", header=T, sep="\t", check.names=F, row.names=1)
rt=as.matrix(rt)
workDir=getwd()

maxK=9
results=ConsensusClusterPlus(rt,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="km",
              distance="euclidean",
              seed=123456,
              plot="png")
calcICL(results, title="consensusScore", plot="png")


Kvec = 2:maxK
x1 = 0.1; x2 = 0.9                     # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="")          # from 2 to maxK
for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}                                 

optK = Kvec[which.min(PAC)]

clusterNum=optK        
cluster=results[[clusterNum]][["consensusClass"]]
write.table(cluster, file="cluster.txt", sep="\t", quote=F, col.names=F)

### by af  ###
##307686155@qq.com###






#install.packages("ggpubr")



library(reshape2)
library(ggpubr)
inputFile="input.txt"      
outFile="vioplot.pdf"     
setwd("C:\\")    


rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)
x=colnames(rt)[1]
colnames(rt)[1]="Type"


data=melt(rt,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")


p=ggviolin(data, x="Gene", y="Expression", color = "Type", 
	     ylab="Gene expression",
	     xlab=x,
	     legend.title=x,
	     add.params = list(fill="white"),
	     palette = c("blue","orange"),
	     width=1, add = "boxplot")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")


pdf(file=outFile, width=6, height=5)
print(p1)
dev.off()


