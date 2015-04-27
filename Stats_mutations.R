# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


rm(list=ls())
setwd("/data/annotate_sophia")

bigtable=read.table("myanno.hg19_multianno.txt",header=1,sep="\t",stringsAsFactors = F)
bigtable$X=NULL
bigtable$Gene.refGene=as.character(unlist(sapply(bigtable$Gene.refGene,function(x){unlist(strsplit(x,":"))[1]})))

library("GenomicRanges")
coordinates=GRanges(seqnames = paste("chr",bigtable$CHR,sep=""),ranges = IRanges(start = bigtable$POS,end = bigtable$POS+nchar(bigtable$REF)-1),strand = "*")
design2=read.table("/data//Gwendal//DESIGN HALOPLEX//AML2nd mai2013//15189-1368439012/15189-1368439012/15189-1368439012_Covered.bed",skip=2)
design2=GRanges(seqnames = design2$V1,ranges = IRanges(start = design2$V2,end = design2$V3),strand="*")
bigtable$design2=(countOverlaps(coordinates,design2)!=0)
design3=read.table("/data//Gwendal//DESIGN HALOPLEX/AML3rdJuly2014_v2/15189-1405499558_Covered.bed",skip=2)
design3=GRanges(seqnames = design3$V1,ranges = IRanges(start = design3$V2,end = design3$V3),strand="*")
bigtable$design3=(countOverlaps(coordinates,design3)!=0)

ADtable1=bigtable[,grep("_AD_alt",colnames(bigtable))]
ADtable2=bigtable[,grep("_AD_ref",colnames(bigtable))]

ADratios=as.matrix(ADtable1)/(as.matrix(ADtable1)+as.matrix(ADtable2))
plot(density(as.numeric(ADratios),na.rm=T),main="Distribution des fractions d'allèle alternatif",xlab="Fraction d'allèle alternatif",ylab="Densité de probabilité")
plot(density(apply(ADratios,1,max,na.rm=T)),main="Distribution des maxima de fréquence d'allèle alternatif",,xlab="Fraction maximale d'allèle alternatif",ylab="Densité de probabilité")

rep=read.csv("Pairs_replicates.tsv",sep=",",colClasses=c("numeric","numeric","character","character"))

MutationsPerSample=bigtable[,grep("GT",colnames(bigtable))]
#MutationsPerSample[ADratios<0.05]="None"

#Get the things that were called in the second version and not the third :
second=apply(MutationsPerSample[,rep[,1]]!="None"&MutationsPerSample[,rep[,2]]=="None",2,which)
#opposite :
third=apply(MutationsPerSample[,rep[,1]]=="None"&MutationsPerSample[,rep[,2]]!="None",2,which)

#We see that many more were wrongly called in the second version :
which(table(unlist(second))>2)
which(table(unlist(third))>2)

#Find out which genes these are in :
bigtable$Gene.refGene[as.numeric(names(which(table(unlist(second))>2)))]
bigtable$Gene.refGene[as.numeric(names(which(table(unlist(third))>2)))]

#Get number of mutations per sample and samples per mutation :
barplot(colSums(MutationsPerSample!="None"),las=2)
hist(colSums(MutationsPerSample!="None"),breaks=50)
barplot(rowSums(MutationsPerSample!="None"),las=2)
hist(rowSums(MutationsPerSample!="None"),breaks=50)

#Get number of mutations per gene and mutated genes per sample :
mutatedpersample=apply(MutationsPerSample,2,function(x){tapply(x!="None",as.factor(bigtable$Gene.refGene),sum)})
rowSums(mutatedpersample)
colSums(mutatedpersample)

write.csv(mutatedpersample,file = "Mutated_gene_per_sample.csv")

MutationsPerSample=cbind(bigtable$Gene.refGene,bigtable$REF,bigtable$ALT,bigtable$AAChange.refGene,bigtable$Func.refGene,bigtable$TYPE,bigtable$ExonicFunc.refGene,MutationsPerSample)

write.csv(MutationsPerSample,file="Mutations_per_sample.csv")

ADrep1=ADratios[,rep[,1]]
ADrep2=ADratios[,rep[,2]]

pdf("Allele_ratios_pairs.pdf")
for(i in 1:ncol(ADrep1))
{
  plot(ADrep1[,i],ADrep2[,i],
       xlab=paste("Sample",rep[i,1],"allele frequency"),ylab=paste("Sample",rep[i,2],"allele frequency"),
       main=paste("cor=",cor(ADrep1[,i],ADrep2[,i],use="pairwise.complete")),
       xlim=c(0,1),ylim=c(0,1))
  abline(0,1)
  #print(table(!is.na(ADrep1[,i]),!is.na(ADrep2[,i]!=0)))
}
dev.off()

ADratios_null=ADratios
ADratios_null[is.na(ADratios)]=0

ADrep1=ADratios_null[,rep[,1]]
ADrep2=ADratios_null[,rep[,2]]

index_changeprot=c(grep("nonsynonymous",bigtable$ExonicFunc.refGene),grep("frameshift",bigtable$ExonicFunc.refGene),grep("gain",bigtable$ExonicFunc.refGene))
isinchangeprot=1:nrow(bigtable)%in%index_changeprot
relevant=bigtable$design2&bigtable$design3

pdf("Allele_ratios_pairs_null.pdf")
for(i in 1:ncol(ADrep1))
{
  isnotnullboth=ADrep1[,i]>0|ADrep2[,i]>0
  plot(ADrep1[isnotnullboth,i],ADrep2[isnotnullboth,i],
       xlab=paste("Sample ",rep[i,1]," (",rep[i,3],")"," allele frequency",sep=""),ylab=paste("Sample ",rep[i,2]," (",rep[i,4],")"," allele frequency",sep=""),
       main=paste("cor=",format(cor(ADrep1[relevant&isnotnullboth,i],ADrep2[relevant&isnotnullboth,i],use="pairwise.complete"),digits = 2)),
       xlim=c(0,1),ylim=c(0,1),col=c(3,1,2)[bigtable$design2-bigtable$design3+2][isnotnullboth],pch=c(1,2)[isinchangeprot+1][isnotnullboth])
  abline(0,1)
  text(ADrep1[isnotnullboth,i],ADrep2[isnotnullboth,i],labels=bigtable$Gene.refGene[isnotnullboth],cex = 0.5,pos = 3)
  #print(table(!is.na(ADrep1[,i]),!is.na(ADrep2[,i]!=0)))
  write.table(bigtable[isnotnullboth,c(1,2,3,4,
                                       which(colnames(bigtable)=="Gene.refGene"),
                                       ncol(bigtable)-1,ncol(bigtable),
                                       as.numeric(rbind(grep(paste("Sample_",rep[i,1],"_",sep=""),colnames(bigtable)),
                                                        grep(paste("Sample_",rep[i,2],"_",sep=""),colnames(bigtable)))))],
              file=paste("comptable",rep[i,1],"-",rep[i,2],".txt",sep=""),
              quote=F,sep="\t",row.names=F)
  legend(x="top",y="top",legend = c("Non-synonymous","Synonymous","Present in both assays","Present only in assay 2","Present only in assay 3"),
         col = c(1,1,1,2,3),pch = c(2,1,NA,NA,NA),lty=c(NA,NA,1,1,1),lwd=c(1,1,1,1,1))
}
dev.off()



#Seems that those that do replicate are replicating quantitatively 
#but there are numerous ones which seem like they're not replicating at all despite high allele counts


smalltable=bigtable[index_changeprot,]

mutatedpersample=apply(MutationsPerSample[index_changeprot,-(1:7)],2,function(x){tapply(x!="None",as.factor(smalltable$Gene.refGene),sum)})
rowSums(mutatedpersample)
colSums(mutatedpersample)

write.csv(mutatedpersample,file = "Mutated_nonsynonymous_gene_per_sample.csv")
write.csv(MutationsPerSample[index_changeprot,],file="Mutations_nonsynonymous_per_sample.csv")


pdf("Allele_ratios_pairs_null_proteinchange.pdf")
for(i in 1:ncol(ADrep1))
{
  plot(ADrep1[index_changeprot,i],ADrep2[index_changeprot,i],
       xlab=paste("Sample",rep[i,1],"allele frequency"),ylab=paste("Sample",rep[i,2],"allele frequency"),
       main=paste("cor=",cor(ADrep1[index_changeprot,i],ADrep2[index_changeprot,i],use="pairwise.complete")),
       xlim=c(0,1),ylim=c(0,1),col=c(3,1,2)[bigtable$design2-bigtable$design3+2][index_changeprot])
  abline(0,1)
  #print(table(!is.na(ADrep1[,i]),!is.na(ADrep2[,i]!=0)))
}
dev.off()


#removing recurrent ones
recurrent=as.numeric(names(which(table(unlist(second))>=5)))
ADrep1=ADratios_null[-recurrent,rep[,1]]
ADrep2=ADratios_null[-recurrent,rep[,2]]

pdf("Allele_ratios_pairs_null_norecurrent.pdf")
for(i in 1:ncol(ADrep1))
{
  plot(ADrep1[,i],ADrep2[,i],
       xlab=paste("Sample",rep[i,1],"allele frequency"),ylab=paste("Sample",rep[i,2],"allele frequency"),
       main=paste("cor=",cor(ADrep1[,i],ADrep2[,i],use="pairwise.complete")),
       xlim=c(0,1),ylim=c(0,1))
  abline(0,1)
  #print(table(!is.na(ADrep1[,i]),!is.na(ADrep2[,i]!=0)))
}
dev.off()

